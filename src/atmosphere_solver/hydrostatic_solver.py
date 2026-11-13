from abc import ABC, abstractmethod
from dataclasses import dataclass
import numpy as np
from periodictable import formula
from scipy.integrate import cumulative_trapezoid
from src.common.phys_constants import kb
from src.common.units import Q_, get_magnitude
from src.parallelization.array_backend import xp
from src.exo_chem.easy_chem_driver import run_easy_chem_full_profile
from src.utilities.logging_module import log
from src.atmosphere_solver.atmosph_data import AtmLayerDyn

@dataclass
class PressureDensitySolverSettings:
    max_iter: int
    abs_tol: float
    rel_tol: float
    logp_tol: float
    damping: float
    min_pressure: float
    anderson_depth: int

@dataclass
class PressureDensitySolveResult:
    layers: AtmLayerDyn
    converged: bool
    iterations: int
    max_abs_pressure_delta: float
    max_rel_pressure_delta: float
    max_log_pressure_delta: float

#
#   hydrostatic main driver
#

def hydrostatic_driver(
        solver,
        *,
        state,
        atmosph_mass,
        planet_radius,
        initialize=False,
        initialization_method=None,
        equilibrium_update=None
):
    """Dispatch initialization or a fixed-composition hydrostatic solve."""
    if initialize and initialization_method == "easychem":
        if equilibrium_update is None:
            log.error("EasyChem initialization requires equilibrium_update")
        return solver.solve_scf(
            altitude=state["altitude"],
            temperature=state["temperature"],
            gravity=state["gravity"],
            mole_frac=state["mole_frac"],
            atmosph_mass=atmosph_mass,
            planet_radius=planet_radius,
            equilibrium_update=equilibrium_update
        )
    # Profile initialization, prescribed initialization, or later evolution:
    # composition is supplied and remains fixed during this solve.
    return solver.solve_single_iter(
        altitude=state["altitude"],
        temperature=state["temperature"],
        gravity=state["gravity"],
        mole_frac=state["mole_frac"],
        atmosph_mass=atmosph_mass,
        planet_radius=planet_radius
    )

#
#   set solver
#

def set_hydro_solver(atmosphere_data, units, **solver_kwargs):
    hydrostatic_data = atmosphere_data.get("hydrostatic_solver", {}) or {}
    solver_type = hydrostatic_data.get("type", "standard")
    settings_data = hydrostatic_data.get("settings", {}) or {}
    #  set up settings
    settings = PressureDensitySolverSettings(
        max_iter=int(
            settings_data["max_iter"]
            if settings_data.get("max_iter") is not None else 50
        ),
        abs_tol=get_magnitude(
            value=settings_data.get("abs_tol"),
            units=units.get("pressure"),
            default_value=Q_(0.0, units.get("pressure"))
        ),
        rel_tol=float(
            settings_data["rel_tol"]
            if settings_data.get("rel_tol") is not None else 1.0e-6
        ),
        logp_tol=float(
            settings_data["logp_tol"] 
            if settings_data.get("logp_tol") is not None
            else 1.0e-6
        ),
        damping=float(
            settings_data["damping"]
            if settings_data.get("damping") is not None else 0.5
        ),
        min_pressure=get_magnitude(
            value=settings_data.get("min_pressure"),
            units=units.get("pressure"),
            default_value=Q_(1e-12, units.get("pressure"))
        ),
        anderson_depth=int(
            settings_data.get("anderson_depth")
            if settings_data.get("anderson_depth") is not None
            else 5
        )
    )
    #  select solver
    solver_classes = {
        "standard": StandardPressureDensitySolver,
        "anderson": AndersonPressureDensitySolver,
    }
    try:
        solver_class = solver_classes[solver_type]
    except KeyError:
        log.error(
            f"unknown hydrostatic solver type '{solver_type}'. "
            f"Valid options: {sorted(solver_classes)}"
        )
    return solver_class(
        units=units, 
        settings=settings,
        **solver_kwargs
    )

#
#   abstract solver class
#

class FixedTemperaturePressureDensitySolver(ABC):
    """
    Base solver for fixed-temperature pressure/density equilibrium.

    This class intentionally knows only about the pressure-density fixed point.
    The parent atmosphere solver should provide callbacks for hydrostatic
    pressure, mean molecular mass, and species number density so this module
    stays separate from radiative transfer and atmosphere setup details.
    """
    solver_name = "fixed-temperature pressure-density"
    def __init__(
        self,
        *,
        units: dict,
        settings: PressureDensitySolverSettings
    ):
        # Internal numerical values use the units defined in this dictionary.
        self._units = units
        self.settings = settings
        self.chem_species = None
        self.chem_masses = None
        self._gas_indices = None
        # Boltzmann constant: numerical magnitude in internal units.
        self._kB = float(
            kb.to(
                f"{self._units['energy']} / {self._units['temperature']}"
            ).magnitude
        )
    # surface pressure
    def _set_base_pressure(
            self,
            *,
            altitude,
            temperature,
            mean_molecular_mass,
            pressure_shape,
            atmosph_mass,
            planet_radius
    ):
        """Determine lowest-grid pressure from the atmospheric mass.

        All values use internal units.
        pressure_shape is dimensionless, with pressure_shape[0] = 1.
        atmospheric_mass is the mass within the supplied altitude domain.
        """
        z = xp.asarray(altitude, dtype=float)
        temperature = xp.asarray(temperature, dtype=float)
        mean_mass = xp.asarray(mean_molecular_mass, dtype=float)
        fz = xp.asarray(pressure_shape, dtype=float)
        # Density per unit boundary pressure: rho(z) / P_boundary
        rho_Pb = (
            fz * mean_mass / (self._kB * temperature)
        )
        # Spherical shell mass per unit boundary pressure
        integrand = (
            4.0 * xp.pi
            * (planet_radius + z) ** 2
            * rho_Pb
        )
        mass_Pb = xp.trapezoid(
            integrand,
            x=z
        )
        return atmosph_mass / float(mass_Pb.item())
    # set chemical species
    def set_chem_species(self, species):
        self.chem_species = tuple(species)
        # precompute masses
        masses = [
            float(
                Q_(
                    formula(
                        entity.symbol
                        if hasattr(entity, "symbol")
                        else entity.sequence
                    ).mass,
                    "amu"
                ).to(self._units["mass"]).magnitude
            )
            for entity in self.chem_species
        ]
        self.chem_masses = xp.asarray(masses, dtype=float)
        self._gas_indices = xp.asarray(
            [
                index
                for index, entity in enumerate(self.chem_species)
                if entity.phase == "gas"
            ],
            dtype=int
        )
    # single iteration hydrostatic solver
    def solve_single_iter(
            self,
            *,
            altitude,
            temperature,
            gravity,
            mole_frac,
            atmosph_mass,
            planet_radius
    ):
        """Solve hydrostatics with fixed temperature and composition.

        All numerical inputs use internal units.
        atmosph_mass is the mass within the supplied altitude domain.
        """
        z = xp.asarray(altitude, dtype=float)
        temperature = xp.asarray(temperature, dtype=float)
        gravity = xp.asarray(gravity, dtype=float)
        mole_frac = xp.asarray(mole_frac, dtype=float)
        # mean gas particle mass
        mole_gas = mole_frac[:, self._gas_indices]
        gas_masses = self.chem_masses[self._gas_indices]
        mean_mass = mole_gas @ gas_masses
        # integrate pressure
        # Integrate d(log P)/dz = -mean_mass * g / (kB * T).
        dlogP_dz = mean_mass * gravity / (self._kB * temperature)
        DlogP = (
            0.5 * (dlogP_dz[1:] + dlogP_dz[:-1]) * xp.diff(z)
        )
        integral = xp.concatenate(
            (xp.zeros(1, dtype=float), xp.cumsum(DlogP))
        )
        pressure_shape = xp.exp(-integral)
        # normalized pressure profile
        # based on total atmospheric mass
        P_0 = self._set_base_pressure(
            altitude=z,
            temperature=temperature,
            mean_molecular_mass=mean_mass,
            pressure_shape=pressure_shape,
            atmosph_mass=atmosph_mass,
            planet_radius=planet_radius
        )
        # set pressure
        pressure = P_0 * pressure_shape
        # number density
        number_density = pressure / (self._kB * temperature)
        # mass density
        mass_density = number_density * mean_mass
        # Preserve complete species ordering; gas concentrations fill gas columns.
        
        exit()



        return AtmLayerDyn(
            altitude=z,
            pressure=pressure,
            temperature=temperature,
            mean_molecular_mass=mean_mass,
            gravity=gravity,
            species_number_density=species_number_density,
            mass_density=mass_density,
            species=self.chem_species,
            mole_frac=mole_frac
        )
    # solve hydrostatic equations
    def solve_scf(self, config: AtmLayerDyn) -> PressureDensitySolveResult:
        """Solve hydrostatics at fixed temperature and species fractions.

        All inputs are numerical magnitudes in internal units.
        species_masses maps species names to particle masses.
        The pressure boundary at the lowest layer is self._P_surf.
        """
        z, temperature, pressure, gravity = self._unpack_config(config)
        log_pressure = np.log(np.maximum(pressure, self.settings.min_pressure))
        # convergence
        converged = False
        iterations = 0
        metrics = {
            "max_abs": np.inf,
            "max_rel": np.inf,
            "max_dlogp": np.inf,
        }
        # starting
        self._start()
        exit()
        for iteration in range(self.settings.max_iter):
            pressure_old = np.exp(log_pressure)
            hydro_pressure, chem_data, mu = self._fixed_point_pressure(
                z=z,
                temperature=temperature,
                gravity=gravity,
                pressure=pressure_old,
            )
            residual = np.log(hydro_pressure) - log_pressure
            metrics = self._pressure_metrics(
                pressure_old=pressure_old,
                pressure_new=hydro_pressure,
                residual=residual,
            )

            iterations = iteration + 1
            self._log_iteration(iterations, metrics)

            if self._has_converged(metrics):
                converged = True
                log_pressure = np.log(hydro_pressure)
                log.info(f"{self.solver_name} converged in {iterations} iterations")
                break

            log_pressure = self._next_log_pressure(
                iteration=iteration,
                log_pressure=log_pressure,
                residual=residual,
            )

        if not converged:
            log.warning(f"{self.solver_name} did not converge after {self.settings.max_iter} iterations")

        pressure = np.maximum(np.exp(log_pressure), self.settings.min_pressure)
        layers = self._final_layers(
            z=z,
            temperature=temperature,
            pressure=pressure,
            gravity=gravity,
        )

        return PressureDensitySolveResult(
            layers=layers,
            converged=converged,
            iterations=iterations,
            max_abs_pressure_delta=metrics["max_abs"],
            max_rel_pressure_delta=metrics["max_rel"],
            max_log_pressure_delta=metrics["max_dlogp"],
        )
    # set initial density
    def _initial_density_profile(self, z):
        # density profile
        density = self._initial_density*np.exp(-z/self._initial_scale_height)
        # compute total mass
        shell_area = 4.0*np.pi*(self.planet_data.planet_radius.to(self._units.get('length')).magnitude + z)**2
        density_mass = np.trapezoid(density*shell_area, z)
        # renormalize mass
        return density*self._atmosph_mass/density_mass
    # mean molecular weight profile
    def _mean_molecular_weight_profile(self, chem_data, z, tol=1.e-5):
        n_layers = len(z)
        # avg. molecular weight
        mu_ofz = Q_(np.zeros(n_layers, dtype=float), "amu")
        moles_fraction_sum = np.zeros_like(mu_ofz)
        # run over species / mole_fraction array
        for species, mole_fraction in chem_data.mole_fraction_profiles.items():
            molecular_weight = self._molecular_weight(species)
            xi_ofz = np.asarray(mole_fraction, dtype=float)
            #print(xi_ofz, molecular_weight)
            mu_ofz += xi_ofz*molecular_weight
            moles_fraction_sum += xi_ofz
        if np.any(moles_fraction_sum < 1.0-tol):
            log.warning("Sum_i x_i(z) < 1")
        return mu_ofz.to(self._units["mass"]).magnitude
    # molecular weight
    def _molecular_weight(self, species):
        try:
            return Q_(formula(species).mass, "amu")
        except Exception as exc:
            log.error(f"cannot compute molecular weight for species {species}: {exc}")
    # species number density
    def _species_number_density(self, density, chem_data, mu):
        # total number density
        total_number_density = density / mu
        # number density of species
        species_number_density = {}
        for species, mole_fraction in chem_data.mole_fraction_profiles.items():
            xi_ofz = np.asarray(mole_fraction, dtype=float)
            species_number_density[species] = xi_ofz*total_number_density
        return species_number_density
    # solve hydrostatic pressure profile
    def _solve_hydrostatic_pressure(self, z, temperature, gravity, density=None, mu=None):
        # if chemical potential provided
        if mu is not None:
            dlogp_dz_abs = mu*gravity/(self._kB*temperature)
            logp_integral = cumulative_trapezoid(dlogp_dz_abs, z, initial=0.0)
            return self._P_surf*np.exp(-logp_integral)
        # else use density
        if density is None:
            log.error("density profile is required when mu is not provided")
        column_weight = cumulative_trapezoid(density*gravity, z, initial=0.0)
        return self._P_surf - column_weight
    # extract data from config
    def _unpack_config(self, config: AtmLayerDyn):
        z = config.altitude
        temperature = config.temperature
        pressure = config.pressure
        gravity = config.gravity
        # check variables state
        if np.any(temperature <= 0.0):
            log.error("temperature must be positive for pressure-density equilibrium solve")
        if np.any(pressure <= 0.0):
            log.error("pressure must be positive for pressure-density equilibrium solve")
        return z, temperature, pressure, gravity
    def _start(self):
        log.warning(
            f"Starting {self.solver_name}: maxiter={self.settings.max_iter}, "
            f"abs_tol={self.settings.abs_tol}, rel_tol={self.settings.rel_tol}, "
            f"logp_tol={self.settings.logp_tol}"
        )
    
    def _fixed_point_pressure(self, *, z, temperature, gravity, pressure):
        # Numerical arrays on the selected NumPy/CuPy backend.
        z = xp.asarray(z, dtype=float)
        temperature = xp.asarray(temperature, dtype=float)
        gravity = xp.asarray(gravity, dtype=float)
        pressure=xp.asarray(pressure, dtype=float)
        # EasyChem needs CPU arrays with units attached.
        chem_data = run_easy_chem_full_profile(
            pressure=Q_(pressure, self.pressure_unit),
            temperature=Q_(temperature, self.temperature_unit),
            atomic_abund=self.atomic_abundances,
            chemical_species=self.chemical_species,
        )
        mu = self._mean_molecular_weight_profile(chem_data, z)
        if np.any(mu <= 0.0):
            log.error("mean molecular mass must be positive for pressure-density equilibrium solve")
        
        hydro_pressure = self._solve_hydrostatic_pressure(
            z=z,
            temperature=temperature,
            gravity=gravity,
            mu=mu,
        )
        hydro_pressure = np.maximum(hydro_pressure, self.settings.min_pressure)
        return hydro_pressure, chem_data, mu
    
    def _pressure_metrics(self, *, pressure_old, pressure_new, residual):
        pressure_delta = np.abs(pressure_new - pressure_old)
        pressure_scale = np.maximum(np.maximum(np.abs(pressure_new), np.abs(pressure_old)), 1e-30)
        return {
            "max_abs": float(np.max(pressure_delta)),
            "max_rel": float(np.max(pressure_delta / pressure_scale)),
            "max_dlogp": float(np.max(np.abs(residual))),
        }
        
    def _has_converged(self, metrics: dict[str, float]) -> bool:
        return (
            metrics["max_dlogp"] <= self.settings.logp_tol
            or metrics["max_abs"] <= self.settings.abs_tol
            or metrics["max_rel"] <= self.settings.rel_tol
        )
        
    def _log_iteration(self, iteration: int, metrics: dict[str, float]):
        log.warning(
            f"{self.solver_name} iter {iteration}: "
            f"max_abs_dP={metrics['max_abs']:.3e} Pa, "
            f"max_rel_dP={metrics['max_rel']:.3e}, "
            f"max_dlogP={metrics['max_dlogp']:.3e}"
        )
    
    def _final_layers(self, *, z, temperature, pressure, gravity) -> AtmLayerDyn:
        chem_data = run_easy_chem_full_profile(
            pressure=Q_(pressure, self._units.get("pressure")),
            temperature=Q_(temperature, self._units.get("temperature")),
            atomic_abund=self.atomic_abundances,
            chemical_species=self.chemical_species,
        )
        mu = self._mean_molecular_weight_profile(chem_data, z)
        density = pressure * mu / (self._kB * temperature)
        species_number_density = self._species_number_density(
            density=density,
            chem_data=chem_data,
            mu=mu,
        )
        return AtmLayerDyn(
            altitude=z,
            pressure=pressure,
            temperature=temperature,
            mean_molecular_mass=mu,
            gravity=gravity,
            species_number_density=species_number_density,
            chemistry=chem_data,
        )
    @abstractmethod
    def _next_log_pressure(self, *, iteration: int, log_pressure: np.ndarray, residual: np.ndarray) -> np.ndarray:
        raise NotImplementedError

#
#    Standard pressure / density solver
#

class StandardPressureDensitySolver(FixedTemperaturePressureDensitySolver):
    solver_name = "standard fixed-temperature pressure-density solver"
    #
    # compute chemical equilibrium + P configuration -> fixed T
    #
    def _compute_equilibrium_atmospheric_configuration(self, config: AtmLayerDyn):
        # set input profiles
        z = np.asarray(config.altitude, dtype=float)
        temperature = np.asarray(config.temperature, dtype=float)
        pressure = np.asarray(config.pressure, dtype=float)
        gravity = np.asarray(config.gravity, dtype=float)
        # solve hydrostatic equilibrium at given T
        pressure = np.maximum(pressure, self._min_pressure.to(self.pressure_unit).magnitude)
        # log pressure
        log_pressure = np.log(pressure)
                
        # Adaptive damping: start conservative, adapt based on convergence
        damping = max(0.3, min(0.9, self._damping_loop))  # Clamp to reasonable range [0.3, 0.9]
        damping_min = 0.1  # Never go below this
        damping_max = 0.95  # Never go above this
        prev_log_pressure_update = None
        oscillation_count = 0
        log.warning(
            "Starting adaptive damped log-pressure equilibrium solver: "
            f"maxiter={self._max_iter_loop}, abs_tol={self._abs_tol}, "
            f"rel_tol={self._rel_tol}, logp_tol={self._logp_tol}, "
            f"initial_damping={damping:.3f}"
        )
        converged = False
        # start loop iterations
        for iteration in range(self._max_iter_loop):
            pressure_old = np.exp(log_pressure)
            # compute chemical structure
            chem_data = run_easy_chem_full_profile(
                pressure=Q_(pressure_old, self.pressure_unit),
                temperature=Q_(temperature, self.temperature_unit),
                atomic_abund=self.chem_input.atomic_abundances,
                chemical_species=self.chem_input.chemical_species
            )
            # average mass
            mu = self._mean_molecular_weight_profile(chem_data, z)
            if np.any(mu <= 0.0):
                log.error("mean molecular mass must be positive for equilibrium atmosphere solve")
            # hydrostatic pressure
            hydro_pressure = self._solve_hydrostatic_pressure(
                z=z,
                temperature=temperature,
                gravity=gravity,
                mu=mu
            )
            hydro_pressure = np.maximum(hydro_pressure, self._min_pressure.to(self.pressure_unit).magnitude)
            # pressure variation
            pressure_delta = np.abs(hydro_pressure - pressure_old)
            # compute pressure scale
            pressure_scale = np.maximum(np.abs(hydro_pressure), np.abs(pressure_old))
            safe_scale = np.maximum(pressure_scale, 1e-30)
            max_abs = np.max(pressure_delta)
            max_rel = np.max(pressure_delta / safe_scale)
            max_dlogp = np.max(np.abs(np.log(hydro_pressure) - log_pressure))
                    
            # Detect oscillations: check if update direction reverses
            log_pressure_update = np.log(hydro_pressure) - log_pressure
            oscillating = False
            if prev_log_pressure_update is not None:
                # Check if signs flip (indicating oscillation)
                sign_flip = np.sum(np.sign(prev_log_pressure_update) != np.sign(log_pressure_update + 1e-14))
                if sign_flip > len(log_pressure_update) * 0.5:  # More than 50% flip
                    oscillating = True
                    oscillation_count += 1
            prev_log_pressure_update = log_pressure_update.copy()
                    
            # Adapt damping based on convergence behavior
            if oscillating and damping < damping_max:
                # Reduce damping if oscillating
                damping = min(damping * 1.5, damping_max)
                log.info(f"  Oscillation detected: increase damping to {damping:.3f}")
            elif not oscillating and oscillation_count == 0 and damping > damping_min:
                # Increase damping slightly if smooth convergence and no prior oscillations
                damping = max(damping * 0.98, damping_min)
                
            # print info
            log.warning(
                f"equilibrium iter {iteration+1}: max_abs_dP={max_abs:.3e} Pa, "
                f"max_rel_dP={max_rel:.3e}, max_dlogP={max_dlogp:.3e}, "
                f"damping={damping:.3f}, osc={'Y' if oscillating else 'N'}, "
                f"abs_tol={self._abs_tol}, rel_tol={self._rel_tol}, logp_tol={self._logp_tol}"
            )
            if max_dlogp <= self._logp_tol:
                converged = True
                pressure = hydro_pressure
                log.info(f"equilibrium atmosphere converged in {iteration + 1} iterations")
                break
            # update log pressure with adaptive damping
            log_pressure = (
                (1.0 - damping)*log_pressure
                + damping*np.log(hydro_pressure)
            )
        if not converged:
            pressure = np.exp(log_pressure)
        if not converged:
            log.warning(f"equilibrium atmosphere did not converge after {self._max_iter_loop} iterations")
        # Recompute chemistry once on the final pressure profile so all returned
        # fields are mutually consistent with the last hydrostatic update.
        chem_data = run_easy_chem_full_profile(
            pressure=Q_(pressure, self.pressure_unit),
            temperature=Q_(temperature, self.temperature_unit),
            atomic_abund=self.chem_input.atomic_abundances,
            chemical_species=self.chem_input.chemical_species
        )
        mu = self._mean_molecular_weight_profile(chem_data, z)
        density = pressure*mu/(self._kB*temperature)
        species_number_density = self._species_number_density(
            density=density,
            chem_data=chem_data,
            mu=mu
        )
        plot_pressure_profile(
            pressure=pressure,
            altitude=z,
            pressure_unit="bar",
            length_unit="km",
            output_file=os.path.join(self.output_dir, "pressure_profile_final.png"),
        )
        plot_chemical_mass_profiles(
            mu=mu,
            chem_data=chem_data,
            altitude=z,
            mass_unit=self.mass_unit,
            length_unit="km",
            output_file=os.path.join(self.output_dir, "chemical_mass_profiles_final.png"),
        )
        return AtmLayerDyn(
            altitude=z,
            pressure=pressure,
            temperature=temperature,
            mean_molecular_mass=mu,
            gravity=gravity,
            species_number_density=species_number_density,
            chemistry=chem_data
        )
    # implement log pressure update
    def _next_log_pressure(self, *, iteration: int, log_pressure: np.ndarray, residual: np.ndarray) -> np.ndarray:
        damping = max(0.0, min(1.0, self.settings.damping))
        return log_pressure + damping * residual

#
#    Anderson pressure / density solver
#

class AndersonPressureDensitySolver(FixedTemperaturePressureDensitySolver):
    solver_name = "Anderson fixed-temperature pressure-density solver"
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._residual_history: list[np.ndarray] = []
    # _start method
    def _start(self):
        self._residual_history.clear()
        super()._start()
        log.info(f"Anderson history depth={self.settings.anderson_depth}")
    # implement log pressure update 
    def _next_log_pressure(self, *, iteration: int, log_pressure: np.ndarray, residual: np.ndarray) -> np.ndarray:
        self._residual_history.append(residual.copy())
        if len(self._residual_history) == 1:
            return log_pressure + self.settings.damping * residual
        depth = min(self.settings.anderson_depth, len(self._residual_history) - 1)
        residual_basis = np.column_stack(self._residual_history[-depth:])
        try:
            alpha, _, _, _ = np.linalg.lstsq(residual_basis, residual, rcond=None)
        except np.linalg.LinAlgError:
            log.warning("Anderson acceleration failed; falling back to damped update")
            return log_pressure + self.settings.damping * residual
        correction = residual_basis @ alpha
        log.info(f"Anderson acceleration applied with depth={depth}")
        return log_pressure + correction
    #
    # Anderson acceleration solver (faster convergence alternative)
    #
    def _compute_equilibrium_atmospheric_configuration_anderson(self, config: AtmLayerDyn, k: int = 5):
        """
        Anderson acceleration for pressure-chemistry equilibrium.
            
        Parameters:
        -----------
        config : AtmLayerDyn
            Initial atmospheric configuration
        k : int
            Number of previous iterations to use for acceleration (default: 5)
                
        Returns:
        --------
        AtmLayerDyn
            Converged atmospheric configuration
        """
        # set input profiles
        z = np.asarray(config.altitude, dtype=float)
        temperature = np.asarray(config.temperature, dtype=float)
        pressure = np.asarray(config.pressure, dtype=float)
        gravity = np.asarray(config.gravity, dtype=float)
        print("BUBU")
        exit()
        # initialize
        converged = False
        pressure = np.maximum(pressure, self._min_pressure.to(self.pressure_unit).magnitude)
        log_pressure = np.log(pressure)
            
        # Storage for Anderson acceleration
        g_history = []  # residuals (log_pressure_new - log_pressure_old)
        x_history = []  # log_pressure values
            
        log.warning(
            "Starting Anderson-accelerated equilibrium solver: "
            f"maxiter={self._max_iter_loop}, abs_tol={self._abs_tol}, "
            f"rel_tol={self._rel_tol}, logp_tol={self._logp_tol}, k={k}"
        )
            
        for iteration in range(self._max_iter_loop):
            pressure_old = np.exp(log_pressure)
                
            # compute chemical structure
            chem_data = run_easy_chem_full_profile(
                pressure=Q_(pressure_old, self.pressure_unit),
                temperature=Q_(temperature, self.temperature_unit),
                atomic_abund=self.chem_input.atomic_abundances,
                chemical_species=self.chem_input.chemical_species
            )
                
            # average mass
            mu = self._mean_molecular_weight_profile(chem_data, z)
            if np.any(mu <= 0.0):
                log.error("mean molecular mass must be positive for equilibrium atmosphere solve")
            
            # hydrostatic pressure
            hydro_pressure = self._solve_hydrostatic_pressure(
                z=z,
                temperature=temperature,
                gravity=gravity,
                mu=mu
            )
            hydro_pressure = np.maximum(hydro_pressure, self._min_pressure.to(self.pressure_unit).magnitude)
                
            # log pressure target
            log_pressure_target = np.log(hydro_pressure)
                
            # compute residual and convergence metrics
            residual = log_pressure_target - log_pressure
            pressure_delta = np.abs(hydro_pressure - pressure_old)
            pressure_scale = np.maximum(np.abs(hydro_pressure), np.abs(pressure_old))
            safe_scale = np.maximum(pressure_scale, 1e-30)
            max_abs = np.max(pressure_delta)
            max_rel = np.max(pressure_delta / safe_scale)
            max_dlogp = np.max(np.abs(residual))
                
            # logging
            log.warning(
                f"anderson iter {iteration+1}: max_abs_dP={max_abs:.3e} Pa, "
                f"max_rel_dP={max_rel:.3e}, max_dlogP={max_dlogp:.3e}, "
                f"abs_tol={self._abs_tol}, rel_tol={self._rel_tol}, logp_tol={self._logp_tol}"
            )
                
            # check convergence
            if max_dlogp <= self._logp_tol:
                converged = True
                pressure = hydro_pressure
                log.info(f"equilibrium atmosphere converged in {iteration + 1} iterations (Anderson)")
                break
            
            # Update history for Anderson acceleration
            x_history.append(log_pressure.copy())
            g_history.append(residual.copy())
                
            # Anderson acceleration: find optimal mixing of past iterations
            if len(g_history) > 1:
                # Use last min(k, len(g_history)-1) residuals
                k_use = min(k, len(g_history) - 1)
                    
                # Stack residuals
                G = np.column_stack(g_history[-k_use:])  # shape: (n_layers, k_use)
                    
                # Solve least squares: minimize ||G @ alpha||^2
                try:
                    alpha, _, _, _ = np.linalg.lstsq(G, residual, rcond=None)
                        
                    # Anderson update: x_new = x_old + sum(alpha_i * g_i)
                    correction = np.zeros_like(log_pressure)
                    for i, g in enumerate(g_history[-k_use:]):
                        correction += alpha[i] * g
                        
                    log_pressure = log_pressure + correction
                    log.info(f"  Anderson acceleration applied with k={k_use}")
                except np.linalg.LinAlgError:
                    # Fall back to standard update if lstsq fails
                    log.warning("  Anderson acceleration failed, using standard update")
                    log_pressure = log_pressure + 0.5 * residual
            else:
                # First iteration: simple update with moderate damping
                log_pressure = log_pressure + 0.5 * residual
            
        if not converged:
            pressure = np.exp(log_pressure)
        if not converged:
            log.warning(f"equilibrium atmosphere did not converge after {self._max_iter_loop} iterations (Anderson)")
        
        # Final chemistry computation
        chem_data = run_easy_chem_full_profile(
            pressure=Q_(pressure, self.pressure_unit),
            temperature=Q_(temperature, self.temperature_unit),
            atomic_abund=self.chem_input.atomic_abundances,
            chemical_species=self.chem_input.chemical_species
        )
        mu = self._mean_molecular_weight_profile(chem_data, z)
        density = pressure*mu/(self._kB*temperature)
        species_number_density = self._species_number_density(
            density=density,
            chem_data=chem_data,
            mu=mu
        )
            
        plot_pressure_profile(
            pressure=pressure,
            altitude=z,
            pressure_unit="bar",
            length_unit="km",
            output_file=os.path.join(self.output_dir, "pressure_profile_final_anderson.png"),
        )
        plot_chemical_mass_profiles(
            mu=mu,
            chem_data=chem_data,
            altitude=z,
            mass_unit=self.mass_unit,
            length_unit="km",
            output_file=os.path.join(self.output_dir, "chemical_mass_profiles_final_anderson.png"),
        )
            
        return AtmLayerDyn(
            altitude=z,
            pressure=pressure,
            temperature=temperature,
            mean_molecular_mass=mu,
            gravity=gravity,
            species_number_density=species_number_density,
            chemistry=chem_data
        )