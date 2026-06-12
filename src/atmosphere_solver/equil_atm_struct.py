import os
import contextlib
import io
import numpy as np
from src.atmosphere_solver.atm_struct_driver import (
    AtmDynResult,
    AtmLayerDyn,
    AtmosphSolver,
)
from src.chemical_env.chem_env_data import ChemEnvResult
from src.common.units import Q_
from src.exo_chem.easy_chem_driver import run_easy_chem_full_profile
from src.utilities.plot_chemical_mass_profiles import plot_chemical_mass_profiles
from src.utilities.plot_pressure_profile import plot_pressure_profile
from src.utilities.logging_module import log

#
#   Layered equilibrium solver
#

class LayeredEquilibriumAtmosphSolver(AtmosphSolver):
    def __init__(self, chem_input, stellar_data, planet_data, atmosphere_data, output_dir=None):
        super().__init__(
            chem_input=chem_input,
            stellar_data=stellar_data,
            planet_data=planet_data,
            atmosphere_data=atmosphere_data,
            output_dir=output_dir,
        )
    # run main driver
    def run(self):
        # 1) set wavelength grids
        self._set_wavelength_grids()
        # 2) set stellar flux
        self._set_stellar_intensity()
        # 3) set base stellar flux
        self._set_base_stellar_flux()
        # 4) set initial atmospheric layer
        initial_atm_config = self._build_initial_state_layers()
        # 5) solve hydrostatic density + equilibrium chemistry at fixed T(z)
        config = self._compute_equilibrium_atmospheric_configuration(initial_atm_config)
        # 6) compute UV and IR optical depths
        tau_z = self._compute_optical_depth(config)
        # iterate until convergence reached with temperature
        exit()
        #
        #  set final results
        #
        result = AtmDynResult(
            atomic_abundances=self.chem_input.atomic_abundances,
            stellar_params=self.stellar_data,
            uv_wavelength_grid=self.UV_wavelength_grid,
            ir_wavelength_grid=self.IR_wavelength_grid,
            stellar_B_lambda=self.stellar_B_lambda,
            bond_albedo=self.bond_albedo,
            spectral_albedo=self.spectral_albedo,
            layers=layers,
        )
        return ChemEnvResult(
            mode="layered_equilibrium",
            layered=result,
        )
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
