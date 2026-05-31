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
        converged = False
        pressure = np.maximum(pressure, self._min_pressure.to(self.pressure_unit).magnitude)
        # log pressure
        log_pressure = np.log(pressure)
        damping = self._damping_loop
        log.warning(
            "Starting damped log-pressure equilibrium solver: "
            f"maxiter={self._max_iter_loop}, abs_tol={self._abs_tol}, "
            f"rel_tol={self._rel_tol}, damping={damping:.3f}"
        )
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
            # print info
            log.warning(
                f"equilibrium iter {iteration+1}: max_abs_dP={max_abs:.3e} Pa, "
                f"max_rel_dP={max_rel:.3e}, max_dlogP={max_dlogp:.3e}, "
                f"damping={damping:.3f}, abs_tol={self._abs_tol}, rel_tol={self._rel_tol}"
            )
            if np.all(pressure_delta <= self._abs_tol + self._rel_tol*safe_scale):
                converged = True
                pressure = hydro_pressure
                log.info(f"equilibrium atmosphere converged in {iteration + 1} iterations")
                break
            # update log pressure
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
