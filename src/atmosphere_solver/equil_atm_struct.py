import os
import numpy as np
from src.atmosphere_solver.atm_struct_driver import AtmosphSolver
from src.atmosphere_solver.atmosph_data import AtmDynResult, AtmLayerDyn
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
        # 1) set base arrays
        z = self._build_altitude_grid()
        g = self._gravity_profile(z)
        # 2) set initial arrays guess
        T0_z = self._initial_temperature_profile(z)
        # 3) set initial variables
        layered_variables = self._hydro_solver.initialize(
            altitude=z,
            temperature=T0_z,
            gravity=g
        )
        plot_chemical_mass_profiles(
            mu=layered_variables.mean_molecular_mass,
            chem_data=layered_variables.chemistry,
            altitude=layered_variables.altitude,
            mass_unit=self._units.get("mass"),
            length_unit="km",
            output_file=os.path.join(self.output_dir, "chemical_mass_profiles_initial_guess.png"),
        )
        # start outer iterations
        converged = False
        for iteration in range(self._radiative_solver.settings.max_iter):
            hydro_result = self._hydro_solver.solve(layered_variables)
            radiative_result = self._radiative_solver.solve(hydro_result.layered_variables)
            # check T convergence
            T0_z = layered_variables.temperature.copy()
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