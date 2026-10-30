import os
from src.atmosphere_solver.atm_struct_driver import AtmosphSolver
from src.utilities.plot_chemical_mass_profiles import plot_chemical_mass_profiles

#
#   Non-equilibrium layered structure
#

class NonEquilLayeredAtmosphSolver(AtmosphSolver):
    def __init__(self, chem_input, chem_system, stellar_data, planet_data, atmosphere_data, output_dir=None):
        super().__init__(
            chem_input=chem_input,
            stellar_data=stellar_data,
            planet_data=planet_data,
            atmosphere_data=atmosphere_data,
            output_dir=output_dir,
        )
        self.chem_system = chem_system
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
