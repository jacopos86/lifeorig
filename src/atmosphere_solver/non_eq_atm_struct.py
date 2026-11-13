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
        # 1) set initial state
        layered_variables = self._build_initial_state_layers()
        exit()
        # 2) set up dynamical solver
        self._set_up_dynamic_solver(layered_variables)
        # 3) evolve
        self._evolve()
        # 4) post process