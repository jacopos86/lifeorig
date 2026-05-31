from dataclasses import dataclass
from src.utilities.logging_module import log
from src.exo_chem.easy_chem_driver import run_easychem_backend
from src.atmosphere_solver.equil_atm_struct import LayeredEquilibriumAtmosphSolver
from src.chemical_env.chem_env_data import ChemEnvInput, ChemEnvResult

#
#   Set input ChemEnv
#

def set_ChemEnvInput(planet_data, env_data) -> ChemEnvInput:
    return ChemEnvInput(
        mode=planet_data["basic_info"].chemical_env.get("mode"),
        atomic_abundances=planet_data["basic_info"].chemical_env.get("atomic_abundances"),
        chemical_species=planet_data["basic_info"].chemical_env.get("chemical_species"),
        pressure=env_data.get("pressure"),
        temperature=env_data.get("temperature")
    )

#
#   Simulate Chemical environment 
#

class SimulateChemEnv:
    def __init__(self, chem_input: ChemEnvInput, planet_data=None, stellar_data=None, atmosphere_data=None, output_dir=None):
        self.chem_input = chem_input
        self.stellar_data = stellar_data
        self.planet_data = planet_data
        self.atmosphere_data = atmosphere_data
        self.output_dir = output_dir
    def run(self) -> ChemEnvResult:
        if self.chem_input.mode == "preset":
            return ChemEnvResult(mode="preset")
        if self.chem_input.mode == "local_equilibrium":
            return self._run_local_equilibrium()
        if self.chem_input.mode == "layered_equilibrium":
            return self._run_layered_equilibrium()
        if self.chem_input.mode == "layered_disequilibrium":
            return self._run_out_equilibrium()
        log.error(f"Unknown chemistry mode: {self.chem_input.mode}")
    def _run_local_equilibrium(self) -> ChemEnvResult:
        if self.chem_input.pressure is None:
            log.error("pressure is required for local_equilibrium")
        if self.chem_input.temperature is None:
            log.error("temperature is required for local_equilibrium")
        result = run_easychem_backend(
            pressure=self.chem_input.pressure, 
            temperature=self.chem_input.temperature, 
            atomic_abund=self.chem_input.atomic_abundances
        )
        return ChemEnvResult(
            mode="local_equilibrium",
            local=result
        )
    def _run_layered_equilibrium(self) -> ChemEnvResult:
        return LayeredEquilibriumAtmosphSolver(
            chem_input=self.chem_input,
            stellar_data=self.stellar_data,
            planet_data=self.planet_data,
            atmosphere_data=self.atmosphere_data,
            output_dir=self.output_dir
        ).run()
    def _run_out_equilibrium(self) -> ChemEnvResult:
        log.error("layered_disequilibrium solver is not implemented")
