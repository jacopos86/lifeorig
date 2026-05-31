from pathlib import Path
from src.geochemical_solver.phreeqc_model import OceanChemState, PHREEQCModel


def solve_ocean_equilibrium(
        temperature: float,
        pressure: float,
        composition: dict[str, float],
        database_path: str | Path = "./external/phreeqc/database/phreeqc.dat",
    ) -> OceanChemState:
    """Solve a first-pass aqueous/ocean equilibrium problem with PHREEQC."""
    model = PHREEQCModel(database_path=database_path)
    return model.solve_ocean_equilibrium(
        temperature=temperature,
        pressure=pressure,
        composition=composition,
    )
