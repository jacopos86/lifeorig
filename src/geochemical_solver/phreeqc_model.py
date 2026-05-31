"""Wrapper for PHREEQC geochemical calculations."""
from dataclasses import dataclass, field
import logging
from pathlib import Path

try:
    from phreeqpy.iphreeqc.phreeqc_dll import IPhreeqc
except ImportError:  # pragma: no cover
    IPhreeqc = None

log = logging.getLogger(__name__)


@dataclass
class OceanChemState:
    """Structured output from an ocean/aqueous PHREEQC calculation."""
    pH: float | None = None
    temperature: float | None = None
    ionic_strength: float | None = None
    totals: dict[str, float] = field(default_factory=dict)
    raw_selected_output: list[dict[str, float]] = field(default_factory=list)

#
#  PHREEQC model
#

class PHREEQCModel:
    """Simple PHREEQC wrapper for aqueous geochemistry."""

    def __init__(self, database_path: str | Path | None = None):
        if IPhreeqc is None:
            raise ImportError("phreeqpy is required for PHREEQCModel")

        self.iphreeqc = IPhreeqc()
        self.database_path = Path(database_path) if database_path is not None else None
        if self.database_path is not None:
            self.load_database(self.database_path)

    def load_database(self, database_path: str | Path) -> None:
        """Load a PHREEQC database file."""
        self.database_path = Path(database_path)
        log.debug("Loading PHREEQC database from %s", self.database_path)
        self.iphreeqc.load_database(str(self.database_path))

    def run_string(self, input_string: str) -> None:
        """Run a PHREEQC input string."""
        log.debug("Running PHREEQC input string")
        self.iphreeqc.run_string(input_string)

    def build_solution_input(
        self,
        solution_id: int,
        temperature: float,
        pressure: float,
        composition: dict[str, float],
        selected_output: list[str] | None = None,
        user_punch: list[str] | None = None,
    ) -> str:
        """Build a PHREEQC solution definition string."""
        lines = [f"SOLUTION {solution_id}", f"    temp {temperature}", f"    pH 7", f"    pressure {pressure}"]
        for species, amount in composition.items():
            lines.append(f"    {species} {amount}")

        if selected_output is not None:
            lines.append("\nSELECTED_OUTPUT")
            for line in selected_output:
                lines.append(f"    {line}")

        if user_punch is not None:
            lines.append("\nUSER_PUNCH")
            for line in user_punch:
                lines.append(f"    {line}")

        lines.append("END")
        return "\n".join(lines)

    def run_solution(
        self,
        solution_id: int,
        temperature: float,
        pressure: float,
        composition: dict[str, float],
        selected_output: list[str] | None = None,
        user_punch: list[str] | None = None,
    ) -> None:
        """Run a PHREEQC solution section and store the result."""
        input_string = self.build_solution_input(
            solution_id=solution_id,
            temperature=temperature,
            pressure=pressure,
            composition=composition,
            selected_output=selected_output,
            user_punch=user_punch,
        )
        self.run_string(input_string)

    def get_selected_output(self):
        """Extract selected output from the last PHREEQC run."""
        return self.iphreeqc.get_selected_output_array()

    def get_selected_output_records(self) -> list[dict[str, float]]:
        """Return selected output as dictionaries keyed by PHREEQC headings."""
        output = self.get_selected_output()
        if not output:
            return []
        headings = [str(item) for item in output[0]]
        records = []
        for row in output[1:]:
            records.append({
                heading: value
                for heading, value in zip(headings, row)
            })
        return records

    def solve_ocean_equilibrium(
        self,
        temperature: float,
        pressure: float,
        composition: dict[str, float],
        solution_id: int = 1,
    ) -> OceanChemState:
        """Run a first-pass ocean equilibrium calculation.

        Parameters are passed directly to PHREEQC:
        - temperature in Celsius
        - pressure in PHREEQC pressure units
        - composition as PHREEQC solution totals, e.g. {"Na": 1.0, "Cl": 1.0}
        """
        selected_output = [
            "-reset false",
            "-pH true",
            "-temperature true",
            "-ionic_strength true",
        ]
        user_punch = []
        if composition:
            headings = " ".join(composition.keys())
            totals = ", ".join(f'TOT("{species}")' for species in composition)
            user_punch = [
                f"-headings {headings}",
                f"PUNCH {totals}",
            ]
        self.run_solution(
            solution_id=solution_id,
            temperature=temperature,
            pressure=pressure,
            composition=composition,
            selected_output=selected_output,
            user_punch=user_punch,
        )
        records = self.get_selected_output_records()
        if not records:
            return OceanChemState()
        latest = records[-1]
        totals = {
            species: latest[species]
            for species in composition
            if species in latest
        }
        return OceanChemState(
            pH=latest.get("pH"),
            temperature=latest.get("temp"),
            ionic_strength=latest.get("mu"),
            totals=totals,
            raw_selected_output=records,
        )
