"""Couple geochemical model outputs to environment and protocell state."""

from __future__ import annotations

from typing import Any


class GeochemicalCoupler:
    """Map geochemical outputs into environment fluxes and protocell inputs."""

    def __init__(self, phreeqc_model: Any, thermal_model: Any):
        self.phreeqc_model = phreeqc_model
        self.thermal_model = thermal_model

    def build_phreeqc_state(
        self,
        fluid_composition: dict[str, float],
        temperature: float,
        pressure: float,
        mineral_phases: dict[str, float] | None = None,
    ) -> str:
        """Build a PHREEQC input string for the current geochemical state."""
        selected_output = ["-reset false", "-pH true", "-temperature true", "-ionic_strength true"]
        user_punch = ["-headings"]
        current_species = " ".join(fluid_composition.keys())
        user_punch.append(f"PUNCH TOT(\"{current_species}\")")
        return self.phreeqc_model.build_solution_input(
            solution_id=1,
            temperature=temperature,
            pressure=pressure,
            composition=fluid_composition,
            selected_output=selected_output,
            user_punch=user_punch,
        )

    def compute_surface_fluxes(
        self,
        geochemical_outputs: dict[str, float],
        flux_scale: float = 1.0,
    ) -> dict[str, float]:
        """Convert geochemical outputs into surface/environment fluxes."""
        return {species: value * flux_scale for species, value in geochemical_outputs.items()}

    def apply_to_environment(self, env: Any, fluxes: dict[str, float], dt: float) -> None:
        """Apply fluxes from geochemistry to an Environment object."""
        if not hasattr(env, "external_field"):
            raise AttributeError("Environment object must expose external_field")

        for species, flux in fluxes.items():
            env.external_field.setdefault(species, 0.0)
            env.external_field[species] += flux * dt
