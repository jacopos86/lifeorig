"""Thermal transport utilities for geochemical systems."""

import numpy as np


class ThermalModel:
    """A simple one-dimensional thermal model for subsurface heat transport."""

    def __init__(
        self,
        thermal_conductivity: float = 2.5,
        heat_capacity: float = 1e3,
        density: float = 3000.0,
        internal_heat_production: float = 1e-6,
    ):
        self.thermal_conductivity = thermal_conductivity
        self.heat_capacity = heat_capacity
        self.density = density
        self.internal_heat_production = internal_heat_production

    def steady_state_profile(
        self,
        depth: np.ndarray,
        surface_temperature: float = 300.0,
        basal_temperature: float | None = None,
    ) -> np.ndarray:
        """Return a simple steady-state temperature profile.

        If a basal temperature is provided, a linear gradient is used.
        Otherwise the model returns a first-order geothermal gradient.
        """
        if basal_temperature is not None:
            gradient = (basal_temperature - surface_temperature) / max(depth[-1], 1.0)
            return surface_temperature + gradient * depth

        # simple geothermal gradient = q / k
        q = self.internal_heat_production * depth[-1]
        gradient = q / max(self.thermal_conductivity, 1e-12)
        return surface_temperature + gradient * depth

    def transient_step(
        self,
        temperature_profile: np.ndarray,
        depth: np.ndarray,
        dt: float,
    ) -> np.ndarray:
        """Advance the temperature profile using explicit 1D conduction."""
        dz = np.diff(depth)
        if np.any(dz <= 0.0):
            raise ValueError("depth must be a strictly increasing array")

        new_temperature = temperature_profile.copy()
        alpha = self.thermal_conductivity / (self.density * self.heat_capacity)
        for i in range(1, len(depth) - 1):
            d2T = (
                temperature_profile[i + 1]
                - 2.0 * temperature_profile[i]
                + temperature_profile[i - 1]
            )
            new_temperature[i] += alpha * dt * d2T / (dz[i - 1] * dz[i])

        return new_temperature
