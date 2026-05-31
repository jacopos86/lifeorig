from dataclasses import dataclass
import numpy as np
from src.common.units import Q_
from src.geochemical_solver.interior_heat_sources import InteriorHeatSources
from src.geochemical_solver.thermal_model import ThermalModel
from src.utilities.logging_module import log


@dataclass
class GeophysicalState:
    depth: np.ndarray
    temperature: np.ndarray
    internal_power: Q_
    internal_heat_flux: Q_
    surface_temperature: Q_
    basal_temperature: Q_ | None = None


@dataclass
class GeophysicalResult:
    planet_params: object
    state: GeophysicalState
    heat_sources: InteriorHeatSources


class GeophysicalDriver:
    """First-pass geophysical/interior driver.

    This is the geophysical analogue of the atmospheric driver: it owns the
    setup of units, heat sources, and a thermal profile. More detailed
    conduction/convection solvers can replace the thermal call later.
    """

    def __init__(
            self,
            planet_data,
            geophysical_data: dict | None = None,
            output_dir=None,
        ):
        self.planet_data = planet_data
        self.geophysical_data = geophysical_data or {}
        self.output_dir = output_dir or "."
        self._set_internal_units()
        self.heat_sources = self._build_heat_sources()
        self.thermal_model = self._build_thermal_model()

    def _set_internal_units(self):
        self.length_unit = "m"
        self.temperature_unit = "K"
        self.power_unit = "W"
        self.heat_flux_unit = "W / m^2"

    def _build_depth_grid(self) -> np.ndarray:
        n_layers = int(self.geophysical_data.get("n_layers", 200))
        max_depth = self.geophysical_data.get("max_depth", self.planet_data.planet_radius)
        return np.linspace(
            0.0,
            max_depth.to(self.length_unit).magnitude,
            n_layers,
        )

    def _build_heat_sources(self) -> InteriorHeatSources:
        heat_sources = self.geophysical_data.get("heat_sources")
        if heat_sources is not None:
            return heat_sources
        return InteriorHeatSources()

    def _build_thermal_model(self) -> ThermalModel:
        thermal_data = self.geophysical_data.get("thermal", {})
        return ThermalModel(
            thermal_conductivity=float(thermal_data.get("thermal_conductivity", 2.5)),
            heat_capacity=float(thermal_data.get("heat_capacity", 1.0e3)),
            density=float(thermal_data.get("density", 3000.0)),
            internal_heat_production=float(thermal_data.get("internal_heat_production", 1.0e-6)),
        )

    def run(self) -> GeophysicalResult:
        time = self.geophysical_data.get("time", Q_(0.0, "year"))
        surface_temperature = self.geophysical_data.get("surface_temperature", Q_(300.0, "K"))
        basal_temperature = self.geophysical_data.get("basal_temperature")
        depth = self._build_depth_grid()
        internal_power = self.heat_sources.total_power(time)
        internal_heat_flux = self.heat_sources.surface_flux(
            time=time,
            planet_radius=self.planet_data.planet_radius,
        )
        temperature = self.thermal_model.steady_state_profile(
            depth=depth,
            surface_temperature=surface_temperature.to(self.temperature_unit).magnitude,
            basal_temperature=(
                basal_temperature.to(self.temperature_unit).magnitude
                if basal_temperature is not None
                else None
            ),
        )
        log.info(f"\t geophysical internal power: {internal_power.to(self.power_unit):.3e}")
        log.info(f"\t geophysical internal heat flux: {internal_heat_flux.to(self.heat_flux_unit):.3e}")
        state = GeophysicalState(
            depth=depth,
            temperature=temperature,
            internal_power=internal_power,
            internal_heat_flux=internal_heat_flux,
            surface_temperature=surface_temperature.to(self.temperature_unit),
            basal_temperature=(
                basal_temperature.to(self.temperature_unit)
                if basal_temperature is not None
                else None
            ),
        )
        return GeophysicalResult(
            planet_params=self.planet_data,
            state=state,
            heat_sources=self.heat_sources,
        )
