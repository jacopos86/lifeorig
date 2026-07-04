from dataclasses import dataclass, field
import numpy as np
from src.environment.environment_base_class import Environment


@dataclass
class PondState:
    """
    Well-mixed shallow pond state.

    This is the simplest surface-pond model: one aqueous reservoir with a
    time-dependent water depth. Drying concentrates solutes; refilling dilutes
    them.
    """
    time: float
    water_depth: float
    temperature: float
    pressure: float
    pH: float
    ionic_strength: float
    uv_flux: float
    sediment_surface_area: float
    volume: float
    concentrations: dict[str, float] = field(default_factory=dict)


class SurfacePond(Environment):
    """
    Basic surface pond environment.

    Geometry:
        one well-mixed shallow water body above a sediment/mineral surface

    Main drivers:
        water-depth cycling, UV exposure, atmospheric input, sediment surface
        chemistry, concentration during drying, dilution during refilling.
    """
    def __init__(self, protocell_list=None, data=None):
        super().__init__(protocell_list or [])
        data = data or {}

        self.area = self._as_float(data.get("area"), default=1.0)
        self.mean_water_depth = self._as_float(data.get("mean_water_depth"), default=1.0e-2)
        self.water_depth_amplitude = self._as_float(data.get("water_depth_amplitude"), default=5.0e-3)
        self.cycle_period = self._as_float(data.get("cycle_period"), default=86400.0)
        self.min_water_depth = self._as_float(data.get("min_water_depth"), default=1.0e-5)

        self.temperature = self._as_float(data.get("temperature"), default=298.15)
        self.pressure = self._as_float(data.get("pressure"), default=1.0e5)
        self.pH = float(data.get("pH", 7.0))
        self.ionic_strength = float(data.get("ionic_strength", 0.01))
        self.uv_flux = float(data.get("uv_flux", 1.0))
        self.sediment_surface_area = self._as_float(
            data.get("sediment_surface_area"),
            default=self.area,
        )

        self.concentrations = dict(data.get("concentrations", {}))
        self.atmospheric_input = dict(data.get("atmospheric_input", {}))
        self.washout_rate = float(data.get("washout_rate", 0.0))
        self.photolysis_rates = dict(data.get("photolysis_rates", {}))

        self.time = 0.0
        self.state = self._build_state()

    def step(self, dt):
        old_depth = self.state.water_depth
        self.time += dt
        new_depth = self._water_depth(self.time)
        self._apply_depth_concentration(old_depth=old_depth, new_depth=new_depth)
        self._apply_sources_and_losses(dt)
        self.state = self._build_state()

    def update_external_conditions(self, dt):
        self.state = self._build_state()

    def apply_transport(self, dt):
        self._apply_sources_and_losses(dt)

    def apply_environmental_effects(self, dt):
        return

    def _build_state(self):
        water_depth = self._water_depth(self.time)
        return PondState(
            time=self.time,
            water_depth=water_depth,
            temperature=self.temperature,
            pressure=self.pressure,
            pH=self.pH,
            ionic_strength=self.ionic_strength,
            uv_flux=self.uv_flux,
            sediment_surface_area=self.sediment_surface_area,
            volume=self.area * water_depth,
            concentrations=dict(self.concentrations),
        )

    def _water_depth(self, time):
        if self.cycle_period <= 0.0:
            return max(self.min_water_depth, self.mean_water_depth)
        phase = 2.0 * np.pi * time / self.cycle_period
        depth = self.mean_water_depth + self.water_depth_amplitude * np.cos(phase)
        return max(self.min_water_depth, depth)

    def _apply_depth_concentration(self, old_depth, new_depth):
        if old_depth <= 0.0 or new_depth <= 0.0:
            return
        factor = old_depth / new_depth
        for species in self.concentrations:
            self.concentrations[species] *= factor

    def _apply_sources_and_losses(self, dt):
        for species, source_rate in self.atmospheric_input.items():
            self.concentrations[species] = self.concentrations.get(species, 0.0) + source_rate * dt
        for species, rate in self.photolysis_rates.items():
            if species in self.concentrations:
                self.concentrations[species] *= np.exp(-rate * self.uv_flux * dt)
        if self.washout_rate > 0.0:
            loss = np.exp(-self.washout_rate * dt)
            for species in self.concentrations:
                self.concentrations[species] *= loss

    def get_environment_state(self):
        self.state = self._build_state()
        return {
            "time": self.state.time,
            "water_depth": self.state.water_depth,
            "temperature": self.state.temperature,
            "pressure": self.state.pressure,
            "pH": self.state.pH,
            "ionic_strength": self.state.ionic_strength,
            "uv_flux": self.state.uv_flux,
            "sediment_surface_area": self.state.sediment_surface_area,
            "volume": self.state.volume,
            "concentrations": self.state.concentrations,
        }

    @staticmethod
    def _as_float(value, default):
        if value is None:
            return float(default)
        if hasattr(value, "to"):
            try:
                return float(value.to_base_units().magnitude)
            except Exception:
                return float(value.magnitude)
        return float(value)
