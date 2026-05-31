from dataclasses import dataclass, field
import numpy as np
from src.common.phys_constants import G
from src.common.units import Q_


@dataclass
class RadiogenicHeating:
    """Radiogenic heating model with optional isotope-like components."""

    initial_power: Q_ = Q_(0.0, "W")
    decay_timescale: Q_ | None = None
    components: dict[str, tuple[Q_, Q_]] = field(default_factory=dict)

    def power(self, time: Q_) -> Q_:
        time_s = time.to("s")
        if self.components:
            total = Q_(0.0, "W")
            for initial_power, half_life in self.components.values():
                decay_constant = np.log(2.0) / half_life.to("s").magnitude
                total += initial_power.to("W") * np.exp(-decay_constant * time_s.magnitude)
            return total.to("W")
        if self.decay_timescale is None:
            return self.initial_power.to("W")
        return (
            self.initial_power.to("W")
            * np.exp(-time_s.magnitude / self.decay_timescale.to("s").magnitude)
        ).to("W")


@dataclass
class TidalHeating:
    """Eccentricity tide heating skeleton."""

    host_mass: Q_ | None = None
    body_radius: Q_ | None = None
    semi_major_axis: Q_ | None = None
    eccentricity: float = 0.0
    love_number_k2: float = 0.0
    tidal_quality_Q: float = 1.0

    def power(self) -> Q_:
        if (
            self.host_mass is None
            or self.body_radius is None
            or self.semi_major_axis is None
            or self.eccentricity <= 0.0
            or self.love_number_k2 <= 0.0
            or self.tidal_quality_Q <= 0.0
        ):
            return Q_(0.0, "W")
        host_mass = self.host_mass.to("kg")
        radius = self.body_radius.to("m")
        semi_major_axis = self.semi_major_axis.to("m")
        grav_const = G.to("m^3 / kg / s^2")
        mean_motion = (grav_const * host_mass / semi_major_axis**3)**0.5
        power = (
            (21.0 / 2.0)
            * (self.love_number_k2 / self.tidal_quality_Q)
            * grav_const
            * host_mass**2
            * radius**5
            * mean_motion
            * self.eccentricity**2
            / semi_major_axis**6
        )
        return power.to("W")


@dataclass
class GravitationalHeating:
    """Optional early-time gravitational/differentiation heat source."""

    initial_power: Q_ = Q_(0.0, "W")
    decay_timescale: Q_ | None = None

    def power(self, time: Q_) -> Q_:
        if self.decay_timescale is None:
            return self.initial_power.to("W")
        return (
            self.initial_power.to("W")
            * np.exp(-time.to("s").magnitude / self.decay_timescale.to("s").magnitude)
        ).to("W")


@dataclass
class InteriorHeatSources:
    """Container for all interior heat sources."""

    radiogenic: RadiogenicHeating = field(default_factory=RadiogenicHeating)
    tidal: TidalHeating = field(default_factory=TidalHeating)
    gravitational: GravitationalHeating = field(default_factory=GravitationalHeating)

    def total_power(self, time: Q_) -> Q_:
        return (
            self.radiogenic.power(time)
            + self.tidal.power()
            + self.gravitational.power(time)
        ).to("W")

    def surface_flux(self, time: Q_, planet_radius: Q_) -> Q_:
        surface_area = 4.0 * np.pi * planet_radius.to("m")**2
        return (self.total_power(time) / surface_area).to("W / m^2")
