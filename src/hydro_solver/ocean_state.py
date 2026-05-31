from dataclasses import dataclass
import numpy as np
from scipy.constants import Boltzmann
from src.common.units import Q_


@dataclass
class OceanPresenceState:
    has_ocean: bool
    reason: str
    surface_temperature: Q_
    surface_pressure: Q_
    boiling_temperature: Q_ | None
    ocean_mass: Q_ | None
    ocean_depth: Q_ | None
    ocean_bottom_pressure: Q_ | None


@dataclass
class CondensedOceanState:
    has_ocean: bool
    reason: str
    surface_temperature: Q_
    surface_pressure: Q_
    surface_h2o_partial_pressure: Q_
    surface_h2o_saturation_pressure: Q_
    atmospheric_h2o_mass: Q_
    saturation_h2o_mass: Q_
    condensed_h2o_mass: Q_
    ocean_depth: Q_
    ocean_bottom_pressure: Q_


def water_boiling_temperature(surface_pressure: Q_) -> Q_:
    """Approximate water boiling temperature from pressure using Antoine's law.

    Valid as a first-pass estimate near ordinary liquid-water pressures.
    """
    pressure_mmhg = surface_pressure.to("mmHg").magnitude
    if pressure_mmhg <= 0.0:
        return Q_(np.nan, "kelvin")
    # Antoine coefficients for water, temperature in Celsius.
    a = 8.07131
    b = 1730.63
    c = 233.426
    temperature_C = b / (a - np.log10(pressure_mmhg)) - c
    return Q_(temperature_C, "degC").to("kelvin")


def water_saturation_vapor_pressure(temperature: Q_) -> Q_:
    """Approximate saturation vapor pressure of water over liquid water."""
    temperature_C = temperature.to("degC").magnitude
    # Magnus formula, pressure in Pa. Good enough for the first condensation check.
    pressure_Pa = 610.94 * np.exp(17.625 * temperature_C / (temperature_C + 243.04))
    return Q_(pressure_Pa, "Pa")


def compute_condensed_ocean_from_atmosphere(
        altitude: Q_,
        temperature: Q_,
        pressure: Q_,
        h2o_number_density,
        planet_radius: Q_,
        gravity: Q_,
        water_density: Q_ = Q_(1000.0, "kg / m^3"),
    ) -> CondensedOceanState:
    """Compute liquid surface water implied by atmospheric H2O condensation.

    The condensed mass is the excess atmospheric H2O above saturation capacity:

        M_liquid = integral max(0, n_H2O - n_sat(T)) m_H2O dV
    """
    z = altitude.to("m").magnitude
    temperature_K = temperature.to("kelvin").magnitude
    pressure_Pa = pressure.to("Pa").magnitude
    n_h2o = np.asarray(h2o_number_density, dtype=float)
    if z.shape != temperature_K.shape or z.shape != pressure_Pa.shape or z.shape != n_h2o.shape:
        raise ValueError("altitude, temperature, pressure, and h2o_number_density must have the same shape")
    p_sat = water_saturation_vapor_pressure(Q_(temperature_K, "kelvin")).to("Pa").magnitude
    n_sat = p_sat / (Boltzmann * temperature_K)
    n_excess = np.maximum(n_h2o - n_sat, 0.0)
    radius_m = planet_radius.to("m").magnitude
    shell_area = 4.0 * np.pi * (radius_m + z)**2
    m_h2o = Q_(18.01528, "amu").to("kg").magnitude
    atmospheric_h2o_mass = Q_(np.trapezoid(n_h2o * m_h2o * shell_area, z), "kg")
    saturation_h2o_mass = Q_(np.trapezoid(n_sat * m_h2o * shell_area, z), "kg")
    condensed_h2o_mass = Q_(np.trapezoid(n_excess * m_h2o * shell_area, z), "kg")
    surface_area = 4.0 * np.pi * radius_m**2
    ocean_depth = (condensed_h2o_mass / (water_density * Q_(surface_area, "m^2"))).to("m")
    ocean_bottom_pressure = (pressure[0].to("Pa") + water_density * gravity * ocean_depth).to("Pa")
    surface_temperature = temperature[0].to("kelvin")
    surface_pressure = pressure[0].to("Pa")
    surface_h2o_partial_pressure = Q_(n_h2o[0] * Boltzmann * temperature_K[0], "Pa")
    surface_h2o_saturation_pressure = Q_(p_sat[0], "Pa")
    has_liquid_conditions = (
        surface_pressure >= Q_(611.657, "Pa")
        and surface_temperature >= Q_(273.15, "kelvin")
        and surface_temperature <= water_boiling_temperature(surface_pressure)
    )
    has_condensed_water = condensed_h2o_mass.to("kg").magnitude > 0.0
    if not has_liquid_conditions:
        reason = "surface pressure-temperature conditions do not allow liquid water"
    elif not has_condensed_water:
        reason = "atmospheric H2O does not exceed saturation capacity"
    else:
        reason = "atmospheric H2O condensation allows liquid surface water"
    return CondensedOceanState(
        has_ocean=has_liquid_conditions and has_condensed_water,
        reason=reason,
        surface_temperature=surface_temperature,
        surface_pressure=surface_pressure,
        surface_h2o_partial_pressure=surface_h2o_partial_pressure,
        surface_h2o_saturation_pressure=surface_h2o_saturation_pressure,
        atmospheric_h2o_mass=atmospheric_h2o_mass,
        saturation_h2o_mass=saturation_h2o_mass,
        condensed_h2o_mass=condensed_h2o_mass,
        ocean_depth=ocean_depth,
        ocean_bottom_pressure=ocean_bottom_pressure,
    )


def compute_ocean_presence(
        surface_temperature: Q_,
        surface_pressure: Q_,
        planet_radius: Q_,
        gravity: Q_,
        water_mass: Q_ | None = None,
        water_mass_fraction: float | None = None,
        planet_mass: Q_ | None = None,
        water_density: Q_ = Q_(1000.0, "kg / m^3"),
    ) -> OceanPresenceState:
    """Check whether a surface ocean is consistent with planetary conditions."""
    temperature_K = surface_temperature.to("kelvin")
    pressure_Pa = surface_pressure.to("Pa")
    if water_mass is None and water_mass_fraction is not None and planet_mass is not None:
        water_mass = water_mass_fraction * planet_mass
    if water_mass is None or water_mass.to("kg").magnitude <= 0.0:
        return OceanPresenceState(
            has_ocean=False,
            reason="no positive water inventory",
            surface_temperature=temperature_K,
            surface_pressure=pressure_Pa,
            boiling_temperature=None,
            ocean_mass=water_mass,
            ocean_depth=None,
            ocean_bottom_pressure=None,
        )
    triple_pressure = Q_(611.657, "Pa")
    freezing_temperature = Q_(273.15, "kelvin")
    if pressure_Pa < triple_pressure:
        return OceanPresenceState(
            has_ocean=False,
            reason="surface pressure below water triple point",
            surface_temperature=temperature_K,
            surface_pressure=pressure_Pa,
            boiling_temperature=None,
            ocean_mass=water_mass,
            ocean_depth=None,
            ocean_bottom_pressure=None,
        )
    boiling_temperature = water_boiling_temperature(pressure_Pa)
    if temperature_K < freezing_temperature:
        return OceanPresenceState(
            has_ocean=False,
            reason="surface temperature below freezing point",
            surface_temperature=temperature_K,
            surface_pressure=pressure_Pa,
            boiling_temperature=boiling_temperature,
            ocean_mass=water_mass,
            ocean_depth=None,
            ocean_bottom_pressure=None,
        )
    if temperature_K > boiling_temperature:
        return OceanPresenceState(
            has_ocean=False,
            reason="surface temperature above boiling point",
            surface_temperature=temperature_K,
            surface_pressure=pressure_Pa,
            boiling_temperature=boiling_temperature,
            ocean_mass=water_mass,
            ocean_depth=None,
            ocean_bottom_pressure=None,
        )
    surface_area = 4.0 * np.pi * planet_radius.to("m").magnitude**2
    ocean_depth = (water_mass / (water_density * Q_(surface_area, "m^2"))).to("m")
    ocean_bottom_pressure = (pressure_Pa + water_density * gravity * ocean_depth).to("Pa")
    return OceanPresenceState(
        has_ocean=True,
        reason="liquid surface ocean allowed",
        surface_temperature=temperature_K,
        surface_pressure=pressure_Pa,
        boiling_temperature=boiling_temperature,
        ocean_mass=water_mass.to("kg"),
        ocean_depth=ocean_depth,
        ocean_bottom_pressure=ocean_bottom_pressure,
    )
