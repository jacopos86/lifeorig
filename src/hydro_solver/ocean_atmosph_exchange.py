from dataclasses import dataclass


@dataclass(frozen=True)
class OceanAtmosphereFluxes:
    """
    Planetary-area averaged atmosphere-ocean fluxes.

    water_mass_flux > 0:
        ocean -> atmosphere

    gas_molar_fluxes[species] > 0:
        atmosphere -> ocean

    sensible_heat_flux > 0:
        ocean -> atmosphere

    latent_heat_flux > 0:
        ocean -> atmosphere
    """
    # change name if necessary
    water_mass_flux: object
    gas_molar_fluxes: dict[str, object]
    sensible_heat_flux: object
    latent_heat_flux: object
    # saturation vapor pressure
    def saturation_vapor_pressure(self, temperature):
        """
        Saturation vapor pressure of water.

        This is a basic Clausius-Clapeyron approximation.
        """
        T = temperature.to("kelvin").magnitude
        T_ref = 273.15
        P_ref = 611.657
        latent_heat = 2.5e6
        Rw = 461.5
        exp = (
            latent_heat / Rw
            * (
                1.0 / T_ref - 1.0 / T
            )
        )
        return Q_(
            P_ref * xp.exp(exp),
            "pascal"
        )
    # specific humidity
    def specific_humidity(self, water_partial_pressure, total_pressure):
        """
        Convert water partial pressure into specific humidity.
        """
        epsilon = 0.622