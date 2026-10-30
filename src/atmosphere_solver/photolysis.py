


@dataclass(frozen=True)
class PhotolysisChannel:
    reaction_id: str
    wavelength: object
    cross_section: object
    quantum_yield: object

#
#   compute photolysis rate here
#

def calculate_photolysis_rate(
        wavelength,
        cross_scetion,
        quantum_yield,
        actinic_flux
):
    integrand = (
        cross_scetion
        * quantum_yield
        * actinic_flux
    )

    return xp.trapezoid(
        integrand,
        x=wavelength,
        axis=-1
    )

def calculate_photolysis_rates(
        channels,
        actinic_flux
):
    return {
        channel.reaction_id: calculate_photolysis_rate(
            wavelength=channel.wavelength,
            cross_scetion=channel.cross_section,
            quantum_yield=channel.quantum_yield,
            actinic_flux=actinic_flux
        )
        for channel in channels
    }