from dataclasses import dataclass
from src.common.units import Q_
from src.stellar_params.black_body_radiation import blackbody_B_lambda
from src.utilities.logging_module import log

#
#   data class for stellar parameters
#

@dataclass
class StellarParams:
    name: str
    spectral_class: str
    effective_temperature: Q_
    radius: Q_
    mass: Q_ | None = None
    luminosity: Q_ | None = None
    # summary
    def log_summary(self):
        log.info("\t --- Stellar Parameters ---")
        log.info(f"\t star_name          : {self.name}")
        log.info(f"\t spectral_class     : {self.spectral_class}")
        log.info(f"\t star_temperature   : {self.effective_temperature}")
        log.info(f"\t star_radius        : {self.radius}")
        if self.mass is not None:
            log.info(f"\t star_mass          : {self.mass}")
        if self.luminosity is not None:
            log.info(f"\t star_luminosity    : {self.luminosity}")
    # B_lambda(T)
    def B_lambda(self, wavelength: Q_) -> Q_:
        return blackbody_B_lambda(
            wavelength=wavelength,
            temperature=self.effective_temperature
        )

#
#   SOLAR PARAMETERS
#

SOLAR_PARAMS = StellarParams(
    name="Sun",
    spectral_class="G2V",
    effective_temperature=Q_(5778.0, "K"),
    radius=Q_(6.957e8, "m"),
    mass=Q_(1.9885e30, "kg"),
    luminosity=None,
)
