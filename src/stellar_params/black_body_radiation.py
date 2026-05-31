import numpy as np
from src.common.units import Q_
from src.common.phys_constants import c, h, kb

#
#   compute black body radiation
#

def blackbody_B_lambda(wavelength: Q_, temperature: Q_) -> Q_:
    """
    Planck spectral radiance B_lambda(T).

    Returns intensity per unit wavelength per steradian.
    """
    # compute B_lambda(T)
    exponent = (h * c / (wavelength * kb * temperature)).to("").magnitude
    return (2.0 * h * c**2 / wavelength**5 / (np.exp(exponent) - 1.0)).to(
        "W / m^3 / steradian"
    )
