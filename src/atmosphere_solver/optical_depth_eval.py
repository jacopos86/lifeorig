from dataclasses import dataclass
import numpy as np
from src.utilities.logging_module import log

#
#   optical depth data class: contains UV / IR
#

@dataclass
class OpticalDepthProfiles:
    tau_UV: np.ndarray
    tau_IR: np.ndarray

def absorption_optical_depth_profile(z, absorption_coefficient):
    #
    #   tau_lambda(z) = int_z^ztop dz' absorption_coefficient_lambda(z')
    #
    z = np.asarray(z, dtype=float)
    absorption_coefficient = np.asarray(absorption_coefficient, dtype=float)
    tau_abs = np.zeros_like(absorption_coefficient)
    for ilayer in range(z.size - 2, -1, -1):
        dz = z[ilayer + 1] - z[ilayer]
        tau_abs[:, ilayer] = (
            tau_abs[:, ilayer + 1]
            + 0.5*(
                absorption_coefficient[:, ilayer]
                + absorption_coefficient[:, ilayer + 1]
            )*dz
        )
    log.info(f"\t max absorption optical depth at surface: {np.max(tau_abs[:, 0]):.3e}")
    return tau_abs
