import numpy as np

from src.common.units import Q_
from src.utilities.logging_module import log


def middle_hill_profile(
        x,
        length=Q_(10.0, "centimeter"),
        hill_height=Q_(1.5, "millimeter"),
        hill_width=Q_(1.0, "centimeter"),
    ):
    """
    1D pool-floor elevation with a smooth hill in the middle.

    The returned value is the local floor elevation above the basin floor.
    A liquid level below the hill separates the pool into two regions; a
    higher liquid level connects them.
    """
    x_cm = x.to("centimeter").magnitude if hasattr(x, "to") else np.asarray(x)
    center_cm = 0.5 * length.to("centimeter").magnitude
    width_cm = hill_width.to("centimeter").magnitude
    height_mm = hill_height.to("millimeter").magnitude
    profile_mm = height_mm * np.exp(-0.5 * ((x_cm - center_cm) / width_cm) ** 2)
    return Q_(profile_mm, "millimeter")


POOL_SPATIAL_PROFILE_FUNCTIONS = {
    "middle_hill_profile": middle_hill_profile,
}


def get_pool_spatial_profile_function(name: str):
    try:
        return POOL_SPATIAL_PROFILE_FUNCTIONS[name]
    except KeyError:
        log.error(f"Unknown pool spatial profile function: {name}")
