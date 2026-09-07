import pint

ureg = pint.UnitRegistry()
ureg.define("AU = astronomical_unit")
Q_ = ureg.Quantity  # type: ignore[misc]

# Define aliases you use often
THz = ureg.terahertz
ps  = ureg.picosecond
ns  = ureg.nanosecond
fs  = ureg.femtosecond
us  = ureg.microsecond
K   = ureg.kelvin
Ang = ureg.angstrom

# Planck constant
h   = ureg.planck_constant

internal_units = {
    "time": "s",
    "energy": "joule",
    "pressure": "Pa",
    "density": "kg / m^3",
    "number_density": "1 / m^3",
    "temperature": "K",
    "length": "m",
    "gravity": "m / s^2",
    "mass": "kg",
    "wavelength": "nanometer",
    "spectral_radiance": "W / m^3 / steradian",
    "spectral_flux": "W / m^3"
}

#
# return magnitude value
#

def get_magnitude(value: Q_, units: str, default_value: Q_):
    if value is None:
        value = default_value
    value = value.to(units)
    return float(value.magnitude)