import pint

ureg = pint.UnitRegistry()
ureg.define("AU = astronomical_unit")
Q_ = ureg.Quantity

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
