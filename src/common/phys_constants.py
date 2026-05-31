import numpy as np
from scipy.constants import physical_constants
from scipy.constants import electron_volt
from scipy.constants import speed_of_light
from src.common.units import h, Q_
#
# physical constants module
#
kb = Q_(physical_constants["Boltzmann constant"][0] / electron_volt, "eV / K")
#
G = Q_(physical_constants["Newtonian constant of gravitation"][0], "m^3 / kg / s^2")
#
c = Q_(speed_of_light, "m / s")
#
hbar = h / (2. * np.pi)
#
