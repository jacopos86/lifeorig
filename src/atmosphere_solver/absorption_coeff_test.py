from hapi import db_begin, fetch

# database
db_begin("hitran_data")

# Example: CO2 around 2300 cm^-1
# HITRAN molecule ID for CO2 = 2, isotope = 1

fetch("CO2", 2, 1, 2200.0, 2400.0)