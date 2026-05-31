import easychem.easychem as ec
import numpy as np
from src.utilities.logging_module import log

# Temperature (K) and pressure (bar)

T = 1500.0
P = 1.0

# -------------------------------
# Initialize solver
# -------------------------------

exo = ec.ExoAtmos()

# -------------------------------
# Set elemental abundances
# -------------------------------

log.info("\t Default atoms: " + str(exo.atoms))

# elemental abundances (solar-like, simplified)

abundances = {
    "H": 1.0,
    "He": 0.1,
    "C": 1e-4,
    "O": 2e-1,
    "N": 1e-5
}

# map to array in correct order

abunds = np.array([abundances.get(a, 0.0) for a in exo.atoms], dtype=float)

exo.updateAtomAbunds(abunds)

# -------------------------------
# Solve equilibrium
# -------------------------------

exo.solve(P, T)

# -------------------------------
# Get results
# -------------------------------

res = exo.result_mol()

log.info("\n\n")
log.info("\t --- Selected species ---")
for sp in ["H2", "H2O", "CO", "CH4", "NH3", "HCN", "O2"]:
    try:
        log.info(f"\t {sp:4s} : " + str(res[sp]))
    except Exception:
        log.error(f"{sp:4s} : not present")