#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$ROOT_DIR/env/bin/python}"
CATALOG_JSON="${CATALOG_JSON:-$ROOT_DIR/src/planet_params/exoplanet_catalog.json}"
EXOPLANET_NAME="${EXOPLANET_NAME:-TRAPPIST-1 e}"
ENVIRONMENT="${ENVIRONMENT:-volcanic_rock}"
OUTPUT_DIR="${OUTPUT_DIR:-$ROOT_DIR/TESTS/EXOPLANETS/${EXOPLANET_NAME// /_}_${ENVIRONMENT}}"
INPUT_JSON="${INPUT_JSON:-$OUTPUT_DIR/input.json}"
CALC_TYPE="${CALC_TYPE:-set_initial_state}"
ATMOSPHERE_MASS_FRACTION="${ATMOSPHERE_MASS_FRACTION:-1e-6}"

mkdir -p "$OUTPUT_DIR"

"$PYTHON_BIN" - "$CATALOG_JSON" "$EXOPLANET_NAME" "$OUTPUT_DIR" "$ENVIRONMENT" "$ATMOSPHERE_MASS_FRACTION" > "$INPUT_JSON" <<'PY'
import json
import sys

catalog_path, planet_name, output_dir, environment, f_atm = sys.argv[1:]

with open(catalog_path) as handle:
    catalog = json.load(handle)

try:
    planet = catalog["planets"][planet_name]
except KeyError:
    known = ", ".join(sorted(catalog["planets"]))
    raise SystemExit(f"Unknown EXOPLANET_NAME={planet_name!r}. Known planets: {known}")

earth_mass_kg = 5.9722e24
earth_radius_m = 6.371e6
solar_radius_m = 6.957e8
solar_mass_kg = 1.9885e30

chemical_species = [
    "H", "H2", "He", "O", "O2", "OH", "H2O", "H2O2",
    "C", "CO", "CO2", "CH4",
    "N", "N2", "NH3", "NO", "NO2", "N2O", "HCN",
]

# The catalog supplies physical planet/star data only. This chemistry block is
# a runnable primitive-atmosphere placeholder, not an observed composition.
primitive_atomic_abundances = {
    "H": 0.20,
    "He": 0.0,
    "C": 0.03,
    "N": 1.0,
    "O": 0.08,
}

input_data = {
    "working_dir": output_dir,
    "planet_model": "custom",
    "environment": environment,
    "catalyst_set": {"distribution": "gaussian", "center": 20, "std": 1.0, "set_size": 10},
    "metabolites_data": {
        "type": "binary",
        "pol_strng_maxsize": 4,
        "metabolites_distr_type": "length_decay",
        "decay_const": 2.0,
        "initial_population_molecules": 1000,
    },
    "QSP_size": 1,
    "protocell_data": {
        "radius": {"value": 1.0, "units": "micrometer"},
        "n_shells": 100,
    },
    "environment_data": {
        "number_pores": 10,
        "pore_radius": {"value": 1.0, "units": "millimeter"},
        "pore_height": {"value": 10.0, "units": "millimeter"},
        "distance_neigh_pores": {"value": 5.0, "units": "millimeter"},
        "temperature": {"value": planet["pl_eqt"], "units": "K"},
        "pressure": {"value": 1.0, "units": "bar"},
        "solvent_data": {
            "name": "H2O",
            "liquid_level_params": {"type": "constant"},
            "density": {"value": 997.0, "units": "kg / m^3"},
            "dynamic_viscosity": {"value": 8.9e-4, "units": "Pa * s"},
            "dielectric_constant": 78.4,
            "diffusion_scale": 1.0,
            "polarity": 1.0,
        },
    },
    "planetary_data": {
        "name": planet_name,
        "planet_radius": {"value": planet["pl_rade"] * earth_radius_m, "units": "m"},
        "planet_mass": {"value": planet["pl_bmasse"] * earth_mass_kg, "units": "kg"},
        "orbital_distance": {"value": planet["pl_orbsmax"], "units": "AU"},
        "rotation_period": {"value": planet["pl_orbper"], "units": "day"},
        "obliquity": 0.0,
        "eccentricity": 0.0,
        "tidal_locked": False,
        "day_night_contrast": 0.0,
        "atmosphere": {
            "n_layers": 250,
            "z_max": {"value": 100.0, "units": "km"},
            "atmosphere_mass_fraction": float(f_atm),
            "max_iter_loop": 100,
            "rel_tol": 1.0e-6,
            "abs_tol": 0.0,
            "damping_loop": 1.0,
        },
        "exo_chemistry": {
            "mode": "layered_equilibrium",
            "chemical_species": chemical_species,
            "atomic_abundances": primitive_atomic_abundances,
        },
    },
    "stellar_data": {
        "name": planet["hostname"],
        "star_temperature": {"value": planet["st_teff"], "units": "K"},
        "star_radius": {"value": planet["st_rad"] * solar_radius_m, "units": "m"},
        "star_mass": {"value": planet["st_mass"] * solar_mass_kg, "units": "kg"},
    },
    "distribution_rates": {"mean": 0.0, "std": 0.05},
    "r_mut": 0.005,
    "evol_params": {"T": 1.0, "dt": 0.001},
}

print(json.dumps(input_data, indent=4))
PY

"$PYTHON_BIN" -m src -ct "$CALC_TYPE" -json_input "$INPUT_JSON"
