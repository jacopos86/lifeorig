#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$ROOT_DIR/.cache/matplotlib}"
mkdir -p "$MPLCONFIGDIR"
PYTHON_BIN="${PYTHON_BIN:-python}"
OUTPUT_DIR="${OUTPUT_DIR:-$ROOT_DIR/TESTS/ATMOSPHERE}"
INPUT_JSON="${INPUT_JSON:-$OUTPUT_DIR/input.json}"
CALC_TYPE="${CALC_TYPE:-env_solver}"

mkdir -p "$OUTPUT_DIR"

cat > "$INPUT_JSON" <<EOF
{
    "working_dir" : "$OUTPUT_DIR",
    "chemical_network": {
        "type": "atmospheric",
        "reaction_files": [
            "ATMOSPHERIC/VPL_ATMOS/PHOTOCHEM/INPUTFILES/TEMPLATES/Archean+haze/reactions.rx"
        ]
    },
    "catalyst_set" : {"distribution": "gaussian", "center": 20, "std": 1.0, "set_size": 10},
    "metabolites_data" : {
        "type": "binary",
        "pol_strng_maxsize" : 4,
        "metabolites_distr_type" : "length_decay",
        "decay_const": 2.0,
        "initial_population_molecules" : 1000
    },
    "QSP_size" : 1,
    "protocell_data" : {
        "radius": {
            "value": 1.0,
            "units": "micrometer"
        },
        "n_shells": 100
    },
    "environment": "volcanic_rock",
    "environment_data" : {
        "number_pores" : 10,
        "pore_radius" : {
            "value": 1.0,
            "units": "millimeter"
        },
        "pore_height" : {
            "value": 10.0,
            "units": "millimeter"
        },
        "distance_neigh_pores" :{
            "value": 5.0,
            "units": "millimeter"
        },
        "solvent_data": {
            "liquid_level_params": {
                "type": "constant"
            }
        }
    },
    "planetary_data" : {
        "planet_radius": { "value": 6371000.0, "units": "m" },
        "planet_mass": { "value": 5.972e24, "units": "kg" },
        "orbital_distance": { "value": 1.0, "units": "AU" },
        "rotation_period": { "value": 24.0, "units": "hour" },
        "obliquity": 0.4091,
        "eccentricity": 0.0167,
        "tidal_locked": false,
        "day_night_contrast": 0.0,
        "atmosphere": {
            "n_layers": 2500,
            "z_max": { "value": 1800.0, "units": "km" },
            "top_pressure": { "value": 1e-5, "units": "Pa"},
            "atmosphere_mass_fraction": 8.6e-7,
            "max_iter_loop": 50,
            "rel_tol": 1e-6,
            "abs_tol": { "value": 0.0, "units": "Pa" },
            "damping_loop": 1.0,
            "radiative_solver": {
                "type": "two_stream",
                "settings": {
                    "max_iter": 100,
                    "abs_tol": 1e-3,
                    "damping": 0.5
                }
            }
        },
        "exo_chemistry": {
            "mode": "${CHEMISTRY_MODE:-layered_equilibrium}",
            "chemical_species": [
                "H",
                "H2",
                "He",
                "O",
                "O2",
                "OH",
                "H2O",
                "H2O2",
                "C",
                "CO",
                "CO2",
                "CH",
                "CH2",
                "CH3",
                "CH4",
                "C2H2",
                "C2H4",
                "C2H6",
                "N",
                "N2",
                "NH",
                "NH2",
                "NH3",
                "NO",
                "NO2",
                "N2O",
                "HCN",
                "HCO",
                "H2CO",
                "CH3OH"
            ],
            "atomic_abundances": {
                "H": 0.20,
                "He": 1.0e-8,
                "C": 0.03,
                "N": 1.0,
                "O": 0.078
            }
        }
    },
    "stellar_data": {
        "spectral_class": "G2V",
        "star_temperature": { "value": 5778, "units": "K" },
        "star_radius": { "value": 6.957e8, "units": "m" },
        "star_mass": { "value": 1.9885e30, "units": "kg" }
    },
    "distribution_rates": {"mean": 0.0, "std": 0.05},
    "r_mut" : 0.005,
    "evol_params": {"T" : 1.0, "dt" : 0.001}
}
EOF

"$PYTHON_BIN" -m src -ct "$CALC_TYPE" -json_input "$INPUT_JSON"
