#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$ROOT_DIR/env/bin/python}"
PLANET_MODEL="${PLANET_MODEL:-Earth}"
ENVIRONMENT="${ENVIRONMENT:-volcanic_rock}"
OUTPUT_DIR="${OUTPUT_DIR:-$ROOT_DIR/TESTS/${PLANET_MODEL}_${ENVIRONMENT}}"
INPUT_JSON="${INPUT_JSON:-$OUTPUT_DIR/input.json}"
CALC_TYPE="${CALC_TYPE:-set_initial_state}"

mkdir -p "$OUTPUT_DIR"

cat > "$INPUT_JSON" <<EOF
{
    "working_dir" : "$OUTPUT_DIR",
    "planet_model": "$PLANET_MODEL",
    "environment": "$ENVIRONMENT",
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
            "name": "H2O",
            "liquid_level_params": {
                "type": "constant"
            },
            "density": {
                "value": 997.0,
                "units": "kg / m^3"
            },
            "dynamic_viscosity": {
                "value": 8.9e-4,
                "units": "Pa * s"
            },
            "dielectric_constant": 78.4,
            "diffusion_scale": 1.0,
            "polarity": 1.0
        }
    },
    "distribution_rates": {"mean": 0.0, "std": 0.05},
    "r_mut" : 0.005,
    "evol_params": {"T" : 1.0, "dt" : 0.001}
}
EOF

"$PYTHON_BIN" -m src -ct "$CALC_TYPE" -json_input "$INPUT_JSON"
