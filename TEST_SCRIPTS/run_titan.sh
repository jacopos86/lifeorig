#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$ROOT_DIR/.cache/matplotlib}"
mkdir -p "$MPLCONFIGDIR"
PYTHON_BIN="${PYTHON_BIN:-python}"
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
    PYTHON_BIN="/home/jacopo/miniforge3/envs/lifeorig_gpu_env/bin/python"
fi
OUTPUT_DIR="${OUTPUT_DIR:-$ROOT_DIR/TESTS/TITAN}"
INPUT_JSON="${INPUT_JSON:-$OUTPUT_DIR/input.json}"
CALC_TYPE="${CALC_TYPE:-set_initial_state}"
CHEM_NETWORK_TYPE="${CHEM_NETWORK_TYPE:-reference_file}"
CHEM_REACTION_FILE="${CHEM_REACTION_FILE:-$ROOT_DIR/reference_reactions/TITAN/titan_reaction_network_v0.txt}"

mkdir -p "$OUTPUT_DIR"

cat > "$INPUT_JSON" <<EOF
{
    "working_dir" : "$OUTPUT_DIR",
    "chemical_network": {
        "type": "$CHEM_NETWORK_TYPE",
        "reaction_file": "$CHEM_REACTION_FILE"
    },
    "metabolites_data": {
        "type": "reference_file",
        "reaction_file": "$CHEM_REACTION_FILE",
        "initial_population_molecules": 1000
    },
    "catalyst_set" : {"distribution": "gaussian", "center": 20, "std": 1.0, "set_size": 10},
    "QSP_size" : 1,
    "planet_model": "Titan",
    "protocell_data" : {
        "n_shells": 100
    },
    "environment_data" : {
        "pool_geometry": {
            "height": {
                "value": 1.0,
                "units": "centimeter"
            },
            "surface_area": {
                "value": 1.0,
                "units": "centimeter ** 2"
            },
            "surface_area_profile": {
                "type": "bottom_weighted",
                "decay_length": {
                    "value": 1.0,
                    "units": "millimeter"
                }
            },
            "reactive_area_to_surface_area_ratio": 0.1
        },
        "solvent_data": {
            "name": "CH4",
            "liquid_level_params": {
                "type": "sinusoidal"
            }
        }
    },
    "distribution_rates": {"mean": 0.0, "std": 0.05},
    "r_mut" : 0.005,
    "evol_params": {"T" : 1.0, "dt" : 0.001}
}
EOF

"$PYTHON_BIN" -m src -ct "$CALC_TYPE" -json_input "$INPUT_JSON"
