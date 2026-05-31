#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$ROOT_DIR/env/bin/python}"
OUTPUT_DIR="${OUTPUT_DIR:-$ROOT_DIR/TESTS/TITAN}"
INPUT_JSON="${INPUT_JSON:-$OUTPUT_DIR/input.json}"
CALC_TYPE="${CALC_TYPE:-set_initial_state}"

mkdir -p "$OUTPUT_DIR"

cat > "$INPUT_JSON" <<EOF
{
    "working_dir" : "$OUTPUT_DIR",
    "catalyst_set" : {"distribution": "gaussian", "center": 20, "std": 1.0, "set_size": 10},
    "QSP_size" : 1,
    "planet_model": "Titan",
    "protocell_data" : {
        "n_shells": 100
    },
    "environment_data" : {
        "pool_spatial_profile": {
            "type": "analytic_1d",
            "domain_length": {
                "value": 1.0,
                "units": "centimeter"
            },
            "bottom_profile": {
                "function": "middle_hill_profile"
            }
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
