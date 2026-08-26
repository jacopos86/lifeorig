#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$ROOT_DIR/.cache/matplotlib}"
export LIFEORIG_REACTION_DB_HOST="${LIFEORIG_REACTION_DB_HOST:-localhost}"
export LIFEORIG_REACTION_DB_NAME="${LIFEORIG_REACTION_DB_NAME:-lifeorig_reactions}"
export LIFEORIG_REACTION_DB_USER="${LIFEORIG_REACTION_DB_USER:-lifeorig}"
if [ -z "${LIFEORIG_REACTION_DB_PASSWORD:-}" ]; then
    echo "ERROR: set LIFEORIG_REACTION_DB_PASSWORD before running Titan network import." >&2
    echo "Example: LIFEORIG_REACTION_DB_PASSWORD='...' $0" >&2
    exit 2
fi
export LIFEORIG_REACTION_DB_PASSWORD
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
        "pol_strng_maxsize": 12,
        "initial_population_molecules": 1000
    },
    "catalyst_set" : {"distribution": "gaussian", "center": 20, "std": 1.0, "set_size": 10},
    "QSP_size" : 1,
    "planet_model": "Titan",
    "environment": "surface_pond",
    "time_grid": {
        "start": {
            "value": 0.0,
            "units": "day"
        },
        "end": {
            "value": 60.0,
            "units": "day"
        },
        "dt": {
            "value": 0.01,
            "units": "day"
        }
    },
    "protocell_data" : {
        "n_shells": 100
    },
    "environment_data" : {
        "pool_geometry": {
            "height": {
                "value": 1.0,
                "units": "centimeter"
            },
            "n_grid_cells": 500,
            "surface_area_z0": {
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
            "composition": {
                "CH4": 0.70,
                "C2H6": 0.20,
                "N2": 0.10
            }
        },
        "liquid_level_params": {
            "base_level": {
                "value": 0.7,
                "units": "centimeter"
            }
        },
        "external_forces": {
            "rainfall": {
                "type": "sinusoidal",
                "base_level": {
                    "value": 0.0,
                    "units": "centimeter / day"
                },
                "amplitude": {
                    "value": 0.05,
                    "units": "centimeter / day"
                },
                "period": {
                    "value": 15.945,
                    "units": "day"
                },
                "phase": 0.0
            },
            "evaporation": {
                "type": "constant",
                "base_level": {
                    "value": 0.01,
                    "units": "centimeter / day"
                }
            }
        }
    },
    "distribution_rates": {"mean": 0.0, "std": 0.05},
    "r_mut" : 0.005
}
EOF

"$PYTHON_BIN" -m src -ct "$CALC_TYPE" -json_input "$INPUT_JSON"
