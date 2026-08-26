#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$ROOT_DIR/.cache/matplotlib}"
export LIFEORIG_REACTION_DB_HOST="${LIFEORIG_REACTION_DB_HOST:-localhost}"
export LIFEORIG_REACTION_DB_NAME="${LIFEORIG_REACTION_DB_NAME:-lifeorig_reactions}"
export LIFEORIG_REACTION_DB_USER="${LIFEORIG_REACTION_DB_USER:-lifeorig}"
if [ -z "${LIFEORIG_REACTION_DB_PASSWORD:-}" ]; then
    echo "ERROR: set LIFEORIG_REACTION_DB_PASSWORD before running Earth-like network import." >&2
    echo "Example: LIFEORIG_REACTION_DB_PASSWORD='...' $0" >&2
    exit 2
fi
export LIFEORIG_REACTION_DB_PASSWORD
mkdir -p "$MPLCONFIGDIR"

PYTHON_BIN="${PYTHON_BIN:-python}"
PLANET_MODEL="${PLANET_MODEL:-Earth}"
ENVIRONMENT="${ENVIRONMENT:-volcanic_rock}"
ENVIRONMENT_SOURCE="${ENVIRONMENT_SOURCE:-explicit}"
OUTPUT_DIR="${OUTPUT_DIR:-$ROOT_DIR/TESTS/${PLANET_MODEL}_${ENVIRONMENT}_network}"
INPUT_JSON="${INPUT_JSON:-$OUTPUT_DIR/input.json}"
CALC_TYPE="${CALC_TYPE:-set_initial_state}"

CHEM_NETWORK_TYPE="${CHEM_NETWORK_TYPE:-reference_file}"
CHEM_REACTION_FILE="${CHEM_REACTION_FILE:-EARTH_LIKE/earth_like_volcanic_rock_combined_network_v1.txt}"
MAX_POLYMER_LENGTH="${MAX_POLYMER_LENGTH:-12}"

mkdir -p "$OUTPUT_DIR"

cat > "$INPUT_JSON" <<EOF
{
    "working_dir": "$OUTPUT_DIR",
    "planet_model": "$PLANET_MODEL",
    "environment": "$ENVIRONMENT",
    "environment_source": "$ENVIRONMENT_SOURCE",

    "chemical_network": {
        "type": "$CHEM_NETWORK_TYPE",
        "reaction_file": "$CHEM_REACTION_FILE"
    },

    "metabolites_data": {
        "type": "reference_file",
        "reaction_file": "$CHEM_REACTION_FILE",
        "pol_strng_maxsize": $MAX_POLYMER_LENGTH,
        "initial_population_molecules": 1000
    },

    "environment_data": {
        "temperature": 288.15,
        "pressure": 101325.0
    },

    "catalyst_set": {
        "distribution": "gaussian",
        "center": 20,
        "std": 1.0,
        "set_size": 10
    },

    "QSP_size": 1,
    "protocell_data": {
        "radius": {
            "value": 1.0,
            "units": "micrometer"
        },
        "n_shells": 100
    },

    "distribution_rates": {
        "mean": 0.0,
        "std": 0.05
    },
    "r_mut": 0.005,
    "evol_params": {
        "T": 1.0,
        "dt": 0.001
    }
}
EOF

"$PYTHON_BIN" -m src -ct "$CALC_TYPE" -json_input "$INPUT_JSON"
