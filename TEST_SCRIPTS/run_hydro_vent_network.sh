#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$ROOT_DIR/env/bin/python}"

ENVIRONMENT="${ENVIRONMENT:-hydro_vent}"
ENVIRONMENT_SOURCE="${ENVIRONMENT_SOURCE:-explicit}"
MOLECULE_SOURCE="${MOLECULE_SOURCE:-chemical_network}"
PLANET_MODEL="${PLANET_MODEL:-Earth}"
OUTPUT_DIR="${OUTPUT_DIR:-$ROOT_DIR/TESTS/Earth_${ENVIRONMENT}}"
INPUT_JSON="${INPUT_JSON:-$OUTPUT_DIR/input.json}"
CALC_TYPE="${CALC_TYPE:-set_initial_state}"

CHEM_NETWORK_TYPE="${CHEM_NETWORK_TYPE:-reference_file}"
CHEM_REACTION_FILE="${CHEM_REACTION_FILE:-$ROOT_DIR/reference_reactions/HYDRO_VENT/hydro_vent_reaction_network_v1.txt}"

mkdir -p "$OUTPUT_DIR"

cat > "$INPUT_JSON" <<EOF
{
    "working_dir": "$OUTPUT_DIR",
    "planet_model": "$PLANET_MODEL",
    "environment": "$ENVIRONMENT",
    "environment_source": "$ENVIRONMENT_SOURCE",
    "molecule_source": "$MOLECULE_SOURCE",

    "chemical_network": {
        "type": "$CHEM_NETWORK_TYPE",
        "reaction_file": "$CHEM_REACTION_FILE"
    },

    "metabolites_data": {
        "type": "reference_file",
        "reaction_file": "$CHEM_REACTION_FILE",
        "initial_population_molecules": 1000
    },

    "environment_data": {
        "n_cells": 20,
        "length": {
            "value": 1.0,
            "units": "meter"
        },
        "cross_section_area": {
            "value": 1.0e-4,
            "units": "meter ** 2"
        },
        "T_hot": {
            "value": 373.15,
            "units": "kelvin"
        },
        "T_cold": {
            "value": 283.15,
            "units": "kelvin"
        },
        "pressure": {
            "value": 3.0e7,
            "units": "pascal"
        },
        "pH_hot": 10.5,
        "pH_cold": 6.0,
        "flow_velocity": {
            "value": 1.0e-5,
            "units": "meter / second"
        },
        "porosity": 0.35,
        "surface_area_per_volume": {
            "value": 1.0e6,
            "units": "1 / meter"
        },
        "vent_concentrations": {
            "H2": 1.0e-3,
            "CO2": 1.0e-3,
            "CH4": 1.0e-5,
            "NH3": 1.0e-5,
            "H2S": 1.0e-4
        },
        "ocean_concentrations": {
            "H2": 1.0e-8,
            "CO2": 1.0e-4,
            "CH4": 1.0e-8,
            "NH3": 1.0e-8,
            "H2S": 1.0e-8
        },
        "mineral_fractions": {
            "FeS": 0.5,
            "NiS": 0.2,
            "clay": 0.3
        },
        "n_pores_per_cell": 0
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
