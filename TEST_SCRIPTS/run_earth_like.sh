#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$ROOT_DIR/env/bin/python}"
INPUT_JSON="${INPUT_JSON:-$ROOT_DIR/TESTS/2/input.json}"
CALC_TYPE="${CALC_TYPE:-set_initial_state}"

"$PYTHON_BIN" -m src -ct "$CALC_TYPE" -json_input "$INPUT_JSON"
