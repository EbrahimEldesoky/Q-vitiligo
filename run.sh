#!/usr/bin/env bash
# ==============================================================================
# QuantumVitiligo Master Execution Runner
# ==============================================================================
# Usage:
#   ./run.sh                  # Run standard Aer simulation (30 generations)
#   ./run.sh --backend aer    # Run Aer local simulator
#   ./run.sh --backend ibm    # Run IBM Quantum cloud hardware
#   ./run.sh test             # Run test suite via pytest
#   ./run.sh --help           # Show all command-line options
# ==============================================================================

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# Activate virtual environment if present
if [ -f "${SCRIPT_DIR}/.venv/bin/activate" ]; then
    source "${SCRIPT_DIR}/.venv/bin/activate"
    PYTHON_EXEC="${SCRIPT_DIR}/.venv/bin/python"
elif command -v python3 &>/dev/null; then
    PYTHON_EXEC="python3"
else
    echo "[ERROR] Python 3 executable not found." >&2
    exit 1
fi

# Subcommand: test
if [ "$1" = "test" ]; then
    shift
    echo "[INFO] Running test suite via pytest..."
    exec "${PYTHON_EXEC}" -m pytest tests/ -v "$@"
fi

# Default execution if no arguments supplied
if [ $# -eq 0 ]; then
    echo "[INFO] No arguments specified. Starting standard discovery run (Aer backend, 30 generations)..."
    exec "${PYTHON_EXEC}" "${SCRIPT_DIR}/run_discovery.py" --backend aer --max-iter 30 --qubits 4
fi

# Forward all arguments directly to the Python runner
exec "${PYTHON_EXEC}" "${SCRIPT_DIR}/run_discovery.py" "$@"
