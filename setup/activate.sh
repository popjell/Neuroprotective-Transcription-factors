# Source this file: source setup/activate.sh
NEURO_PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export R_LIBS_USER="$NEURO_PROJECT_ROOT/.Rlibs"
export RETICULATE_PYTHON="$NEURO_PROJECT_ROOT/.venv/bin/python"
export PATH="$NEURO_PROJECT_ROOT/.venv/bin:$PATH"
