#!/usr/bin/env bash
# Native make check helper: independently compare the portable C Schur core
# with a pure Python dense reference. No external numerical library required.
set -euo pipefail
if ! command -v python3 >/dev/null 2>&1; then
  echo 'SKIP: stdlib differential check requires Python 3'
  exit 77
fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PYTHONPATH= python3 "$SCRIPT_DIR/../tools/differential_schur_stdlib.py"
