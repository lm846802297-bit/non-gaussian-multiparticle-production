#!/usr/bin/env bash
set -euo pipefail
: "${ARCHIVED_SRC:?Set ARCHIVED_SRC to the extracted baseline src directory}"
cd "$(dirname "$(realpath "$0")")"
export PYTHONPATH="$ARCHIVED_SRC" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p out
for args in '5 6 20 0.1' '6 6 20 0.1' '6 6 24 0.1' '5 8 20 0.1' '6 8 20 0.1' '6 8 24 0.1' '5 6 20 0.025' '5 6 20 0.0125' '5 6 20 0' '5 7 20 0.1';do
  set -- $args
  key="k$1-n$2-L$3-lam$4"
  python run_out_modes.py --K "$1" --N "$2" --L "$3" --lam "$4" --out "out/$key.json" > "out/$key.log" 2>&1
done
python run_out_modes.py --K 5 --N 6 --L 20 --lam .1 --step-scale .5 --rtol 2e-11 --atol 2e-13 --out out/k5-n6-L20-lam0.1-tight.json > out/k5-n6-L20-lam0.1-tight.log 2>&1
python validate_out_modes.py > summary.json
sha256sum README.md replay.sh run_out_modes.py validate_out_modes.py summary.json out/* > SHA256SUMS.txt
sha256sum -c SHA256SUMS.txt
