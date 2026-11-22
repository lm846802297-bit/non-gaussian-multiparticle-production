From the repository root:
```sh
repo="$PWD"
stage="$repo/research/ctp-autonomous-research/stage-01-k5-fock"
export PYTHONPATH="$repo/src:$stage/scripts"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
work="$(mktemp -d)"
mkdir -p "$work/run"
cp "$stage"/scripts/*.py "$work/"
cp "$stage"/results/*.json "$work/run/"
cd "$work"
python -m unittest -v test_scan test_independent_source
python verify_runs.py
python scan.py 5 6 run/replay-k5-n6.json
python scan_segmented.py 5 7 run/replay-k5-n7.json
python independent_source.py --N 7 --output run/replay-source-K5N7.json
```
Replays write new outputs. Wall-clock metadata is not expected to be identical.
The software versions and CPU class are recorded in results/hardware.txt.
