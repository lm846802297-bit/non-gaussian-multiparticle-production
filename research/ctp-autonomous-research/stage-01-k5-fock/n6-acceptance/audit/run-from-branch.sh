#!/usr/bin/env bash
set -euo pipefail
: "${1:?Pass the exact Phase VII .tar.xz archive}"
archive="$(realpath "$1")"
out="$(realpath -m "${2:-gate-verify-replay.json}")"
here="$(cd "$(dirname "$0")" && pwd)"
gate="$(dirname "$here")"
legacy="${GATE_DOP_STAGE:-$(dirname "$gate")}"
opt="${GATE_PHASE7_SRC:-$legacy/phase7-reconciled/src}"
pinned="${GATE_PINNED_SRC:-$(git -C "$here" rev-parse --show-toplevel)/src}"
w="$(mktemp -d)"
trap 'rm -rf "$w"' EXIT
mkdir "$w/input"
cp "$here"/{gate_audit.py,test_gate_audit.py,replay.sh} "$w/"
ln -s "$archive" "$w/input/phase7-cutoff-scan-001.tar.xz"
for f in "$legacy"/results/k5-n6{,-segmented}.json "$legacy/scripts/scan.py" \
         "$opt/phase7_fast_ops.py" "$pinned/sd_residual.py" \
         "$gate/cf4/"{gate_cf4.py,k5-n6-cf4.json} \
         "$gate/direct/"{direct_moments.py,k5-n6-direct-moments.json} \
         "$gate/context/n5-to-n6-fock-context-A.json"; do
 ln -s "$f" "$w/input/$(basename "$f")"
done
cd "$w"
python -m unittest -v test_gate_audit
GATE_INPUT_DIR="$w/input" bash replay.sh "$out"
