#!/usr/bin/env bash
set -euo pipefail
: "${GATE_INPUT_DIR:?Set GATE_INPUT_DIR to the prepared gate native-input directory}"
d="$GATE_INPUT_DIR"
python "$(dirname "$(realpath "$0")")/gate_audit.py" \
  --cf "$d/k5-n6-cf4.json" --direct "$d/k5-n6-direct-moments.json" \
  --dop "$d/k5-n6.json" --direct-dop-reference "$d/k5-n6-segmented.json" \
  --cf-runner "$d/gate_cf4.py" --direct-runner "$d/direct_moments.py" \
  --dop-runner "$d/scan.py" --pinned-module "$d/sd_residual.py" \
  --restored-fast-ops "$d/phase7_fast_ops.py" --restored-archive "$d/phase7-cutoff-scan-001.tar.xz" \
  --fock-context "$d/n5-to-n6-fock-context-A.json" \
  --out "${1:-gate_validation.json}"
