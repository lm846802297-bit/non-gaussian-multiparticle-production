#!/usr/bin/env bash
set -euo pipefail
stage="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
research="$(dirname "$stage")"
if [ "$#" -ne 1 ]; then echo "Usage: bash replay-audits.sh NEW_OUTPUT_DIRECTORY" >&2; exit 2; fi
out="$1"
mkdir -p "$out"
out="$(cd "$out" && pwd)"
test ! -e "$out/high_verification.json"
test ! -e "$out/high_signed_audit.json"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
mkdir "$tmp/run"
ln -s "$stage/next_scan.py" "$tmp/next_scan.py"
for name in k5-n6-C k5-n7-C k6-n6-A k6-n6-A-L24; do
  path="$research/stage-02-ctp-regulators/scan/run/$name.json"
  test -s "$path"
  ln -s "$path" "$tmp/run/$name.json"
done
for name in k6-n7-A k6-n7-A-L24; do
  path="$stage/run/$name.json"
  test -s "$path"
  ln -s "$path" "$tmp/run/$name.json"
done
base="$research/stage-01-k5-fock/results"
test -s "$base/k5-n6-segmented.json"
test -s "$base/k5-n7.json"
(cd "$tmp" && CTP_BASE_RESULTS="$base" python "$stage/verify_high.py" > "$out/high_verification.json")
python "$stage/audit_high.py" --base "$base" --main "$tmp/run" --out "$out/high_signed_audit.json"
cmp "$out/high_verification.json" "$stage/high_verification.json"
cmp "$out/high_signed_audit.json" "$stage/high_signed_audit.json"
