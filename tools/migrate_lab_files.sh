#!/usr/bin/env bash
# Retire the old copies of what now lives in the numbered modules. Dry run unless --execute. Run it only when NO download or alignment job is running on this machine (it checks),
# and only after `python tools/run_tests.py` is green here and every other machine has pulled this layout (see MIGRATE.md).
#
#   tools/migrate_lab_files.sh check [--tests]
#   tools/migrate_lab_files.sh retire-old-download [--execute]      the old download modules, scripts and tests of cichlid_sr_sequencing/ (now in 02_download_reads/)
#   tools/migrate_lab_files.sh retire-yhsd [--execute]              cichlid_sr_sequencing/YH_SD_Paper -> legacy/YH_SD_Paper_v34 (now in common/ and modules 01, 04-07, 10)
#   tools/migrate_lab_files.sh retire-old-scripts [--execute]       downloadShortReadData.py and grabENA.py -> legacy/
#   tools/migrate_lab_files.sh upgrade-file-manager [--execute]     make the FileManager methods import the manifest from 02_download_reads (keeps the old place as the fallback)
# Module 03 (the alignment pipeline) is NOT moved by this script: see MIGRATE.md.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cmd="${1:-help}"; EXEC=0; [ "${2:-}" = "--execute" ] && EXEC=1
OLD=cichlid_sr_sequencing

OLD_DOWNLOAD=(
  $OLD/helper_modules/ena_fetch.py $OLD/helper_modules/sra_fetch.py $OLD/helper_modules/ena_scheduler.py $OLD/helper_modules/grab_ena_v2.py $OLD/helper_modules/download_manifest.py
  $OLD/downloadShortReadData_v2.py $OLD/downloadManifest.py $OLD/unit_scripts/ena_bench.py $OLD/unit_scripts/patch_file_manager.py
  $OLD/tests/test_ena_v2.py $OLD/tests/test_manifest_mode.py $OLD/tests/test_sra.py $OLD/tests/test_download_manifest.py $OLD/tests/MalinskyPaperData.csv
)
OLD_SCRIPTS=( $OLD/downloadShortReadData.py $OLD/grabENA.py )

act() { if [ "$EXEC" = 1 ]; then echo "+ $*"; "$@"; else echo "would: $*"; fi; }
busy() { [ "${MIGRATE_NO_BUSY_CHECK:-0}" = 1 ] && return 0; pgrep -af 'downloadShortReadData|download_reads|grab_ena|alignFastQReads|align_pipeline|ena_bench' | grep -v "$$" | grep -v 'migrate_lab_files' || true; }

if [ "$EXEC" = 1 ]; then
  b="$(busy)"; if [ -n "$b" ]; then echo "refusing: download or alignment processes are running on this machine:"; echo "$b"; exit 3; fi
  if [ -n "$(git -C "$REPO" status --porcelain --untracked-files=no)" ]; then echo "refusing: the working tree has uncommitted changes; commit or stash them first"; exit 4; fi
fi

case "$cmd" in
  check)
    if [ "${2:-}" = "--tests" ]; then echo "== tests"; python3 "$REPO/tools/run_tests.py" | tail -12; fi
    echo "== running jobs"; b="$(busy)"; [ -n "$b" ] && echo "$b" || echo "none"
    echo "== old files still present"; for f in "${OLD_DOWNLOAD[@]}" "${OLD_SCRIPTS[@]}"; do [ -e "$REPO/$f" ] && echo "  $f"; done; [ -d "$REPO/$OLD/YH_SD_Paper" ] && echo "  $OLD/YH_SD_Paper/"; true ;;
  retire-old-download)
    for f in "${OLD_DOWNLOAD[@]}"; do [ -e "$REPO/$f" ] && act git -C "$REPO" rm -q -- "$f"; done; true ;;
  retire-yhsd)
    if [ -d "$REPO/$OLD/YH_SD_Paper" ]; then act mkdir -p "$REPO/legacy"; act git -C "$REPO" mv "$OLD/YH_SD_Paper" legacy/YH_SD_Paper_v34; else echo "no $OLD/YH_SD_Paper"; fi ;;
  retire-old-scripts)
    act mkdir -p "$REPO/legacy"; for f in "${OLD_SCRIPTS[@]}"; do [ -e "$REPO/$f" ] && act git -C "$REPO" mv "$f" "legacy/$(basename "$f")"; done; true ;;
  upgrade-file-manager)
    export PYTHONPATH="$REPO/common:$REPO/02_download_reads${PYTHONPATH:+:$PYTHONPATH}"
    if [ "$EXEC" = 1 ]; then python3 -m download_reads patch-file-manager; else python3 -m download_reads patch-file-manager --dry-run; fi ;;
  *) sed -n 2,12p "${BASH_SOURCE[0]}" ;;
esac
