#!/usr/bin/env bash
# Run this on your LAPTOP, from anywhere inside your clone, after you have committed and pushed.
# Copies the committed YH_SD_Paper folder to utaka (no git is needed on utaka) and runs the self-test there.
#   bash YH_SD_Paper/scripts/deploy_to_utaka.sh                 # defaults below
#   ALLOW_DIRTY=1 bash YH_SD_Paper/scripts/deploy_to_utaka.sh   # deploy even with uncommitted changes (stamped as dirty)
#   REMOTE=user@host DEST=relative/path/in/home/ bash ...       # override
set -euo pipefail
REMOTE="${REMOTE:-mcgrath-lab@utaka.biosci.gatech.edu}"
DEST="${DEST:-CichlidSRSequencing/cichlid_sr_sequencing/YH_SD_Paper/}"       # relative to your home directory on utaka
CONDA_SH="${CONDA_SH:-miniconda3/etc/profile.d/conda.sh}"                    # relative to home on utaka
ENV_NAME="${ENV_NAME:-yh_sd_paper}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"                      # .../YH_SD_Paper
COMMIT="$(git -C "$HERE" rev-parse --short HEAD 2>/dev/null || echo nogit)"
DIRTY="$(git -C "$HERE" status --porcelain -- "$HERE" 2>/dev/null | wc -l | tr -d ' ')"
BRANCH="$(git -C "$HERE" rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
echo "Package: $HERE"
echo "Branch $BRANCH, commit $COMMIT, uncommitted changes in package: $DIRTY"
if [ "$DIRTY" != "0" ] && [ "${ALLOW_DIRTY:-0}" != "1" ]; then
  echo "ERROR: uncommitted changes. Commit (and push) first, or run with ALLOW_DIRTY=1."; exit 1
fi
FLAG=0; [ "$DIRTY" != "0" ] && FLAG=1
printf 'commit: %s\ndirty: %s\nbranch: %s\ndeployed: %s\nfrom: %s\n' "$COMMIT" "$FLAG" "$BRANCH" "$(date '+%Y-%m-%dT%H:%M:%S')" "$(hostname)" > "$HERE/DEPLOYED_FROM.txt"

if command -v rsync >/dev/null 2>&1; then
  rsync -az --delete --exclude '__pycache__' --exclude '*.pyc' --exclude '.git' "$HERE/" "$REMOTE:$DEST"
else
  echo "rsync not found; falling back to tar over ssh (does not delete files removed since the last deploy)"
  tar -C "$HERE" --exclude '__pycache__' --exclude '*.pyc' -cf - . | ssh "$REMOTE" "mkdir -p '$DEST' && tar -C '$DEST' -xf -"
fi
echo "Copied to $REMOTE:~/$DEST"
ssh "$REMOTE" "bash -lc 'source ~/$CONDA_SH && conda activate $ENV_NAME && cd ~/$DEST && cat DEPLOYED_FROM.txt | head -2 && python -m yhsd.cli selftest 2>&1 | tail -1'"
