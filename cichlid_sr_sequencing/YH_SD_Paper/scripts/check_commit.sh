#!/usr/bin/env bash
# Run on the LAPTOP before `git add`/commit: lists package files that git would silently skip because an ignore rule matches them.
#   bash YH_SD_Paper/scripts/check_commit.sh
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IGN="$(git -C "$HERE" status --porcelain --ignored -- "$HERE" 2>/dev/null | grep '^!!' | grep -v -E '__pycache__|DEPLOYED_FROM\.txt|\.pyc$' || true)"
MISSING=""
for f in yhsd/cli.py yhsd/common.py config/paper.yaml; do
  [ -f "$HERE/$f" ] || MISSING="$MISSING $f"
done
if [ -n "$MISSING" ]; then echo "ERROR: missing from the folder:$MISSING"; exit 2; fi
if [ -n "$IGN" ]; then
  echo "WARNING: these package files are IGNORED by git and would NOT be committed:"; echo "$IGN"
  echo "Find the rule:   git check-ignore -v <file>"
  echo "Force-add:       git add -f <file>      (or fix the ignore rule)"; exit 1
fi
echo "OK: no package files are being ignored by git"
