#!/bin/bash
# Checks that the active environment can run modules 02 and 04 (the panel measure/call/retro/parity, the pool analyses, the downloads).
#   conda activate enadl && source tools/env.sh && bash tools/check_env.sh
bad=0
echo "tools:"
for t in samtools minimap2 mosdepth gatk aria2c fasterq-dump vdb-validate curl rclone; do
  if p=$(command -v $t); then printf "  %-14s %s\n" $t "$p"; else printf "  %-14s MISSING\n" $t; bad=1; fi
done
echo "samtools: $(samtools --version 2>/dev/null | head -1)   rclone: $(rclone version 2>/dev/null | head -1)"
if rclone listremotes 2>/dev/null | grep -q '^ptm_dropbox:'; then echo "  rclone remote ptm_dropbox: found"; else echo "  rclone remote ptm_dropbox: NOT FOUND (the panel does not need it, retro, the pools and the FileManager do)"; bad=1; fi
echo "python $(python -V 2>&1 | cut -d' ' -f2):"
python - <<'PY' || bad=1
import importlib, sys
need, optional = ('numpy', 'pandas', 'yaml', 'scipy', 'matplotlib', 'openpyxl'), ('parasail', 'gspread')
fail = 0
for m in need + optional:
    try:
        mod = importlib.import_module(m); print(f"  {m:11s} {getattr(mod, '__version__', 'ok')}" + ("  (optional)" if m in optional else ""))
    except ImportError:
        hint = {'parasail': '  (optional: ~100x faster transposon measurement)', 'gspread': '  (optional: only the sheet readers need it)'}.get(m, '')
        print(f"  {m:11s} MISSING" + hint); fail += m in need
try:
    from sd_genotypes.panel import measure as pm
    print("  sd_genotypes panel modules:", ", ".join(sorted(pm.MODULES)))
except Exception as e:
    print("  sd_genotypes NOT importable:", e, "(run `source tools/env.sh` from the repo root first)"); fail += 1
sys.exit(1 if fail else 0)
PY
[ $bad = 0 ] && echo "environment OK" || echo "environment INCOMPLETE (see MISSING above)"
exit $bad
