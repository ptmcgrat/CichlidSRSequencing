"""Wrapper under the old name: `python downloadManifest.py ...` = `python -m download_reads manifest ...` (kept for one release; see MIGRATE.md)."""
import os, sys
sys.path[:0] = [os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'common'), os.path.dirname(os.path.abspath(__file__))]
from download_reads.manifest_tool import main
main()
