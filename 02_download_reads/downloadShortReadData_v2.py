"""Wrapper under the old name: `python downloadShortReadData_v2.py ...` = `python -m download_reads run ...` (kept for one release; see MIGRATE.md)."""
import os, sys
sys.path[:0] = [os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'common'), os.path.dirname(os.path.abspath(__file__))]
from download_reads.run import main
sys.exit(main())
