"""Read-only access to the lab's cloud data: Dropbox (through rclone) and the Google sheets.

This packages the parts of the lab FileManager that sexloci needs and nothing else. Sheets are read-only: Google is asked for
read-only scopes unless the object is created with read_only=False, which only `gsdf-genotype --write-sample-database` does (see
samplesheet.py, which can add one column and nothing else).

  * `list_files` / `download`     rclone lsf / rclone copy, mapping <local master dir> <-> <cloud master dir>
  * `read_sample_tables`          SampleDatabase and DNAReadsDatabase sheets
  * `read_alignment_table`        the alignment worksheet named after the genome version
  * `ensure_gff`                  genome annotation (downloaded on demand)

Authentication is lazy (first sheet read), so listing and downloading work without gspread installed.
`rclone_runner` and `_read_sheet` can be replaced, which is how the tests exercise this code without a network.
"""
from __future__ import annotations

import os
import platform
import subprocess
import time
from typing import Callable, List, Optional, Tuple

import pandas as pd

KNOWN_GENOMES = ["Mzebra_UMD2a", "Mzebra_GT3", "Mzebra_GT3_NCBI", "Mconophoros_GT1"]
# genome version -> annotation files inside <master>/Genomes/<genome version>/ (as defined by the lab FileManager)
GENOME_ANNOTATION = {"Mzebra_GT3_NCBI": {"gff": "genomic.gff", "gtf": "genomic.gtf"}}
CREDENTIAL_FILE = "cichlidsrsequencing_api_creds.json"
_BIG_FILE_RCLONE_FLAGS = ["--multi-thread-streams", "96", "--multi-thread-cutoff", "100Mi"]


def default_local_master_dir() -> str:
    """Same rule as the lab FileManager (lab servers use /Data/<user>/Temp/, everyone else ~/Temp/)."""
    if os.getenv("SEXLOCI_LOCAL_MASTER"):
        return os.environ["SEXLOCI_LOCAL_MASTER"]
    node = platform.node().lower()
    if node in ("ebb-utaka.biosci.gatech.edu", "utaka.biosci.gatech.edu") or any(
            k in node for k in ("utaka", "yellowhead", "burtoni")):
        return "/Data/" + os.getenv("USER", "") + "/Temp/CichlidSequencingData/"
    return os.getenv("HOME", os.path.expanduser("~")).rstrip("/") + "/Temp/CichlidSequencingData/"


def _subprocess_rclone(args: List[str]) -> Tuple[int, str, str]:
    try:
        p = subprocess.run(["rclone", *args], capture_output=True, encoding="utf-8")
    except FileNotFoundError:
        raise RuntimeError("rclone is not installed or not on PATH") from None
    return p.returncode, p.stdout, p.stderr


class LabData:
    def __init__(self, rclone_remote="ptm_dropbox:/", cloud_master_dir="COS/BioSci/BioSci-McGrath/Apps/CichlidSequencingData/",
                 local_master_dir: Optional[str] = None, sample_sheet_id: str = "", alignment_sheet_id: str = "",
                 read_only: bool = True, rclone_runner: Optional[Callable] = None, log: Callable = print):
        self.local_master_dir = (local_master_dir or default_local_master_dir()).rstrip("/") + "/"
        self.cloud_master_dir = rclone_remote + cloud_master_dir
        self.sample_sheet_id, self.alignment_sheet_id = sample_sheet_id, alignment_sheet_id
        self.read_only = read_only
        self._rclone = rclone_runner or _subprocess_rclone
        self.log = log
        self._gc = None
        self._sample_cache = None
        self.credential_file = self.local_master_dir + CREDENTIAL_FILE

    @classmethod
    def from_config(cls, cfg, **kw):
        return cls(cfg.rclone_remote, cfg.cloud_master_dir, cfg.local_master_dir, cfg.sample_sheet_id,
                   cfg.alignment_sheet_id, cfg.sheets_read_only, **kw)

    # ------------------------------------------------------------------ cloud files (rclone)
    def _cloud_path(self, local_path: str) -> str:
        if not local_path.startswith(self.local_master_dir):
            raise ValueError(f"{local_path} is not inside the local master directory {self.local_master_dir}")
        return local_path.replace(self.local_master_dir, self.cloud_master_dir, 1)

    def _run(self, args):
        rc, out, err = self._rclone(args)
        if rc != 0:
            raise RuntimeError(f"rclone {' '.join(args[:1])} failed (exit {rc}): {err.strip()[:300]}")
        return out

    def list_files(self, local_dir: str) -> List[str]:
        out = self._run(["lsf", self._cloud_path(local_dir.rstrip("/") + "/")])
        return [x for x in out.splitlines() if x and not x.endswith("/")]

    def download(self, local_path: str, fast: bool = False) -> str:
        """Fetch one cloud file to the mirrored local path. Raises FileNotFoundError / RuntimeError on failure."""
        name = os.path.basename(local_path)
        parent = os.path.dirname(local_path).rstrip("/") + "/"
        cloud_parent = self._cloud_path(parent)
        if name not in self.list_files(parent):
            raise FileNotFoundError(f"cannot find {cloud_parent}{name}")
        os.makedirs(parent, exist_ok=True)
        self._run(["copy"] + (_BIG_FILE_RCLONE_FLAGS if fast else []) + [cloud_parent + name, parent])
        if not os.path.exists(local_path):
            raise FileNotFoundError(f"download failed: {local_path}")
        return local_path

    def sync(self, local_path: str, fast: bool = False, force: bool = False) -> str:
        """Make the local copy match the cloud copy and say what happened: 'downloaded', 'updated' or 'unchanged'.

        `rclone copy` itself skips a file only when size and modification time agree with the cloud copy, so a file replaced in the
        cloud under the same name is fetched again. `force` deletes the local copy first."""
        stat = lambda p: (os.path.getsize(p), int(os.path.getmtime(p))) if os.path.exists(p) else None
        before = stat(local_path)
        if force and before:
            os.remove(local_path)
        self.download(local_path, fast)
        if before is None or force:
            return "downloaded"
        return "unchanged" if stat(local_path) == before else "updated"

    # ------------------------------------------------------------------ genome annotation
    def genome_dir(self, genome_version: str) -> str:
        return self.local_master_dir + "Genomes/" + genome_version + "/"

    def ensure_gff(self, genome_version: str) -> Optional[str]:
        files = GENOME_ANNOTATION.get(genome_version)
        if not files:
            return None
        path = self.genome_dir(genome_version) + files["gff"]
        if not os.path.exists(path):
            try:
                self.download(path)
            except Exception as e:
                self.log(f"[sexloci] annotation not available for {genome_version} ({e}); regions will not list genes")
                return None
        return path

    # ------------------------------------------------------------------ Google sheets (read only)
    def _authenticate(self):
        if self._gc is not None:
            return
        try:
            import gspread
            from google.oauth2.service_account import Credentials
        except ImportError as e:
            raise RuntimeError("Reading the lab sheets needs gspread, gspread-dataframe and google-auth "
                               "(pip install gspread gspread-dataframe google-auth)") from e
        if not os.path.exists(self.credential_file):
            self.download(self.credential_file)
        base = "https://www.googleapis.com/auth/"
        scopes = [base + ("spreadsheets.readonly" if self.read_only else "spreadsheets")]
        if not self.read_only:
            scopes.append(base + "drive")
        self._gc = gspread.authorize(Credentials.from_service_account_file(self.credential_file, scopes=scopes))

    def _read_sheet(self, spreadsheet_id: str, worksheet: str) -> pd.DataFrame:
        from gspread_dataframe import get_as_dataframe
        self._authenticate()
        last = None
        for attempt in range(3):                       # transient API / rate-limit errors
            try:
                ws = self._gc.open_by_key(spreadsheet_id).worksheet(worksheet)
                return get_as_dataframe(ws, evaluate_formulas=True)
            except Exception as e:
                if type(e).__name__ == "WorksheetNotFound":
                    raise
                last = e
                time.sleep(2 * (attempt + 1))
        raise RuntimeError(f"could not read worksheet '{worksheet}': {last}") from last

    def open_sample_worksheet(self):
        """The SampleDatabase worksheet as a gspread object, for writing. Refused unless this object was created with read_only=False."""
        if self.read_only:
            raise PermissionError("this LabData is read-only; writing the Sample Database needs read_only=False")
        self._authenticate()
        return self._gc.open_by_key(self.sample_sheet_id).worksheet("SampleDatabase")

    def read_sample_tables(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        if self._sample_cache is None:
            self._sample_cache = (self._read_sheet(self.sample_sheet_id, "SampleDatabase"),
                                  self._read_sheet(self.sample_sheet_id, "DNAReadsDatabase"))
        return self._sample_cache

    def read_alignment_table(self, genome_version: str) -> pd.DataFrame:
        try:
            return self._read_sheet(self.alignment_sheet_id, genome_version)
        except Exception as e:
            if type(e).__name__ == "WorksheetNotFound":
                raise ValueError(f"No alignment worksheet named '{genome_version}'. Known genome versions: "
                                 + ", ".join(KNOWN_GENOMES)) from e
            raise
