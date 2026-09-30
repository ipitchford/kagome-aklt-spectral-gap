"""Portable output paths. Generated caches are never proof inputs on a fresh replay."""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
EVIDENCE=Path(os.environ.get("AKLT_EVIDENCE",str(ROOT/"evidence"/"new")))
WORK=Path(os.environ.get("AKLT_WORK",str(ROOT/".work")))
EVIDENCE.mkdir(parents=True,exist_ok=True)
WORK.mkdir(parents=True,exist_ok=True)
