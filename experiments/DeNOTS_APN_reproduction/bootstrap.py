"""Resolve data and unpack the required author sources into a temporary cache."""
import hashlib
import os
import tempfile
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE / "data"
os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
archive = HERE / "author_sources.zip"
VENDOR = Path(tempfile.gettempdir()) / "denots_apn_sources" / hashlib.sha256(archive.read_bytes()).hexdigest()[:16]
with zipfile.ZipFile(archive) as sources:
    for member in sources.infolist():
        destination = (VENDOR / member.filename).resolve()
        if VENDOR.resolve() not in destination.parents:
            raise ValueError("Invalid author source archive path")
        if not destination.exists():
            sources.extract(member, VENDOR)
