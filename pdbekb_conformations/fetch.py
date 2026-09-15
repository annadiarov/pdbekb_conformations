"""
fetch.py -- networked access to PDBe-KB superposition data.

Two sources, exactly as described by the PDBe team:

1. The "superposition" Graph-API endpoint. This is the ONLY thing you need to
   query to discover the segments for an accession -- you never have to guess
   the `segment_<start>_<end>` path fragment yourself (that was the original
   pain point). The JSON lists, per segment, the clusters and their members.

       https://www.ebi.ac.uk/pdbe/graph-api/uniprot/superposition/{accession}

2. The FTP area, which holds the precomputed matrices per segment. Paths are
   fully determined once the API has given you (start, end):

       .../superposition/{L}/{ACC}/segment_{start}_{end}/{ACC}_{start}_{end}_score_data.npz
       .../superposition/{L}/{ACC}/segment_{start}_{end}/{ACC}_{start}_{end}_linkage_data.npz
       .../superposition/{L}/{ACC}/{ACC}.json          (viewer transformation matrices)

   where L = ACC[0].

Everything is cached on disk so re-runs are cheap and offline-friendly.
"""

from __future__ import annotations

import io
import json
import logging
import time
from pathlib import Path
from typing import Optional

import numpy as np
import requests

log = logging.getLogger(__name__)

GRAPH_API = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/superposition/{acc}"
FTP_ROOT = "https://ftp.ebi.ac.uk/pub/databases/pdbe-kb/superposition"

# Structure files for the RMSD step come from the main PDBe / RCSB servers.
MMCIF_URL = "https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb}_updated.cif"
SIFTS_MAP_API = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/{pdb}"

DEFAULT_HEADERS = {"User-Agent": "bindome-conformations/0.1 (research use)"}


class Fetcher:
    """Thin caching HTTP client. One instance per run."""

    def __init__(
        self,
        cache_dir: str | Path = "cache",
        pause_s: float = 0.2,
        max_retries: int = 4,
        timeout: int = 60,
    ):
        self.cache = Path(cache_dir)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.pause_s = pause_s
        self.max_retries = max_retries
        self.timeout = timeout
        self._session = requests.Session()
        self._session.headers.update(DEFAULT_HEADERS)

    # -- low level -----------------------------------------------------------
    def _get(self, url: str) -> Optional[bytes]:
        """GET with retries. Returns None on a clean 404 (missing data)."""
        for attempt in range(1, self.max_retries + 1):
            try:
                r = self._session.get(url, timeout=self.timeout)
                if r.status_code == 404:
                    log.debug("404 (absent) %s", url)
                    return None
                r.raise_for_status()
                time.sleep(self.pause_s)  # be polite to the EBI servers
                return r.content
            except requests.RequestException as exc:
                wait = min(2 ** attempt, 30)
                log.warning("GET failed (%s/%s) %s -- %s; retry in %ss",
                            attempt, self.max_retries, url, exc, wait)
                time.sleep(wait)
        log.error("giving up on %s", url)
        return None

    def _cached(self, key: str, url: str) -> Optional[bytes]:
        f = self.cache / key
        if f.exists():
            return f.read_bytes()
        data = self._get(url)
        if data is not None:
            f.parent.mkdir(parents=True, exist_ok=True)
            f.write_bytes(data)
        return data

    # -- superposition API ---------------------------------------------------
    def superposition(self, acc: str) -> Optional[dict]:
        """Return the parsed superposition JSON for a UniProt accession.

        Shape (as returned by the Graph-API):
            { "<acc>": [ {segment}, {segment}, ... ] }
        where each segment is
            { "segment_start": int, "segment_end": int,
              "clusters": [ [ {member}, ... ], ... ] }
        and each member carries pdb_id / auth_asym_id / struct_asym_id /
        is_representative. Key names are handled defensively in parse.py.
        """
        raw = self._cached(f"api/{acc}.json", GRAPH_API.format(acc=acc))
        if raw is None:
            return None
        return json.loads(raw)

    # -- FTP matrices --------------------------------------------------------
    def _seg_dir_url(self, acc: str, start: int, end: int) -> str:
        return f"{FTP_ROOT}/{acc[0]}/{acc}/segment_{start}_{end}"

    def score_npz(self, acc: str, start: int, end: int) -> Optional[dict]:
        url = f"{self._seg_dir_url(acc, start, end)}/{acc}_{start}_{end}_score_data.npz"
        return self._load_npz(f"npz/{acc}_{start}_{end}_score.npz", url)

    def linkage_npz(self, acc: str, start: int, end: int) -> Optional[dict]:
        url = f"{self._seg_dir_url(acc, start, end)}/{acc}_{start}_{end}_linkage_data.npz"
        return self._load_npz(f"npz/{acc}_{start}_{end}_linkage.npz", url)

    def transforms(self, acc: str) -> Optional[dict]:
        """Per-accession viewer transformation matrices (GESAMT output)."""
        raw = self._cached(f"transforms/{acc}.json", f"{FTP_ROOT}/{acc[0]}/{acc}/{acc}.json")
        return json.loads(raw) if raw is not None else None

    def _load_npz(self, key: str, url: str) -> Optional[dict]:
        raw = self._cached(key, url)
        if raw is None:
            return None
        with np.load(io.BytesIO(raw), allow_pickle=True) as z:
            return {k: z[k] for k in z.files}

    # -- structures (for the optional RMSD step) -----------------------------
    def mmcif(self, pdb: str) -> Optional[bytes]:
        pdb = pdb.lower()
        return self._cached(f"mmcif/{pdb}.cif", MMCIF_URL.format(pdb=pdb))

    def sifts_uniprot_segments(self, pdb: str) -> Optional[dict]:
        raw = self._cached(f"sifts/{pdb}.json", SIFTS_MAP_API.format(pdb=pdb.lower()))
        return json.loads(raw) if raw is not None else None
