"""
Real PDB Fetcher - Queries RCSB PDB and AlphaFold for protein structures.
Used by DockingService to get real target structures for DiffDock.
"""
import logging
import os
import requests
from typing import Optional

logger = logging.getLogger(__name__)

_RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
_RCSB_DOWNLOAD = "https://files.rcsb.org/download/{pdb_id}.pdb"
_ALPHAFOLD_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
_UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search?query={gene}&format=json&size=1&fields=accession"

# Well-known targets for quick lookup (gene symbol → PDB ID)
_KNOWN_PDB = {
    "BACE1": "2QMG", "BACE-1": "2QMG",
    "APP":   "1AAP",
    "MAPT":  "2N3X", "TAU": "2N3X",
    "ACHE":  "4EY7", "ACHE1": "4EY7",
    "PPARG": "2PRG",
    "LRRK2": "7LI3",
    "SNCA":  "1XQ8", "ALPHA-SYNUCLEIN": "1XQ8",
    "MAO-B": "2V5Z", "MAOB": "2V5Z",
    "GSK3B": "1Q5K", "GSK3": "1Q5K",
    "EGFR":  "1M17",   # kinase domain with erlotinib bound (correct ATP pocket)
    "ACE":   "1O86",
    "HMGCR": "1HWK",
    "DPP4":  "1NNY",
    "VEGFR": "2OH4",
    "BCR-ABL": "2HYY",
    "HER2":  "3PP0",
    "CDK5":  "1H4L",
    "HDAC1": "4BKX",
    "SIRT1": "4ZZH",
    "MTOR":  "4JSN",
    "PI3K":  "2RD0",
}


class RealPDBFetcher:
    """Fetches protein structures from RCSB PDB or AlphaFold."""

    def __init__(self):
        self._cache: dict = {}

    def fetch_protein_structure(self, target_name: str) -> str:
        """
        Fetch PDB string for a target protein.
        Strategy:
          1. Fast lookup in _KNOWN_PDB map
          2. RCSB free-text search
          3. AlphaFold via UniProt lookup
        Returns PDB string or "" on failure.
        """
        key = target_name.strip().upper()
        if key in self._cache:
            return self._cache[key]

        # 1. Known PDB map
        pdb_id = _KNOWN_PDB.get(key) or _KNOWN_PDB.get(key.replace(" ", "_"))
        if pdb_id:
            pdb = self._download_pdb(pdb_id)
            if pdb:
                self._cache[key] = pdb
                return pdb

        # 2. RCSB free-text search
        pdb_id = self._rcsb_search(target_name)
        if pdb_id:
            pdb = self._download_pdb(pdb_id)
            if pdb:
                self._cache[key] = pdb
                return pdb

        # 3. AlphaFold via UniProt
        uniprot_id = self._uniprot_lookup(target_name)
        if uniprot_id:
            pdb = self._alphafold_fetch(uniprot_id)
            if pdb:
                self._cache[key] = pdb
                return pdb

        logger.warning(f"No structure found for {target_name}")
        return ""

    def _download_pdb(self, pdb_id: str) -> str:
        try:
            r = requests.get(_RCSB_DOWNLOAD.format(pdb_id=pdb_id), timeout=8)
            if r.status_code == 200:
                logger.info(f"PDB downloaded: {pdb_id}")
                return r.text
        except Exception as e:
            logger.debug(f"PDB download failed {pdb_id}: {e}")
        return ""

    def _rcsb_search(self, name: str) -> Optional[str]:
        try:
            payload = {
                "query": {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {"value": name},
                },
                "return_type": "entry",
                "request_options": {
                    "paginate": {"start": 0, "rows": 1},
                    "results_content_type": ["experimental"],
                    "sort": [{"sort_by": "score", "direction": "desc"}],
                    "scoring_strategy": "combined",
                },
            }
            r = requests.post(_RCSB_SEARCH, json=payload, timeout=6)
            if r.status_code == 200:
                hits = r.json().get("result_set", [])
                if hits:
                    return hits[0]["identifier"]
        except Exception as e:
            logger.debug(f"RCSB search failed for {name}: {e}")
        return None

    def _uniprot_lookup(self, gene: str) -> Optional[str]:
        try:
            r = requests.get(
                _UNIPROT_SEARCH.format(gene=gene),
                timeout=5,
                headers={"Accept": "application/json"},
            )
            if r.status_code == 200:
                results = r.json().get("results", [])
                if results:
                    return results[0]["primaryAccession"]
        except Exception as e:
            logger.debug(f"UniProt lookup failed for {gene}: {e}")
        return None

    def _alphafold_fetch(self, uniprot_id: str) -> str:
        try:
            r = requests.get(_ALPHAFOLD_URL.format(uniprot_id=uniprot_id), timeout=8)
            if r.status_code == 200:
                logger.info(f"AlphaFold structure fetched for {uniprot_id}")
                return r.text
        except Exception as e:
            logger.debug(f"AlphaFold fetch failed {uniprot_id}: {e}")
        return ""
