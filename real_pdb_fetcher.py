"""
Real PDB Fetcher — resolve a gene/target to a real 3-D protein structure.

Selection strategy (most trustworthy first):
  1. Curated gene→PDB map (hand-verified well-known targets).
  2. UniProt IDENTITY path — resolve gene→UniProt accession (gene-exact, human,
     reviewed), then:
       a. query RCSB for experimental structures that reference THAT accession,
          ranked by resolution → guarantees the structure actually contains the
          target protein (no free-text mismatch);
       b. else AlphaFold predicted model for that accession, pLDDT-filtered to
          drop disordered low-confidence regions (which otherwise create phantom
          pockets).
  3. RCSB free-text search (last resort only).

Modern/large RCSB entries are often mmCIF-only (no legacy .pdb); those are fetched
as .cif and converted to PDB via Biopython. Results persist to a disk cache so we
do not re-hit the network every restart.
"""
import logging
from typing import Optional

from services import http_client
from services.diskcache import DiskCache

logger = logging.getLogger(__name__)

_RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
_RCSB_DOWNLOAD = "https://files.rcsb.org/download/{pdb_id}.pdb"
_RCSB_CIF = "https://files.rcsb.org/download/{pdb_id}.cif"
_ALPHAFOLD_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
_UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"

# AlphaFold per-residue confidence (pLDDT, stored in the B-factor column). Atoms
# below this are disordered/unreliable and are stripped before docking.
_PLDDT_MIN = 50.0

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
    """Fetches real protein structures, selecting by UniProt identity for accuracy."""

    def __init__(self):
        self._cache: dict = {}                       # per-process memo
        # v2: holo-preferring, fragment-rejecting selection (was picking isolated
        # domains like the ABL SH3 instead of the drug-binding kinase domain).
        self._disk = DiskCache("protein_pdb_v2", ttl=30 * 86400)  # 30-day disk cache
        # v2: schema now includes `ligands` (nonpolymer count) for holo preference.
        self._meta_cache = DiskCache("rcsb_meta_v2", ttl=30 * 86400)  # entry metadata

    def fetch_protein_structure(self, target_name: str) -> str:
        """Return a real PDB string for a target gene/protein, or "" if none found."""
        key = (target_name or "").strip().upper()
        if not key:
            return ""
        if key in self._cache:
            return self._cache[key]
        disk = self._disk.get(key)
        if disk is not None:
            self._cache[key] = disk
            return disk

        pdb = self._resolve(target_name, key)
        if pdb and len(pdb) > 200:
            self._cache[key] = pdb
            self._disk.set(key, pdb)
            return pdb
        logger.warning(f"No real structure found for {target_name}")
        return ""

    def _resolve(self, target_name: str, key: str) -> str:
        # 1. Curated map (trusted)
        pdb_id = _KNOWN_PDB.get(key) or _KNOWN_PDB.get(key.replace(" ", "_"))
        if pdb_id:
            pdb = self._download_pdb(pdb_id)
            if pdb:
                logger.info(f"{target_name}: curated PDB {pdb_id}")
                return pdb

        # 2. UniProt identity path
        uniprot_id = self._uniprot_lookup(target_name)
        if uniprot_id:
            pdb_id = self._rcsb_by_uniprot(uniprot_id)
            if pdb_id:
                pdb = self._download_pdb(pdb_id)
                if pdb:
                    logger.info(f"{target_name}: RCSB {pdb_id} via UniProt {uniprot_id}")
                    return pdb
            af = self._alphafold_fetch(uniprot_id)
            if af:
                logger.info(f"{target_name}: AlphaFold model for {uniprot_id} (pLDDT-filtered)")
                return af

        # 3. Free-text search (last resort)
        pdb_id = self._rcsb_search(target_name)
        if pdb_id:
            pdb = self._download_pdb(pdb_id)
            if pdb:
                logger.info(f"{target_name}: RCSB {pdb_id} via free-text (unverified)")
                return pdb
        return ""

    # ── Downloads ───────────────────────────────────────────────────────────────
    def _download_pdb(self, pdb_id: str) -> str:
        # 1. Legacy .pdb format (fast, directly usable)
        text = http_client.get_text(_RCSB_DOWNLOAD.format(pdb_id=pdb_id), timeout=10)
        if text and text.lstrip().startswith(("HEADER", "ATOM", "REMARK", "CRYST", "TITLE")):
            return text
        # 2. Many modern / large entries are mmCIF-only — fetch CIF and convert
        cif = http_client.get_text(_RCSB_CIF.format(pdb_id=pdb_id), timeout=15)
        if cif:
            pdb = self._cif_to_pdb(cif, pdb_id)
            if pdb:
                logger.info(f"PDB via mmCIF conversion: {pdb_id}")
                return pdb
        return ""

    @staticmethod
    def _cif_to_pdb(cif_text: str, pdb_id: str = "X") -> str:
        """Convert an mmCIF string to legacy PDB text (first model only) via Biopython."""
        try:
            import io
            from Bio.PDB import MMCIFParser, PDBIO, Select

            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure(pdb_id, io.StringIO(cif_text))
            first_model_id = list(structure.get_models())[0].id

            class _FirstModel(Select):
                def accept_model(self, model):
                    return model.id == first_model_id

            out = io.StringIO()
            pio = PDBIO()
            pio.set_structure(structure)
            pio.save(out, select=_FirstModel())
            text = out.getvalue()
            return text if text and text.count("ATOM ") > 20 else ""
        except Exception as e:
            logger.debug(f"CIF→PDB conversion error for {pdb_id}: {e}")
            return ""

    def _alphafold_fetch(self, uniprot_id: str) -> str:
        pdb = http_client.get_text(_ALPHAFOLD_URL.format(uniprot_id=uniprot_id), timeout=12)
        if not pdb:
            return ""
        return self._filter_plddt(pdb)

    @staticmethod
    def _filter_plddt(pdb_text: str, cutoff: float = _PLDDT_MIN) -> str:
        """Drop AlphaFold atoms whose pLDDT (B-factor column) is below cutoff.

        Removes disordered tails/loops that would otherwise be detected as fake
        binding pockets. Falls back to the unfiltered model if filtering would
        discard nearly everything.
        """
        kept, total = [], 0
        for ln in pdb_text.splitlines():
            if ln.startswith(("ATOM", "HETATM")):
                total += 1
                try:
                    if float(ln[60:66]) < cutoff:
                        continue
                except ValueError:
                    pass
                kept.append(ln)
            elif ln.startswith(("HEADER", "TITLE", "CRYST", "TER", "END")):
                kept.append(ln)
        atom_lines = sum(1 for ln in kept if ln.startswith("ATOM"))
        if total and atom_lines < max(20, 0.2 * total):
            return pdb_text  # too aggressive — keep the full model
        return "\n".join(kept) + "\n"

    # ── Identity resolution ───────────────────────────────────────────────────────
    def _uniprot_lookup(self, gene: str) -> Optional[str]:
        """Resolve a gene symbol to a human, reviewed UniProt accession (exact first)."""
        queries = [
            f"(gene_exact:{gene}) AND (organism_id:9606) AND (reviewed:true)",
            f"(gene:{gene}) AND (organism_id:9606)",
            gene,
        ]
        for q in queries:
            data = http_client.get_json(
                _UNIPROT_SEARCH,
                params={"query": q, "format": "json", "size": 1, "fields": "accession"},
                headers={"Accept": "application/json"}, timeout=8, default=None)
            if data:
                results = data.get("results", [])
                if results:
                    return results[0]["primaryAccession"]
        return None

    def _rcsb_by_uniprot(self, uniprot_id: str) -> Optional[str]:
        """Best experimental RCSB structure that references this UniProt accession.

        Guarantees the structure contains the target protein (identity), then ranks
        candidates to prefer a COMPACT, well-resolved structure (few chains, modest
        atom count) over a giant multi-chain assembly — a large complex often hides
        the orthosteric pocket from cavity detection.
        """
        attr = "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers"
        payload = {
            "query": {
                "type": "group", "logical_operator": "and",
                "nodes": [
                    {"type": "terminal", "service": "text", "parameters": {
                        "attribute": f"{attr}.database_accession",
                        "operator": "exact_match", "value": uniprot_id}},
                    {"type": "terminal", "service": "text", "parameters": {
                        "attribute": f"{attr}.database_name",
                        "operator": "exact_match", "value": "UniProt"}},
                ],
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": 0, "rows": 15},
                "results_content_type": ["experimental"],
                "sort": [{"sort_by": "rcsb_entry_info.resolution_combined",
                          "direction": "asc"}],
            },
        }
        data = http_client.post_json(_RCSB_SEARCH, json=payload, timeout=10, default=None)
        if not data:
            return None
        ids = [h["identifier"] for h in data.get("result_set", [])]
        if not ids:
            return None

        # Pull lightweight metadata for the top candidates and score them. Pull a
        # wider pool than before: sorting by resolution floats tiny high-res domain
        # FRAGMENTS to the top (e.g. an isolated SH3/SH2 domain), so a holo catalytic-
        # domain structure can sit lower in the list — we need it in the scored pool.
        cands = []
        for pid in ids[:20]:
            meta = self._rcsb_entry_meta(pid)
            if meta:
                cands.append((pid, meta))
        if not cands:
            return ids[0]  # metadata unavailable — fall back to best-resolution hit

        def _score(meta: dict) -> float:
            # Lower is better. Priorities, in order:
            #   1. HOLO — a bound non-polymer ligand marks the REAL druggable pocket
            #      (dock box centred on it). This is what a small molecule actually
            #      binds, so it dominates.
            #   2. Not a FRAGMENT — reject tiny isolated domains (SH3/SH2/etc., a few
            #      hundred atoms per chain) that lack the catalytic/orthosteric pocket.
            #      This is the bug that put imatinib into the ABL SH3 domain (3EG3)
            #      instead of the kinase domain (1IEP).
            #   3. Sane assembly size — gently avoid giant multi-chain assemblies.
            #   4. Resolution — a mild tie-breaker only.
            chains = meta.get("chains") or 1
            atoms  = meta.get("atoms") or 0
            res    = meta.get("resolution")
            res    = res if res is not None else 3.5
            has_ligand = (meta.get("ligands") or 0) > 0
            holo_bonus = -5.0 if has_ligand else 0.0
            per_chain  = atoms / max(1, chains)
            frag_penalty = 5.0 if per_chain < 1200 else 0.0     # <~150 residues = fragment
            size_penalty = atoms / 20000.0                       # gentle large-assembly penalty
            return chains * 1.5 + size_penalty + res * 0.5 + holo_bonus + frag_penalty

        best = min(cands, key=lambda c: _score(c[1]))
        return best[0]

    def _rcsb_entry_meta(self, pdb_id: str) -> Optional[dict]:
        """Resolution / chain count / atom count for an RCSB entry (cached)."""
        cached = self._meta_cache.get(pdb_id)
        if cached is not None:
            return cached or None
        data = http_client.get_json(
            f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}", timeout=8, default=None)
        meta = {}
        if data:
            info = data.get("rcsb_entry_info", {})
            res_list = info.get("resolution_combined") or []
            meta = {
                "resolution": (res_list[0] if res_list else None),
                "chains": info.get("deposited_polymer_entity_instance_count"),
                "atoms": info.get("deposited_atom_count"),
                # Bound non-polymer entities (ligands/cofactors; excludes water). Used
                # to prefer HOLO structures — see _rcsb_by_uniprot._score.
                "ligands": info.get("nonpolymer_entity_count"),
            }
        self._meta_cache.set(pdb_id, meta)
        return meta or None

    def _rcsb_search(self, name: str) -> Optional[str]:
        payload = {
            "query": {"type": "terminal", "service": "full_text",
                      "parameters": {"value": name}},
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": 0, "rows": 1},
                "results_content_type": ["experimental"],
                "sort": [{"sort_by": "score", "direction": "desc"}],
                "scoring_strategy": "combined",
            },
        }
        data = http_client.post_json(_RCSB_SEARCH, json=payload, timeout=8, default=None)
        if data:
            hits = data.get("result_set", [])
            if hits:
                return hits[0]["identifier"]
        return None
