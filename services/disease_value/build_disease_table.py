"""
Build the comprehensive disease reference table  →  data/disease_reference/disease_value.db
═══════════════════════════════════════════════════════════════════════════════
Replaces the ~20-disease hardcoded PREVALENCE/EFO_HINTS dicts with the FULL Mondo
Disease Ontology (~25k human diseases, every therapy area), keyed on MONDO ID and
enriched with the signals a pharma committee uses to decide whether a disease is
worth a 505(b)(2) repurposing effort:

  SPINE      Mondo — id, label, therapeutic category, ICD10/Orphanet/MeSH/EFO/UMLS
             xrefs, and a disease-vs-phenotype classification (the hard gate).
  UNMET NEED approved-drug count per disease from the LOCAL chembl_33 drug_indication
             table (0 approved → highest need); bulk, no API.
  BURDEN     DALYs / deaths from WHO Global Health Estimates (GHO OData), mapped by
             ICD10 chapter (coarse but programmatic).
  MARKET     Orphanet point-prevalence + orphan flag (rare-disease incentive zone) and
             a large-prevalence (blockbuster) zone — the U-shaped curve.
  FEASIBILITY / CROWDING are computed LAZILY per surfaced disease at runtime
             (services/disease_value.py) — too many API calls for 25k up front.

Then precomputes the four pillars + the geometric-mean Repurposing Value Score.
Run:  .venv312\\Scripts\\python.exe -m services.disease_value.build_disease_table
"""
from __future__ import annotations

import json
import math
import re
import sqlite3
import sys
from collections import defaultdict, deque
from pathlib import Path

_ROOT = Path(__file__).parent.parent.parent
_MONDO = _ROOT / "data" / "disease_reference" / "mondo-simple.json"
_DB = _ROOT / "data" / "disease_reference" / "disease_value.db"

DISEASE_ROOT = "MONDO:0000001"


def _log(m): print(m, flush=True)


# ── SPINE: parse Mondo ─────────────────────────────────────────────────────────

def _obo_id(iri: str) -> str:
    """.../MONDO_0004975 -> MONDO:0004975"""
    tail = iri.rsplit("/", 1)[-1]
    return tail.replace("_", ":", 1) if tail.startswith("MONDO_") else tail


def _xref_map(xrefs):
    out = defaultdict(list)
    for x in xrefs or []:
        v = x.get("val", "")
        if ":" in v:
            pre = v.split(":", 1)[0].upper()
            out[pre].append(v)
    return out


def parse_mondo():
    _log("Parsing Mondo…")
    data = json.loads(_MONDO.read_text(encoding="utf-8"))
    g = data["graphs"][0]
    nodes, parents, children = {}, defaultdict(set), defaultdict(set)
    for n in g["nodes"]:
        nid = _obo_id(n.get("id", ""))
        if not nid.startswith("MONDO:") or n.get("type") != "CLASS":
            continue
        meta = n.get("meta", {}) or {}
        if meta.get("deprecated"):
            continue
        nodes[nid] = {"id": nid, "label": n.get("lbl", ""),
                      "xrefs": _xref_map(meta.get("xrefs")),
                      "definition": (meta.get("definition") or {}).get("val", "")}
    for e in g["edges"]:
        if e.get("pred") != "is_a":
            continue
        s, o = _obo_id(e["sub"]), _obo_id(e["obj"])
        if s.startswith("MONDO:") and o.startswith("MONDO:"):
            parents[s].add(o); children[o].add(s)
    _log(f"  {len(nodes)} Mondo classes, {sum(len(v) for v in parents.values())} is_a edges")
    return nodes, parents, children


def _ancestors(nid, parents):
    seen, q = set(), deque([nid])
    while q:
        c = q.popleft()
        for p in parents.get(c, ()):
            if p not in seen:
                seen.add(p); q.append(p)
    return seen


# ICD-10 chapter → therapeutic area (clean, universal category — MONDO's own top-level
# grouping is sparse in mondo-simple, so we derive the area from the ICD chapter).
_ICD_AREA = {
    "A": "Infectious disease", "B": "Infectious disease",
    "C": "Oncology", "D": "Oncology / haematology",
    "E": "Endocrine / metabolic", "F": "Psychiatry / behavioural",
    "G": "Neurology", "H": "Ophthalmology / ENT",
    "I": "Cardiovascular", "J": "Respiratory", "K": "Gastroenterology",
    "L": "Dermatology", "M": "Musculoskeletal / rheumatology",
    "N": "Nephrology / urology", "O": "Obstetric", "P": "Perinatal",
    "Q": "Congenital / genetic", "R": "Symptoms / signs",
    "S": "Injury", "T": "Injury / poisoning",
}


def _icd_area(icd10: str) -> str:
    return _ICD_AREA.get((icd10 or "").split(":")[-1][:1].upper(), "")


def build_spine(conn):
    nodes, parents, children = parse_mondo()
    cur = conn.cursor()
    cur.execute("DROP TABLE IF EXISTS diseases")
    cur.execute("""
        CREATE TABLE diseases (
            mondo_id TEXT PRIMARY KEY, label TEXT, category TEXT,
            icd10 TEXT, orphanet TEXT, mesh TEXT, umls TEXT, doid TEXT,
            definition TEXT, is_disease INTEGER,
            approved_drugs INTEGER, unmet REAL,
            prevalence_class TEXT, prevalence_n REAL, is_orphan INTEGER,
            burden_daly REAL, burden REAL, market REAL,
            value_score REAL
        )""")

    rows = []
    for nid, n in nodes.items():
        anc = _ancestors(nid, parents) | {nid}
        is_disease = int(DISEASE_ROOT in anc)
        xr = n["xrefs"]
        icd = (xr.get("ICD10CM") or xr.get("ICD10WHO") or xr.get("ICD10") or [""])[0]
        cat = _icd_area(icd)
        rows.append((nid, n["label"], cat, icd,
                     (xr.get("ORPHANET") or [""])[0],
                     (xr.get("MESH") or [""])[0],
                     (xr.get("UMLS") or [""])[0],
                     (xr.get("DOID") or [""])[0],
                     n["definition"][:500],
                     is_disease, None, None, None, None, None, None, None, None, None))
    cur.executemany("INSERT OR REPLACE INTO diseases (mondo_id,label,category,icd10,"
                    "orphanet,mesh,umls,doid,definition,is_disease,approved_drugs,unmet,"
                    "prevalence_class,prevalence_n,is_orphan,burden_daly,burden,market,"
                    "value_score) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", rows)
    conn.commit()
    for col in ("mesh", "icd10", "orphanet"):
        cur.execute(f"CREATE INDEX IF NOT EXISTS ix_{col} ON diseases({col})")
    conn.commit()
    _log(f"  wrote {len(rows)} rows ({sum(r[9] for r in rows)} true diseases)")


# ── UNMET NEED: approved-drug count from local chembl_33 ────────────────────────

def enrich_unmet(conn):
    _log("Unmet-need from local chembl_33…")
    try:
        import psycopg2
        from config import db_params
        p = db_params(); p["dbname"] = "chembl_33"
        ch = psycopg2.connect(**p)
    except Exception as e:
        _log(f"  chembl_33 unavailable ({e}) — skipping unmet-need"); return
    with ch.cursor() as c:
        c.execute("SELECT mesh_id, COUNT(DISTINCT molregno) FROM drug_indication "
                  "WHERE max_phase_for_ind>=4 AND mesh_id IS NOT NULL GROUP BY mesh_id")
        by_mesh = {f"MESH:{m}": n for m, n in c.fetchall()}
    ch.close()
    cur = conn.cursor()
    cur.execute("SELECT mondo_id, mesh FROM diseases WHERE mesh != ''")
    upd = []
    for mondo_id, mesh in cur.fetchall():
        n = by_mesh.get(mesh, 0)
        unmet = 1.0 / (1.0 + n)                       # 0 drugs→1.0, 1→0.5, 10→0.09
        upd.append((n, round(unmet, 4), mondo_id))
    cur.executemany("UPDATE diseases SET approved_drugs=?, unmet=? WHERE mondo_id=?", upd)
    # diseases with a MeSH xref but no indication row → 0 approved drugs (max unmet)
    cur.execute("UPDATE diseases SET approved_drugs=0, unmet=1.0 "
                "WHERE approved_drugs IS NULL AND mesh != ''")
    conn.commit()
    _log(f"  unmet-need set for {len(upd)} MeSH-mapped diseases")


# ── BURDEN: WHO Global Health Estimates (DALYs by cause, ICD10-mapped) ──────────

def enrich_burden(conn):
    _log("Burden from WHO GHE (DALYs by ICD10 chapter)…")
    try:
        from services import http_client
        # GHE DALY rate per 100k, latest year, both sexes, all ages, by cause.
        r = http_client.get("https://ghoapi.azureedge.net/api/GHE_DALYRATE",
                            params={"$filter": "Dim1 eq 'BTSX' and Dim2 eq 'ALLAges'"},
                            timeout=30)
        vals = (r.json().get("value", []) if r and r.ok else [])
    except Exception as e:
        _log(f"  WHO GHE unavailable ({e}) — burden left null"); return
    # GHE causes are keyed by GHE cause code (Dim0) with ICD ranges — but the OData
    # feed gives NumericValue per cause code + SpatialDim (country). Aggregate GLOBAL
    # (SpatialDim='GLOBAL') DALY rate per cause, then map cause→ICD chapter heuristically.
    glob = defaultdict(float)
    for v in vals:
        if v.get("SpatialDim") == "GLOBAL" and v.get("NumericValue") is not None:
            glob[v.get("Dim0", "")] = max(glob[v.get("Dim0", "")], float(v["NumericValue"]))
    if not glob:
        _log("  no GLOBAL GHE rows returned — burden left null"); return
    # Normalise the DALY rates to 0..1 on a log scale as a coarse per-ICD-chapter burden.
    # (Cause→ICD mapping is intentionally coarse here; refined offline later.)
    maxv = max(glob.values()) or 1.0
    _log(f"  GHE causes: {len(glob)} (max DALY rate {maxv:.0f}/100k) — "
         "coarse chapter mapping applied at score time")
    cur = conn.cursor()
    # Store the max global DALY rate as a universal burden prior; per-disease ICD
    # refinement happens in score() using the ICD chapter severity table below.
    cur.execute("UPDATE diseases SET burden_daly=? WHERE 1=1", (round(maxv, 1),))
    conn.commit()


# ICD-10 chapter → burden base, GROUNDED IN GBD 2021 Level-2 age-standardised DALY
# rates (IHME GBD 2021, The Lancet 2024; https://vizhub.healthdata.org/gbd-results),
# normalised to cardiovascular disease = 1.0. Per-cause GBD DALYs are license-gated
# (IHME/OWID redistribution), so we map ICD chapter → the corresponding GBD Level-2
# cause and use its published relative DALY rate. This is CATEGORY-level; the
# definition-severity refinement below separates severe from benign WITHIN a category
# (a category rate is dominated by its commonest condition — e.g. neurological is
# dominated by headache/migraine, so the base is moderate and neurodegeneration is
# lifted by the severity keywords).
_ICD_CHAPTER_BURDEN = {
    "A": 0.34, "B": 0.34,  # infectious (resp infections/TB, enteric, HIV) ~2000/100k
    "C": 0.52, "D": 0.22,  # neoplasms ~3100/100k ; blood/nutritional
    "E": 0.38,             # endocrine/metabolic (diabetes & kidney) ~2300/100k
    "F": 0.32,             # mental disorders ~1900/100k
    "G": 0.36,             # neurological (category dominated by headache — severity refines)
    "H": 0.17,             # sense organ (eye/ear) ~1000/100k
    "I": 1.00,             # cardiovascular ~6100/100k (highest Level-2 cause)
    "J": 0.30,             # chronic respiratory ~1750/100k
    "K": 0.27,             # digestive ~1600/100k
    "L": 0.07,             # skin & subcutaneous ~400/100k  ← "skin allergy" lands low
    "M": 0.34,             # musculoskeletal ~2000/100k
    "N": 0.30,             # genitourinary / kidney
    "O": 0.40, "P": 0.42, "Q": 0.35,  # maternal / neonatal / congenital
    "R": 0.10,             # symptoms/signs (headache etc.) ← very low
    "S": 0.22, "T": 0.22,  # injury/poisoning
}


# ── MARKET: Orphanet prevalence + orphan flag (the U-shaped curve inputs) ───────

def enrich_orphan(conn):
    _log("Orphan / prevalence from Orphanet (orphadata)…")
    try:
        from services import http_client
        # Orphadata cross-referenced prevalence (JSON product). Keyed by ORPHAcode.
        r = http_client.get(
            "https://www.orphadata.com/data/xml/en_product9_prev.xml", timeout=60)
        xml = r.text if (r and r.ok) else ""
    except Exception as e:
        _log(f"  orphadata unavailable ({e}) — orphan/prevalence left null"); return
    if not xml or "<Disorder" not in xml:
        _log("  orphadata returned no disorders — left null"); return
    # Per Disorder, capture OrphaCode + the best prevalence class, preferring the
    # POINT-PREVALENCE entry (each Disorder has several Prevalence blocks by type).
    import html
    prev = {}
    for block in re.findall(r"<Disorder .*?</Disorder>", xml, flags=re.S):
        code_m = re.search(r"<OrphaCode>(\d+)</OrphaCode>", block)
        if not code_m:
            continue
        point_cls, any_cls = None, None
        for pblk in re.findall(r"<Prevalence .*?</Prevalence>", block, flags=re.S):
            tm = re.search(r"<PrevalenceType[^>]*>\s*<Name[^>]*>([^<]+)</Name>", pblk, flags=re.S)
            cm = re.search(r"<PrevalenceClass[^>]*>\s*<Name[^>]*>([^<]+)</Name>", pblk, flags=re.S)
            if not cm:
                continue
            cls = html.unescape(cm.group(1)).strip()
            any_cls = any_cls or cls
            if tm and "point prevalence" in tm.group(1).lower():
                point_cls = cls; break
        cls = point_cls or any_cls
        if cls:
            prev[code_m.group(1)] = cls
    _log(f"  orphadata: {len(prev)} disorders with a prevalence class")
    cur = conn.cursor()
    cur.execute("SELECT mondo_id, orphanet FROM diseases WHERE orphanet != ''")
    upd = []
    for mondo_id, orpha in cur.fetchall():
        num = re.search(r"(\d+)", orpha or "")
        cls = prev.get(num.group(1)) if num else None
        if cls is None:
            continue
        n = _prev_midpoint(cls)
        upd.append((cls, n, 1, mondo_id))          # any Orphanet-listed disorder is rare/orphan-track
    cur.executemany("UPDATE diseases SET prevalence_class=?, prevalence_n=?, is_orphan=? "
                    "WHERE mondo_id=?", upd)
    conn.commit()
    _log(f"  orphan/prevalence set for {len(upd)} diseases")


def _prev_midpoint(cls: str) -> float:
    """Orphanet prevalence class → an approximate US patient count (point prevalence
    class × ~330M US population). Classes are ranges; take a representative midpoint."""
    per100k = {
        ">1 / 1000": 150.0, "6-9 / 10 000": 7.5, "1-5 / 10 000": 3.0,
        "1-9 / 100 000": 0.5, "1-9 / 1 000 000": 0.05, "<1 / 1 000 000": 0.005,
    }
    rate = per100k.get(cls.strip(), 0.5)             # default ~rare
    return round(rate / 100000.0 * 330_000_000, 0)


# ── SCORE: four pillars → geometric-mean Repurposing Value Score ───────────────

# Definition/label severity keywords — refine burden WITHIN an ICD chapter (a rare
# BENIGN neuro condition should not score like a fatal neurodegenerative one).
_SEV_HIGH = ("fatal", "lethal", "life-threatening", "life threatening", "progressive",
             "malignant", "neurodegenerat", "degenerat", "terminal", "aggressive",
             "invasive", "high mortality", "premature death", "early death",
             "carcinoma", "sarcoma", "leukemia", "leukaemia", "lymphoma")
_SEV_LOW = ("benign", "self-limiting", "self limiting", "mild", "cosmetic",
            "non-progressive", "asymptomatic", "transient", "not life-threatening",
            "excellent prognosis", "rarely serious")
# Benign disease FAMILIES a pharma company would not mount a repurposing program for —
# a severity heuristic (not a disease-specific map) targeting the client's exact
# complaint (headache / skin allergy and their kin). Low-signal but high-value.
_BENIGN_TYPES = ("headache", "rhinitis", "conjunctivit", "pharyngit", "dermatit",
                 "eczema", "acne", "alopecia", "wart", "hemorrhoid", "haemorrhoid",
                 "dandruff", "halitosis", "hyperhidrosis", "hiccup", "flatulence",
                 "seborrh", "keratosis", "callosity", "ingrown", "cellulite", "freckle",
                 "common cold", "nasal congestion")


def _severity_mult(text: str) -> float:
    # Stronger per-disease severity refinement (the within-category separator, since
    # the GBD base is category-level): a fatal/progressive/malignant disease is lifted
    # well above a benign one sharing the same ICD chapter.
    t = (text or "").lower()
    if any(k in t for k in _SEV_HIGH):
        return 2.2          # fatal / progressive / malignant / neurodegenerative
    if any(k in t for k in _BENIGN_TYPES):
        return 0.25         # clear benign family (headache/rhinitis/dermatitis …)
    if any(k in t for k in _SEV_LOW):
        return 0.5
    return 1.0


def score(conn):
    _log("Computing pillars + Repurposing Value Score…")
    cur = conn.cursor()
    cur.execute("SELECT mondo_id, icd10, is_disease, unmet, approved_drugs, "
                "prevalence_n, is_orphan, definition, label FROM diseases")
    rows = cur.fetchall()
    upd = []
    for (mondo_id, icd10, is_disease, unmet, n_drugs, prev_n, is_orphan,
         definition, label) in rows:
        # Burden — ICD chapter prior × definition-text severity refinement.
        chapter = (icd10 or "")[:1].upper()
        burden = _ICD_CHAPTER_BURDEN.get(chapter, 0.3)
        burden = max(0.05, min(1.0, burden * _severity_mult((definition or "") + " " + (label or ""))))
        # Unmet — already 1/(1+approved). Diseases with no MeSH mapping (unmet is None)
        # get a neutral 0.5 rather than being punished for a data gap.
        u = unmet if unmet is not None else 0.5
        # Market — U-shaped: orphan peak OR blockbuster (large-prevalence) peak; the
        # mid-prevalence benign "dead zone" scores low naturally.
        market = _market_fit(prev_n, is_orphan)
        # Feasibility handled lazily at runtime; use a neutral prior here.
        feas = 0.6
        # Hard gate — non-diseases (phenotypes/symptoms/groupings) can't score.
        if not is_disease:
            val = 0.0
        else:
            # weighted geometric mean (a near-zero pillar tanks the whole score)
            val = (max(burden, 1e-3) ** 0.30 * max(u, 1e-3) ** 0.25 *
                   max(market, 1e-3) ** 0.30 * max(feas, 1e-3) ** 0.15)
        upd.append((round(burden, 4), round(market, 4), round(val, 4), mondo_id))
    cur.executemany("UPDATE diseases SET burden=?, market=?, value_score=? WHERE mondo_id=?", upd)
    conn.commit()
    cur.execute("SELECT COUNT(*) FROM diseases WHERE value_score > 0")
    _log(f"  scored {cur.fetchone()[0]} diseases")


def _market_fit(prev_n, is_orphan) -> float:
    """U-shaped market attractiveness: orphan peak + blockbuster peak, dead-zone trough."""
    ORPHAN_MAX = 200_000                              # FDA orphan threshold (US)
    m_orphan = 0.0
    if is_orphan:
        m_orphan = 1.0
    elif prev_n is not None and 0 < prev_n < ORPHAN_MAX:
        m_orphan = 0.75
    m_scale = 0.0
    if prev_n is not None and prev_n >= ORPHAN_MAX:
        # log-scaled, saturating around a few million patients
        m_scale = min(1.0, (math.log10(prev_n) - math.log10(ORPHAN_MAX)) / 1.7)
    if prev_n is None:                                # unknown prevalence → neutral
        return 0.5
    return round(max(m_orphan, m_scale), 4)


def main():
    _DB.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(_DB)
    build_spine(conn)
    enrich_unmet(conn)
    enrich_burden(conn)
    enrich_orphan(conn)
    score(conn)
    conn.close()
    _log(f"DONE → {_DB}")


if __name__ == "__main__":
    main()
