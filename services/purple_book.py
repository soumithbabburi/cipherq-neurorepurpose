"""
FDA Purple Book — authoritative status of licensed biologics & biosimilars.

The Purple Book is the biologics counterpart to the Orange Book: it lists every
FDA-licensed biological product, whether it is a reference product (BLA 351(a)) or
a biosimilar/interchangeable (BLA 351(k)), the reference product each biosimilar
copies, and the reference-product exclusivity dates. It is the authoritative source
for "which biologics are open (or opening) to biosimilar competition, and how
crowded is that space already".

This module powers the Biosimilar Opportunity Radar: it groups reference biologics,
computes their exclusivity cliff, counts existing biosimilars/interchangeables, and
ranks each reference product by how attractive it is as a biosimilar target
(cliff timing x whitespace x franchise size).

Data: the FDA "Purple Book Monthly Data" CSV, downloaded + cached locally and
refreshed monthly. See https://purplebooksearch.fda.gov/downloads
"""
import calendar
import csv
import json
import logging
import re
import time
from datetime import date, datetime
from pathlib import Path
from typing import Dict, List, Optional

import requests

logger = logging.getLogger(__name__)

PB_PAGE = "https://purplebooksearch.fda.gov/downloads"
# accessdata hosts the actual files, e.g.
#   /drugsatfda_docs/PurpleBook/2026/purplebook-search-June-data-download.csv
PB_HOST = "https://www.accessdata.fda.gov/drugsatfda_docs/PurpleBook"
PB_DIR = Path(__file__).parent.parent / "data" / "purple_book"
_MAX_AGE = 30 * 86400          # refresh monthly
_UA = {"User-Agent": "Mozilla/5.0"}

# BPCIA grants a reference biologic 12 years of exclusivity from first licensure.
_REF_EXCLUSIVITY_YEARS = 12

# Repurposing Value Score per biologic (built offline by build_value_cache()), so the
# radar can rank commercial attractiveness without hitting the ChEMBL DB per request.
_VALUE_CACHE = PB_DIR / "value_cache.json"

_index: Optional[Dict] = None
_value_map: Optional[Dict[str, Dict]] = None


# ────────────────────────────────────────────────────────────────────────────
# download / cache
# ────────────────────────────────────────────────────────────────────────────
def _candidate_urls() -> List[str]:
    """Newest-first list of monthly CSV URLs to try (this month back 18 months)."""
    urls, y, m = [], date.today().year, date.today().month
    for _ in range(18):
        month_name = calendar.month_name[m]              # e.g. "June"
        for mn in (month_name, month_name.lower()):      # FDA is inconsistent on case
            urls.append(f"{PB_HOST}/{y}/purplebook-search-{mn}-data-download.csv")
        m -= 1
        if m == 0:
            m, y = 12, y - 1
    return urls


def _scrape_urls() -> List[str]:
    """Fallback: scrape the downloads page for accessdata CSV links (newest first)."""
    try:
        r = requests.get(PB_PAGE, timeout=25, headers=_UA)
        if r.ok:
            links = re.findall(r'href=["\']([^"\']*PurpleBook/[^"\']*data-download\.csv)["\']',
                               r.text, re.I)
            return [(u if u.startswith("http") else "https://www.accessdata.fda.gov" + u)
                    for u in links]
    except Exception as e:
        logger.debug(f"Purple Book page scrape failed: {e}")
    return []


def _download() -> bool:
    PB_DIR.mkdir(parents=True, exist_ok=True)
    for url in _candidate_urls() + _scrape_urls():
        try:
            r = requests.get(url, timeout=90, headers=_UA)
            # CSV starts with a BOM + the "Purple Book ... Report" preamble line
            if r.ok and b"," in r.content[:200] and b"Purple Book" in r.content[:200]:
                (PB_DIR / "purplebook.csv").write_bytes(r.content)
                logger.info(f"Purple Book refreshed from {url}")
                return True
        except Exception as e:
            logger.debug(f"Purple Book download {url} failed: {e}")
    return False


def _ensure() -> bool:
    csv_path = PB_DIR / "purplebook.csv"
    if csv_path.exists() and (time.time() - csv_path.stat().st_mtime) < _MAX_AGE:
        return True
    if _download():
        return True
    return csv_path.exists()      # use stale copy if a refresh fails


# ────────────────────────────────────────────────────────────────────────────
# parse
# ────────────────────────────────────────────────────────────────────────────
def _parse_date(s: str) -> Optional[date]:
    s = (s or "").strip()
    if not s or s.upper() in ("N/A", "NA", "NONE"):
        return None
    for fmt in ("%d-%b-%y", "%d-%b-%Y", "%m/%d/%Y", "%Y-%m-%d", "%b %d, %Y"):
        try:
            return datetime.strptime(s, fmt).date()
        except ValueError:
            continue
    return None


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", (s or "").lower()).strip()


def _rows():
    """Yield product dicts from the cached CSV, skipping the report preamble."""
    path = PB_DIR / "purplebook.csv"
    if not path.exists():
        return
    with open(path, encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f)
        header = None
        for row in reader:
            if header is None:
                # the real header row begins with the "N/R/U" change-flag column
                if row and row[0].strip().upper() in ("N/R/U", "N/R/U "):
                    header = [h.strip() for h in row]
                continue
            if not any(c.strip() for c in row):
                continue
            yield dict(zip(header, row))


def _load_index() -> Dict:
    """Build reference-biologic records with their biosimilar competitors."""
    global _index
    if _index is not None:
        return _index
    if not _ensure():
        _index = {"refs": {}, "as_of": None}
        return _index

    refs: Dict[str, Dict] = {}          # normalized proper name -> reference record

    def _ref(proper: str) -> Dict:
        key = _norm(proper)
        if key not in refs:
            refs[key] = {
                "proper_name": proper.strip(),
                "brands": set(), "applicants": set(), "blas": set(),
                "routes": set(), "n_presentations": 0,
                "centers": set(),
                "first_approval": None, "ref_exclusivity_exp": None,
                "patent_list": False,
                "biosimilars": [],           # list of {name, brand, applicant, bla, approval, interchangeable}
                "_bs_blas": set(),
            }
        return refs[key]

    biosimilar_rows = []
    for r in _rows():
        proper = (r.get("Proper Name") or "").strip()
        ltype = (r.get("License Type") or "").strip()
        if not proper or proper.upper() == "N/A":
            continue

        if "351(a)" in ltype:                       # reference / original biologic
            rec = _ref(proper)
            rec["n_presentations"] += 1
            brand = (r.get("Proprietary Name") or "").strip()
            if brand and brand.upper() != "N/A":
                rec["brands"].add(brand)
            appl = (r.get("Applicant") or "").strip()
            if appl:
                rec["applicants"].add(appl)
            bla = (r.get("BLA Number") or "").strip()
            if bla:
                rec["blas"].add(bla)
            route = (r.get("Route of Administration") or "").strip()
            if route:
                rec["routes"].add(route)
            center = (r.get("Center") or "").strip().upper()
            if center:
                rec["centers"].add(center)
            appr = _parse_date(r.get("Approval Date", "")) or \
                _parse_date(r.get("Date of First Licensure", ""))
            if appr and (rec["first_approval"] is None or appr < rec["first_approval"]):
                rec["first_approval"] = appr
            rex = _parse_date(r.get("Ref. Product Exclusivity Exp. Date", "")) or \
                _parse_date(r.get("Exclusivity Expiration Date", ""))
            if rex and (rec["ref_exclusivity_exp"] is None or rex > rec["ref_exclusivity_exp"]):
                rec["ref_exclusivity_exp"] = rex
            if (r.get("Patent List Provided") or "").strip().upper() == "YES":
                rec["patent_list"] = True

        elif "351(k)" in ltype:                     # biosimilar / interchangeable
            biosimilar_rows.append(r)

    # attach biosimilars to their reference product
    for r in biosimilar_rows:
        ref_proper = (r.get("Ref. Product Proper Name") or "").strip()
        if not ref_proper or ref_proper.upper() == "N/A":
            # some rows name the ref only by proper name of the biosimilar's molecule
            ref_proper = (r.get("Proper Name") or "").strip()
        rec = _ref(ref_proper)
        bla = (r.get("BLA Number") or "").strip()
        brand = (r.get("Proprietary Name") or "").strip()
        # interchangeability is signalled by an "Inter. Approval Date" (the date FDA
        # granted interchangeable status) or an interchangeable supplement number —
        # the "Licensure" column just reads "Licensed" for every biosimilar.
        interchangeable = bool(_parse_date(r.get("Inter. Approval Date", ""))) or \
            (r.get("Inter. Supplement Number") or "").strip() not in ("", "0")
        # one product spans many presentation rows — collapse on brand (fallback BLA)
        dedup_key = _norm(brand) or bla
        if dedup_key and dedup_key in rec["_bs_blas"]:
            if interchangeable:
                for b in rec["biosimilars"]:
                    if (_norm(b["brand"]) or b["bla"]) == dedup_key:
                        b["interchangeable"] = True
            continue
        if dedup_key:
            rec["_bs_blas"].add(dedup_key)
        rec["biosimilars"].append({
            "name": (r.get("Proper Name") or "").strip(),
            "brand": brand,
            "applicant": (r.get("Applicant") or "").strip(),
            "bla": bla,
            "approval": (r.get("Approval Date") or "").strip(),
            "interchangeable": interchangeable,
        })

    _index = {"refs": refs, "as_of": _as_of()}
    return _index


def _as_of() -> Optional[str]:
    path = PB_DIR / "purplebook.csv"
    if not path.exists():
        return None
    try:
        with open(path, encoding="utf-8-sig") as f:
            first = f.readline()
        m = re.search(r"Report\s*-\s*(.+)$", first.strip())
        if m:
            return m.group(1).strip()
    except Exception:
        pass
    return datetime.fromtimestamp(path.stat().st_mtime).strftime("%B %Y")


# ────────────────────────────────────────────────────────────────────────────
# market value — each biologic's most valuable indication (Repurposing Value Score)
# ────────────────────────────────────────────────────────────────────────────
def build_value_cache() -> Dict[str, Dict]:
    """Resolve every reference biologic → its indications → best disease value score,
    and persist {normalized proper name: {value_score, tier, best_indication}} to disk.
    Reuses the platform's local-first ChEMBL resolution + disease_value table, so it
    needs the chembl_33 DB and the disease_value reference table. Run offline / on
    refresh; the radar reads the cached JSON at request time. Fail-soft per drug."""
    try:
        from services.reverse_repurposing import _local_chembl_id_for_name, _local_indications
        from services.disease_value import value_for
    except Exception as e:
        logger.warning("build_value_cache: dependencies unavailable: %s", e)
        return {}

    idx = _load_index()
    out: Dict[str, Dict] = {}
    for key, rec in idx.get("refs", {}).items():
        # try the proper (INN) name first, then any brand
        names = [rec["proper_name"]] + sorted(rec["brands"])
        cid = ""
        for nm in names:
            cid = _local_chembl_id_for_name(nm)
            if cid:
                break
        if not cid:
            continue
        best = None
        for ind in _local_indications([cid]):
            dv = value_for(ind.get("name", ""), ind.get("efo_id", ""))
            if dv and (best is None or dv["value_score"] > best["value_score"]):
                best = {"value_score": round(dv["value_score"], 3), "tier": dv["tier"],
                        "best_indication": ind.get("name", "")}
        if best:
            out[key] = best
    try:
        PB_DIR.mkdir(parents=True, exist_ok=True)
        _VALUE_CACHE.write_text(json.dumps(out), encoding="utf-8")
        logger.info("build_value_cache: %d biologics scored", len(out))
    except Exception as e:
        logger.warning("build_value_cache: write failed: %s", e)
    return out


def _values() -> Dict[str, Dict]:
    global _value_map
    if _value_map is not None:
        return _value_map
    if _VALUE_CACHE.exists():
        try:
            _value_map = json.loads(_VALUE_CACHE.read_text(encoding="utf-8"))
        except Exception:
            _value_map = {}
    else:
        _value_map = {}
    return _value_map


# ────────────────────────────────────────────────────────────────────────────
# development activity — biosimilar clinical trials & literature per molecule
# ────────────────────────────────────────────────────────────────────────────
_ACTIVITY_CACHE = PB_DIR / "activity_cache.json"
_activity_map: Optional[Dict[str, Dict]] = None

_CT_API = "https://clinicaltrials.gov/api/v2/studies"
_NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _ct_term(molecule: str) -> str:
    # ClinicalTrials.gov Essie syntax: require both the (quoted) molecule and
    # "biosimilar" so an unmatched name can't broaden to the whole corpus.
    return f'"{molecule.strip()}" AND biosimilar'


def _pm_term(molecule: str, recent: bool = False) -> str:
    # Field-tag + quote the molecule. Without [tiab], PubMed auto-term-mapping drops
    # an unrecognised drug name and collapses the query to "biosimilar" alone —
    # returning the entire biosimilar literature (the bulevirtide=3525 bug).
    t = f'"{molecule.strip()}"[tiab] AND biosimilar[tiab]'
    return t + " AND 2021:3000[dp]" if recent else t


def _ct_count(molecule: str) -> int:
    """Total ClinicalTrials.gov studies for a biosimilar of this molecule."""
    from services import http_client
    d = http_client.get_json(_CT_API, default={},
                             params={"query.term": _ct_term(molecule),
                                     "countTotal": "true", "pageSize": 1,
                                     "format": "json"})
    try:
        return int(d.get("totalCount") or 0)
    except (TypeError, ValueError):
        return 0


def _pubmed_count(term: str) -> int:
    """PubMed hit count for a fully-formed search term."""
    from services import http_client
    d = http_client.get_json(f"{_NCBI}/esearch.fcgi", default={},
                             params={"db": "pubmed", "term": term, "retmax": 0,
                                     "retmode": "json"})
    try:
        return int(d.get("esearchresult", {}).get("count") or 0)
    except (TypeError, ValueError):
        return 0


def _activity_counts(molecule: str) -> Dict:
    """Biosimilar-development footprint for one reference molecule: how many
    clinical trials and (recent) publications concern a biosimilar of it."""
    n_trials = _ct_count(molecule)
    n_papers = _pubmed_count(_pm_term(molecule))
    n_recent = _pubmed_count(_pm_term(molecule, recent=True))
    idx = n_trials * 3 + n_recent          # trials weigh more than a paper
    return {"n_trials": n_trials, "n_papers": n_papers,
            "n_recent_papers": n_recent, "activity_index": idx}


def build_activity_cache(center: str = "CDER", sleep: float = 0.34) -> Dict[str, Dict]:
    """Precompute biosimilar clinical-trial + literature counts per reference biologic
    and cache to disk, so the radar can rank/badge by development activity without
    hitting live APIs per request. Defaults to the CDER universe (the biosimilar-
    relevant molecules). Run offline / on refresh."""
    idx = _load_index()
    want = (center or "").strip().upper()
    out: Dict[str, Dict] = {}
    names = []
    for key, rec in idx.get("refs", {}).items():
        if rec["n_presentations"] == 0 and not rec["biosimilars"]:
            continue
        if want and rec["centers"] and want not in rec["centers"]:
            continue
        names.append((key, rec["proper_name"]))
    for i, (key, proper) in enumerate(names):
        try:
            out[key] = _activity_counts(proper)
        except Exception as e:
            logger.debug("activity %s failed: %s", proper, e)
        if sleep:
            time.sleep(sleep)
        if (i + 1) % 25 == 0:
            logger.info("build_activity_cache: %d/%d", i + 1, len(names))
    try:
        PB_DIR.mkdir(parents=True, exist_ok=True)
        _ACTIVITY_CACHE.write_text(json.dumps(out), encoding="utf-8")
        logger.info("build_activity_cache: %d biologics", len(out))
    except Exception as e:
        logger.warning("build_activity_cache: write failed: %s", e)
    return out


def _activity() -> Dict[str, Dict]:
    global _activity_map
    if _activity_map is not None:
        return _activity_map
    if _ACTIVITY_CACHE.exists():
        try:
            _activity_map = json.loads(_ACTIVITY_CACHE.read_text(encoding="utf-8"))
        except Exception:
            _activity_map = {}
    else:
        _activity_map = {}
    return _activity_map


def biosimilar_evidence(molecule: str, max_items: int = 8) -> Dict:
    """Live detail for one biologic: its biosimilar clinical trials + key papers.
    Used by the per-card 'Evidence & activity' panel (lazy, on demand)."""
    from services import http_client
    q = (molecule or "").strip()
    if not q:
        return {"trials": [], "papers": []}

    trials = []
    d = http_client.get_json(_CT_API, default={},
                             params={"query.term": _ct_term(q), "countTotal": "true",
                                     "pageSize": max_items, "format": "json"})
    n_trials = int(d.get("totalCount") or 0) if isinstance(d, dict) else 0
    for s in (d.get("studies", []) if isinstance(d, dict) else []):
        ps = s.get("protocolSection", {})
        nct = ps.get("identificationModule", {}).get("nctId", "")
        trials.append({
            "nct": nct,
            "title": ps.get("identificationModule", {}).get("briefTitle", "")[:110],
            "status": ps.get("statusModule", {}).get("overallStatus", ""),
            "phase": ", ".join(ps.get("designModule", {}).get("phases", [])) or "N/A",
            "url": f"https://clinicaltrials.gov/study/{nct}" if nct else "",
        })

    papers = []
    sr = http_client.get_json(f"{_NCBI}/esearch.fcgi", default={},
                              params={"db": "pubmed", "term": _pm_term(q),
                                      "retmax": max_items, "retmode": "json",
                                      "sort": "date"})
    ids = sr.get("esearchresult", {}).get("idlist", []) if isinstance(sr, dict) else []
    n_papers = 0
    try:
        n_papers = int(sr.get("esearchresult", {}).get("count") or 0)
    except (TypeError, ValueError):
        pass
    if ids:
        smr = http_client.get_json(f"{_NCBI}/esummary.fcgi", default={},
                                   params={"db": "pubmed", "id": ",".join(ids),
                                           "retmode": "json"})
        res = smr.get("result", {}) if isinstance(smr, dict) else {}
        for p in ids:
            if p in res:
                papers.append({"pmid": p, "title": res[p].get("title", "")[:110],
                               "journal": res[p].get("source", ""),
                               "year": (res[p].get("pubdate", "") or "")[:4],
                               "url": f"https://pubmed.ncbi.nlm.nih.gov/{p}/"})
    return {"trials": trials, "papers": papers,
            "n_trials": n_trials, "n_papers": n_papers}


# ────────────────────────────────────────────────────────────────────────────
# opportunity scoring
# ────────────────────────────────────────────────────────────────────────────
def _cliff(rec: Dict) -> Optional[date]:
    """Best estimate of when reference-product exclusivity lifts (competition opens)."""
    if rec["ref_exclusivity_exp"]:
        return rec["ref_exclusivity_exp"]
    if rec["first_approval"]:
        try:
            fa = rec["first_approval"]
            return fa.replace(year=fa.year + _REF_EXCLUSIVITY_YEARS)
        except ValueError:                       # Feb-29 edge
            return rec["first_approval"].replace(year=rec["first_approval"].year + _REF_EXCLUSIVITY_YEARS, day=28)
    return None


def _score(rec: Dict, cliff: Optional[date], today: date) -> Dict:
    """
    0-100 opportunity score from Purple Book facts alone:
      * timing   — is the exclusivity cliff already open, imminent, or far off
      * whitespace — how few biosimilars have entered (crowded space is worth less)
      * franchise  — bigger reference product (more presentations) = bigger prize
      * friction   — an active patent list signals a thicket (small penalty)
    """
    n_bs = len(rec["biosimilars"])

    # timing: peaks when the cliff is open but recent, or opening soon
    if cliff is None:
        timing, timing_note = 30, "exclusivity date unknown"
    else:
        months = (cliff.year - today.year) * 12 + (cliff.month - today.month)
        if months <= 0:
            years_open = -months // 12
            if n_bs > 0:
                timing, timing_note = 100, "market open — biosimilars present"
            elif years_open <= 8:
                timing, timing_note = 100, "exclusivity expired — market open"
            else:
                # open for many years with zero entrants: likely not a viable target
                timing, timing_note = 45, f"open {years_open} yr, no entrants"
        elif months <= 24:
            timing, timing_note = 85, f"cliff in ~{months} mo"
        elif months <= 60:
            timing, timing_note = 55, f"cliff in ~{months // 12} yr"
        else:
            timing, timing_note = 20, f"protected ~{months // 12} yr"

    # whitespace: 0 biosimilars = wide open; falls off as the field fills
    whitespace = max(0, 100 - n_bs * 22)

    # franchise size proxy: number of distinct product presentations
    franchise = min(100, 30 + rec["n_presentations"] * 8)

    friction = 8 if rec["patent_list"] else 0

    score = round(0.45 * timing + 0.35 * whitespace + 0.20 * franchise - friction, 1)
    score = max(0.0, min(100.0, score))

    if score >= 70:
        tier = "High opportunity"
    elif score >= 50:
        tier = "Attractive"
    elif score >= 32:
        tier = "Watch"
    else:
        tier = "Low / protected"

    return {"score": score, "tier": tier, "timing_note": timing_note,
            "n_biosimilars": n_bs, "whitespace": whitespace}


def opportunity_radar(limit: int = 200, only_open: bool = False,
                      query: str = "", center: str = "CDER",
                      sort: str = "opportunity") -> Dict:
    """
    Ranked reference biologics as biosimilar targets.

    only_open — restrict to reference products whose exclusivity has already lifted.
    query     — case-insensitive filter on proper/brand name.
    center    — regulatory center; defaults to CDER (therapeutic mAbs & proteins, the
                real biosimilar universe). CBER products (vaccines, plasma, allergenics)
                are excluded unless center is set to "" (all).
    sort      — "opportunity" (default), "value" (market attractiveness),
                "combined" (0.6 x opportunity + 0.4 x market value), or
                "activity" (biosimilar clinical-trial + literature footprint).
    """
    idx = _load_index()
    if not idx.get("refs"):
        return {"available": False, "as_of": None, "rows": [], "counts": {}}

    today = date.today()
    q = _norm(query)
    want_center = (center or "").strip().upper()
    vmap = _values()
    amap = _activity()
    rows = []
    for rec in idx["refs"].values():
        # skip records that are only referenced (e.g. a biosimilar whose ref product
        # is not itself CDER-listed) and have no franchise footprint
        if rec["n_presentations"] == 0 and not rec["biosimilars"]:
            continue
        # restrict to the biosimilar-relevant regulatory center (CDER by default)
        if want_center and rec["centers"] and want_center not in rec["centers"]:
            continue
        if q and q not in _norm(rec["proper_name"]) and \
           not any(q in _norm(b) for b in rec["brands"]):
            continue
        cliff = _cliff(rec)
        sc = _score(rec, cliff, today)
        is_open = bool(cliff and cliff <= today)
        if only_open and not is_open:
            continue
        n_inter = sum(1 for b in rec["biosimilars"] if b["interchangeable"])
        val = vmap.get(_norm(rec["proper_name"]))
        act = amap.get(_norm(rec["proper_name"]))
        # combined priority: opportunity blended with market value when we have it
        combined = round(0.6 * sc["score"] + 0.4 * val["value_score"] * 100, 1) \
            if val else sc["score"]
        rows.append({
            "proper_name": rec["proper_name"],
            "brands": sorted(rec["brands"]),
            "applicants": sorted(rec["applicants"]),
            "routes": sorted(rec["routes"]),
            "first_approval": rec["first_approval"].isoformat() if rec["first_approval"] else None,
            "cliff": cliff.isoformat() if cliff else None,
            "cliff_open": is_open,
            "n_biosimilars": sc["n_biosimilars"],
            "n_interchangeable": n_inter,
            "patent_list": rec["patent_list"],
            "score": sc["score"], "tier": sc["tier"],
            "timing_note": sc["timing_note"],
            "value": val,                      # {value_score, tier, best_indication} or None
            "activity": act,                   # {n_trials, n_papers, n_recent_papers, activity_index} or None
            "combined_score": combined,
            "biosimilars": sorted(rec["biosimilars"],
                                  key=lambda b: b["brand"] or b["name"]),
        })

    if sort == "value":
        rows.sort(key=lambda r: (r["value"]["value_score"] if r["value"] else -1), reverse=True)
    elif sort == "combined":
        rows.sort(key=lambda r: r["combined_score"], reverse=True)
    elif sort == "activity":
        rows.sort(key=lambda r: (r["activity"]["activity_index"] if r["activity"] else -1), reverse=True)
    else:
        rows.sort(key=lambda r: r["score"], reverse=True)
    counts = {
        "reference_products": len(rows),
        "cliff_open": sum(1 for r in rows if r["cliff_open"]),
        "whitespace_open": sum(1 for r in rows if r["cliff_open"] and r["n_biosimilars"] == 0),
    }
    return {"available": True, "as_of": idx.get("as_of"),
            "has_value": bool(vmap), "has_activity": bool(amap), "sort": sort,
            "rows": rows[:limit], "counts": counts}


def biologic_status(name: str) -> Dict:
    """Purple Book record for a single molecule (reference or biosimilar), by name."""
    idx = _load_index()
    key = _norm(name)
    rec = idx.get("refs", {}).get(key)
    if not rec:
        # maybe the name is a biosimilar brand — search
        for r in idx.get("refs", {}).values():
            if any(key == _norm(b["brand"]) or key == _norm(b["name"]) for b in r["biosimilars"]):
                rec = r
                break
    if not rec:
        return {"available": True, "found": False}
    cliff = _cliff(rec)
    return {
        "available": True, "found": True,
        "proper_name": rec["proper_name"], "brands": sorted(rec["brands"]),
        "cliff": cliff.isoformat() if cliff else None,
        "cliff_open": bool(cliff and cliff <= date.today()),
        "n_biosimilars": len(rec["biosimilars"]),
        "n_interchangeable": sum(1 for b in rec["biosimilars"] if b["interchangeable"]),
        "biosimilars": rec["biosimilars"],
    }


if __name__ == "__main__":
    r = opportunity_radar(limit=15)
    print(f"Purple Book as of {r['as_of']} — {r['counts']}")
    for row in r["rows"][:15]:
        print(f"  {row['score']:5.1f}  {row['tier']:16}  {row['proper_name'][:28]:28}  "
              f"biosims={row['n_biosimilars']} inter={row['n_interchangeable']}  "
              f"cliff={row['cliff']}  {row['timing_note']}")
