"""
FDA Orange Book — authoritative drug patent & exclusivity status.

The Orange Book is THE source for which patents and regulatory exclusivities
protect an approved drug, with real expiry dates (so active vs expired is a fact,
not an estimate). Used for the Patents tab and the 505(b)(2) dossier.

Data: the FDA "Orange Book Data Files" zip (patent.txt, products.txt,
exclusivity.txt), downloaded + cached locally and refreshed monthly.
"""
import csv
import io
import logging
import re
import time
import zipfile
from datetime import date, datetime
from pathlib import Path
from typing import Dict, List, Optional

import requests

logger = logging.getLogger(__name__)

OB_PAGE = "https://www.fda.gov/drugs/drug-approvals-and-databases/orange-book-data-files"
OB_MEDIA = "https://www.fda.gov/media/76860/download?attachment"
OB_DIR = Path(__file__).parent.parent / "data" / "orange_book"
_MAX_AGE = 30 * 86400          # refresh monthly
_UA = {"User-Agent": "Mozilla/5.0"}

_index: Optional[Dict] = None


def _resolve_zip_url() -> str:
    """The media id occasionally changes — fall back to scraping the data-files page."""
    try:
        r = requests.get(OB_PAGE, timeout=20, headers=_UA)
        if r.ok:
            for m in re.findall(r'href=["\']([^"\']*media/\d+/download[^"\']*)', r.text, re.I):
                # the data zip is the only application/zip on the page; try it first
                url = m if m.startswith("http") else "https://www.fda.gov" + m
                if "attachment" not in url:
                    url += ("&" if "?" in url else "?") + "attachment"
                try:
                    h = requests.head(url, timeout=15, headers=_UA, allow_redirects=True)
                    if "zip" in h.headers.get("Content-Type", "").lower():
                        return url
                except Exception:
                    continue
    except Exception as e:
        logger.debug(f"OB page scrape failed: {e}")
    return OB_MEDIA


def _download() -> bool:
    OB_DIR.mkdir(parents=True, exist_ok=True)
    for url in (OB_MEDIA, _resolve_zip_url()):
        try:
            r = requests.get(url, timeout=90, headers=_UA)
            if r.ok and r.content[:2] == b"PK":
                zipfile.ZipFile(io.BytesIO(r.content)).extractall(OB_DIR)
                logger.info(f"Orange Book refreshed from {url}")
                return True
        except Exception as e:
            logger.debug(f"OB download {url} failed: {e}")
    return False


def _ensure() -> bool:
    pat, prod = OB_DIR / "patent.txt", OB_DIR / "products.txt"
    if pat.exists() and prod.exists() and (time.time() - pat.stat().st_mtime) < _MAX_AGE:
        return True
    if _download():
        return True
    return pat.exists() and prod.exists()      # use stale copy if a refresh fails


def _rows(fname: str):
    path = OB_DIR / fname
    if not path.exists():
        return
    with open(path, encoding="latin-1", newline="") as f:
        yield from csv.DictReader(f, delimiter="~")


def _load_index() -> Dict:
    global _index
    if _index is not None:
        return _index
    if not _ensure():
        _index = {}
        return _index
    ing_to_appl: Dict[str, set] = {}
    appl_info: Dict[str, Dict] = {}
    for row in _rows("products.txt"):
        appl = (row.get("Appl_No") or "").strip()
        ing = (row.get("Ingredient") or "").strip()
        if not appl or not ing:
            continue
        for token in ing.split(";"):                       # combination products
            ing_to_appl.setdefault(token.strip().upper(), set()).add(appl)
        appl_info[appl] = {
            "trade": (row.get("Trade_Name") or "").strip(),
            "applicant": (row.get("Applicant_Full_Name") or row.get("Applicant") or "").strip(),
        }
    pats_by_appl: Dict[str, List[Dict]] = {}
    for row in _rows("patent.txt"):
        appl = (row.get("Appl_No") or "").strip()
        if appl:
            pats_by_appl.setdefault(appl, []).append(row)
    excl_by_appl: Dict[str, List[Dict]] = {}
    for row in _rows("exclusivity.txt"):
        appl = (row.get("Appl_No") or "").strip()
        if appl:
            excl_by_appl.setdefault(appl, []).append(row)
    _index = {"ing": ing_to_appl, "pat": pats_by_appl, "excl": excl_by_appl, "info": appl_info}
    return _index


def _parse_date(s: str) -> Optional[date]:
    s = (s or "").strip()
    for fmt in ("%b %d, %Y", "%Y-%m-%d", "%m/%d/%Y"):
        try:
            return datetime.strptime(s, fmt).date()
        except ValueError:
            continue
    return None


_EXCL_CODES = {
    "NCE": "New Chemical Entity (5 yr)", "ODE": "Orphan Drug (7 yr)",
    "NP": "New Product (3 yr)", "NPP": "New Patient Population (3 yr)",
    "NDF": "New Dosage Form (3 yr)", "NS": "New Strength (3 yr)",
    "NC": "New Combination (3 yr)", "PED": "Pediatric (6 mo)",
    "I": "New Indication (3 yr)", "M": "Misc. (3 yr)", "GAIN": "Antibiotic (5 yr)",
}


def orange_book_protection(drug_name: str) -> Dict:
    """All Orange Book patents + exclusivities for a drug, with active/expired status."""
    idx = _load_index()
    if not idx:
        return {"available": False, "patents": [], "exclusivities": []}
    name = (drug_name or "").strip().upper()
    if not name:
        return {"available": True, "patents": [], "exclusivities": []}

    appls = set()
    for ing, a in idx["ing"].items():
        if name == ing or (len(name) > 3 and (name in ing or ing in name)):
            appls |= a

    today = date.today()
    patents, seen = [], set()
    for appl in appls:
        info = idx["info"].get(appl, {})
        for row in idx["pat"].get(appl, []):
            pno = (row.get("Patent_No") or "").strip()
            exp_txt = (row.get("Patent_Expire_Date_Text") or "").strip()
            if not pno or (pno, exp_txt) in seen:
                continue
            seen.add((pno, exp_txt))
            expd = _parse_date(exp_txt)
            kinds = []
            if (row.get("Drug_Substance_Flag") or "").strip().upper() == "Y":
                kinds.append("drug substance")
            if (row.get("Drug_Product_Flag") or "").strip().upper() == "Y":
                kinds.append("drug product")
            if (row.get("Patent_Use_Code") or "").strip():
                kinds.append("method-of-use")
            pid_url = (f"https://patents.google.com/patent/US{pno}/en"
                       if pno[:1].isdigit() else f"https://patents.google.com/patent/{pno}/en")
            patents.append({
                "id": pno, "expiry": exp_txt, "expiry_iso": expd.isoformat() if expd else None,
                "status": "active" if (expd and expd >= today) else "expired",
                "type": ", ".join(kinds) or "—",
                "use_code": (row.get("Patent_Use_Code") or "").strip(),
                "trade": info.get("trade", ""), "applicant": info.get("applicant", ""),
                "url": pid_url, "source": "FDA Orange Book",
            })

    exclusivities = []
    for appl in appls:
        info = idx["info"].get(appl, {})
        for row in idx["excl"].get(appl, []):
            code = (row.get("Exclusivity_Code") or "").strip()
            exp_txt = (row.get("Exclusivity_Date") or "").strip()
            expd = _parse_date(exp_txt)
            exclusivities.append({
                "code": code, "label": _EXCL_CODES.get(code, _EXCL_CODES.get(code.split("-")[0], code)),
                "expiry": exp_txt, "expiry_iso": expd.isoformat() if expd else None,
                "status": "active" if (expd and expd >= today) else "expired",
                "trade": info.get("trade", ""),
            })

    patents.sort(key=lambda p: (p["status"] != "active", p.get("expiry_iso") or "9999"))
    # de-dupe exclusivities (code+date)
    ex_seen, ex_uniq = set(), []
    for e in sorted(exclusivities, key=lambda x: (x["status"] != "active", x.get("expiry_iso") or "9999")):
        k = (e["code"], e["expiry"])
        if k not in ex_seen:
            ex_seen.add(k); ex_uniq.append(e)
    return {"available": True, "patents": patents, "exclusivities": ex_uniq,
            "n_applications": len(appls)}


if __name__ == "__main__":
    import json
    for d in ["sildenafil", "tadalafil", "apixaban"]:
        r = orange_book_protection(d)
        print(d, "->", len(r["patents"]), "patents,", len(r["exclusivities"]), "exclusivities")
        for p in r["patents"][:4]:
            print("   ", p["id"], p["status"], p["expiry"], p["type"])
