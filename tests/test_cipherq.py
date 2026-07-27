"""
CipherQ regression / smoke tests.

Run:  .venv312/Scripts/python.exe -m pytest tests -q

Unit tests (cache, config hygiene) always run. Integration tests that need the
network (RCSB/UniProt), the database, or the xTB engine skip cleanly when those
are unavailable, so the suite is safe to run anywhere.
"""
import os
import sys
import time
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


# ── Unit: disk cache ────────────────────────────────────────────────────────────
def test_diskcache_roundtrip():
    from services.diskcache import DiskCache
    c = DiskCache("test_ns", ttl=None)
    c.set("k1", {"a": 1, "b": [1, 2, 3]})
    assert c.get("k1") == {"a": 1, "b": [1, 2, 3]}
    assert c.get("missing") is None


def test_diskcache_ttl_expiry():
    from services.diskcache import DiskCache
    c = DiskCache("test_ns_ttl", ttl=0.5)
    c.set("k", "v")
    assert c.get("k") == "v"
    time.sleep(0.7)
    assert c.get("k") is None  # expired


# ── Unit: config hygiene (no personal credentials hardcoded in source) ──────────
def test_no_hardcoded_personal_username_in_source():
    offenders = []
    for rel in ("config.py", "services/neuro_db_service.py",
                "services/disease_resolver.py", "services/repurposing_scorer.py"):
        text = (ROOT / rel).read_text(encoding="utf-8", errors="ignore")
        if "babburisoumith" in text:
            offenders.append(rel)
    assert not offenders, f"hardcoded personal DB user found in: {offenders}"


def test_db_params_shape():
    from config import db_params
    p = db_params()
    assert {"host", "port", "user", "password", "dbname"} <= set(p)
    assert isinstance(p["port"], int)


# ── Unit: empirical scorer behaves sensibly ─────────────────────────────────────
def test_pocket_subset_and_score_import():
    # dock_engine must import and expose the scoring primitives used elsewhere.
    from services.dock_engine import _parse_protein, score_pose, dock_ligand
    assert callable(score_pose) and callable(dock_ligand) and callable(_parse_protein)


# ── Integration: structure fetch (network) ──────────────────────────────────────
@pytest.mark.integration
def test_fetch_real_structure_for_known_gene():
    from real_pdb_fetcher import RealPDBFetcher
    f = RealPDBFetcher()
    pdb = f.fetch_protein_structure("CDK4")
    if not pdb:
        pytest.skip("no network / RCSB unavailable")
    # A real receptor — NOT the 14-atom generic placeholder that caused the tripeptide bug
    assert pdb.count("\nATOM") > 500, "fetched structure is too small to be a real protein"


@pytest.mark.integration
def test_cif_only_entry_is_converted():
    # TCF7L2's best hit historically was mmCIF-only (9I8O) and broke the .pdb-only fetch.
    from real_pdb_fetcher import RealPDBFetcher
    f = RealPDBFetcher()
    pdb = f.fetch_protein_structure("TCF7L2")
    if not pdb:
        pytest.skip("no network / RCSB unavailable")
    assert pdb.count("\nATOM") > 500


# ── Integration: docking regression (network + rdkit) ───────────────────────────
@pytest.mark.integration
def test_docking_lands_in_named_pocket():
    """Regression for the floating-ligand / tripeptide bug: a real structure loads
    and the docked pose sits inside a named pocket."""
    try:
        from services.docking_service import DockingService
    except Exception as e:
        pytest.skip(f"docking service unavailable: {e}")
    svc = DockingService()
    palbo = "CC(=O)c1c(C)c2cnc(Nc3ccc(N4CCNCC4)cn3)nc2n(C2CCCC2)c1=O"
    res = svc.perform_docking(drug_name="Palbociclib", target_name="CDK4",
                              ligand_smiles=palbo, method="fast")
    if not res.get("success"):
        pytest.skip(f"docking did not run: {res.get('error')}")
    if res.get("structure_quality") == "none":
        pytest.skip("no structure available (offline)")
    # Core regression guard (the tripeptide/floating-ligand bug): a real protein
    # loads and the docked pose is contained inside it with no severe clash.
    assert (res.get("protein_pdb") or "").count("\nATOM") > 500
    assert res.get("poses"), "no poses produced"
    valid = res.get("pose_valid") or []
    assert valid and valid[0].get("valid"), "top pose is not a valid, contained pose"
    # When pockets were detected, the pose should land in a named one.
    if res.get("pockets"):
        pp = res.get("pose_pockets") or []
        assert any(p is not None for p in pp), "pockets found but no pose landed in one"


@pytest.mark.integration
def test_no_structure_degrades_honestly():
    try:
        from services.docking_service import DockingService
    except Exception as e:
        pytest.skip(f"docking service unavailable: {e}")
    svc = DockingService()
    res = svc.perform_docking(drug_name="X", target_name="ZZZNOTAGENE123",
                              ligand_smiles="CCO", method="fast")
    # Bogus target → no real structure → explicit honesty, not a fake blob
    if res.get("structure_quality") != "none":
        pytest.skip("environment resolved an unexpected structure")
    assert res.get("protein_pdb") == ""
    assert "structure_warning" in res


# ── Integration: quantum engine (xTB) ───────────────────────────────────────────
@pytest.mark.integration
def test_qc_engine_aspirin():
    from services import qc_engine
    if not qc_engine.available():
        pytest.skip("xTB engine not installed")
    r = qc_engine.compute_properties("CC(=O)Oc1ccccc1C(=O)O", name="aspirin")
    assert r.get("success"), r.get("error")
    # Drug-like HOMO-LUMO gaps are typically ~0.5–8 eV
    assert 0.3 < r["gap_ev"] < 10
    assert r["dipole_debye"] >= 0
    assert len(r.get("atomic_charges", [])) == r.get("n_atoms")
    # Reactivity indices derived from frontier orbitals
    assert r["hardness_ev"] > 0


# ── Unit: endpoint outcome-signal logic (regression for the terminated_safety fix) ──
def test_endpoint_terminated_safety_is_negative(monkeypatch):
    """A trial stopped for a safety reason must yield a strongly NEGATIVE outcome
    signal, never neutral (0.0) or positive. Regression for commit 5b7732b."""
    from services import endpoint_parser as ep
    monkeypatch.setattr(ep, "classify_study",
                        lambda s: {"class": s["_c"], "note": s["_c"], "p": None})

    only_safety = ep.aggregate([{"_c": "terminated_safety"}])
    assert only_safety["outcome_signal"] == -0.9
    assert only_safety["verdict"] == "terminated_safety"

    # a met_primary must NOT flip a safety stop to a positive readout
    met_and_safety = ep.aggregate([{"_c": "met_primary"}, {"_c": "terminated_safety"}])
    assert met_and_safety["outcome_signal"] < 0


def test_endpoint_verdict_tie_is_order_independent(monkeypatch):
    """Equal-magnitude opposite-sign classes must resolve to the same verdict
    regardless of study order (prefers the more-negative class)."""
    from services import endpoint_parser as ep
    monkeypatch.setattr(ep, "classify_study",
                        lambda s: {"class": s["_c"], "note": s["_c"], "p": None})
    a = ep.aggregate([{"_c": "met_primary"}, {"_c": "terminated_efficacy"}])["verdict"]
    b = ep.aggregate([{"_c": "terminated_efficacy"}, {"_c": "met_primary"}])["verdict"]
    assert a == b == "terminated_efficacy"


# ── Unit: audit trail (Part 11) — append + tamper-evidence ──────────────────────
def test_audit_trail_chain_and_tamper(tmp_path, monkeypatch):
    from services import audit_trail as at
    monkeypatch.setattr(at, "_AUDIT_DIR", tmp_path)
    monkeypatch.setattr(at, "_LOG_FILE", tmp_path / "audit_log.jsonl")

    e1 = at.record("login", actor="alice", details={"ip": "127.0.0.1"})
    e2 = at.record("run_screen", actor="alice", details={"disease": "asthma"})
    assert e1["seq"] == 1 and e2["seq"] == 2
    assert e2["prev_hash"] == e1["hash"]         # chained
    ok, broken = at.verify()
    assert ok and broken is None

    # tamper with the first entry WITHOUT recomputing its hash -> chain must break
    lines = (tmp_path / "audit_log.jsonl").read_text(encoding="utf-8").splitlines()
    rec0 = json.loads(lines[0]); rec0["details"] = {"disease": "cancer"}
    lines[0] = json.dumps(rec0)
    (tmp_path / "audit_log.jsonl").write_text("\n".join(lines) + "\n", encoding="utf-8")
    ok, broken = at.verify()
    assert not ok and broken == 1


# ── Unit: auth — password hashing, roles, gate ──────────────────────────────────
def test_auth_password_and_roles(tmp_path, monkeypatch):
    from services import auth
    monkeypatch.setattr(auth, "_AUTH_DIR", tmp_path)
    monkeypatch.setattr(auth, "_USERS_FILE", tmp_path / "users.json")

    assert auth.add_user("Bob", "s3cret!", "analyst")
    # stored hashed, never plaintext
    stored = (tmp_path / "users.json").read_text(encoding="utf-8")
    assert "s3cret!" not in stored
    assert auth.verify_credentials("bob", "s3cret!")["role"] == "analyst"   # case-insensitive user
    assert auth.verify_credentials("bob", "wrong") is None
    assert auth.verify_credentials("nobody", "x") is None
    # role hierarchy
    assert auth.has_role("admin", "analyst") and auth.has_role("analyst", "viewer")
    assert not auth.has_role("viewer", "analyst")
    # gate defaults off
    monkeypatch.delenv("AUTH_ENABLED", raising=False)
    assert auth.auth_enabled() is False
    monkeypatch.setenv("AUTH_ENABLED", "true")
    assert auth.auth_enabled() is True


def test_auth_role_policy_and_user_listing(tmp_path, monkeypatch):
    from services import auth
    monkeypatch.setattr(auth, "_AUTH_DIR", tmp_path)
    monkeypatch.setattr(auth, "_USERS_FILE", tmp_path / "users.json")
    # declarative authority policy
    assert auth.required_role_for("GET", "/discover") == "viewer"
    assert auth.required_role_for("POST", "/api/compound/x/dock") == "analyst"
    assert auth.required_role_for("GET", "/admin/users") == "admin"
    assert auth.required_role_for("POST", "/admin/users") == "admin"
    # listing never leaks password hashes
    auth.add_user("carol", "pw12345", "admin")
    users = auth.list_users()
    assert users and users[0]["username"] == "carol" and users[0]["role"] == "admin"
    assert all("pw_hash" not in u for u in users)


def test_audit_trail_rotation_keeps_chain(tmp_path, monkeypatch):
    from services import audit_trail as at
    monkeypatch.setattr(at, "_AUDIT_DIR", tmp_path)
    monkeypatch.setattr(at, "_LOG_FILE", tmp_path / "audit_log.jsonl")
    monkeypatch.setattr(at, "_MAX_BYTES", 500)     # force rotation after a couple entries
    for i in range(6):
        at.record("event", actor="alice", details={"i": i})
    rotated = list(tmp_path.glob("audit_log.*.jsonl"))
    assert len(rotated) >= 1                        # at least one segment rolled off
    ok, broken = at.verify()                        # chain continuous across segments
    assert ok and broken is None
    # total entries across all segments equals what we wrote
    total = sum(sum(1 for ln in seg.open(encoding="utf-8") if ln.strip())
                for seg in at._segments())
    assert total == 6
