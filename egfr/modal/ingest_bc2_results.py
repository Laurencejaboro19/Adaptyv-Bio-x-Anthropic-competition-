"""
MODULE: scripts/ingest_bc2_results.py
Pull a BindCraft2 Modal campaign off the results volume into runs/ and register every
design in the DesignDB, so nothing useful depends on the Modal volume or a session.

WHAT IS KEPT. For each shard under <tag>/ on the `raptorbind-bc2-results` volume:
  - campaign.json, result.json, records.json, campaign.log, timed.log
  - the whole BindCraft2 project tree (1_Trajectories, 2_Refolded, 3_Ranked)
Locally under runs/bc2_<tag>/<shard>/.

WHAT GOES IN THE DB (runs/designs.db, src/design_db.py):
  run      bc2_<tag>, config = every shard's campaign JSON, target = human 305-515 seq
  project  egfr_conditional_h409 (registered once)
  designs  one row per MPNN candidate BindCraft2 refolded (2_Refolded/!_Refolded.csv),
           lot = bc<NNN>.<letter>: backbone number per trajectory, letter per candidate.
           Accepted candidates (3_Ranked) get stage "generate"; the rest
           "generate_candidate", so a later screen can recover near-misses.
  metrics  every numeric column of the candidate's row, stage "generate"; strings such
           as Interface_Target_Residues go to runs/bc2_<tag>/designs.json
  artifacts the per-shard campaign.json (kind dataset) and designs.json (kind table)
Idempotent: existing (run, design_id) rows keep their lot; new backbones continue the
numbering.

INDEX:
  pull_tag -- modal volume get <tag>/ into runs/bc2_<tag>/ (skips if already local)
  load_shards -- read every shard's campaign/result/records + refolded + ranked tables
  assign_lots -- deterministic bc lots per trajectory/candidate, reusing DB rows
  ingest -- register run/project/designs/metrics/artifacts, write designs.json
  ingest_refold_metrics -- helper the refold stage calls to add refold_<species> metrics
  main -- CLI

USAGE:
  python scripts/ingest_bc2_results.py --tag egfr_h409_ablate            # pull + ingest
  python scripts/ingest_bc2_results.py --tag egfr_h409_ablate --no-pull  # local tree only
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
from design_db import DesignDB  # noqa: E402
from lot import assign_backbone_lot  # noqa: E402


def append_variant(backbone_lot: str, index: int) -> str:
    """bc004.a .. bc004.z, then bc004.aa, bc004.ab, ... -- sequence_candidates can be 30+,
    which lot.py's single-letter variant_letter (0-25) cannot express."""
    letters = ""
    i = index
    while True:
        letters = chr(ord("a") + i % 26) + letters
        i = i // 26 - 1
        if i < 0:
            break
    return f"{backbone_lot}.{letters}"

VOLUME = "raptorbind-bc2-results"
DB_PATH = REPO / "runs" / "designs.db"
PROJECT_ID = "egfr_conditional_h409"
GENERATOR = "bindcraft2"
LOT_GENERATOR = "bindcraft"          # lot.py prefix "bc"; the designs row says bindcraft2
TARGET_FASTA = REPO / "inputs" / "targets" / "Adaptyv EGFR" / "egfr_d3ext_human.pdb"
THREE_TO_ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E',
                'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
                'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}


def _target_seq() -> str:
    return "".join(THREE_TO_ONE[l[17:20]] for l in TARGET_FASTA.read_text().splitlines()
                   if l.startswith("ATOM") and l[12:16].strip() == "CA")


def pull_tag(tag: str, dest_root: Path, force: bool = False) -> Path:
    """Download <tag>/ from the volume into runs/bc2_<tag>/ (contents at top level)."""
    dest = dest_root / f"bc2_{tag}"
    if dest.exists() and any(dest.iterdir()) and not force:
        print(f"already local: {dest} (use --force-pull to refresh)")
        return dest
    tmp = dest_root / f".pull_{tag}"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    cmd = [sys.executable, "-m", "modal", "volume", "get", "--force", VOLUME, tag, str(tmp)]
    print("pulling:", " ".join(cmd[2:]))
    subprocess.run(cmd, check=True)
    src = tmp / tag if (tmp / tag).exists() else tmp
    dest.mkdir(parents=True, exist_ok=True)
    for child in src.iterdir():
        target = dest / child.name
        if target.exists():
            shutil.rmtree(target) if target.is_dir() else target.unlink()
        shutil.move(str(child), str(target))
    shutil.rmtree(tmp, ignore_errors=True)
    print(f"pulled into {dest}")
    return dest


def _read_csv(path: Path) -> list:
    return list(csv.DictReader(path.open(newline=""))) if path.exists() else []


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def load_shards(local: Path) -> list:
    """One dict per shard directory with everything the ingest needs."""
    shards = []
    dirs = sorted(p for p in local.iterdir() if p.is_dir() and "seed" in p.name)
    if not dirs and (local / "project").exists():
        dirs = [local]           # an unsharded single-container campaign (bindcraft2_campaign)
    for d in dirs:
        campaign = json.loads((d / "campaign.json").read_text()) if (d / "campaign.json").exists() else {}
        result = json.loads((d / "result.json").read_text()) if (d / "result.json").exists() else {}
        status = json.loads((d / "status.json").read_text()) if (d / "status.json").exists() else {}
        project = d / "project"
        trajectories = []
        refolded = []
        ranked = []
        for csvp in project.rglob("!_Trajectories.csv"):
            trajectories += _read_csv(csvp)
        for csvp in project.rglob("!_Refolded.csv"):
            refolded += _read_csv(csvp)
        for csvp in project.rglob("!_Ranked.csv"):
            ranked += _read_csv(csvp)
        structures = {p.stem: p for p in project.rglob("*.cif")} if project.exists() else {}
        structures.update({p.stem: p for p in project.rglob("*.pdb")} if project.exists() else {})
        shards.append({"dir": d, "shard": d.name, "campaign": campaign, "result": result,
                       "status": status, "trajectories": trajectories, "refolded": refolded,
                       "ranked": ranked, "structures": structures,
                       "label": status.get("label") or ("dual" if status.get("dual") else "single"),
                       "seed": campaign.get("campaign_seed"), "gpu": status.get("gpu")})
    return shards


def _structure_for(structures: dict, design: str) -> Optional[Path]:
    if design in structures:
        return structures[design]
    hits = sorted(p for stem, p in structures.items() if design in stem)
    # prefer the ranked/accepted complex, then a refolded complex, over a monomer
    hits.sort(key=lambda p: (0 if "3_Ranked" in str(p) else 1 if "Complexes" in str(p) else 2, str(p)))
    return hits[0] if hits else None


def assign_lots(db: DesignDB, run_id: str, shards: list) -> dict:
    """{design_id: (lot, backbone_lot)}; existing rows keep their lot."""
    existing = {r["design_id"]: r["lot"] for r in db.conn.execute(
        "SELECT design_id, lot FROM designs WHERE run_id = ?", (run_id,))}
    used = set()
    for lot in existing.values():
        head = lot.split(".")[0]
        if head.startswith("bc") and head[2:].isdigit():
            used.add(int(head[2:]))
    next_n = (max(used) + 1) if used else 1
    out, backbone_of = {}, {}
    for sh in shards:
        traj_names = sorted((t.get("design") or "").strip() for t in sh["trajectories"])
        traj_names = [t for t in traj_names if t]
        groups: dict = {}
        for row in sh["refolded"]:
            did = (row.get("design") or "").strip()
            if not did:
                continue
            backbone = next((t for t in sorted(traj_names, key=len, reverse=True) if did.startswith(t)), did)
            groups.setdefault(backbone, []).append(did)
        for backbone, members in sorted(groups.items()):
            known = [existing[m].split(".")[0] for m in members if m in existing]
            if known:
                blot = known[0]
            else:
                blot = assign_backbone_lot(LOT_GENERATOR, next_n)
                next_n += 1
            backbone_of[backbone] = blot
            taken = {existing[m] for m in members if m in existing}
            idx = 0
            for m in members:
                if m in existing:
                    out[m] = (existing[m], blot)
                    continue
                while append_variant(blot, idx) in taken:
                    idx += 1
                out[m] = (append_variant(blot, idx), blot)
                taken.add(out[m][0])
                idx += 1
    return out


def ingest(tag: str, local: Path, db_path: Path = DB_PATH, project_id: str = PROJECT_ID) -> dict:
    shards = load_shards(local)
    run_id = f"bc2_{tag}"
    target_seq = _target_seq()
    designs_json = []
    n_acc = n_cand = 0
    with DesignDB(db_path) as db:
        db.register_project(
            project_id, name="EGFR domain III conditional binder (pH 6.4 on / 7.4 off, human+mouse)",
            target="EGFR_d3ext 305-515 (human; mouse threaded)",
            goal="De novo 60-120 aa binder on the H409 face, >=2 His in the interface, "
                 "binds at pH 6.4 and releases at 7.4, cross-reactive human/mouse",
            issue_ref="PDFA / Adaptyv x Anthropic Challenge 1")
        done = all(sh["status"].get("state") in ("done", "cancelled") or sh["result"] for sh in shards)
        db.register_run(run_id, config={"tag": tag, "volume": VOLUME, "local": str(local),
                                        "shards": {sh["shard"]: sh["campaign"] for sh in shards}},
                        target=target_seq, status="completed" if done else "partial")
        lots = assign_lots(db, run_id, shards)
        for sh in shards:
            # 3_Ranked names the file <design>_seq0 while 2_Refolded says <design>_candidateN,
            # so an accepted design is recognised by its SEQUENCE, not its name.
            accepted_seqs = {(r.get("Binder_Sequence") or "").strip() for r in sh["ranked"]} - {""}
            ranked_by_seq = {(r.get("Binder_Sequence") or "").strip(): r for r in sh["ranked"]}
            for row in sh["refolded"]:
                did = (row.get("design") or "").strip()
                seq = (row.get("Binder_Sequence") or "").strip()
                if not did or not seq or did not in lots:
                    continue
                lot, blot = lots[did]
                is_acc = seq in accepted_seqs
                if is_acc:
                    # the ranked row carries rank and the final-model metrics; merge it in
                    row = {**row, **{k: v for k, v in ranked_by_seq[seq].items() if v not in (None, "")}}
                stage = "generate" if is_acc else "generate_candidate"
                structure = _structure_for(sh["structures"], did)
                if is_acc:
                    ranked_name = (ranked_by_seq[seq].get("design") or "").strip()
                    structure = _structure_for(sh["structures"], ranked_name) or structure
                db_id = db.insert_design(run_id=run_id, lot=lot, design_id=did, generator=GENERATOR,
                                         stage=stage, binder_seq=seq, parent_lot=blot,
                                         structure_path=str(structure) if structure else None)
                if structure:
                    db.update_structure_path(db_id, str(structure))
                db.link_design_to_project(db_id, project_id)
                metrics = {k: _num(v) for k, v in row.items() if _num(v) is not None and k != "design"}
                metrics["bc2_accepted"] = 1.0 if is_acc else 0.0
                metrics["binder_len"] = float(len(seq))
                if sh["seed"] is not None:
                    metrics["campaign_seed"] = float(sh["seed"])
                db.insert_metrics(db_id, "generate", metrics)
                strings = {k: v for k, v in row.items() if _num(v) is None and v not in (None, "")}
                designs_json.append({"lot": lot, "db_id": db_id, "design_id": did, "shard": sh["shard"],
                                     "variant": sh["label"], "seed": sh["seed"], "gpu": sh["gpu"],
                                     "accepted": is_acc, "binder_seq": seq,
                                     "structure_path": str(structure) if structure else None,
                                     "metrics": metrics, "strings": strings})
                n_acc += int(is_acc); n_cand += int(not is_acc)
            if (sh["dir"] / "campaign.json").exists():
                db.insert_artifact(project_id, kind="dataset", path=str(sh["dir"] / "campaign.json"),
                                   caption=f"BindCraft2 campaign JSON, {tag}/{sh['shard']} [{sh['label']}]",
                                   provenance={"run_id": run_id, "shard": sh["shard"], "seed": sh["seed"],
                                               "gpu": sh["gpu"], "label": sh["label"],
                                               "minutes": sh["result"].get("minutes"),
                                               "trajectories": sh["result"].get("trajectories_started"),
                                               "accepted": sh["result"].get("accepted")},
                                   tags=["bindcraft2", "egfr", "campaign-config", sh["label"]])
        out_json = local / "designs.json"
        out_json.write_text(json.dumps(designs_json, indent=2))
        db.insert_artifact(project_id, kind="table", path=str(out_json),
                           caption=f"Every BindCraft2 candidate under {tag}: lot, sequence, structure, metrics, interface residue strings",
                           provenance={"run_id": run_id, "n_accepted": n_acc, "n_candidates": n_cand},
                           tags=["bindcraft2", "egfr", "designs"])
    print(f"run {run_id}: {n_acc} accepted design(s) [stage generate], {n_cand} refolded candidate(s) "
          f"[stage generate_candidate], {len(shards)} shard(s); designs.json -> {out_json}")
    return {"run_id": run_id, "accepted": n_acc, "candidates": n_cand, "designs_json": str(out_json)}


def ingest_refold_metrics(db_path: Path, run_id: str, design_id: str, species: str, metrics: dict,
                          structure_path: Optional[str] = None) -> Optional[int]:
    """Attach Boltz-2 refold metrics (stage refold_<species>) to an ingested design."""
    with DesignDB(db_path) as db:
        row = db.lookup_design_id(run_id, design_id)
        if not row:
            return None
        db.insert_metrics(row["id"], f"refold_{species}", {k: _num(v) for k, v in metrics.items()})
        if structure_path and species == "human":
            db.update_structure_path(row["id"], structure_path)
        return row["id"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("INDEX:")[0])
    ap.add_argument("--tag", required=True)
    ap.add_argument("--no-pull", action="store_true", help="use runs/bc2_<tag>/ as is")
    ap.add_argument("--force-pull", action="store_true")
    ap.add_argument("--db", default=str(DB_PATH))
    args = ap.parse_args(argv)
    local = REPO / "runs" / f"bc2_{args.tag}"
    if not args.no_pull:
        local = pull_tag(args.tag, REPO / "runs", force=args.force_pull)
    if not local.exists():
        sys.exit(f"nothing at {local}")
    ingest(args.tag, local, Path(args.db))


if __name__ == "__main__":
    main()
