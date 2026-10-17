"""
MODULE: scripts/modal_bc2_egfr_refold.py
Refold the EGFR H409 campaign's BindCraft2 designs with stock Boltz-2 against BOTH the
human and the mouse domain III, bring every structure home, and score each complex for
confidence, ipSAE, footprint agreement with the design-time contacts, and the pH switch.

Reuses `refold_score` / `fold_target` from scripts/modal_binder_campaign.py (same image,
same YAML, same ipSAE call, L4) so the numbers are comparable with every prior campaign in
this repo, and `bc2_results` from scripts/modal_uplift_kits.py for the design records.
Target MSAs are built once per species by fold_target and cached in the MSA volume.

Scoring is local (scripts/score_ph_switch.py) on the returned CIFs, so it needs
biopython, numpy and propka in the local interpreter.

INDEX:
  HUMAN_SEQ / MOUSE_SEQ -- the 305-515 constructs, read from the target files at import
  load_designs -- accepted designs (records.json) from every shard under a tag
  validate -- entrypoint: refold x2 species, save files, score, write validation.csv

USAGE:
  modal run --detach scripts/modal_bc2_egfr_refold.py --tag egfr_h409_l4 --out runs/egfr_h409_l4_refold
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import modal  # noqa: E402
from modal_binder_campaign import app, fold_target, refold_score  # noqa: E402
from modal_uplift_kits import bc2_results  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
TARGET_DIR = REPO / "inputs" / "targets" / "Adaptyv EGFR"
TARGET_OFFSET = 304          # file residue 1 of the refold == mature 305
ANCHOR = 409                 # the pH anchor His, mature numbering
THREE_TO_ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E',
                'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
                'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}


def _seq(pdb: Path) -> str:
    return "".join(THREE_TO_ONE[l[17:20]] for l in pdb.read_text().splitlines()
                   if l.startswith("ATOM") and l[12:16].strip() == "CA")


HUMAN_SEQ = _seq(TARGET_DIR / "egfr_d3ext_human.pdb")
MOUSE_SEQ = _seq(TARGET_DIR / "egfr_d3ext_mouse.pdb")
assert len(HUMAN_SEQ) == len(MOUSE_SEQ) == 211, (len(HUMAN_SEQ), len(MOUSE_SEQ))

_img = modal.Image.debian_slim().add_local_python_source(
    "modal_binder_campaign", "modal_uplift_kits", "modal_bc2_egfr_refold")


@app.function(image=_img, timeout=60 * 10, volumes={"/results": bc2_results})
def load_designs(tag: str, include_rejected: bool = False) -> list:
    """Accepted designs from every shard; optionally every MPNN candidate BindCraft2 refolded."""
    bc2_results.reload()
    root = Path("/results") / tag
    out = []
    for d in sorted(root.glob("*seed*")):      # seed4417 and base_seed5000 alike
        rj = d / "records.json"
        if rj.exists():
            for r in json.loads(rj.read_text()):
                out.append({"design_id": r["design_id"], "binder_seq": r["binder_seq"],
                            "seed": d.name, "accepted": True, "scores": r.get("scores", {})})
        if include_rejected:
            for csvp in (d / "project").rglob("!_Refolded.csv"):
                import io
                for row in csv.DictReader(io.StringIO(csvp.read_text())):
                    seq = (row.get("Binder_Sequence") or "").strip()
                    did = (row.get("design") or "").strip()
                    if seq and did and not any(o["design_id"] == did for o in out):
                        out.append({"design_id": did, "binder_seq": seq, "seed": d.name,
                                    "accepted": False, "scores": {
                                        k: row.get(k) for k in ("i_pTM", "i_pAE", "pLDDT",
                                        "Hotspot_Contact_Fraction", "Interface_Target_Residues",
                                        "Interface_H_Count", "outcome")}})
    print(f"{len(out)} design(s) under {tag} ({sum(o['accepted'] for o in out)} accepted)")
    return out


def _save_files(files: dict, dest: Path) -> dict:
    dest.mkdir(parents=True, exist_ok=True)
    paths = {}
    for rel, data in (files or {}).items():
        p = dest / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(data)
        paths[rel] = p
    return paths


@app.local_entrypoint()
def validate(tag: str = "egfr_h409_l4", out: str = "runs/egfr_h409_l4_refold",
             include_rejected: bool = False, use_msa: bool = True, select: str = "volume",
             skip_done: bool = False):
    """select: 'volume' (accepted designs read off the Modal volume), 'screened'
    (runs/bc2_<tag>/to_refold.json written by scripts/screen_bc2_candidates.py after the
    ingest), or 'accepted' (runs/bc2_<tag>/designs.json rows BindCraft2 accepted)."""
    from score_ph_switch import score_complex   # local scoring

    if select == "volume":
        designs = load_designs.remote(tag, include_rejected)
    else:
        local = REPO / "runs" / f"bc2_{tag}"
        if select == "screened":
            designs = json.loads((local / "to_refold.json").read_text())
        else:
            designs = [{"design_id": d["design_id"], "binder_seq": d["binder_seq"], "seed": d["shard"],
                        "accepted": True, "scores": {**d["metrics"], **d["strings"]}}
                       for d in json.loads((local / "designs.json").read_text()) if d["accepted"]]
        for d in designs:
            d.setdefault("seed", d.get("shard")); d.setdefault("accepted", False)
    out_dir = Path(out); out_dir.mkdir(parents=True, exist_ok=True)
    if skip_done:
        before = len(designs)
        designs = [d for d in designs if not all(
            (out_dir / d["design_id"] / sp / "score.json").exists() for sp in ("human", "mouse"))]
        print(f"skip_done: {before - len(designs)} design(s) already refolded under {out_dir}")
    if not designs:
        print("no designs to refold yet")
        return
    (out_dir / "designs.json").write_text(json.dumps(designs, indent=2))

    species = {"human": HUMAN_SEQ, "mouse": MOUSE_SEQ}
    if use_msa:
        # One MSA per species, cached in the volume for every refold that follows.
        for name, seq in species.items():
            fold_target.remote(seq, use_msa=True)
            print(f"MSA cached for {name}")

    jobs = [(d, sp) for d in designs for sp in species]
    payload = [{"design_id": f"{d['design_id']}__{sp}", "binder_seq": d["binder_seq"]} for d, sp in jobs]
    # One (design, target_seq, use_msa) tuple per job; starmap keeps the order so rows
    # line up with `jobs`.
    rows = list(refold_score.starmap(
        [(p, species[sp], use_msa) for p, (_, sp) in zip(payload, jobs)], order_outputs=True))

    table = []
    for (d, sp), r in zip(jobs, rows):
        row = {"design_id": d["design_id"], "seed": d["seed"], "accepted": d["accepted"],
               "species": sp, "binder_len": len(d["binder_seq"]),
               "bc2_iptm": d["scores"].get("i_pTM"), "bc2_ipae": d["scores"].get("i_pAE"),
               "bc2_hotspot_frac": d["scores"].get("Hotspot_Contact_Fraction"),
               "bc2_interface_H": d["scores"].get("Interface_H_Count")}
        if not isinstance(r, dict) or r.get("status") != "ok":
            row["status"] = "failed"; row["error"] = (r or {}).get("error") if isinstance(r, dict) else repr(r)
            table.append(row); continue
        row.update({k: r.get(k) for k in ("iptm", "ptm", "plddt", "ipae", "ipsae", "ipsae_AB", "ipsae_BA", "pdockq")})
        paths = _save_files(r.get("files"), out_dir / d["design_id"] / sp)
        cif = next((p for rel, p in paths.items() if rel.endswith("_model_0.cif")), None) or \
              next((p for rel, p in paths.items() if rel.endswith(".cif")), None)
        pae = next((p for rel, p in paths.items() if rel.endswith("_model_0.npz") and "pae_" in rel), None)
        row["cif"] = str(cif) if cif else None
        if cif:
            sc = score_complex(str(cif), "A", "B", TARGET_OFFSET, ANCHOR,
                               str(d["scores"].get("Interface_Target_Residues") or ""),
                               str(pae) if pae else None)
            row.update({
                "n_target_if": sc["n_target_if"], "n_binder_if": sc["n_binder_if"],
                "footprint": ",".join(map(str, sc["target_footprint"])),
                "footprint_hc": ",".join(map(str, sc["target_footprint_hc"])),
                "mean_contact_pae": sc.get("mean_contact_pae"),
                "anchor_contacted": sc.get("anchor_contacted"),
                "anchor_contacted_hc": sc.get("anchor_contacted_hc"),
                "jaccard_vs_bc2": (sc.get("footprint_vs_bc2") or {}).get("jaccard"),
                "jaccard_hc_vs_bc2": (sc.get("footprint_hc_vs_bc2") or {}).get("jaccard"),
                "n_if_his": sc["n_interface_his"], "n_if_his_good": sc["n_interface_his_good"],
                "n_if_his_bad": sc["n_interface_his_bad"],
                "anchor_d_carboxylate": (sc.get("anchor_geometry") or {}).get("d_carboxylate"),
                "anchor_carboxylate": (sc.get("anchor_geometry") or {}).get("carboxylate"),
                "anchor_d_cation": (sc.get("anchor_geometry") or {}).get("d_cation"),
                "anchor_burial": (sc.get("anchor_geometry") or {}).get("burial_on_binding"),
                "anchor_verdict": (sc.get("anchor_geometry") or {}).get("verdict"),
                "anchor_area_share": sc.get("anchor_area_share"),
                "anchor_pka_complex": sc.get("anchor_pka_complex"), "anchor_pka_free": sc.get("anchor_pka_free"),
                "log10_Ka_ratio_all": (sc.get("linkage_all_groups") or {}).get("log10_Ka_ratio_6.4_over_7.4"),
                "log10_Ka_ratio_his": (sc.get("linkage_interface_his") or {}).get("log10_Ka_ratio_6.4_over_7.4"),
                "his_detail": json.dumps(sc["interface_his"]),
                "propka_error": sc.get("propka_error"),
            })
            (out_dir / d["design_id"] / sp / "score.json").write_text(json.dumps(sc, indent=2))
        row["status"] = "ok"
        table.append(row)
        # Persist to the DesignDB alongside the generate-stage row the ingest created
        # (scripts/ingest_bc2_results.py), so the refold verdict is recoverable later.
        try:
            from ingest_bc2_results import ingest_refold_metrics, DB_PATH
            db_metrics = {k: row.get(k) for k in (
                "iptm", "ptm", "plddt", "ipae", "ipsae", "ipsae_AB", "ipsae_BA", "pdockq",
                "n_target_if", "n_binder_if", "mean_contact_pae", "jaccard_vs_bc2",
                "jaccard_hc_vs_bc2", "n_if_his", "n_if_his_good", "n_if_his_bad",
                "anchor_d_carboxylate", "anchor_d_cation", "anchor_burial", "anchor_area_share",
                "anchor_pka_complex", "anchor_pka_free", "log10_Ka_ratio_all", "log10_Ka_ratio_his")}
            db_metrics["anchor_contacted"] = 1.0 if row.get("anchor_contacted") else 0.0
            if ingest_refold_metrics(DB_PATH, f"bc2_{tag}", d["design_id"], sp, db_metrics,
                                     structure_path=row.get("cif")) is None:
                print(f"  (not in DesignDB yet: run `python scripts/ingest_bc2_results.py --tag {tag}` first)")
        except Exception as e:  # noqa: BLE001
            print(f"  DesignDB write skipped for {d['design_id']}/{sp}: {e}")

    vpath = out_dir / "validation.csv"
    if skip_done and vpath.exists():
        old_rows = list(csv.DictReader(vpath.open()))
        seen = {(t["design_id"], t["species"]) for t in table}
        table = [r for r in old_rows if (r["design_id"], r["species"]) not in seen] + table
    cols = sorted({k for t in table for k in t}, key=lambda k: (k not in ("design_id", "species", "accepted"), k))
    with vpath.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(table)
    print(f"\n{'design':<36}{'sp':>6}{'len':>5}{'iptm':>7}{'ipsae':>7}{'ipae':>7}{'H409':>6}{'jacc':>6}"
          f"{'nHis':>5}{'good':>5}{'bad':>4}{'dCOO':>6}{'dpKa':>6}{'logK':>6}")
    for t in table:
        if t.get("status") != "ok":
            print(f"{t['design_id'][:35]:<36}{t['species']:>6}  FAILED {t.get('error')}"); continue
        dpka = (t.get("anchor_pka_complex") or 0) - (t.get("anchor_pka_free") or 0) if t.get("anchor_pka_complex") else None
        print(f"{t['design_id'][:35]:<36}{t['species']:>6}{t['binder_len']:>5}"
              f"{(t.get('iptm') or float('nan')):>7.3f}{(t.get('ipsae') or float('nan')):>7.3f}"
              f"{(t.get('ipae') or float('nan')):>7.3f}{str(t.get('anchor_contacted')):>6}"
              f"{(t.get('jaccard_vs_bc2') if t.get('jaccard_vs_bc2') is not None else float('nan')):>6.2f}"
              f"{t.get('n_if_his', 0):>5}{t.get('n_if_his_good', 0):>5}{t.get('n_if_his_bad', 0):>4}"
              f"{(t.get('anchor_d_carboxylate') or float('nan')):>6.1f}"
              f"{(dpka if dpka is not None else float('nan')):>6.2f}"
              f"{(t.get('log10_Ka_ratio_all') if t.get('log10_Ka_ratio_all') is not None else float('nan')):>6.2f}")
    print(f"\nwrote {out_dir / 'validation.csv'}; structures under {out_dir}/<design>/<species>/")
