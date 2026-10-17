"""
MODULE: scripts/modal_bc2_egfr.py
EGFR domain III conditional (pH 6.4-on / 7.4-off), human+mouse cross-reactive binder campaign
on BindCraft2, fanned out across cheap GPUs: N containers x 1 accepted design each.

EPITOPE (analysis of 2026-10-01, agrees with the PDFA collaborator document): the H409 face
of domain III -- hotspots T406/H409/G410/Q411, mature numbering, identical in human and
mouse, free of glycan shadow, accessible in the tethered ectodomain (1NQL), and the direct
EGF / cetuximab / panitumumab contact surface. H409 is the only conserved, exposed His on a
bindable face; its protonation at pH 6.4 is the switch, so the binder must present
carboxylates to it and keep Lys/Arg away.

HOW THE INTERFACE-HISTIDINE REQUEST IS MET. BindCraft2 has no pH setting (filters.py only
reports Binder_Net_Charge at pH 7.4). Two real knobs are used instead:
  aa_bias            -- a GLOBAL propensity (1 neutral, 2 favoured, 0.4 disfavoured, 0 off)
                        applied in hallucination and MPNN redesign; H favoured, D/E nudged
                        up, K/R nudged down.
  filters.Interface_H_Count -- an ACCEPTANCE gate on the number of contacting BINDER His
                        (settings/core default.json has no such entry; outputs.md documents
                        Interface_<X>_Count for all 20 residues). With target H409 in the
                        footprint this gives >= 2 His in the interface.
Nothing in BindCraft2 can place a His at a chosen interface position; the geometry (which
target residue each His faces) is scored afterwards from the Boltz-2 refold.

WHY FAN OUT ON L4. Measured 2026-09-24: BindCraft2 fits L4 and A100-40GB; memory was never
the constraint. Each shard carries its own campaign_seed (identical seeds would return the
same design N times) and a max_trajectories cap, so wall-clock to K accepted designs is the
time for the K fastest shards, not N designs in series. The JAX compile cache is shared via
the weights volume so only the first shard per GPU type pays the compile.

INDEX:
  HOTSPOTS / COLDSPOTS / CAMPAIGN -- the EGFR H409 campaign, as BindCraft2 keys
  campaign_json -- adapter JSON + campaign keys (+ optional mouse second target)
  _run_shard -- one shard: stream + timestamp the log, persist records and project tree
  SHARDS -- Modal functions, one per GPU type
  launch -- entrypoint: spawn N single-target (+ M dual-target) shards and exit
  VARIANTS / ablate -- one-change-at-a-time shards to find what stops docking
  status -- entrypoint: per-shard progress read back from the results volume
  spend -- entrypoint: estimated $ per shard from container ages / run time; --record logs it

USAGE:
  export RB_ADAPTER_DIR=<dir with generators/ + egfr_d3ext_{human,mouse}.pdb>
  modal run --detach scripts/modal_bc2_egfr.py::launch --n 4 --dual-n 1 --tag egfr_h409_l4
  modal run scripts/modal_bc2_egfr.py::status --tag egfr_h409_l4
  modal volume get raptorbind-bc2-results egfr_h409_l4 <local_dir>
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

import modal  # noqa: E402
from modal_uplift_kits import (  # noqa: E402
    app, bindcraft2_image, bindcraft2_weights, bc2_results,
)

REPO = Path(__file__).resolve().parent.parent
TARGET_DIR = REPO / "inputs" / "targets" / "Adaptyv EGFR"


def _default_adapter_dir() -> str:
    """Build the /adapter payload from the repo itself: the two generator modules the
    container imports plus both target files. A fixed path under runs/ (not committed)
    rather than a session scratchpad, so the image is reproducible from a checkout."""
    import shutil
    d = REPO / "runs" / ".bc2_adapter_egfr"
    (d / "generators").mkdir(parents=True, exist_ok=True)
    (d / "generators" / "__init__.py").write_text(
        "from .base import Generator\nfrom .bindcraft2 import BindCraft2Generator\n")
    for m in ("base.py", "bindcraft2.py"):
        shutil.copy(REPO / "src" / "generators" / m, d / "generators" / m)
    for pdb in (HUMAN_PDB, MOUSE_PDB):
        shutil.copy(TARGET_DIR / pdb, d / pdb)
    return str(d)


HUMAN_PDB = "egfr_d3ext_human.pdb"
MOUSE_PDB = "egfr_d3ext_mouse.pdb"
# Built only on the client: the container imports this module too, and there the repo
# does not exist (the payload is already mounted at /adapter).
ADAPTER_DIR = os.environ.get("RB_ADAPTER_DIR") or (
    _default_adapter_dir() if modal.is_local() else "/adapter")

# Mature EGFR numbering, chain A in both files. The adapter strips the chain prefix.
HOTSPOTS = ["A406", "A409", "A410", "A411"]
# Species-divergent exposed surface within reach of the footprint (10 from the collaborator's
# list + S460P and G461N, the two substitutions nearest the hotspots at 10.2 and 8.0 A, which
# that list omitted) plus the five N-glycan sequons.
COLDSPOTS = "324,328,337,353,388,389,418,420,443,460,461,467,468,471,473,504"

# Campaign-level keys layered on BindCraft2's core -> modality(binder) -> property -> target.
# Every key here is in bindcraft/settings.py CAMPAIGN_SETTING_NAMES or a documented filter.
CAMPAIGN = {
    # Lysine trick: exposed target residues outside a 10 A shell around the hotspots become
    # Lys during screen/refine/anneal, so the binder has nowhere else to sit; the real
    # sequence is restored before hardening and validation. Property default floor is 0.5
    # of hotspots contacted; 0.75 = 3 of 4, so H409 is very likely in the footprint.
    "forced_targeting": True,
    "min_hotspot_contact_final": 0.75,
    # Both chain ends away from the target: the assay construct carries a tag.
    "termini_accessible": True,
    # Composition prior for the pH mechanism (see module docstring). C stays excluded.
    "aa_bias": {"C": 0, "H": 2.0, "D": 1.3, "E": 1.3, "K": 0.8, "R": 0.8},
    "filters": {
        # >= 1 contacting binder His; with target H409 that is >= 2 His in the interface.
        "Interface_H_Count": {"threshold": 1, "higher": True},
        # 12 contacting binder residues (core default 7): a His-tag-decoy guard, since a
        # binder gripping one His cannot reach 12 interface residues.
        "Interface_Residues": {"threshold": 12, "higher": True},
    },
    # Core defaults kept on purpose and recorded here so the run is self-describing:
    # min_iptm_final 0.7, max_ipae_final 0.35 (filters.i_pAE), pTM 0.55,
    # Unbound_Binder_pLDDT 0.8 (binder preset 0.7 via min_monomer_plddt_final),
    # Binder_RMSD <= 3.5, Backbone_Clashes 0, max_coldspot_contact_final 0.05 (set
    # automatically once coldspots exist), redesign_interface false, mpnn_variant negative,
    # desperation true (irrelevant under a 12-trajectory cap: it needs 750 dry
    # trajectories), autotune true.
    "workers_per_gpu": 1,          # one worker; the fan-out is across containers
}


def campaign_json(seed: int, project: Path, dual: bool, max_trajectories: int,
                  binder_lengths: list, n_designs: int = 1,
                  adapter_root: str = "/adapter", overrides: Optional[dict] = None) -> dict:
    """Assemble the campaign through the pipeline's own adapter, then layer the EGFR keys.

    `overrides` (a dict) is applied last: top-level keys replace CAMPAIGN's, `filters`
    merges by metric, and the pseudo-keys `coldspots` / `hotspots` rewrite every
    target's residue string (so an ablation can drop or thin the coldspots)."""
    import types
    sys.path.insert(0, adapter_root)
    stub = types.ModuleType("runner")
    stub.run_in_environment = lambda **kw: None
    sys.modules.setdefault("runner", stub)
    from generators.bindcraft2 import BindCraft2Generator

    gen = BindCraft2Generator()
    cfg = {"modality": "binder", "number_of_final_designs": n_designs,
           "campaign_name": "EGFR_d3_H409", "target_chain": "A",
           "binder_lengths": list(binder_lengths)}
    campaign = gen._build_campaign_config(
        Path(adapter_root) / HUMAN_PDB, HOTSPOTS, project, cfg)
    human = campaign["targets"][0]
    human["name"] = "EGFR_d3_human"
    human["coldspots"] = COLDSPOTS
    if dual:
        # Same epitope, same numbering, equal weight: BindCraft2's ortholog-pair mode
        # (examples/pdl1_ortholog_pair.json) rotates the design through both folds.
        campaign["targets"].append({
            "name": "EGFR_d3_mouse",
            "target_path": str(Path(adapter_root) / MOUSE_PDB),
            "chains": "A", "hotspots": human["hotspots"], "coldspots": COLDSPOTS,
            "objective": "target", "weight": 1.0,
        })
    campaign.update(json.loads(json.dumps(CAMPAIGN)))
    for key, value in (overrides or {}).items():
        if key in ("coldspots", "hotspots"):
            for t in campaign["targets"]:
                t[key] = value
        elif key == "filters":
            campaign.setdefault("filters", {}).update(value)
        else:
            campaign[key] = value
    campaign["campaign_seed"] = int(seed)
    campaign["max_trajectories"] = int(max_trajectories)
    return campaign


def _run_shard(seed: int, tag: str, dual: bool, max_trajectories: int,
               binder_lengths: str, n_designs: int = 1, overrides: str = "",
               workers: int = 1, label: str = "") -> dict:
    import shutil
    import subprocess
    import threading
    import time

    lengths = [int(x) for x in binder_lengths.split(",")]
    work = Path(f"/tmp/bc2_{tag}_{seed}")
    project = work / "project"
    project.mkdir(parents=True, exist_ok=True)

    gpu_name = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total",
                               "--format=csv,noheader"], capture_output=True,
                              text=True).stdout.strip()

    # Shared JAX compile cache: BindCraft2 writes <project>/compile_cache/<GPU name>/, so a
    # symlink to the weights volume lets every later shard on the same GPU type skip the
    # compile. Committed at the end of the run.
    shared_cache = Path("/weights/compile_cache")
    shared_cache.mkdir(parents=True, exist_ok=True)
    link = project / "compile_cache"
    if not link.exists():
        link.symlink_to(shared_cache)
    warm = sorted(p.name for p in shared_cache.iterdir()) if shared_cache.exists() else []

    ov = json.loads(overrides) if overrides else {}
    campaign = campaign_json(seed, project, dual, max_trajectories, lengths, n_designs,
                             overrides=ov)
    campaign["workers_per_gpu"] = int(workers)
    cpath = work / "campaign.json"
    cpath.write_text(json.dumps(campaign, indent=2))

    res_dir = Path("/results") / tag / (f"{label}_seed{seed}" if label else f"seed{seed}")
    res_dir.mkdir(parents=True, exist_ok=True)
    (res_dir / "campaign.json").write_text(cpath.read_text())
    (res_dir / "status.json").write_text(json.dumps(
        {"state": "running", "gpu": gpu_name, "dual": dual, "warm_cache": warm,
         "label": label, "workers": workers, "overrides": ov,
         "started_epoch": time.time()}))      # so `spend` needs no app mapping
    bc2_results.commit()

    peak = {"mib": 0}
    stop = threading.Event()

    def sample():
        while not stop.is_set():
            out = subprocess.run(["nvidia-smi", "--query-gpu=memory.used",
                                  "--format=csv,noheader,nounits"],
                                 capture_output=True, text=True).stdout.split()
            for v in out:
                try:
                    peak["mib"] = max(peak["mib"], int(v))
                except ValueError:
                    pass
            time.sleep(2.0)

    th = threading.Thread(target=sample, daemon=True)
    th.start()
    t0 = time.time()
    markers = []          # (elapsed_s, trajectory line)
    stages = []           # (elapsed_s, stage readout line)
    accepted_at = []      # elapsed_s when an acceptance was logged
    lines = []
    env = dict(os.environ, BINDCRAFT_WORKERS_PER_GPU=str(workers), BINDCRAFT_WEIGHTS="/weights")
    # Stream BindCraft2's stdout: plain log for the adapter, a timestamped copy for timing.
    with (res_dir / "campaign.log").open("w") as log, \
            (res_dir / "timed.log").open("w") as tlog:
        proc = subprocess.Popen(["bindcraft", "design", str(cpath)], cwd=str(work),
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, bufsize=1, env=env)
        for line in proc.stdout:
            el = time.time() - t0
            lines.append(line)
            log.write(line)
            tlog.write(f"{el:9.1f}s  {line}")
            if line.startswith("=== trajectory "):
                markers.append((round(el, 1), line.strip()))
                print(f"[seed {seed}] {el/60:6.1f} min  {line.strip()}", flush=True)
            if " design stage " in line or "trajectory rejected" in line:
                stages.append((round(el, 1), line.strip()))
                print(f"[seed {seed}] {el/60:6.1f} min  {line.strip()}", flush=True)
                log.flush(); tlog.flush()
                bc2_results.commit()           # every stage readout lands promptly
            if "accepted" in line and ("design accepted" in line or " accepted  " in line):
                accepted_at.append(round(el, 1))
            if len(lines) % 40 == 0:
                log.flush(); tlog.flush()
                bc2_results.commit()
        proc.wait()
    elapsed = time.time() - t0
    stop.set()
    th.join(timeout=3)

    blob = "".join(lines)
    oom = any(k in blob for k in ("RESOURCE_EXHAUSTED", "Out of memory", "out of memory"))

    # Parse with OUR adapter (3_Ranked/!_Ranked.csv -> candidates.jsonl).
    import types
    sys.path.insert(0, "/adapter")
    sys.modules.setdefault("runner", types.ModuleType("runner"))
    sys.modules["runner"].run_in_environment = lambda **kw: None
    from generators.bindcraft2 import BindCraft2Generator
    normalized = BindCraft2Generator().normalize_outputs(work, HOTSPOTS)
    records = []
    cj = Path(normalized["candidates_jsonl"])
    if cj.exists():
        records = [json.loads(ln) for ln in cj.read_text().splitlines() if ln.strip()]

    # Persist the whole project tree (minus the compile cache) so the design-time
    # complexes are available for footprint + pH scoring without re-running the GPU.
    dest = res_dir / "project"
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(project, dest, symlinks=False,
                    ignore=shutil.ignore_patterns("compile_cache"))
    for r in records:
        # Re-point the design-time pose at its copy inside the volume.
        p = Path(r["init_pose_path"])
        try:
            r["init_pose_path"] = str(dest / p.relative_to(project))
        except ValueError:
            pass

    traj_times = [round(b[0] - a[0], 1) for a, b in zip(markers, markers[1:])]
    result = {
        "seed": seed, "tag": tag, "label": label, "workers": workers, "overrides": ov,
        "gpu": gpu_name, "dual": dual, "rc": proc.returncode,
        "stages": stages,
        "minutes": round(elapsed / 60, 1), "peak_mib": peak["mib"], "oom": oom,
        "warm_cache_on_entry": warm,
        "compile_to_first_trajectory_s": markers[0][0] if markers else None,
        "trajectories_started": len(markers),
        "trajectory_seconds": traj_times,
        "accepted_at_s": accepted_at,
        "accepted": len(records),
        "records": records,
    }
    (res_dir / "records.json").write_text(json.dumps(records, indent=2))
    (res_dir / "result.json").write_text(json.dumps(result, indent=2))
    (res_dir / "status.json").write_text(json.dumps(
        {"state": "done", **{k: v for k, v in result.items() if k != "records"}}, indent=2))
    bc2_results.commit()
    bindcraft2_weights.commit()   # the shared compile cache
    print(f"[seed {seed}] DONE rc={proc.returncode} {elapsed/60:.1f} min peak={peak['mib']}MiB "
          f"trajectories={len(markers)} accepted={len(records)} "
          f"first_traj_at={result['compile_to_first_trajectory_s']}s "
          f"per_traj={traj_times}{'  OOM!' if oom else ''}", flush=True)
    if proc.returncode != 0:
        print(blob[-3000:], flush=True)
    return result


# Modal 1.x does not auto-include local Python modules; both this module and the one it
# imports must be in the image (see agent/memory/modal_detach_client_disconnect.md).
_image = (bindcraft2_image
          .add_local_dir(ADAPTER_DIR, "/adapter")
          .add_local_python_source("modal_uplift_kits", "modal_bc2_egfr"))


def _mk(gpu: str, name: str):
    fn = app.function(
        image=_image, gpu=gpu, timeout=60 * 180, name=name,
        volumes={"/weights": bindcraft2_weights, "/results": bc2_results},
        env={"BINDCRAFT_WEIGHTS": "/weights", "BINDCRAFT_WORKERS_PER_GPU": "1"},
        max_containers=16, retries=0,
    )
    return fn(_run_shard)


SHARDS = {g: _mk(g, f"bc2_egfr_{g.replace(':', '_').replace('-', '_')}")
          for g in ("L4", "A10", "A100-40GB", "H100")}


# The reader container imports THIS module (which imports modal_uplift_kits), so it needs
# the local sources too -- a bare debian_slim dies with ModuleNotFoundError.
_reader_image = modal.Image.debian_slim().add_local_python_source(
    "modal_uplift_kits", "modal_bc2_egfr")


@app.function(image=_reader_image, timeout=60 * 5,
              volumes={"/results": bc2_results})
def _read_status(tag: str) -> list:
    bc2_results.reload()
    root = Path("/results") / tag
    out = []
    if not root.exists():
        return out
    for d in sorted(root.glob("*seed*")):     # seed4417 and base_seed5000 alike
        st = {}
        if (d / "status.json").exists():
            st = json.loads((d / "status.json").read_text())
        tl = (d / "timed.log").read_text().splitlines() if (d / "timed.log").exists() else []
        marks = [ln for ln in tl if "=== trajectory " in ln]
        stage_lines = [ln for ln in tl if " design stage " in ln or "trajectory rejected" in ln]
        out.append({
            "shard": d.name, "state": st.get("state", "?"), "gpu": st.get("gpu"),
            "dual": st.get("dual"), "label": st.get("label"), "workers": st.get("workers"),
            "stages": [ln.strip()[:150] for ln in stage_lines],
            "trajectories_started": len(marks),
            "last_marker": marks[-1].strip() if marks else None,
            "last_elapsed_min": round(float(tl[-1].split("s")[0]) / 60, 1) if tl else None,
            "accepted": st.get("accepted"), "minutes": st.get("minutes"),
            "trajectory_seconds": st.get("trajectory_seconds"),
            "compile_s": st.get("compile_to_first_trajectory_s"),
            "tail": "\n".join(ln[12:] for ln in tl[-6:]),
        })
    return out


@app.local_entrypoint()
def launch(n: int = 4, dual_n: int = 0, gpu: str = "H100", tag: str = "egfr_h409",
           max_trajectories: int = 12, binder_lengths: str = "60,120", seed0: int = 4417,
           variant: str = "", workers: int = 3, n_designs: int = 1):
    """Spawn the shards and exit; results land in the volume under <tag>/[variant_]seed*/.

    `variant` names an entry of VARIANTS (overrides on top of CAMPAIGN); empty = CAMPAIGN
    as is. `workers` = BindCraft2 workers per container (3 fits an H100/H200 at ~9 GB
    each; 1 on an L4). .spawn() rather than .map(): nothing local waits, so a sleeping
    laptop costs nothing (run with --detach so the app itself outlives this process).
    """
    fn = SHARDS[gpu]
    ov = json.dumps(VARIANTS[variant]) if variant else ""
    handles = []
    for i in range(n):
        h = fn.spawn(seed=seed0 + i, tag=tag, dual=False, max_trajectories=max_trajectories,
                     binder_lengths=binder_lengths, n_designs=n_designs, overrides=ov,
                     workers=workers, label=variant)
        handles.append((variant or "single", seed0 + i, h.object_id))
    for j in range(dual_n):
        h = fn.spawn(seed=seed0 + 100 + j, tag=tag, dual=True,
                     max_trajectories=max_trajectories, binder_lengths=binder_lengths,
                     n_designs=n_designs, overrides=ov, workers=workers,
                     label=(variant + "_dual") if variant else "dual")
        handles.append(("dual", seed0 + 100 + j, h.object_id))
    for kind, seed, oid in handles:
        print(f"SPAWNED {kind:<6} seed={seed}  {oid}  -> {tag}/seed{seed}/")
    print(f"\n{len(handles)} shard(s) on {gpu} x {workers} worker(s); binder {binder_lengths} aa; "
          f"max {max_trajectories} trajectories each; variant {variant or 'CAMPAIGN'}; "
          f"hotspots {','.join(HOTSPOTS)}")
    print(f"follow with: modal run scripts/modal_bc2_egfr.py::status --tag {tag}")


# Which element of the campaign stops the binder docking? Each variant is one change
# from the baseline that produced i_pTM ~0.8 at screen on eGFP/PLA2G2A (hotspots only).
COLD_FAR = "328,337,353,388,389,420,443,468,471,473,504"   # drops 324,418,460,461,467 (<= 13 A)
BASE = {"forced_targeting": False, "min_hotspot_contact_final": 0.5,
        "termini_accessible": False, "aa_bias": {"C": 0}, "coldspots": ""}
VARIANTS = {
    "base":     dict(BASE),
    "forced":   {**BASE, "forced_targeting": True, "min_hotspot_contact_final": 0.75},
    "cold16":   {**BASE, "coldspots": COLDSPOTS},
    "coldfar":  {**BASE, "coldspots": COLD_FAR},
    "aabias":   {**BASE, "aa_bias": CAMPAIGN["aa_bias"]},
    "termini":  {**BASE, "termini_accessible": True},
    # candidate production config: everything except the near coldspots
    "fixed":    {"coldspots": COLD_FAR},
    # Rate-measurement candidates (2026-10-01): hallucination is the cost (~24 min/L4
    # trajectory) and MPNN+validation is cheap (~7 s per candidate on L4), so give each
    # hallucination 30 shots at the acceptance gates instead of 10.
    # The His-bias trajectory's 10 redesigns validated at i_pTM 0.20-0.29 with 3-6 interface
    # residues and 0-1 interface His, so the hard Interface_H_Count gate costs acceptances
    # without delivering His: it becomes a recorded metric (threshold 0 -- a NULL threshold
    # disables the check AND drops the measurement from the tables) and the 2-4 His
    # requirement is enforced at screening / MPNN lead-opt instead. 80-120 aa buries more
    # interface than 60-80 (eGFP 120-140 aa redesigns had 22-31 interface residues vs 3-16
    # here). desperation_trajectories 40 lets BindCraft2's own escalation ladder engage
    # after 40 dry trajectories; the Boltz-2 refold is the independent check on its
    # higher false-positive rate.
    "prodA":    {"coldspots": COLD_FAR, "sequence_candidates": 30, "binder_lengths": [80, 120],
                 "desperation_trajectories": 40,
                 "filters": {"Interface_H_Count": {"threshold": 0, "higher": True, "mandatory": False}}},  # forced + termini + aabias
    "prodB":    {"coldspots": COLD_FAR, "sequence_candidates": 30, "binder_lengths": [80, 120],
                 "desperation_trajectories": 40, "termini_accessible": False,
                 "filters": {"Interface_H_Count": {"threshold": 0, "higher": True, "mandatory": False}}},  # forced + aabias
    "prodC":    {**BASE, "coldspots": COLD_FAR, "aa_bias": CAMPAIGN["aa_bias"],
                 "sequence_candidates": 30, "binder_lengths": [80, 120],
                 "desperation_trajectories": 40,
                 "filters": {"Interface_H_Count": {"threshold": 0, "higher": True, "mandatory": False}}},  # hotspots + aabias, no lysine trick
    # PRODUCTION (settled 2026-10-01 from the rate run): prodC plus a one-coldspot
    # tolerance, so BindCraft2's own acceptance matches the post-hoc screen (its default
    # 0.05 ceiling rejected 8 of 10 i_pTM 0.80-0.84 redesigns for one far-coldspot contact).
    # Measured on one H100 x 3 workers: 21 trajectories/h, ~38% of trajectories give a
    # screened backbone, ~1 in 5 a both-species backbone after the Boltz-2 refold.
    "prod":     {**BASE, "coldspots": COLD_FAR, "aa_bias": CAMPAIGN["aa_bias"],
                 "sequence_candidates": 30, "binder_lengths": [80, 120],
                 "desperation_trajectories": 40, "max_coldspot_contact_final": 0.10,
                 "filters": {"Interface_H_Count": {"threshold": 0, "higher": True, "mandatory": False}}},
}


@app.local_entrypoint()
def ablate(tag: str = "egfr_h409_ablate", gpu: str = "L4", max_trajectories: int = 3,
           binder_lengths: str = "60,120", seed0: int = 5000,
           variants: str = ",".join(VARIANTS), h100_fixed: bool = True):
    """One shard per variant on `gpu`; optionally one `fixed` shard on H100 with three
    workers, to measure the per-trajectory cost on the fast tier for the same complex."""
    fn = SHARDS[gpu]
    for i, name in enumerate([v.strip() for v in variants.split(",") if v.strip()]):
        ov = VARIANTS[name]
        h = fn.spawn(seed=seed0 + i, tag=tag, dual=False, max_trajectories=max_trajectories,
                     binder_lengths=binder_lengths, overrides=json.dumps(ov), label=name)
        print(f"SPAWNED {name:<8} seed={seed0 + i} {gpu}  {h.object_id}")
    if h100_fixed:
        h = SHARDS["H100"].spawn(seed=seed0 + 50, tag=tag, dual=False,
                                 max_trajectories=max_trajectories * 2,
                                 binder_lengths=binder_lengths,
                                 overrides=json.dumps(VARIANTS["fixed"]), workers=3,
                                 label="fixed_h100w3")
        print(f"SPAWNED fixed_h100w3 seed={seed0 + 50} H100 x3 workers  {h.object_id}")
    print(f"follow with: modal run scripts/modal_bc2_egfr.py::status --tag {tag}")


@app.function(image=_reader_image, timeout=60 * 5, volumes={"/results": bc2_results})
def _read_shards(tags: list) -> list:
    """Per shard: state, gpu, label, minutes if done, last logged second if not."""
    bc2_results.reload()
    out = []
    for tag in tags:
        root = Path("/results") / tag
        if not root.exists():
            continue
        for d in sorted(root.glob("*seed*")):
            st = json.loads((d / "status.json").read_text()) if (d / "status.json").exists() else {}
            tl = (d / "timed.log").read_text().splitlines() if (d / "timed.log").exists() else []
            last_s = float(tl[-1].split("s")[0]) if tl else 0.0
            cancelled = (json.loads((d / "cancelled.json").read_text())
                         if (d / "cancelled.json").exists() else None)
            if cancelled:
                st["state"] = "cancelled"
            out.append({"tag": tag, "shard": d.name, "state": st.get("state", "?"),
                        "cancelled": cancelled, "started_epoch": st.get("started_epoch"),
                        "gpu": st.get("gpu", ""), "label": st.get("label") or ("dual" if st.get("dual") else "single"),
                        "minutes": st.get("minutes"), "last_logged_s": last_s,
                        "trajectories": st.get("trajectories_started"), "accepted": st.get("accepted")})
    return out


def _gpu_key(gpu_name: str) -> str:
    """nvidia-smi name -> scripts/modal_cost.py rate key."""
    n = (gpu_name or "").upper()
    for key in ("H200", "H100", "A100-80GB", "A100-40GB", "L40S", "A10", "L4", "T4"):
        if key.replace("-", " ") in n.replace("-", " ") or key in n:
            if key == "A100-40GB" and "80GB" in n:
                continue
            return key
    return "L4"


@app.local_entrypoint()
def spend(tags: str = "egfr_h409_l4,egfr_h409_ablate", record: bool = False,
          note: str = ""):
    """Estimated GPU spend per shard, from container ages (live) or run time (done).

    Modal exposes no billing endpoint, so this is rate x container-seconds, the same
    estimate scripts/modal_cost.py keeps for every other campaign. A cancelled shard
    (state still 'running' but no live container) is charged its last logged second,
    which is a lower bound. `--record` appends one record per shard to the cost log.
    """
    import subprocess
    from datetime import datetime
    from modal_cost import GPU_HOURLY_USD, record_usage

    now = datetime.now().astimezone()
    live = json.loads(subprocess.run([sys.executable, "-m", "modal", "container", "list", "--json"],
                                     capture_output=True, text=True).stdout or "[]")
    by_app = {}
    for c in live:
        by_app.setdefault(c["app_id"], []).append(
            datetime.fromisoformat(c["start_time"].replace(" ", "T")))
    # tag -> app id: the first launch and the ablation are different apps; map by which
    # app has as many live containers as the tag has running shards.
    shards = _read_shards.remote([t.strip() for t in tags.split(",") if t.strip()])
    running_by_tag = {}
    for sh in shards:
        if sh["state"] == "running":
            running_by_tag[sh["tag"]] = running_by_tag.get(sh["tag"], 0) + 1
    app_for_tag = {}
    for tag, n_run in sorted(running_by_tag.items(), key=lambda kv: -kv[1]):
        # exact container-count match first, then the largest app that fits
        # closest container count wins; a just-finished container can still be listed,
        # so the app may briefly show MORE containers than the tag has running shards
        cands = [(abs(len(st) - n_run), -len(st), aid) for aid, st in by_app.items()
                 if aid not in app_for_tag.values()]
        if cands:
            app_for_tag[tag] = sorted(cands)[0][2]
    rows, total = [], 0.0
    for sh in shards:
        key = _gpu_key(sh["gpu"])
        if key == "H200" and "h100" in sh["shard"].lower():
            key = "H100"          # requested H100, served by an H200, billed as H100
        rate = GPU_HOURLY_USD[key]
        if sh["state"] == "done" and sh["minutes"]:
            secs, how = sh["minutes"] * 60 + 60, "run time + 1 min start"
        elif sh.get("cancelled"):
            secs, how = sh["cancelled"]["container_seconds"], "cancelled: " + sh["cancelled"].get("note", "")[:60]
        elif sh.get("started_epoch"):
            secs, how = now.timestamp() - sh["started_epoch"] + 45, "own start time + 45 s pull"
        else:
            app_id = app_for_tag.get(sh["tag"])
            starts = by_app.get(app_id, [])
            if starts:
                secs, how = (now - min(starts)).total_seconds(), "live container age"
            else:
                secs, how = sh["last_logged_s"], "cancelled: last logged second (lower bound)"
        cost = rate * secs / 3600.0
        total += cost
        rows.append((sh["tag"], sh["shard"], sh["label"], key, sh["state"], secs / 60, cost, how,
                     sh["trajectories"], sh["accepted"]))
        if record:
            record_usage("modal_bc2_egfr.py", key, secs, 1, entrypoint=sh["tag"],
                         note=f"{sh['shard']} [{sh['label']}] {sh['state']}; {how}; "
                              f"trajectories={sh['trajectories']} accepted={sh['accepted']}" + (f"; {note}" if note else ""))
    print(f"{'tag':<18}{'shard':<22}{'variant':<13}{'gpu':<10}{'state':<9}{'min':>7}{'$':>7}  basis")
    for r in rows:
        print(f"{r[0]:<18}{r[1]:<22}{r[2]:<13}{r[3]:<10}{r[4]:<9}{r[5]:>7.1f}{r[6]:>7.2f}  {r[7]}")
    print(f"\nTOTAL estimated GPU spend: ${total:.2f} as of {now:%Y-%m-%d %H:%M %Z}"
          + ("  (recorded to runs/modal_cost_log.jsonl)" if record else ""))


@app.local_entrypoint()
def status(tag: str = "egfr_h409_l4"):
    rows = _read_status.remote(tag)
    if not rows:
        print(f"nothing under {tag}/ yet")
        return
    for r in rows:
        print(f"\n=== {r['shard']} [{'dual' if r['dual'] else 'single'}"
              f"{(' ' + r['label']) if r.get('label') else ''}] state={r['state']} "
              f"gpu={r['gpu']} elapsed={r['last_elapsed_min']} min "
              f"trajectories={r['trajectories_started']} accepted={r['accepted']} "
              f"compile_s={r['compile_s']} per_traj_s={r['trajectory_seconds']}")
        for ln in r.get("stages", []):
            print("   ", ln)


if __name__ == "__main__":
    # Print the single-target campaign JSON for inspection without touching Modal.
    print(json.dumps(campaign_json(4417, Path("<project>"), False, 12, [60, 120],
                                   adapter_root=ADAPTER_DIR), indent=2))
