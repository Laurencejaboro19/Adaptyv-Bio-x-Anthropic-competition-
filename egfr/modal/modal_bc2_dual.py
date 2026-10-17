#!/usr/bin/env python3
"""
PDFA - Adaptyv x Anthropic Challenge 1
Conditional EGFR binder: bind at pH 6.4, release at 7.4, human AND mouse.

10 accepted designs, both species as co-targets, BindCraft2 on Modal.

This supersedes modal_bc2_egfr.py. The image build and the campaign settings
below come from a campaign that actually ran (44 trajectories -> 14 screened
backbones -> 3 winners), not from reading documentation.

Run (PowerShell):
    python -m modal run modal_bc2_dual.py::selfcheck            # ~$0.30, proves the image
    python -m modal run --detach modal_bc2_dual.py::launch      # the campaign
    python -m modal run modal_bc2_dual.py::status
    python -m modal volume get pdfa-bc2-results egfr_dual ./results

--detach matters: Modal tears down an ephemeral app when the local client
disconnects, which kills a running campaign. Results are persisted to a Volume
rather than returned, so closing PowerShell is safe.
"""
from __future__ import annotations

import json
import os
import pathlib

import modal

# ===========================================================================
# campaign constants
# ===========================================================================
HOTSPOTS = "406,409,410,411"

# 11 coldspots. 324, 418 and 467 are deliberately ABSENT: measured at 13.0,
# 13.3 and 10.6 A from the nearest hotspot atom. Coldspots closer than ~15 A
# leave the binder nowhere to land -- an ablation with them in place collapsed
# i_pTM to 0.08 against 0.5-0.9 without. Protection against those three comes
# from the mouse co-target instead.
COLDSPOTS = "328,337,353,388,389,420,443,468,471,473,504"

BINDER_LENGTHS = [80, 120]
SHARDS = 10              # 10 shards x 1 accepted design = 10 designs
MAX_TRAJ_PER_SHARD = 30  # spend cap: 300 trajectories total
WORKERS_PER_GPU = 3      # ~11x the throughput of 1 worker, measured
BASE_SEED = 4417
GPU = "H100"

HUMAN = "egfr_d3ext_human.pdb"
MOUSE = "egfr_d3ext_mouse.pdb"

HERE = pathlib.Path(__file__).resolve().parent

# ===========================================================================
app = modal.App("pdfa-egfr-bc2-dual")

weights = modal.Volume.from_name("pdfa-bc2-weights", create_if_missing=True)
results = modal.Volume.from_name("pdfa-bc2-results", create_if_missing=True)

BINDCRAFT2_COMMIT = "301efdd1937fb963cc40a5b0ecc1bc2f9b2b2d15"

image = (
    # Plain Ubuntu, NOT a CUDA base and NOT debian_slim: BindCraft2's jax cuda13
    # wheels carry their own CUDA and cuDNN.
    modal.Image.from_registry("ubuntu:24.04", add_python=None)
    .apt_install("python3", "python3-pip", "python3-dev", "build-essential",
                 "ca-certificates", "git")
    .env({"PIP_BREAK_SYSTEM_PACKAGES": "1", "PIP_NO_CACHE_DIR": "1",
          "PYTHONUNBUFFERED": "1"})
    # ubuntu:24.04 has python3 but no `python`, and Modal cannot introspect the
    # interpreter version without one. A symlink, not add_python=, which would
    # install a SECOND interpreter and strand the editable install below.
    .run_commands("ln -sf /usr/bin/python3 /usr/local/bin/python")
    .run_commands(
        "git clone https://github.com/PacesaLab/BindCraft2.git /opt/bindcraft && "
        f"git -C /opt/bindcraft checkout {BINDCRAFT2_COMMIT}")
    # Editable on purpose: settings/ and scaffolds/ sit at the repo root and the
    # runtime resolves them relative to it. A plain install orphans both trees.
    .run_commands("python3 -m pip install -e '/opt/bindcraft[cuda13]'")
    .run_commands("python3 -m pip install 'gemmi==0.6.5' 'ipsae==1.0.1'")
    # THE IMPORTANT ONE. The cuda13 wheels hide their libraries under
    # site-packages/nvidia/*/lib, where the dynamic loader does not look.
    # Without this, jax warns once and then runs on CPU about 100x slower --
    # silently. The grep asserts the fix worked, so the build fails loudly
    # instead of the campaign quietly costing 100x.
    .run_commands(
        "python3 -c \"import nvidia, pathlib; print('\\n'.join(sorted(str(p) "
        "for r in nvidia.__path__ for p in pathlib.Path(r).glob('*/lib'))))\" "
        "> /etc/ld.so.conf.d/bindcraft-cuda.conf && ldconfig && "
        "ldconfig -p | grep -q libcupti")
    .run_commands("cd / && python3 -m bindcraft.selfcheck cuda13 --shipped-only")
    .env({"NVIDIA_VISIBLE_DEVICES": "all",
          "NVIDIA_DRIVER_CAPABILITIES": "compute,utility"})
    .add_local_file(HERE / HUMAN, f"/adapter/{HUMAN}", copy=True)
    .add_local_file(HERE / MOUSE, f"/adapter/{MOUSE}", copy=True)
)


def campaign_json(seed: int, project: str) -> dict:
    def target(name, pdb):
        return {
            "name": name,
            "target_path": f"/adapter/{pdb}",
            "chains": "A",
            "hotspots": HOTSPOTS,
            "coldspots": COLDSPOTS,
            "objective": "target",
            "weight": 1,
        }

    return {
        "campaign_name": "EGFR_d3_H409_dual",
        "modality": "binder",
        "project_folder": project,
        "campaign_seed": seed,
        "number_of_final_designs": 1,
        "max_trajectories": MAX_TRAJ_PER_SHARD,
        "desperation_trajectories": 40,
        "sequence_candidates": 30,

        # both species co-optimised during design, rather than designing on
        # human and filtering on mouse afterwards
        "targets": [target("EGFR_d3_human", HUMAN),
                    target("EGFR_d3_mouse", MOUSE)],

        "binder_lengths": BINDER_LENGTHS,

        # forced_targeting OFF; hotspot contact is enforced at acceptance
        # instead, which is what the successful run did.
        "forced_targeting": False,
        "min_hotspot_contact_final": 0.5,
        "max_coldspot_contact_final": 0.1,
        "termini_accessible": False,

        # H up: interface histidines are the switch. D/E up: they are the
        # carboxylate partners that raise the His pKa on binding. K/R down:
        # a cation beside a protonated His lowers the pKa and kills the
        # mechanism. C off: free cysteines scramble during expression.
        "aa_bias": {"C": 0, "H": 2.0, "D": 1.3, "E": 1.3, "K": 0.8, "R": 0.8},

        "filters": {
            "Interface_H_Count": {"threshold": 0, "higher": True,
                                  "mandatory": False},
            "Interface_Residues": {"threshold": 12, "higher": True},
        },
        "workers_per_gpu": WORKERS_PER_GPU,
    }


# ===========================================================================
@app.function(image=image, gpu="L40S", timeout=1800, volumes={"/weights": weights})
def selfcheck():
    """Proves the image before the campaign spends anything."""
    import subprocess
    out = {}
    out["gpu"] = subprocess.run(
        ["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader"],
        capture_output=True, text=True).stdout.strip()
    for cmd in (["python3", "-m", "bindcraft.selfcheck", "cuda13"],
                ["bindcraft", "design", "--help"]):
        r = subprocess.run(cmd, capture_output=True, text=True)
        out[" ".join(cmd)] = (r.stdout + r.stderr)[-1200:]
    # the thing that silently costs 100x if wrong
    r = subprocess.run(
        ["python3", "-c",
         "import jax; print('jax devices:', jax.devices())"],
        capture_output=True, text=True)
    out["jax_devices"] = (r.stdout + r.stderr)[-600:]
    # does our dual-target config parse?
    cfg = campaign_json(BASE_SEED, "/tmp/probe")
    cfg["max_trajectories"] = 0
    pathlib.Path("/tmp/probe.json").write_text(json.dumps(cfg, indent=1))
    r = subprocess.run(["bindcraft", "design", "/tmp/probe.json"],
                       capture_output=True, text=True, timeout=900)
    out["dry_run_rc"] = r.returncode
    out["dry_run_tail"] = (r.stdout + r.stderr)[-2000:]
    for k, v in out.items():
        print(f"\n===== {k} =====\n{v}")
    return out


@app.function(image=image, gpu=GPU, timeout=60 * 180,
              volumes={"/weights": weights, "/results": results})
def shard(i: int, tag: str) -> dict:
    import subprocess, time, threading

    seed = BASE_SEED + i
    work = pathlib.Path(f"/tmp/bc2_{tag}_{seed}")
    project = work / "project"
    project.mkdir(parents=True, exist_ok=True)

    # Shared JAX compile cache on the weights volume: BindCraft2 writes
    # <project>/compile_cache/<GPU name>/, so only the first shard on a given
    # GPU type pays the compile.
    cache = pathlib.Path("/weights/compile_cache")
    cache.mkdir(parents=True, exist_ok=True)
    link = project / "compile_cache"
    if not link.exists():
        link.symlink_to(cache)

    cfg = campaign_json(seed, str(project))
    cpath = work / "campaign.json"
    cpath.write_text(json.dumps(cfg, indent=2))

    out = pathlib.Path("/results") / tag / f"shard{i:02d}_seed{seed}"
    out.mkdir(parents=True, exist_ok=True)
    (out / "campaign.json").write_text(cpath.read_text())
    (out / "status.json").write_text(json.dumps(
        {"state": "running", "shard": i, "seed": seed, "started": time.time()}))
    results.commit()

    env = dict(os.environ, BINDCRAFT_WORKERS_PER_GPU=str(WORKERS_PER_GPU),
               BINDCRAFT_WEIGHTS="/weights")
    t0 = time.time()
    accepted = 0
    with (out / "campaign.log").open("w") as log:
        p = subprocess.Popen(["bindcraft", "design", str(cpath)], cwd=str(work),
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             text=True, bufsize=1, env=env)
        for line in p.stdout:
            log.write(line)
            if "accept" in line.lower():
                accepted += 1
            print(f"[{i:02d}] {line}", end="")
        rc = p.wait()
    dt = (time.time() - t0) / 60

    # keep everything the scoring step needs
    import shutil
    for sub in ("Accepted", "Trajectory", "final_design_stats.csv",
                "mpnn_design_stats.csv", "designs"):
        src = project / sub
        if src.exists():
            dst = out / sub
            if src.is_dir():
                shutil.copytree(src, dst, dirs_exist_ok=True)
            else:
                shutil.copy(src, dst)
    (out / "status.json").write_text(json.dumps(
        {"state": "done", "shard": i, "seed": seed, "rc": rc,
         "minutes": round(dt, 1), "accept_lines": accepted}))
    results.commit()
    weights.commit()
    print(f"[{i:02d}] rc={rc} {dt:.1f} min")
    return {"shard": i, "rc": rc, "minutes": round(dt, 1)}


@app.local_entrypoint()
def launch(tag: str = "egfr_dual", n: int = SHARDS):
    rows = list(shard.map(range(n), [tag] * n))
    mins = [r["minutes"] for r in rows]
    ok = sum(1 for r in rows if r["rc"] == 0)
    print("\n" + "=" * 60)
    print(f"tag                  : {tag}")
    print(f"shards clean         : {ok}/{len(rows)}")
    if mins:
        print(f"per shard            : {min(mins):.1f} - {max(mins):.1f} min")
        print(f"GPU-hours            : {sum(mins)/60:.1f}")
    print(f"\npull results:  python -m modal volume get pdfa-bc2-results {tag} ./results")


@app.local_entrypoint()
def status(tag: str = "egfr_dual"):
    print(read_status.remote(tag))


@app.function(image=modal.Image.debian_slim(), volumes={"/results": results})
def read_status(tag: str) -> str:
    import glob
    results.reload()
    out = []
    for f in sorted(glob.glob(f"/results/{tag}/*/status.json")):
        out.append(f"  {f.split('/')[-2]:28} {open(f).read()}")
    return "\n".join(out) or f"nothing under /results/{tag} yet"
