"""
MODULE: scripts/modal_uplift_kits.py
Modal packaging for the uplifting inference kits (#043) — Boltz-2 and BoltzGen.

WHY THIS EXISTS: the kits ship GPU config cards for a100/h100/h200 only. No consumer
card is supported, so the kits are a CLOUD-ONLY proposition for this project. All GPU work
here runs on Modal, so the kits are packaged as Modal images rather than local conda envs.

The images clone the kit repo at build time (sparse, blob-filtered) and replay the kit's
own environment/Dockerfile steps, so the pinned stack stays the kit's, byte for byte. We
deliberately do NOT re-resolve their dependencies.

INDEX:
  boltz2_kit_image -- Modal image for the Boltz-2 kit
  boltzgen_kit_image -- Modal image for the BoltzGen kit
  kit_check -- run `run.sh check` for a mode and report the exit code (the #043 WI-3 gate)
  kit_predict -- run one Boltz-2 complex through a given mode, return metrics
  main -- local entrypoint: gate check, then the stock-vs-exact-vs-fast A/B

USAGE:
  modal run scripts/modal_uplift_kits.py                      # check + A/B on the default case
  modal run scripts/modal_uplift_kits.py --only-check         # just the gate
  modal run scripts/modal_uplift_kits.py --gpu A100-80GB      # the other supported card
"""
from __future__ import annotations

import json
import os

import modal

app = modal.App("raptorbind-uplift-kits")

KIT_REPO = "https://github.com/anthropics/uplifting-biomolecular-modeling.git"
# Weight caches, one per kit, so a fresh container does not re-download.
boltz2_weights = modal.Volume.from_name("raptorbind-kit-boltz2-weights", create_if_missing=True)
boltzgen_weights = modal.Volume.from_name("raptorbind-kit-boltzgen-weights", create_if_missing=True)

# Modal GPU string. The kits' cards are a100/h100/h200; H100 is the pinned target and
# avoids the A100-only `exact` refusal bands documented in boltz2/STOCK.md.
GPU = os.environ.get("RB_KIT_GPU", "H100")
# Modal GPU string -> the kit's config card name.
CARD = {"H100": "h100", "H200": "h200", "A100": "a100", "A100-80GB": "a100", "A100-40GB": "a100"}


def _clone(kit: str) -> str:
    """Sparse-clone one kit plus the shared core it pins."""
    return (
        f"git clone --filter=blob:none --sparse --depth 1 {KIT_REPO} /kit && "
        f"cd /kit && git sparse-checkout set {kit} common"
    )


def _install_lock(kit: str) -> str:
    """Replay the kit Dockerfile's pip steps against its own requirements.lock.

    --no-deps throughout: the lock is complete, and a transitive resolve here could lift
    a package off the version the stack was tested on.
    """
    return (
        f"cd /kit/{kit} && "
        "grep -E '^(pip|wheel)==' environment/requirements.lock > /tmp/pip.txt && "
        "python -m pip install --no-cache-dir --no-deps -r /tmp/pip.txt && "
        "printf 'setuptools==82.0.1\\n' > /tmp/bc.txt && "
        "grep -v -E '^(#|pip==|wheel==|boltz==|boltzgen==)' environment/requirements.lock > /tmp/stack.txt && "
        "python -m pip install --no-cache-dir --no-deps --build-constraint /tmp/bc.txt -r /tmp/stack.txt"
    )


# Both kits are assembled on the interpreter their own Dockerfile names (3.11), with the
# C toolchain Triton needs to build kernels at first use.
_base = (
    modal.Image.from_registry("python:3.11.12-slim-bookworm", add_python=None)
    .apt_install("git", "build-essential", "gfortran", "curl")
)

boltz2_kit_image = (
    _base
    .run_commands(_clone("boltz2"))
    .run_commands(_install_lock("boltz2"))
    # Stock boltz from the wheel the kit vendors, not from PyPI.
    .run_commands(
        "python -m pip install --no-cache-dir --no-deps "
        "/kit/boltz2/stock/boltz-2.2.1-py3-none-any.whl"
    )
    # Pipeline-wide contract (CLAUDE.md): ipsae + gemmi everywhere. gemmi is already in
    # their lock at 0.6.5, our exact pin, so only ipsae is added.
    .run_commands("python -m pip install --no-cache-dir --no-deps ipsae==1.0.1")
    .run_commands("cd /kit/boltz2 && bash run.sh install")
    .env({"PYTHONHASHSEED": "0", "CFLAGS": "-g0",
          "LD_LIBRARY_PATH": "/usr/local/lib/python3.11/site-packages/nvidia/cu13/lib"})
)

boltzgen_kit_image = (
    _base
    .run_commands(_clone("boltzgen"))
    .run_commands(_install_lock("boltzgen"))
    .run_commands(
        "cd /kit/boltzgen && ls stock/*.whl >/dev/null 2>&1 && "
        "python -m pip install --no-cache-dir --no-deps stock/*.whl || true"
    )
    .run_commands("python -m pip install --no-cache-dir --no-deps ipsae==1.0.1")
    .run_commands("cd /kit/boltzgen && bash run.sh install")
    .env({"PYTHONHASHSEED": "0", "CFLAGS": "-g0"})
)


# BindCraft2 — not a kit; its own repo, packaged the same way for the same reason (all
# GPU work is on Modal). LICENSE: source-available, hosting-restricted, explicitly not OSI
# open source; RUO-only here, never served to third parties. See docker/bindcraft2-cloud/.
BINDCRAFT2_COMMIT = "301efdd1937fb963cc40a5b0ecc1bc2f9b2b2d15"
bindcraft2_weights = modal.Volume.from_name("raptorbind-bindcraft2-weights", create_if_missing=True)
# A campaign outlives the local client, so its result is PERSISTED rather than returned.
# Modal tears down an ephemeral app when the client disconnects ("Stopping app - local
# client disconnected"), which killed the first attempt mid-run; launch with
# `modal run --detach` and read the outcome back from here.
bc2_results = modal.Volume.from_name("raptorbind-bc2-results", create_if_missing=True)

bindcraft2_image = (
    # Plain Ubuntu, not a CUDA base: BindCraft2's jax cuda13 wheels carry their own CUDA
    # and cuDNN. This follows their containers/Dockerfile rather than our house base image.
    modal.Image.from_registry("ubuntu:24.04", add_python=None)
    .apt_install("python3", "python3-pip", "python3-dev", "build-essential",
                 "ca-certificates", "git")
    .env({"PIP_BREAK_SYSTEM_PACKAGES": "1", "PIP_NO_CACHE_DIR": "1", "PYTHONUNBUFFERED": "1"})
    # ubuntu:24.04 ships `python3` with no `python` on PATH, and Modal's image
    # introspection cannot determine the interpreter version without one ("We were unable
    # to determine the version of Python installed in the Image"). A symlink is the right
    # fix here: add_python= would install a SECOND interpreter, and the editable
    # BindCraft2 install below would land in the other one.
    .run_commands("ln -sf /usr/bin/python3 /usr/local/bin/python")
    .run_commands(
        "git clone https://github.com/PacesaLab/BindCraft2.git /opt/bindcraft && "
        f"git -C /opt/bindcraft checkout {BINDCRAFT2_COMMIT}")
    # Editable on purpose: settings/ and scaffolds/ live at the repo root and the runtime
    # resolves them relative to it; a plain install orphans both trees.
    .run_commands("python3 -m pip install -e '/opt/bindcraft[cuda13]'")
    .run_commands("python3 -m pip install 'gemmi==0.6.5' 'ipsae==1.0.1'")
    # The cuda13 wheels hide their libraries under site-packages/nvidia/*/lib where the
    # loader does not look. Without this, jax warns once and runs on CPU ~100x slower --
    # silently. The grep asserts it worked.
    .run_commands(
        "python3 -c \"import nvidia, pathlib; print('\\n'.join(sorted(str(p) "
        "for r in nvidia.__path__ for p in pathlib.Path(r).glob('*/lib'))))\" "
        "> /etc/ld.so.conf.d/bindcraft-cuda.conf && ldconfig && ldconfig -p | grep -q libcupti")
    .run_commands("cd / && python3 -m bindcraft.selfcheck cuda13 --shipped-only")
    .env({"NVIDIA_VISIBLE_DEVICES": "all", "NVIDIA_DRIVER_CAPABILITIES": "compute,utility"})
)


BC2_GPU = os.environ.get("RB_BC2_GPU", "L4")   # BindCraft2 does not need an H100


def _bc2_shard(gpu: str):
    """One campaign shard per container: N containers x 1 accepted design each.

    Fan-out beats one big container for this workload. BindCraft2's own parallelism is
    workers-per-GPU on a single host, which caps at one machine's memory and makes the
    whole campaign one failure domain; independent shards scale to Modal's container
    limit, let a cheap GPU be the unit, and lose only one design if a shard dies. Each
    shard gets its own `campaign_seed` -- without that they would draw identical
    trajectories and return the same design N times.
    """
    return app.function(
        image=bindcraft2_test_image, gpu=gpu, timeout=60 * 90,
        volumes={"/weights": bindcraft2_weights, "/results": bc2_results},
        env={"BINDCRAFT_WEIGHTS": "/weights", "BINDCRAFT_WORKERS_PER_GPU": "1"},
        max_containers=32, retries=1,
    )


# The real BindCraft2 campaign test (#043 WI-9). Ships our own adapter into the container
# so the run exercises src/generators/bindcraft2.py end to end -- building the campaign JSON
# and, more importantly, PARSING REAL 3_Ranked/!_Ranked.csv output, which it has only ever
# seen as synthetic fixtures.
_ADAPTER = os.environ.get(
    "RB_ADAPTER_DIR",
    "/private/tmp/claude-501/-Users-mike-Code-Raptorbind/"
    "28fe5a35-1395-4e57-9734-0c12caf013f3/scratchpad/adapterpkg")
bindcraft2_test_image = (
    bindcraft2_image
    .add_local_dir(_ADAPTER, "/adapter")
)


@app.function(image=bindcraft2_test_image, gpu=GPU, timeout=60 * 150,
              volumes={"/weights": bindcraft2_weights, "/results": bc2_results},
              env={"BINDCRAFT_WEIGHTS": "/weights"})
def bindcraft2_campaign(modality: str = "binder", n_designs: int = 1,
                        overrides: str = "", hotspots: str = "A65,A67,A79,A82,A89",
                        target_pdb: str = "target.pdb", campaign_name: str = "PLA2G2A_P02_c84",
                        binder_lengths: str = "", tag: str = "") -> dict:
    """Run one real BindCraft2 campaign and parse it with OUR adapter.

    Defaults to PLA2G2A P02_c84 (the epitope the 2026-08-31 deep run used). Pass
    target_pdb/hotspots/campaign_name to point it elsewhere -- e.g. eGFP 2Y0G at the
    E2_minimizer epitope, the only eGFP epitope that produced yield last time.

    binder_lengths is "min,max" and overrides the modality default (binder.json: 60-180).
    """
    import subprocess, sys, threading, time
    from pathlib import Path
    sys.path.insert(0, "/adapter")

    # Stub the pipeline's runner: generate() would shell out through it, but in here we
    # already ARE the container, so the campaign is launched directly below.
    import types
    stub = types.ModuleType("runner")
    stub.run_in_environment = lambda **kw: None
    sys.modules.setdefault("runner", stub)
    from generators.bindcraft2 import BindCraft2Generator

    gen = BindCraft2Generator()
    work = Path("/tmp/bc2run"); work.mkdir(parents=True, exist_ok=True)
    project = work / "project"

    # THE REAL CODE PATH: our adapter builds BindCraft2's campaign JSON.
    gen_cfg = {"modality": modality, "number_of_final_designs": n_designs,
               "campaign_name": campaign_name, "target_chain": "A"}
    if binder_lengths:
        gen_cfg["binder_lengths"] = [int(x) for x in binder_lengths.split(",")]
    campaign = gen._build_campaign_config(
        Path("/adapter") / target_pdb, hotspots.split(","), project, gen_cfg)
    campaign_path = work / "campaign.json"
    campaign_path.write_text(json.dumps(campaign, indent=2))
    print("campaign.json our adapter produced:\n" + campaign_path.read_text())

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
            time.sleep(1.0)

    cmd = ["bindcraft", "design", str(campaign_path)]
    for kv in filter(None, overrides.split(";")):
        cmd += ["--set", kv]
    print("running:", " ".join(cmd), flush=True)

    tag = tag or f"{modality}_{n_designs}"
    res_dir = Path("/results") / tag
    res_dir.mkdir(parents=True, exist_ok=True)
    (res_dir / "campaign.json").write_text(campaign_path.read_text())
    (res_dir / "status.json").write_text(json.dumps({"state": "running"}))
    bc2_results.commit()

    th = threading.Thread(target=sample, daemon=True); th.start()
    t0 = time.time()
    # Stream to a log file AND the container stdout, so a detached run is followable
    # with `modal app logs` instead of being a black box for an hour.
    log_path = res_dir / "campaign.log"
    with log_path.open("w") as log:
        proc_p = subprocess.Popen(cmd, cwd=str(work), stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT, text=True, bufsize=1)
        lines = []
        for line in proc_p.stdout:
            lines.append(line)
            log.write(line); log.flush()
            print(line.rstrip(), flush=True)
            if len(lines) % 40 == 0:
                bc2_results.commit()
        proc_p.wait()
    elapsed = time.time() - t0
    stop.set(); th.join(timeout=3)

    class _P:
        returncode = proc_p.returncode
        stdout = "".join(lines)
        stderr = ""
    proc = _P()
    print(f"bindcraft exited {proc.returncode} after {elapsed/60:.1f} min, "
          f"peak GPU {peak['mib']} MiB", flush=True)

    tree = sorted(str(f.relative_to(work)) for f in work.rglob("*")
                  if f.is_file())[:60]
    # THE OTHER REAL CODE PATH: our adapter parses whatever BindCraft2 actually wrote.
    normalized = gen.normalize_outputs(work, hotspots.split(","))
    records = []
    cj = Path(normalized["candidates_jsonl"])
    if cj.exists():
        records = [json.loads(ln) for ln in cj.read_text().splitlines() if ln.strip()]

    result = {
        "rc": proc.returncode, "minutes": round(elapsed / 60, 1),
        "peak_mib": peak["mib"], "modality": modality,
        "accepted_count": normalized["accepted_count"],
        "records": records[:10],
        "tree": tree,
        "stdout": proc.stdout[-6000:], "stderr": proc.stderr[-6000:],
    }
    (res_dir / "status.json").write_text(json.dumps({"state": "done", **{
        k: v for k, v in result.items() if k not in ("stdout", "stderr")}}, indent=2))
    (res_dir / "result.json").write_text(json.dumps(result, indent=2))
    # Keep whatever BindCraft2 actually wrote, so the adapter can be re-tested offline
    # against real output instead of re-running an hour-long campaign.
    if project.exists():
        subprocess.run(["cp", "-r", str(project), str(res_dir / "project")], check=False)
    bc2_results.commit()
    return result


@app.function(image=modal.Image.debian_slim(), timeout=60 * 5,
              volumes={"/results": bc2_results})
def bc2_fetch(tag: str = "binder_1", tail: int = 60) -> dict:
    """Read a detached campaign's state back out of the results volume."""
    from pathlib import Path
    d = Path("/results") / tag
    bc2_results.reload()
    if not d.exists():
        return {"state": "absent", "have": sorted(p.name for p in Path("/results").iterdir())}
    status = json.loads((d / "status.json").read_text()) if (d / "status.json").exists() else {}
    log = (d / "campaign.log").read_text() if (d / "campaign.log").exists() else ""
    files = sorted(str(f.relative_to(d)) for f in d.rglob("*") if f.is_file())
    print(f"state={status.get('state')} files={len(files)} loglines={len(log.splitlines())}")
    print("\n".join(log.splitlines()[-tail:]))
    return {"status": status, "files": files[:80], "log_tail": "\n".join(log.splitlines()[-tail:])}


def _run(cmd: list[str], cwd: str) -> dict:
    """Run a kit command, capture rc/stdout/stderr without raising."""
    import subprocess
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    return {"rc": proc.returncode, "stdout": proc.stdout[-6000:], "stderr": proc.stderr[-6000:]}


@app.function(image=boltz2_kit_image, timeout=60 * 60,
              volumes={"/weights": boltz2_weights},
              env={"BOLTZ_CACHE": "/weights"})
def fetch_weights() -> dict:
    """Fill the Boltz-2 cache Volume once. No GPU needed — this is a download."""
    from pathlib import Path
    if (Path("/weights") / "mols").exists():
        return {"status": "already present",
                "entries": sorted(p.name for p in Path("/weights").iterdir())}
    out = _run(["bash", "run.sh", "install", "--weights", "/weights"], "/kit/boltz2")
    boltz2_weights.commit()
    out["entries"] = sorted(p.name for p in Path("/weights").iterdir())
    return out


@app.function(image=boltz2_kit_image, gpu=GPU, timeout=60 * 30,
              volumes={"/weights": boltz2_weights},
              # Without BOLTZ_CACHE every mode reports "NOT ACTIVE: BOLTZ_CACHE is not
              # set" (rc 3) for a reason that has nothing to do with the GPU.
              env={"BOLTZ_CACHE": "/weights"})
def kit_check(mode: str, card: str) -> dict:
    """`run.sh check` — a dry run that resolves and gates a mode on this GPU.

    This is #043 WI-3, the check that gates every other kit work item: it says whether the
    kit's optimizations actually engage on the GPU we can rent.
    """
    import subprocess
    gpu_name = subprocess.run(
        ["nvidia-smi", "--query-gpu=name,compute_cap,memory.total", "--format=csv,noheader"],
        capture_output=True, text=True).stdout.strip()
    out = _run(["bash", "run.sh", "check", "--config", card, "--mode", mode], "/kit/boltz2")
    out.update({"mode": mode, "card": card, "gpu": gpu_name})
    return out


@app.function(image=boltz2_kit_image, gpu=GPU, timeout=60 * 60,
              volumes={"/weights": boltz2_weights},
              env={"BOLTZ_CACHE": "/weights"})
def kit_warm(mode: str, card: str) -> dict:
    """`run.sh warm` — compile the mode's Triton/NVRTC kernels.

    Without this, the first timed run in a mode measures kernel COMPILATION, not
    inference, and a kit mode can look slower than stock. (It did, on the first pass:
    `fast` came in at 66s against stock's 57s, entirely from a cold cache.)
    """
    import time
    t0 = time.time()
    out = _run(["bash", "run.sh", "warm", "--config", card, "--mode", mode, "--out", "/tmp/warm"],
               "/kit/boltz2")
    out.update({"mode": mode, "seconds": round(time.time() - t0, 1)})
    return out


@app.function(image=boltz2_kit_image, gpu=GPU, timeout=60 * 60,
              volumes={"/weights": boltz2_weights},
              env={"BOLTZ_CACHE": "/weights"})
def kit_predict(mode: str, card: str, target_seq: str, binder_seq: str,
                design_id: str, seed: int = 0) -> dict:
    """Fold one target+binder complex in `mode` and return Boltz-2's own confidences.

    The point is the A/B: `exact` claims identical outputs and `fast` claims small numeric
    differences. Since the refolded structure is this project's single source of truth,
    that claim gets measured on real metrics rather than trusted.
    """
    import subprocess, time
    from pathlib import Path

    # Fetch the Boltz-2 cache into the mounted Volume once; later runs reuse it.
    if not (Path("/weights") / "mols").exists():
        print("[kit] fetching Boltz-2 weights into the volume (first run only)")
        w = _run(["bash", "run.sh", "install", "--weights", "/weights"], "/kit/boltz2")
        print(w["stdout"][-2000:], w["stderr"][-2000:])
        boltz2_weights.commit()

    work = Path(f"/tmp/{design_id}_{mode}")
    work.mkdir(parents=True, exist_ok=True)
    # Same shape src/predictors/boltz2.py writes: target = chain A, binder = chain B,
    # both single-sequence (a de novo binder has no homologs by construction).
    (work / "complex.yaml").write_text(
        "version: 1\nsequences:\n"
        f"  - protein:\n      id: A\n      sequence: {target_seq}\n      msa: empty\n"
        f"  - protein:\n      id: B\n      sequence: {binder_seq}\n      msa: empty\n"
    )
    out_dir = work / "out"
    env = dict(os.environ, BOLTZ_CACHE="/weights")
    cmd = ["bash", "run.sh", "pred", "--config", card, "--mode", mode,
           "--input", str(work / "complex.yaml"), "--out_dir", str(out_dir),
           "--seed", str(seed), "--recycling_steps", "3", "--sampling_steps", "200",
           "--diffusion_samples", "1", "--output_format", "mmcif"]
    t0 = time.time()
    proc = subprocess.run(cmd, cwd="/kit/boltz2", capture_output=True, text=True, env=env)
    elapsed = time.time() - t0

    # Kit modes nest stock's file set under by_seed/<name>/s<seed>/; `off` does not.
    # This is the same dual layout src/predictors/boltz2.py:_find_boltz2_outputs walks.
    confidences = sorted(out_dir.rglob("confidence_*.json"))
    metrics = json.loads(confidences[0].read_text()) if confidences else {}
    return {
        "mode": mode, "rc": proc.returncode, "seconds": round(elapsed, 1),
        "confidence_path": str(confidences[0].relative_to(out_dir)) if confidences else None,
        "iptm": metrics.get("iptm"), "ptm": metrics.get("ptm"),
        "complex_plddt": metrics.get("complex_plddt"),
        "confidence_score": metrics.get("confidence_score"),
        "stdout": proc.stdout[-4000:], "stderr": proc.stderr[-4000:],
    }


@app.function(image=boltzgen_kit_image, gpu=GPU, timeout=60 * 30,
              volumes={"/weights": boltzgen_weights},
              env={"BOLTZGEN_CACHE": "/weights"})
def boltzgen_kit_check(mode: str = "exact", card: str = "h100") -> dict:
    """`run.sh check` for the BoltzGen kit, and the installed upstream version.

    Unlike Boltz-2 this is also a VERSION BUMP: the kit pins boltzgen 0.3.2 against our
    0.2.0, so the version it reports here is part of the answer.
    """
    import subprocess
    nvsmi = subprocess.run(
        ["nvidia-smi", "--query-gpu=name,compute_cap,memory.total", "--format=csv,noheader"],
        capture_output=True, text=True).stdout.strip()
    # boltzgen exposes no __version__ attribute; the installed distribution does.
    ver = subprocess.run(
        ["python", "-c", "from importlib.metadata import version; print(version('boltzgen'))"],
        capture_output=True, text=True)
    out = _run(["bash", "run.sh", "check", "--config", card, "--mode", mode], "/kit/boltzgen")
    print(f"gpu             : {nvsmi}")
    print(f"boltzgen version: {(ver.stdout or ver.stderr).strip()}")
    print(f"check rc        : {out['rc']}")
    for line in (out["stdout"] or out["stderr"]).strip().splitlines()[-8:]:
        print(f"  {line}")
    out.update({"gpu": nvsmi, "version": (ver.stdout or ver.stderr).strip(), "mode": mode})
    return out


@app.function(image=bindcraft2_image, gpu=GPU, timeout=60 * 30,
              volumes={"/weights": bindcraft2_weights},
              env={"BINDCRAFT_WEIGHTS": "/weights"})
def bindcraft2_check() -> dict:
    """Prove BindCraft2 is installed AND that jax reached the GPU.

    The failure this guards against is silent: if the cuda13 wheels' libraries are not on
    the loader path, jax warns once and runs on the CPU ~100x slower without erroring.
    """
    import subprocess
    import jax
    devices = [str(d) for d in jax.devices()]
    backend = jax.default_backend()
    listed = subprocess.run(["bindcraft", "design", "--list-targets"],
                            capture_output=True, text=True)
    from pathlib import Path
    gpu_ok = backend == "gpu" and any("cuda" in d.lower() or "gpu" in d.lower()
                                      for d in devices)
    modalities = sorted(m.stem for m in Path("/opt/bindcraft/settings/modality").glob("*.json"))
    nvsmi = subprocess.run(
        ["nvidia-smi", "--query-gpu=name,compute_cap,memory.total", "--format=csv,noheader"],
        capture_output=True, text=True).stdout.strip()
    result = {
        "gpu": nvsmi, "jax_version": jax.__version__, "jax_backend": backend,
        "jax_devices": devices, "gpu_ok": gpu_ok, "modalities": modalities,
        "list_targets_rc": listed.returncode,
        "list_targets": (listed.stdout or listed.stderr)[-1200:],
    }
    print(f"gpu            : {nvsmi}")
    print(f"jax            : {jax.__version__}  backend={backend}  devices={devices}")
    print(f"GPU REACHED    : {gpu_ok}   <-- False means jax silently fell back to CPU")
    print(f"modalities ({len(modalities)}) : {', '.join(modalities)}")
    print(f"list-targets rc: {listed.returncode}")
    print((listed.stdout or listed.stderr)[-800:])
    return result


@app.function(image=boltz2_kit_image, gpu=GPU, timeout=60 * 60,
              volumes={"/weights": boltz2_weights},
              env={"BOLTZ_CACHE": "/weights"})
def kit_memory(mode: str, card: str, binder_seq: str, design_id: str,
               n_gpu: int = 1, repeats: int = 1) -> dict:
    """Measure PEAK GPU memory for a mode by sampling nvidia-smi during the run.

    `big` is the low-memory mode and `--n_gpu P` shards its pair representation across P
    GPUs -- both exist for inputs that will not otherwise fit. This measures whether either
    does anything at OUR complex sizes, and how much headroom we actually have.

    `repeats` tiles the binder to synthesise a larger complex, so the memory curve can be
    sampled well past our real ~170-230 token range without inventing a fake target.
    """
    import subprocess, threading, time
    from pathlib import Path

    peak = {"mib": 0}
    stop = threading.Event()

    def sample():
        while not stop.is_set():
            out = subprocess.run(
                ["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"],
                capture_output=True, text=True).stdout.split()
            for value in out:
                try:
                    peak["mib"] = max(peak["mib"], int(value))
                except ValueError:
                    pass
            time.sleep(0.25)

    binder = binder_seq * repeats
    tokens = len(PLA2G2A) + len(binder)
    work = Path(f"/tmp/mem_{mode}_{n_gpu}_{repeats}")
    work.mkdir(parents=True, exist_ok=True)
    (work / "complex.yaml").write_text(
        "version: 1\nsequences:\n"
        f"  - protein:\n      id: A\n      sequence: {PLA2G2A}\n      msa: empty\n"
        f"  - protein:\n      id: B\n      sequence: {binder}\n      msa: empty\n")

    cmd = ["bash", "run.sh", "pred", "--config", card, "--mode", mode,
           "--input", str(work / "complex.yaml"), "--out_dir", str(work / "out"),
           "--seed", "0", "--recycling_steps", "3", "--sampling_steps", "200",
           "--diffusion_samples", "1", "--output_format", "mmcif"]
    if n_gpu > 1:
        cmd += ["--n_gpu", str(n_gpu)]

    thread = threading.Thread(target=sample, daemon=True)
    thread.start()
    t0 = time.time()
    proc = subprocess.run(cmd, cwd="/kit/boltz2", capture_output=True, text=True)
    elapsed = time.time() - t0
    stop.set()
    thread.join(timeout=2)

    blob = proc.stdout + proc.stderr
    shard_line = next((ln for ln in blob.splitlines() if "sharding=" in ln), "")
    total = subprocess.run(
        ["nvidia-smi", "--query-gpu=memory.total,count", "--format=csv,noheader,nounits"],
        capture_output=True, text=True).stdout.strip()
    result = {"mode": mode, "n_gpu": n_gpu, "tokens": tokens, "rc": proc.returncode,
              "seconds": round(elapsed, 1), "peak_mib": peak["mib"], "gpu_total": total,
              "sharding": shard_line[shard_line.find("sharding="):][:24] if shard_line else "",
              "tail": blob[-2000:]}
    print(f"[{mode} n_gpu={n_gpu} tokens={tokens}] rc={proc.returncode} "
          f"{elapsed:.1f}s peak={peak['mib']}MiB {result['sharding']}")
    return result


@app.function(image=boltz2_kit_image, gpu="H100:2", timeout=60 * 60,
              volumes={"/weights": boltz2_weights},
              env={"BOLTZ_CACHE": "/weights"})
def kit_shard_check(repeats: int = 20) -> dict:
    """Does `--n_gpu 2` actually engage sharding, and is it refused outside `big`?

    Two GPUs on one host. Confirms the escape hatch works if a target ever outgrows one
    card -- not something our ~172-token complexes need, but worth knowing is there.
    """
    import subprocess
    from pathlib import Path
    binder = "SAEKEGNHVRVSGSYSAEEAAEILRKLTGASDAEVEREGDDLVATAKT" * repeats
    work = Path("/tmp/shard"); work.mkdir(parents=True, exist_ok=True)
    (work / "complex.yaml").write_text(
        "version: 1\nsequences:\n"
        f"  - protein:\n      id: A\n      sequence: {PLA2G2A}\n      msa: empty\n"
        f"  - protein:\n      id: B\n      sequence: {binder}\n      msa: empty\n")
    print(subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total",
                          "--format=csv,noheader"], capture_output=True, text=True).stdout)

    out = {}
    for label, mode, n_gpu in [("big n_gpu=1", "big", 1), ("big n_gpu=2", "big", 2),
                               ("exact n_gpu=2 (expect refusal)", "exact", 2)]:
        cmd = ["bash", "run.sh", "check", "--config", "h100", "--mode", mode]
        if n_gpu > 1:
            cmd += ["--n_gpu", str(n_gpu)]
        r = subprocess.run(cmd, cwd="/kit/boltz2", capture_output=True, text=True)
        blob = r.stdout + r.stderr
        shard = next((ln for ln in blob.splitlines() if "sharding=" in ln), "")
        frag = shard[shard.find("n_gpu="):][:40] if shard else blob.strip().splitlines()[-1][:180]
        print(f"  {label:<32} rc={r.returncode}  {frag}")
        out[label] = {"rc": r.returncode, "frag": frag, "tail": blob[-700:]}
    return out


# A batch of designs from the same run, for the throughput test below.
BATCH = [
    ("bg_000", "MKGIARAIRDYEATEPHELSFKKGDLLEVLGTTGGDWVLARNLETGEVGLVPAEHLELIE"),
    ("bg_001", "NDVSGVVNGKPAVSVNGEVYVYDEEKGEWVKDDSVTPEEVEANKVA"),
    ("bg_002", "AIFKAKTKEEAIKKALEAGLEVKEVFEIHGKYYALAE"),
    ("bg_003", "MLEVGDIVVALRDYEAKTPDELSFKKGDRLLVLEVIDEEWVRAASLNTGETGKIPKEYVEIVGRS"),
    ("bg_004", "SARSVQIFVKTLTGKTITLEVETSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRLD"),
    ("bg_005", "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRPL"),
    ("bg_006", "AVAVSSKGVVATGEDELAAIEEVKKIAGEDAGATTFS"),
    ("bg_007", "GLTVEEAAELLSRNFGVSLETARRAVELAIELLEQGKSEEEVAEEISRRLGISREAAEAAVLLARMRL"),
    ("bg_081", "SAEKEGNHVRVSGSYSAEEAAEILRKLTGASDAEVEREGDDLVATAKT"),
    ("bg_488", "SFEELKKKYKELTGADVSDLTEEEIKELLAKLV"),
    ("mpnn_494", "GELYRALVDYTATRDDELSFKEGELLELLEEVDDVLWHARSLVTGKEGLVPRFYVERV"),
    ("bg_222", "MVGDLVRALKDYEAQTPDELSFKKGDILKVLEVIDDDWIKAQSLSTGEVGEVPVSYVEVIEE"),
]


@app.function(image=boltz2_kit_image, gpu=GPU, timeout=60 * 90,
              volumes={"/weights": boltz2_weights},
              env={"BOLTZ_CACHE": "/weights"})
def kit_batch(mode: str, card: str, n: int, reps: int = 1) -> dict:
    """Fold N complexes in ONE call, the way the kit is meant to be used.

    The per-design test was the kit's worst case: `route=worker` means the kit runs a
    PERSISTENT worker, and `--input` accepts a DIRECTORY, so model load and graph capture
    are paid once and amortized across the batch. One design per container measures mostly
    startup — which is why `exact` barely beat stock and `fast` lost to it.

    A campaign refolds hundreds of designs, so the batch number is the one that matters.
    Returns total wall-clock and per-design metrics; `reps` re-runs the whole batch inside
    the same container so the first (cold) pass can be separated from steady state.
    """
    import subprocess, time
    from pathlib import Path

    work = Path(f"/tmp/batch_{mode}")
    inp = work / "inputs"
    inp.mkdir(parents=True, exist_ok=True)
    designs = BATCH[:n]
    for did, seq in designs:
        (inp / f"{did}.yaml").write_text(
            "version: 1\nsequences:\n"
            f"  - protein:\n      id: A\n      sequence: {PLA2G2A}\n      msa: empty\n"
            f"  - protein:\n      id: B\n      sequence: {seq}\n      msa: empty\n"
        )

    passes = []
    for rep in range(reps):
        out_dir = work / f"out{rep}"
        cmd = ["bash", "run.sh", "pred", "--config", card, "--mode", mode,
               "--input", str(inp), "--out_dir", str(out_dir), "--seed", "0",
               "--recycling_steps", "3", "--sampling_steps", "200",
               "--diffusion_samples", "1", "--output_format", "mmcif"]
        t0 = time.time()
        proc = subprocess.run(cmd, cwd="/kit/boltz2", capture_output=True, text=True)
        elapsed = time.time() - t0
        metrics = {}
        for conf in out_dir.rglob("confidence_*.json"):
            # .../<design>/confidence_<design>_model_0.json -- recover the design name
            stem = conf.stem.replace("confidence_", "").rsplit("_model", 1)[0]
            metrics[stem] = json.loads(conf.read_text()).get("iptm")
        passes.append({"rep": rep, "rc": proc.returncode, "seconds": round(elapsed, 1),
                       "n_out": len(metrics), "iptm": metrics,
                       "tail": (proc.stderr or proc.stdout)[-2500:]})
        print(f"[{mode}] rep{rep}: rc={proc.returncode} {elapsed:.1f}s "
              f"for {len(metrics)}/{len(designs)} designs "
              f"({elapsed/max(len(designs),1):.1f}s per design)")
    return {"mode": mode, "n": len(designs), "passes": passes}


# Test case: PLA2G2A + design bg_081 from
# runs/pla_PLA2G2A_deep_P02_c84_20260831_211845/refold/bg_081/.
# 124 + 48 = 172 tokens, outside every A100 `exact` refusal band.
PLA2G2A = ("NLVNFHRMIKLTTGKEAALSYGFYGCHCGVGGRGSPKDATDRCCVTHDCCYKRLEKRGCGTKFLSYKF"
           "SNSGSRITCAKQDSCRSQLCECDKAAATCFARNKTTYNKKYQYYSNKHCRGSTPRC")
# Four designs from that run, each with the iptm the 2026-08-31 local run recorded.
# Chosen to span the interface-quality range (0.33 -> 0.73) rather than one lucky case.
DESIGNS = [
    ("bg_081",   "SAEKEGNHVRVSGSYSAEEAAEILRKLTGASDAEVEREGDDLVATAKT", 0.327487),
    ("bg_488",   "SFEELKKKYKELTGADVSDLTEEEIKELLAKLV", 0.712719),
    ("mpnn_494", "GELYRALVDYTATRDDELSFKEGELLELLEEVDDVLWHARSLVTGKEGLVPRFYVERV", 0.734295),
    ("bg_222",   "MVGDLVRALKDYEAQTPDELSFKKGDILKVLEVIDDDWIKAQSLSTGEVGEVPVSYVEVIEE", 0.382839),
]


@app.local_entrypoint()
def main(only_check: bool = False, gpu: str = "", modes: str = "exact,fast"):
    card = CARD.get(gpu or GPU, "h100")
    print("=== fetching the Boltz-2 cache into the volume (once) ===")
    w = fetch_weights.remote()
    print(f"    {w.get('status', 'fetched')}: {w.get('entries')}\n")
    print(f"=== #043 WI-3: does the kit engage on Modal's {gpu or GPU} (card={card})? ===\n")
    for result in kit_check.starmap([(m, card) for m in ["exact", "fast", "big"]]):
        tail = (result["stdout"] or result["stderr"]).strip().splitlines()
        print(f"--- mode={result['mode']}  rc={result['rc']}  gpu={result['gpu']}")
        for line in tail[-6:]:
            print(f"    {line}")
        print()
    if only_check:
        return

    print("=== warming kernel caches (so timings measure inference, not compilation) ===")
    for w in kit_warm.starmap([(m, card) for m in modes.split(",") if m != "off"]):
        print(f"    mode={w['mode']:<6} rc={w['rc']}  {w['seconds']}s")

    print("\n=== #043 WI-5: stock-vs-kit A/B, 4 designs x modes, H100, seed 0 ===")
    jobs = [(m, card, PLA2G2A, seq, did, 0)
            for did, seq, _ in DESIGNS for m in modes.split(",")]
    rows = list(kit_predict.starmap(jobs))
    by_design: dict = {}
    for job, r in zip(jobs, rows):
        by_design.setdefault(job[4], {})[r["mode"]] = r

    mode_list = modes.split(",")
    print(f"\n{'design':<10}{'local_iptm':>12}" + "".join(f"{m:>12}" for m in mode_list)
          + f"{'exact-off':>12}{'fast-off':>12}")
    drifts = []
    for did, _, local_iptm in DESIGNS:
        got = by_design.get(did, {})
        cells = "".join(
            f"{got[m]['iptm']:>12.6f}" if got.get(m, {}).get("iptm") is not None else f"{'ERR':>12}"
            for m in mode_list)
        off_v = got.get("off", {}).get("iptm")
        ex_v = got.get("exact", {}).get("iptm")
        fa_v = got.get("fast", {}).get("iptm")
        d_ex = f"{ex_v - off_v:>+12.6f}" if (off_v is not None and ex_v is not None) else f"{'-':>12}"
        d_fa = f"{fa_v - off_v:>+12.6f}" if (off_v is not None and fa_v is not None) else f"{'-':>12}"
        if off_v is not None and fa_v is not None:
            drifts.append(abs(fa_v - off_v))
        print(f"{did:<10}{local_iptm:>12.6f}{cells}{d_ex}{d_fa}")

    print(f"\n{'design':<10}" + "".join(f"{m + '_sec':>12}" for m in mode_list))
    for did, _, _ in DESIGNS:
        got = by_design.get(did, {})
        print(f"{did:<10}" + "".join(f"{got[m]['seconds']:>12.1f}" if m in got else f"{'-':>12}"
                                     for m in mode_list))

    if drifts:
        print(f"\nfast-vs-off |d_iptm|: max={max(drifts):.6f} mean={sum(drifts)/len(drifts):.6f}")
    print("\nipae/ipSAE gate at 0.45 ipae; an iptm shift of this size moves designs across "
          "the screening gates, so `fast` is not a drop-in for a source-of-truth refold.")
