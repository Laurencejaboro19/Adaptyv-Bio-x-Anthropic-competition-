# Block 1: write rfd4_app.py to Downloads\tfr1_rfd4, build the cloud image and download the RFD4-Proteina weights (one paste, creates its own folder)
$dir = "$env:USERPROFILE\Downloads\tfr1_rfd4"
New-Item -ItemType Directory -Force -Path $dir | Out-Null
$code = @'
"""RFD4-Proteina on Modal.  Entrypoints: setup, smoke, design, fetch."""
import io, os, subprocess, tarfile, textwrap, time
import modal

COMMIT = "aaccdd2ad0720c180468df5f9cf037a713d403e5"       # RosettaCommons/RFD4-Proteina HEAD I read
HF_REV = "a06e2509c7033b010fc43614e9fd840303588502"       # nvidia/RFD4-Proteina pinned weights
ACCEL = os.environ.get("RFD4_ACCEL", "cu12")              # cu12 (safe) or cu13
GPU = os.environ.get("RFD4_GPU", "L40S")                  # tested by the authors: A100 H100 L40S
EXTRA, BACKEND = {"cu13": ("all-cuda13", "cu130"), "cu12": ("all-cuda12", "cu126")}[ACCEL]
PIP_FIXED = ("pip @ https://files.pythonhosted.org/packages/3f/06/3e4420302363fbb71be3c7f35ce042f54bfca330e3e66b3a2e3c4cca579f/"
             "pip-26.2.2-py3-none-any.whl#sha256=fd7f4dfa8190e014013b9f08f3e973f4e4c26ffd179a204895533a7254813a70")

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("git", "ca-certificates", "build-essential", "libxrender1", "libx11-6", "libxext6")
    .run_commands(f'python -m pip install --no-cache-dir "{PIP_FIXED}"', 'python -m pip install --no-cache-dir "uv>=0.13.4"')
    .run_commands(
        f"git clone https://github.com/RosettaCommons/RFD4-Proteina /opt/rfd4 && cd /opt/rfd4 && git checkout {COMMIT}",
        "printf 'certifi>=2026.10.19\\n' > /tmp/c.txt",
        f"cd /opt/rfd4 && uv pip install --system --no-cache --no-sources -c /tmp/c.txt "
        f"--only-binary tmol,tmol-ops-cu12,tmol-ops-cu13 --extra {EXTRA} --torch-backend {BACKEND} -r pyproject.toml",
        "cd /opt/rfd4 && uv pip install --system --no-cache --no-sources --no-deps . && python -m pip check",
        "uv pip install --system --no-cache gemmi",
    )
    .env({"DATA_PATH": "/data", "HF_HOME": "/cache/hf", "PYTHONNOUSERSITE": "1", "TMOL_DISABLE_WHEEL_FETCH": "1",
          "TMOL_JIT_FALLBACK": "0", "TMOL_USE_JIT": "0", "TORCHINDUCTOR_CACHE_DIR": "/cache/inductor"})
)
vol = modal.Volume.from_name("rfd4-cache", create_if_missing=True)
app = modal.App("rfd4-proteina", image=image)
CK = [f"hf_ckpt_path=nvidia/RFD4-Proteina", f"hf_ckpt_revision={HF_REV}"]


def _run(cmd):
    print("$", " ".join(cmd), flush=True)
    r = subprocess.run(cmd, cwd="/opt/rfd4", text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(r.stdout[-6000:], flush=True)
    return r.returncode


@app.function(gpu=GPU, volumes={"/cache": vol}, timeout=3600)
def setup():
    os.makedirs("/data", exist_ok=True)
    import torch
    print("torch", torch.__version__, "cuda ok:", torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else "")
    from huggingface_hub import snapshot_download
    p = snapshot_download("nvidia/RFD4-Proteina", revision=HF_REV)
    print("weights at", p)
    for root, _, fs in os.walk(p):
        for f in fs:
            fp = os.path.join(root, f); print(f"  {os.path.getsize(fp)/1e9:6.2f} GB  {f}")
    rc = _run(["rfd4-generate", "--help"])
    vol.commit()
    return {"help_rc": rc}


@app.function(gpu=GPU, volumes={"/cache": vol}, timeout=3 * 3600)
def smoke(n: int = 2):
    spec = "/opt/rfd4/docs/how-to/input_prep/examples/condition_spec_ppi.yaml"   # the authors' PD-L1 example
    out = f"/cache/runs/smoke_{int(time.time())}"
    rc = _run(["rfd4-generate", f"inputs={spec}", f"num_replicates={n}", "batch_size=1", f"out_dir={out}",
               "~metrics", "compile_lengths=null", *CK])
    vol.commit()
    files = [os.path.join(r, f) for r, _, fs in os.walk(out) for f in fs] if os.path.isdir(out) else []
    print("OUTPUT FILES:", len(files)); [print(" ", f) for f in files[:20]]
    return {"rc": rc, "out": out, "n_files": len(files)}


@app.function(gpu=GPU, volumes={"/cache": vol}, timeout=4 * 3600)
def design(tag: str, pdb_text: str, chain: str, hot: str, contig_binder: str, target_range: str, n: int, batch: int):
    import gemmi
    d = f"/cache/jobs/{tag}"; os.makedirs(d, exist_ok=True)
    open(f"{d}/t.pdb", "w").write(pdb_text)
    s = gemmi.read_structure(f"{d}/t.pdb"); s.setup_entities(); s.assign_label_seq_id(force=True)
    doc = s.make_mmcif_document(); b = doc.sole_block()
    t = b.find('_atom_site.', ['label_asym_id', 'auth_asym_id', 'label_seq_id', 'auth_seq_id'])
    for r in t:  # keep chain id and residue numbers as in the PDB; RFD4 matches contigs against them
        r['_atom_site.label_asym_id'] = r['_atom_site.auth_asym_id']
        r['_atom_site.label_seq_id'] = r['_atom_site.auth_seq_id']
    doc.write_file(f"{d}/target.cif")
    hot_sel = ",".join(x.strip() for x in hot.split(","))
    yml = textwrap.dedent(f'''\
        {tag}:
          input: target.cif
          contig: "{contig_binder}, /0, {chain}{target_range}"
          conditions:
            C_CRDGRP: {{1: [{{select: "sel('{chain}/*/*/[N,CA,C,O]') and not has_nan_coord()"}}]}}
            C_SEQ: {{true: [{{select: "sel('{chain}')"}}]}}
            C_HOT: {{true: [{{select: "sel('{chain}/*/[{hot_sel}]')"}}]}}
        ''')
    open(f"{d}/spec.yaml", "w").write(yml); print(yml)
    out = f"/cache/runs/{tag}"
    rc = _run(["rfd4-generate", f"inputs={d}/spec.yaml", f"num_replicates={n}", f"batch_size={batch}",
               f"out_dir={out}", "~metrics", "compile_lengths=null", *CK])
    vol.commit()
    files = [os.path.join(r, f) for r, _, fs in os.walk(out) for f in fs] if os.path.isdir(out) else []
    print("OUTPUT FILES:", len(files))
    return {"rc": rc, "out": out, "n_files": len(files)}


@app.function(volumes={"/cache": vol}, timeout=1800)
def pack(path: str) -> bytes:
    vol.reload()
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as t:
        t.add(path, arcname=os.path.basename(path))
    return buf.getvalue()


@app.local_entrypoint()
def run_setup():
    print(setup.remote())


@app.local_entrypoint()
def run_smoke(n: int = 2):
    print(smoke.remote(n))


@app.local_entrypoint()
def run_design(tag: str, pdb: str, hot: str, contig_binder: str = "60-100", target_range: str = "191-383",
               chain: str = "A", n: int = 20, batch: int = 1):
    print(design.remote(tag, open(pdb).read(), chain, hot, contig_binder, target_range, n, batch))


@app.local_entrypoint()
def fetch(tag: str, outdir: str = "."):
    data = pack.remote(f"/cache/runs/{tag}")
    os.makedirs(outdir, exist_ok=True)
    p = os.path.join(outdir, f"{tag}.tar.gz"); open(p, "wb").write(data)
    print("saved", p, len(data), "bytes")
'@
[System.IO.File]::WriteAllText("$dir\rfd4_app.py", $code, (New-Object System.Text.UTF8Encoding($false)))
Set-Location $dir
$env:PYTHONUTF8="1"; $env:PYTHONIOENCODING="utf-8"
python -m modal run rfd4_app.py::run_setup
