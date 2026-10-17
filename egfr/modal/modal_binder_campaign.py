"""
MODULE: scripts/modal_binder_smoke.py
Modal GPU smoke test for the Raptorbind binder-design pipeline: generate de-novo
binders toward a target sequence with BoltzGen, then refold + score each with
Boltz-2 -- all on Modal's on-demand GPUs, fanning the per-design refolds out
across many containers.

This exercises the whole cloud-GPU path (image build from the repo's existing
docker/ recipes, model-weight caching in Modal Volumes, GPU generation +
prediction, parallel fan-out) end-to-end, so it's the go/no-go test for using
Modal as the compute backend.

Flow (local_entrypoint `main`):
  0. fold_target  -- Boltz-2 folds the target sequence -> target.pdb (chain A),
                     which also warms the Boltz-2 weight cache Volume.
  1. generate     -- BoltzGen designs N binders (length range) toward target.pdb;
                     returns each binder's sequence (parsed from the output CIFs).
  2. refold_score -- Boltz-2 refolds each target+binder complex and extracts
                     iptm / pLDDT / ipae / ipSAE / pDockQ; fanned via .map().
  3. rank + save the results locally.

Run (from repo root, with the persistent Modal venv).
ALWAYS pass --detach for anything longer than a smoke test: an ephemeral Modal app is
tethered to the launching client's heartbeat, so a slept laptop, closed lid or dropped
network stops the app and CANCELS every in-flight input -- GPU work discarded. Verify with
`modal app list` (want "ephemeral (detached)", not "stopped"). Note --detach keeps the
remote work alive but does NOT save results, because the local entrypoint doing the
collection still dies -- see docs/operations/binder_campaign_playbook.md "Mitigations,
weakest to strongest", and use --stage refold_resume to recover.

  ~/modal-venv/bin/modal run scripts/modal_binder_smoke.py --n-designs 4   # quick validate
  ~/modal-venv/bin/modal run --detach scripts/modal_binder_smoke.py        # full: 100 designs
  ~/modal-venv/bin/modal run --detach scripts/modal_binder_smoke.py --stage fold   # target fold only
  ~/modal-venv/bin/modal run --detach scripts/modal_binder_campaign.py --gen-mode structure \
      --hotspots "A68,A72,A75"   # constrain binders to contact target residues 68/72/75
                                  # (hotspots require --gen-mode structure; ignored in cofold)

Single-sequence mode (msa: empty) throughout -- a 6-residue target and de-novo
binders have no meaningful MSA, and it avoids the external MSA-server dependency.
"""
import json
import os

import modal

app = modal.App("raptorbind-binder-campaign")

_CU128 = "https://download.pytorch.org/whl/cu128"

# Images built natively by Modal (no Dockerfile / no local Docker daemon) --
# same package recipe as docker/{boltz2,boltzgen}-cloud/Dockerfile, but via
# from_registry(CUDA base) + pip_install instead of replaying the miniforge/mamba
# bootstrap. Modal layer-caches each step, so rebuilds are fast.
boltz2_image = (
    modal.Image.from_registry(
        "nvidia/cuda:12.8.0-cudnn-runtime-ubuntu22.04", add_python="3.11")
    .apt_install("git")
    .pip_install("numpy==1.26.4")
    .pip_install("torch==2.10.0", index_url=_CU128)
    .pip_install(
        "triton==3.6.0",
        "cuequivariance==0.10.0", "cuequivariance-ops-cu12==0.10.0",
        "cuequivariance-ops-torch-cu12==0.10.0", "cuequivariance-torch==0.10.0",
        "boltz==2.2.1", "gemmi==0.6.5", "ipsae==1.0.1")
)
boltzgen_image = (
    modal.Image.from_registry(
        "nvidia/cuda:12.8.0-cudnn-runtime-ubuntu22.04", add_python="3.12")
    .apt_install("git", "build-essential")
    .pip_install("torch==2.9.0", index_url=_CU128)
    .pip_install(
        "triton>=3.4.0",
        "cuequivariance==0.8.1", "cuequivariance-ops-cu12==0.8.1",
        "cuequivariance-ops-torch-cu12==0.8.1", "cuequivariance-torch==0.8.1",
        "gemmi==0.6.5", "ipsae==1.0.1")
    .pip_install("git+https://github.com/HannesStark/boltzgen.git")
)
# soluble ProteinMPNN (dauparas/ProteinMPNN, MIT) — soluble weights are bundled
# in the repo clone, so this stage has NO HuggingFace dependency. Redesigns the
# binder chain on a fixed target backbone.
mpnn_image = (
    modal.Image.from_registry(
        "nvidia/cuda:12.8.0-cudnn-runtime-ubuntu22.04", add_python="3.10")
    .apt_install("git")
    .pip_install("torch==2.8.0", index_url=_CU128)
    .pip_install("biopython==1.86", "gemmi==0.7.4", "ipsae==1.0.1", "numpy")
    .run_commands(
        # Pinned, not floating (#043): this is upstream HEAD as of 2026-09-24 (a
        # 2023-06-27 commit) and exactly what the uplifting proteinmpnn kit pins, so
        # our image and the kit run identical code. An unpinned clone meant a rebuild
        # could silently change the inverse-folding model between campaigns.
        "git clone https://github.com/dauparas/ProteinMPNN.git /opt/ProteinMPNN && "
        "git -C /opt/ProteinMPNN checkout 8907e6671bfbfc92303b5f79c4b5e6ce47cdef57")
)

# Persistent weight caches so a fresh container doesn't re-download (~6GB for Boltz-2).
boltz_cache = modal.Volume.from_name("raptorbind-boltz-cache", create_if_missing=True)
gen_cache = modal.Volume.from_name("raptorbind-boltzgen-cache", create_if_missing=True)
# Per-target MSA cache, shared across every container in a campaign. Without it
# each refold container queries the ColabFold server independently -- ~1,600
# identical queries for one mCherry campaign -- which is both slow and the reason
# refold fan-out had to stay throttled. Keyed by a hash of the target sequence,
# so it is correct by construction: a different target cannot collide with it.
msa_cache = modal.Volume.from_name("raptorbind-msa-cache", create_if_missing=True)
MSA_CACHE_DIR = "/msa"


def _msa_key(seq: str) -> str:
    import hashlib
    return hashlib.sha256(seq.strip().upper().encode()).hexdigest()[:32]

# BoltzGen weights: point at a public mirror. The official `boltzgen/*` HF repos
# were 503-ing (HF-side outage) while this namespace serves the identical
# checkpoints (HTTP 200). BoltzGen's checkpoint flags accept huggingface:<repo>:<file>.
BG_MIRROR = "moritztng/boltzgen"

GPU = "L4"                       # 24GB; cheapest tier that fits these small jobs
BOLTZ_CACHE_DIR = "/cache/boltz"  # passed to `boltz --cache`
# Fan-out width. The account GPU-concurrency limit was raised to 50 on
# 2026-08-28 (it was 10; an earlier attempt at 20 hit the old wall).
#
# GENERATION and REFOLD are capped separately, because they are limited by
# different things. Generation is pure GPU and scales with the account limit.
# Refold with --use-msa additionally queries the ColabFold MSA server once per
# container -- there is still no MSA cache (see the note in refold_score: an
# a3m-in-Volume attempt was abandoned when boltz silently dropped every record
# that used it). At 10 containers that query rate was self-limiting; at 50 it is
# 5x higher against a shared public service, so the MSA path stays lower until
# caching lands. Override either with an env var.
#
# UPDATE 2026-08-28: MSAs are now cached per target in a shared Volume and built
# once by fold_target, so a campaign issues roughly ONE server query per target
# instead of one per design. The MSA path is no longer the bottleneck it was;
# VERIFIED 2026-08-28: a refold consuming the cached CSV produces structures with
# no fallback to the server, so the MSA path is capped at 50 like the others. The
# earlier "boltz silently drops cached MSAs" conclusion was wrong in its cause --
# the cached file had been the ColabFold intermediate .a3m rather than the merged
# msa/*.csv that boltz actually consumes.
GEN_MAX_CONTAINERS = int(os.environ.get("GEN_MAX_CONTAINERS", "50"))
REFOLD_MAX_CONTAINERS = int(os.environ.get("REFOLD_MAX_CONTAINERS", "50"))
REFOLD_MSA_MAX_CONTAINERS = int(os.environ.get("REFOLD_MSA_MAX_CONTAINERS", "50"))
FUNC_TIMEOUT = 60 * 60            # 1h per container


# ---------------------------------------------------------------------------
# helpers that run INSIDE the containers
# ---------------------------------------------------------------------------

def _run(cmd, cwd=None, env=None):
    """Run a subprocess, streaming output; raise with tail on failure."""
    import subprocess
    import os
    full_env = {**os.environ, **(env or {})}
    print(f"$ {' '.join(cmd)}", flush=True)
    p = subprocess.run(cmd, cwd=cwd, env=full_env, capture_output=True, text=True)
    if p.returncode != 0:
        tail = (p.stdout or "")[-2000:] + "\n---STDERR---\n" + (p.stderr or "")[-3000:]
        raise RuntimeError(f"command failed (rc={p.returncode}): {' '.join(cmd)}\n{tail}")
    return p.stdout


_AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}


def _parse_hotspots(hotspots: str, default_chain: str = "A") -> str:
    """Parse a hotspot string like "A68,A72,A75" or "68,72,75" into a
    comma-separated residue-number string for BoltzGen's `binding_types`
    (only tokens matching default_chain are kept; bare numbers assume it).
    Mirrors src/generators/boltzgen.py's `_parse_hotspots` convention so the
    same hotspot notation works in both the local pipeline and this cloud path.
    """
    if not hotspots:
        return ""
    residues = []
    for hs in hotspots.replace(" ", ",").split(","):
        hs = hs.strip()
        if not hs:
            continue
        if hs[0].isalpha() and len(hs) > 1:
            chain, res_num = hs[0].upper(), hs[1:]
            if chain == default_chain.upper():
                residues.append(res_num)
        else:
            residues.append(hs)
    return ",".join(residues)


# Natural proteins BoltzGen has been observed to emit verbatim as "de novo"
# designs. Measured across the eGFP and mCherry campaigns: ~9% of BoltzGen output
# is ubiquitin, uniformly 8.3-8.8% per epitope, and one identical 75-aa sequence
# appeared under BOTH targets -- conclusive memorisation rather than convergent
# design, since a target-specific design cannot be identical for two different
# targets. Add entries here as more are identified.
NATURAL_SEQUENCES = {
    "ubiquitin": ("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRT"
                  "LSDYNIQKESTLHLVLRLRGG"),
}

# Threshold set empirically over 1636 real BoltzGen designs, not guessed. The
# 6-mer containment distribution against ubiquitin is sharply trimodal:
#
#   0.00          1421 designs (86.9%)  genuinely de novo
#   0.01 - 0.6      26 designs ( 1.6%)  ubiquitin FOLD, sequence substantially
#                                       diverged -- e.g. MQIFVKLQDGQVITLEVSPIDT...
#                                       This is what a successful MPNN rescue
#                                       looks like, and it already occurs at low
#                                       rate on its own.
#   >= 0.6         191 designs (11.7%)  verbatim or near-verbatim copies,
#                                       including N-terminally extended ones
#                                       (GPMQIFVKTLTGK..., TPMQIFVKTLTGK...)
#
# 0.6 therefore separates "this IS the natural protein" from "this borrows its
# fold". Only the former needs rescuing. A simple prefix test misses 84 of the
# 191 because of those N-terminal extensions, which is why the earlier ~9%
# estimate understated the contamination.
NATURAL_KMER_THRESHOLD = 0.6
_NATURAL_K = 6


def _natural_containment(seq: str, ref: str) -> float:
    """Fraction of seq's k-mers that also occur in ref. Order-independent, so it
    is unaffected by N-terminal extensions and tolerates point mutations."""
    if not seq or len(seq) < _NATURAL_K:
        return 0.0
    kmers = {seq[i:i + _NATURAL_K] for i in range(len(seq) - _NATURAL_K + 1)}
    if not kmers:
        return 0.0
    ref_kmers = {ref[i:i + _NATURAL_K] for i in range(len(ref) - _NATURAL_K + 1)}
    return len(kmers & ref_kmers) / len(kmers)


def _natural_scan(seq: str):
    """Return (name, containment) for the closest natural sequence -- name is ""
    below threshold. Containment is returned regardless so a rescue can be
    measured: an MPNN variant of a flagged backbone should show it falling."""
    best_name, best = "", 0.0
    for name, ref in NATURAL_SEQUENCES.items():
        c = _natural_containment(seq, ref)
        if c > best:
            best_name, best = name, c
    return (best_name if best >= NATURAL_KMER_THRESHOLD else ""), round(best, 3)


def _natural_match(seq: str) -> str:
    """Name of the natural protein this sequence reproduces, or ""."""
    return _natural_scan(seq)[0]


def _chain_sequences(struct_path):
    """Return {chain_id: one_letter_seq} for a PDB/CIF via gemmi (first model)."""
    import gemmi
    st = gemmi.read_structure(str(struct_path))
    st.setup_entities()
    out = {}
    if len(st) == 0:
        return out
    for chain in st[0]:
        seq = "".join(_AA3TO1.get(r.name, "") for r in chain if r.name in _AA3TO1)
        if seq:
            out[chain.name] = out.get(chain.name, "") + seq
    return out


def _collect_files(root, max_bytes=25_000_000):
    """Collect every file under `root` into {relpath: bytes} for return to the
    local entrypoint (so structures + run files land on the local machine).
    Skips files larger than max_bytes (logged)."""
    from pathlib import Path
    root = Path(root)
    files = {}
    for p in sorted(root.rglob("*")):
        if p.is_file():
            sz = p.stat().st_size
            if sz > max_bytes:
                print(f"  (skipping large file {p.relative_to(root)}: {sz} bytes)")
                continue
            files[str(p.relative_to(root))] = p.read_bytes()
    return files


# ---------------------------------------------------------------------------
# Stage 0: fold the target sequence -> target.pdb (Boltz-2)
# ---------------------------------------------------------------------------

@app.function(image=boltz2_image, gpu=GPU,
              volumes={"/cache": boltz_cache, MSA_CACHE_DIR: msa_cache},
              timeout=FUNC_TIMEOUT, retries=2)
def fold_target(target_seq: str, use_msa: bool = False) -> str:
    """Boltz-2 monomer fold of the target sequence; returns the PDB text (chain A).

    use_msa=True drops `msa: empty` and passes --use_msa_server, so Boltz-2
    builds a real MSA (ColabFold server) instead of folding single-sequence.
    Single-sequence is fine for some targets and catastrophic for others: eGFP
    folds to 0.40 A of its crystal without an MSA, while mCherry -- same fold
    family, same size -- comes out compact but topologically wrong (6.3 A mean
    internal CA-CA distance error). Always verify the target fold against a
    reference structure before designing against it.
    """
    from pathlib import Path
    work = Path("/tmp/foldtarget")
    work.mkdir(parents=True, exist_ok=True)
    yaml_path = work / "target.yaml"
    msa_field = "" if use_msa else ", msa: empty"
    yaml_path.write_text(
        "version: 1\nsequences:\n"
        f"  - protein: {{id: A, sequence: {target_seq}{msa_field}}}\n"
    )
    out = work / "out"
    cmd = ["boltz", "predict", str(yaml_path), "--out_dir", str(out),
           "--model", "boltz2", "--accelerator", "gpu", "--devices", "1",
           "--recycling_steps", "1", "--sampling_steps", "50",
           "--diffusion_samples", "1", "--output_format", "pdb",
           "--cache", BOLTZ_CACHE_DIR, "--no_kernels"]
    if use_msa:
        cmd.append("--use_msa_server")
    _run(cmd)
    boltz_cache.commit()  # persist downloaded weights

    # Capture the MSA this fold just paid for, so every later refold of the same
    # target reuses it instead of re-querying the server. The fold is the natural
    # place to do it: it runs exactly once per target, before any refold.
    if use_msa:
        # Cache the MERGED MSA that boltz itself produced (msa/<record>_<n>.csv),
        # NOT the .a3m files under msa/*_tmp_env/. Those are ColabFold pipeline
        # intermediates (uniref.a3m, bfd.mgnify30...a3m); boltz declines them,
        # which is why the first attempt at MSA caching looked like boltz
        # "silently dropping" cached records. The CSV is the documented input
        # format for the yaml `msa:` field and round-trips.
        cands = [c for c in sorted(out.rglob("msa/*.csv")) if "_tmp_env" not in str(c)]
        if cands:
            dest = Path(MSA_CACHE_DIR) / f"{_msa_key(target_seq)}.csv"
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(cands[0].read_bytes())
            msa_cache.commit()
            print(f"cached MSA {dest.name} ({dest.stat().st_size} bytes) "
                  f"from {cands[0].relative_to(out)}")
        else:
            print("WARNING: --use_msa_server ran but produced no merged msa/*.csv; "
                  "refolds will fall back to querying the server")

    pdbs = sorted(out.rglob("*_model_0.pdb")) or sorted(out.rglob("*.pdb"))
    if not pdbs:
        raise RuntimeError(f"no PDB produced; tree: {[str(p) for p in out.rglob('*')][:40]}")
    text = pdbs[0].read_text()
    print(f"folded target -> {pdbs[0].name} ({len(text)} bytes)")
    return text


# ---------------------------------------------------------------------------
# Stage 1: BoltzGen generation
# ---------------------------------------------------------------------------

@app.function(image=boltzgen_image, gpu=GPU, volumes={"/root/.cache": gen_cache},
              timeout=FUNC_TIMEOUT, retries=3)
def generate(n_designs: int, target_ref: str, gen_mode: str = "cofold",
             length_min: int = 30, length_max: int = 80, omit_aas: str = "C",
             steps: str = "design inverse_folding", hotspots: str = "") -> dict:
    """Run BoltzGen to design n_designs binders toward the target; return
    {designs, files}. Binder sequences are parsed from the generated complex
    structures (the chain that isn't the target).

    gen_mode:
      "cofold"    -- target passed as a SEQUENCE (target_ref = seq); BoltzGen
                     co-folds target+binder jointly (right for a flexible/IDR
                     target that shouldn't be pinned to one predicted pose).
                     hotspots is ignored here -- there's no fixed structure/
                     chain numbering to anchor residue positions to.
      "structure" -- target passed as a fixed PDB (target_ref = pdb text);
                     BoltzGen designs against that rigid conformation.
                     hotspots (e.g. "A68,A72,A75" or "68,72,75") constrains
                     the binder to contact those target residues via
                     BoltzGen's `binding_types` block.

    Steps are just `design inverse_folding` -- BoltzGen's own
    design_folding/analysis/filtering are redundant since stage 2 refolds +
    scores every design with Boltz-2."""
    from pathlib import Path
    work = Path("/tmp/bg")
    work.mkdir(parents=True, exist_ok=True)
    out_dir = work / "bg_out"
    spec = work / "design_spec.yaml"

    if gen_mode == "cofold":
        target_seqs = {target_ref}
        spec.write_text(
            "entities:\n"
            f"  - protein: {{id: B, sequence: {length_min}..{length_max}}}\n"
            f"  - protein: {{id: A, sequence: {target_ref}}}\n"
        )
        print(f"co-fold mode: target sequence '{target_ref}' co-folded with the binder")
        if hotspots:
            print(f"WARNING: hotspots={hotspots!r} ignored in cofold mode "
                  f"(no fixed target structure to anchor residue numbers to); "
                  f"pass --gen-mode structure to apply hotspots.")
    else:  # structure mode: fixed predicted target
        target_pdb = work / "target.pdb"
        target_pdb.write_text(target_ref)
        target_seqs = set(_chain_sequences(target_pdb).values())
        binding_residues = _parse_hotspots(hotspots, default_chain="A")
        lines = [
            "entities:",
            f"  - protein: {{id: B, sequence: {length_min}..{length_max}}}",
            "  - file:",
            f"      path: {target_pdb}",
            "      include:",
            "        - chain: {id: A}",
        ]
        if binding_residues:
            lines += [
                "      binding_types:",
                "        - chain:",
                "            id: A",
                f"            binding: {binding_residues}",
            ]
        lines.append('      structure_groups: "all"')
        spec.write_text("\n".join(lines) + "\n")
        print("structure mode: designing against the fixed target.pdb"
              + (f" (hotspots: {binding_residues})" if binding_residues else ""))
    # Checkpoints from the mirror (official boltzgen/* HF repos were 503-ing).
    # `design inverse_folding` needs only the design + inverse-fold checkpoints.
    # `--` before the positional guards against any greedy nargs='+' option
    # (e.g. --steps / --design_checkpoints) swallowing the design_spec path.
    cmd = ["boltzgen", "run", "--protocol", "protein-anything",
           "--steps"] + steps.split() + [
        "--num_designs", str(n_designs), "--output", str(out_dir),
        "--cache", "/root/.cache/boltzgen",
        # NB: leave --moldir at its default (boltzgen/inference-data, a *dataset*
        # repo that isn't affected by the outage and is already cached). The
        # mirror is a *model* repo, so pointing moldir at it 401s (looks under
        # /datasets/). Only the checkpoints (model repos) come from the mirror.
        # BoltzGen's pipeline constructor eagerly resolves ALL checkpoint kinds
        # up front (even ones these steps won't use), so every default that
        # points at the down official org must be redirected to the mirror.
        "--inverse_fold_checkpoint", f"huggingface:{BG_MIRROR}:boltzgen1_ifold.ckpt",
        "--folding_checkpoint", f"huggingface:{BG_MIRROR}:boltz2_conf_final.ckpt",
        "--affinity_checkpoint", f"huggingface:{BG_MIRROR}:boltz2_aff.ckpt",
        "--design_checkpoints", f"huggingface:{BG_MIRROR}:boltzgen1_diverse.ckpt",
    ]
    if omit_aas:  # forbid these residues in the inverse-folded sequences (e.g. "C")
        cmd += ["--inverse_fold_avoid", omit_aas]
        print(f"forbidding residues in generation (inverse_fold_avoid): {omit_aas}")
    cmd += ["--", str(spec)]
    _run(cmd)
    gen_cache.commit()

    # Log the output tree once, so if extraction misses we can see the layout.
    all_files = sorted(out_dir.rglob("*"))
    print(f"boltzgen output tree ({len(all_files)} entries), sample:")
    for p in all_files[:60]:
        if p.is_file():
            print("   ", p.relative_to(out_dir))

    # Prefer the INVERSE-FOLDED designs (the final assigned binder sequences);
    # the raw intermediate_designs/ are pre-inverse-folding, so taking both
    # double-counts each design. Fall back to raw, then anywhere in the tree.
    structs = []
    for sub in ("intermediate_designs_inverse_folded", "intermediate_designs"):
        d = out_dir / sub
        if d.is_dir():
            structs = sorted(d.glob("*.cif")) + sorted(d.glob("*.mmcif")) + sorted(d.glob("*.pdb"))
            if structs:
                print(f"using designs from {sub}/ ({len(structs)} structures)")
                break
    if not structs:
        for ext in ("*.cif", "*.mmcif", "*.pdb"):
            structs += sorted(out_dir.rglob(ext))

    designs, seen = [], set()
    for s in structs:
        chains = _chain_sequences(s)
        for cid, seq in chains.items():
            if seq in target_seqs or len(seq) < max(10, length_min - 10):
                continue  # skip the target chain / too-short fragments
            if seq in seen:
                continue
            seen.add(seq)
            designs.append({"design_id": f"design_{len(designs):03d}",
                            "binder_seq": seq, "binder_chain": cid,
                            "source_file": str(s.relative_to(out_dir))})
    # Return the full BoltzGen output tree too, so the local entrypoint can
    # persist the generated structures + run files locally.
    gen_files = _collect_files(out_dir)
    # Fallback: parse binder sequences from any FASTA if no structures yielded them.
    if not designs:
        for fa in sorted(out_dir.rglob("*.fasta")) + sorted(out_dir.rglob("*.fa")):
            seq_lines, cur = [], []
            for line in fa.read_text().splitlines():
                if line.startswith(">"):
                    if cur:
                        seq_lines.append("".join(cur)); cur = []
                else:
                    cur.append(line.strip())
            if cur:
                seq_lines.append("".join(cur))
            for seq in seq_lines:
                if seq and seq not in target_seqs and seq not in seen and len(seq) >= 10:
                    seen.add(seq)
                    designs.append({"design_id": f"design_{len(designs):03d}",
                                    "binder_seq": seq, "binder_chain": "B",
                                    "source_file": str(fa.relative_to(out_dir))})
    # Attach each design's complex structure (target A + binder B) so the MPNN
    # stage can redesign chain B on that backbone.
    for d in designs:
        d["generation_mode"] = gen_mode
        try:
            d["complex_cif"] = (out_dir / d["source_file"]).read_text()
        except OSError:
            d["complex_cif"] = None
    print(f"extracted {len(designs)} unique binder sequences (mode={gen_mode})")
    return {"designs": designs, "files": gen_files}


# ---------------------------------------------------------------------------
# Stage 1b: soluble ProteinMPNN — redesign the binder chain (variants)
# ---------------------------------------------------------------------------

@app.function(image=mpnn_image, gpu=GPU, timeout=FUNC_TIMEOUT, retries=2)
def mpnn_variants(designs_batch: list, num_per: int = 1, temperature: float = 0.15,
                  design_chain: str | None = None, omit_aas: str = "C") -> dict:
    """Soluble-MPNN redesign of the binder chain (design_chain) on each design's
    complex backbone, keeping the target chain fixed. Returns
    [{parent_id, binder_seq, mpnn_score, binder_chain}] (num_per best per input).

    Follows the repo's contract: parse_multiple_chains -> assign_fixed_chains
    --chain_list <binder> -> protein_mpnn_run --use_soluble_model. Weights are
    bundled in the clone, so no network fetch."""
    import re
    from pathlib import Path
    import gemmi

    # The binder chain to redesign = whatever `generate` identified (BoltzGen
    # labels the binder chain A and the target chain B in its output CIFs -- the
    # opposite of the design_spec ids). Derive it from the batch rather than
    # hardcoding, so we never redesign the target by mistake.
    if design_chain is None:
        chains = [d.get("binder_chain") for d in designs_batch if d.get("binder_chain")]
        design_chain = max(set(chains), key=chains.count) if chains else "A"
    print(f"mpnn: redesigning chain {design_chain} (target chain fixed)")

    MP = "/opt/ProteinMPNN"
    work = Path("/tmp/mpnn")
    inp = work / "inputs"
    inp.mkdir(parents=True, exist_ok=True)
    parent_by_stem = {}
    for d in designs_batch:
        if not d.get("complex_cif"):
            continue
        stem = d["design_id"]
        cif = work / f"{stem}.cif"
        cif.write_text(d["complex_cif"])
        st = gemmi.read_structure(str(cif))
        st.setup_entities()
        st.write_pdb(str(inp / f"{stem}.pdb"))  # ProteinMPNN needs PDB
        parent_by_stem[stem] = d
    if not parent_by_stem:
        return {"variants": [], "files": {}}

    out = work / "out"
    _run(["python", f"{MP}/helper_scripts/parse_multiple_chains.py",
          "--input_path", str(inp), "--output_path", str(work / "parsed.jsonl")])
    _run(["python", f"{MP}/helper_scripts/assign_fixed_chains.py",
          "--input_path", str(work / "parsed.jsonl"),
          "--output_path", str(work / "assigned.jsonl"), "--chain_list", design_chain])
    mpnn_cmd = ["python", f"{MP}/protein_mpnn_run.py",
                "--jsonl_path", str(work / "parsed.jsonl"),
                "--chain_id_jsonl", str(work / "assigned.jsonl"),
                "--out_folder", str(out), "--sampling_temp", str(temperature),
                "--num_seq_per_target", str(num_per), "--use_soluble_model",
                "--save_score", "1", "--batch_size", "1"]
    if omit_aas:  # forbid these residues in the redesigned sequences (e.g. "C")
        mpnn_cmd += ["--omit_AAs", omit_aas]
        print(f"forbidding residues in MPNN (omit_AAs): {omit_aas}")
    _run(mpnn_cmd)

    variants = []
    for fa in sorted((out / "seqs").glob("*.fa")):
        stem = fa.stem
        if stem not in parent_by_stem:
            continue
        recs, hdr = [], None
        for ln in fa.read_text().splitlines():
            if ln.startswith(">"):
                hdr = ln
            elif ln.strip() and hdr and "sample=" in hdr:  # design records only (native has no sample=)
                m = re.search(r"score=([-0-9.]+)", hdr)
                sc = float(m.group(1)) if m else 999.0
                seq = ln.strip()
                if "/" in seq:
                    seq = max(seq.split("/"), key=len)  # binder is the long chain
                recs.append((sc, seq))
                hdr = None
        recs.sort(key=lambda x: x[0])  # lower ProteinMPNN score = better
        for sc, seq in recs[:num_per]:
            variants.append({"parent_id": stem, "binder_seq": seq,
                             "mpnn_score": round(sc, 4), "binder_chain": design_chain})
    print(f"mpnn: {len(variants)} variants from {len(parent_by_stem)} designs")
    return {"variants": variants, "files": _collect_files(out)}


# ---------------------------------------------------------------------------
# Stage 2: Boltz-2 refold + score one target+binder complex
# ---------------------------------------------------------------------------

@app.function(image=boltz2_image, gpu=GPU,
              volumes={"/cache": boltz_cache, MSA_CACHE_DIR: msa_cache},
              timeout=FUNC_TIMEOUT, max_containers=REFOLD_MAX_CONTAINERS, retries=2)
def refold_score(design: dict, target_seq: str, use_msa: bool = False) -> dict:
    """Refold target(A)+binder(B) with Boltz-2, extract iptm/plddt/ipae/ipsae/pdockq.

    use_msa=True builds a real MSA for the TARGET only; the de-novo binder keeps
    `msa: empty`, having no homologs by construction. Worth the extra runtime
    whenever the target does not fold correctly single-sequence -- the eGFP
    campaign found that only ~50% of single-sequence co-folds reproduced the
    target within 3 A of its own monomer fold, and pLDDT did not flag the
    failures.
    """
    from pathlib import Path
    did = design["design_id"]
    binder_seq = design["binder_seq"]
    work = Path(f"/tmp/refold_{did}")
    work.mkdir(parents=True, exist_ok=True)
    yaml_path = work / "complex.yaml"
    # MSA for the TARGET only; the de-novo binder has no homologs by
    # construction.
    #
    # Prefer the a3m that fold_target already cached for this exact sequence, so
    # a campaign makes ONE MSA query instead of one per design. An earlier
    # attempt at this was abandoned because boltz silently dropped every record
    # that referenced a cached path -- exit 0, no prediction. That failure is now
    # detectable: the artefact check below turns "no structure" into an explicit
    # failure, and `cached_msa` is retried once against the live server before
    # the design is given up on, so a bad cache costs one retry rather than a
    # whole campaign of empty results.
    cached = Path(MSA_CACHE_DIR) / f"{_msa_key(target_seq)}.csv"
    used_cache = use_msa and cached.exists() and cached.stat().st_size > 0
    if not use_msa:
        target_msa = ", msa: empty"
    elif used_cache:
        target_msa = f", msa: {cached}"
    else:
        target_msa = ""
    yaml_path.write_text(
        "version: 1\nsequences:\n"
        f"  - protein: {{id: A, sequence: {target_seq}{target_msa}}}\n"
        f"  - protein: {{id: B, sequence: {binder_seq}, msa: empty}}\n"
    )
    out = work / "out"
    try:
        cmd = ["boltz", "predict", str(yaml_path), "--out_dir", str(out),
               "--model", "boltz2", "--accelerator", "gpu", "--devices", "1",
               "--recycling_steps", "3", "--sampling_steps", "200",
               "--diffusion_samples", "1", "--output_format", "mmcif",
               "--cache", BOLTZ_CACHE_DIR, "--no_kernels"]
        if use_msa and not used_cache:
            cmd.append("--use_msa_server")
        _run(cmd)
        # A cached a3m that boltz declines leaves no structure and exits 0 --
        # the exact silent failure that made the first attempt at MSA caching
        # look like it worked. Detect it here and redo this one design against
        # the live server, so a bad cache costs a retry, not the campaign.
        if used_cache and not (sorted(out.rglob("*_model_0.cif"))
                               or sorted(out.rglob("*.cif"))):
            print(f"{did}: cached MSA yielded no structure -- retrying via MSA server")
            yaml_path.write_text(
                "version: 1\nsequences:\n"
                f"  - protein: {{id: A, sequence: {target_seq}}}\n"
                f"  - protein: {{id: B, sequence: {binder_seq}, msa: empty}}\n"
            )
            import shutil as _sh
            _sh.rmtree(out, ignore_errors=True)
            _run(cmd + ["--use_msa_server"])
            used_cache = False
    except RuntimeError as e:
        return {"design_id": did, "binder_seq": binder_seq, "status": "failed",
                "error": str(e)[:400]}

    metrics = {"design_id": did, "binder_seq": binder_seq,
               "binder_len": len(binder_seq), "status": "ok"}
    # confidence JSON
    conf = sorted(out.rglob("confidence_*_model_0.json")) or sorted(out.rglob("confidence_*.json"))
    if conf:
        c = json.loads(conf[0].read_text())
        metrics["iptm"] = c.get("iptm")
        metrics["ptm"] = c.get("ptm")
        metrics["plddt"] = c.get("complex_plddt")
        ipde = c.get("complex_ipde")
        metrics["ipae"] = round(min(ipde / 30.0, 1.0), 4) if isinstance(ipde, (int, float)) else None
        metrics["confidence_score"] = c.get("confidence_score")
    # ipSAE / pDockQ from the PAE npz + cif
    try:
        import ipsae
        paes = sorted(out.rglob("pae_*_model_0.npz")) or sorted(out.rglob("pae_*.npz"))
        if paes:
            pae = paes[0]
            name = pae.stem[len("pae_"):]
            cif = pae.parent / f"{name}.cif"
            if not cif.exists():
                cifs = sorted(pae.parent.glob("*.cif"))
                cif = cifs[0] if cifs else None
            if cif:
                r = ipsae.calculate_ipsae(str(pae), str(cif), 15.0, 10.0, output_dir=None)
                if isinstance(r, dict):
                    metrics.update(_extract_ipsae_metrics(r))
    except Exception as e:  # noqa: BLE001
        metrics["ipsae_error"] = str(e)[:200]
    # Persist the refolded complex structure + all run files back to the caller.
    metrics["files"] = _collect_files(out)
    # Identify the predicted complex structure file (for structure_path in the DB).
    cifs = sorted(out.rglob("*_model_0.cif")) or sorted(out.rglob("*.cif"))
    metrics["structure_relpath"] = str(cifs[0].relative_to(out)) if cifs else None
    if metrics["structure_relpath"] is None:
        # boltz can exit 0 having silently dropped the record (a malformed MSA
        # reference leaves only processed/manifest.json). Without this the run
        # reports "N/N refolded OK" with every metric None, which reads as
        # success and hides the failure until someone inspects the metrics.
        metrics["status"] = "failed"
        metrics["error"] = ("boltz exited 0 but produced no structure; "
                            f"output tree: {[str(p.relative_to(out)) for p in out.rglob('*')][:20]}")
        print(f"{did}: FAILED -- boltz produced no structure")
        return metrics
    print(f"{did}: iptm={metrics.get('iptm')} plddt={metrics.get('plddt')} "
          f"ipae={metrics.get('ipae')} ipsae={metrics.get('ipsae')} "
          f"(AB={metrics.get('ipsae_AB')} BA={metrics.get('ipsae_BA')}) "
          f"pdockq={metrics.get('pdockq')} pdockq2={metrics.get('pdockq2')} "
          f"files={len(metrics['files'])}")
    return metrics


# ---------------------------------------------------------------------------
# orchestration (runs locally; dispatches to the GPU functions)
# ---------------------------------------------------------------------------

@app.local_entrypoint()
def main(target_seq: str = "GGGGSHHHHHH", n_boltzgen: int = 100, mpnn_per: int = 1,
         length_min: int = 30, length_max: int = 80, gen_mode: str = "cofold",
         omit_aas: str = "__cfg__", hotspots: str = "__cfg__",
         config: str = "configs/binder_campaign.yaml",
         stage: str = "all", run_dir: str = "", run_label: str = "",
         use_msa: bool = False, natural_mpnn_per: int = 8):
    """Campaign: BoltzGen n_boltzgen designs (gen_mode cofold|structure) ->
    soluble-MPNN mpnn_per variant(s) per design -> Boltz-2 refold+score ALL ->
    persist locally + DesignDB. generation_mode is persisted per design.

    omit_aas: residues forbidden in binder sequences (BoltzGen inverse_fold_avoid
    + MPNN omit_AAs). Default sentinel '__cfg__' reads binder_campaign.omit_aas
    from `config` (default "C" -- cysteine forbidden). Pass e.g. --omit-aas "" to
    allow all residues.

    hotspots: target residues the binder should contact, e.g. "A68,A72,A75"
    (chain-prefixed) or "68,72,75" (bare numbers, chain A assumed). Default
    sentinel '__cfg__' reads binder_campaign.hotspots from `config` (default
    "" -- unconstrained); pass --hotspots explicitly to override per-run.
    Only honored in --gen-mode structure (needs a fixed target pose to number
    residues against); ignored with a warning in cofold mode. Threaded into
    BoltzGen's design_spec.yaml as a `binding_types` block -- same notation
    as the local pipeline's hotspot config (src/generators/boltzgen.py).

    stage='refold_resume' + run_dir=<existing>: re-refold only the designs in
    that run's designs.json that don't yet have a structure on disk (recovers a
    run whose tail refolds were cancelled), merging into its results.json + DB.

    natural_mpnn_per: MPNN variants generated per backbone whose sequence
    reproduces a known natural protein (see NATURAL_SEQUENCES), instead of
    mpnn_per. ~9% of BoltzGen output is memorised ubiquitin; those backbones
    still fold and dock, so rather than discarding them this drives their
    SEQUENCE away from the natural one while keeping the pose. Set 0 to give
    them no special treatment.

    use_msa: build a real MSA for the TARGET (ColabFold server) instead of
    folding it single-sequence, in both the target fold and every refold. Verify
    per target before designing: eGFP folds to 0.40 A of its crystal without an
    MSA, mCherry does not fold correctly at all without one.

    run_label: name for the run directory under runs/. Defaults to a truncated,
    hashed form of the target sequence; set it explicitly when several runs
    share a target and differ only in --hotspots, since the directory name is
    otherwise the only thing telling them apart."""
    if stage == "refold_resume":
        _refold_resume(target_seq, run_dir, use_msa=use_msa)
        return

    omit_aas = _resolve_omit_aas(omit_aas, config)
    hotspots = _resolve_hotspots(hotspots, config)
    if hotspots and gen_mode != "structure":
        print(f"WARNING: --hotspots {hotspots!r} requires --gen-mode structure "
              f"(currently {gen_mode!r}); hotspots will be ignored.")
    run_dir, run_label, db_path = _local_run_dir(f"{target_seq}_{gen_mode}", run_label)
    print(f"=== CAMPAIGN target={target_seq} mode={gen_mode} boltzgen={n_boltzgen} "
          f"mpnn_per={mpnn_per} binder_len={length_min}-{length_max} "
          f"forbid_residues={omit_aas or '(none)'} hotspots={hotspots or '(none)'} "
          f"natural_mpnn_per={natural_mpnn_per} ===")
    print(f"local run dir: {run_dir}\nDesignDB: {db_path}")

    # Target reference for generation: co-fold uses the raw sequence (no fold --
    # an IDR has no meaningful single structure to pin); structure mode folds it.
    if gen_mode == "cofold":
        target_ref = target_seq
        (run_dir / "target.fasta").write_text(f">target\n{target_seq}\n")
        print("Stage 0: co-fold mode -> target passed as sequence (no pre-fold)")
    else:
        print("Stage 0: structure mode -> folding target ...")
        target_ref = fold_target.remote(target_seq, use_msa=use_msa)
        (run_dir / "target.pdb").write_text(target_ref)
        print(f"  wrote target.pdb ({len(target_ref)} bytes)")
    if stage == "fold":
        return

    # Stage 1: BoltzGen generation (fanned). Constants via kwargs; vary per-worker count.
    n_workers = min(GEN_MAX_CONTAINERS, max(1, n_boltzgen))
    per = -(-n_boltzgen // n_workers)
    print(f"Stage 1: BoltzGen generating {n_boltzgen} designs ({n_workers} workers x {per}) ...")
    gen_results = list(generate.map(
        [per] * n_workers,
        kwargs={"target_ref": target_ref, "gen_mode": gen_mode,
                "length_min": length_min, "length_max": length_max,
                "omit_aas": omit_aas, "hotspots": hotspots},
        return_exceptions=True))
    bg_designs, seen = [], set()
    for wi, res in enumerate(gen_results):
        if not isinstance(res, dict) or "designs" not in res:
            print(f"  (gen worker {wi} failed: {str(res)[:150]})")
            continue
        _write_files(run_dir / "generate" / f"worker_{wi}", res.get("files", {}))
        for d in res["designs"]:
            if d["binder_seq"] in seen:
                continue
            seen.add(d["binder_seq"])
            d["design_id"] = f"bg_{len(bg_designs):03d}"
            d["generator"] = "boltzgen"
            d["parent"] = None
            bg_designs.append(d)
    print(f"  {len(bg_designs)} unique BoltzGen designs")
    if not bg_designs:
        print("no designs produced; stopping")
        return

    # Flag memorised natural proteins BEFORE the MPNN stage, so the extra
    # sampling happens where it is needed rather than being discovered after the
    # refolds are already paid for.
    #
    # These backbones are not discarded: they fold and they dock (4 of mCherry's
    # 47 trusted designs sat on ubiquitin backbones), so the pose is worth
    # keeping. What is wrong is the SEQUENCE -- it is not de novo, with the
    # novelty/IP consequences that implies. MPNN redesigns sequence on a fixed
    # backbone, which is exactly the right instrument: give these backbones far
    # more variants so the sequence is driven away from the natural one while the
    # pose is retained.
    for d in bg_designs:
        d["natural_match"], d["natural_containment"] = _natural_scan(d["binder_seq"])
    flagged = [d for d in bg_designs if d["natural_match"]]
    if flagged:
        by_name = {}
        for d in flagged:
            by_name[d["natural_match"]] = by_name.get(d["natural_match"], 0) + 1
        print(f"  memorised natural sequences: {len(flagged)}/{len(bg_designs)} "
              f"({len(flagged)/len(bg_designs):.1%}) {by_name} "
              f"-> {natural_mpnn_per} MPNN variants each instead of {mpnn_per}")

    # Stage 1b: soluble-MPNN variants on the BoltzGen backbones (fanned).
    variant_designs = []
    if mpnn_per > 0:
        raw = []
        # Two passes so the flagged backbones can get a different variant count;
        # mpnn_variants applies one num_per per call.
        passes = [(("clean", [d for d in bg_designs if not d["natural_match"]]), mpnn_per)]
        if flagged and natural_mpnn_per > 0:
            passes.append((("natural", flagged), natural_mpnn_per))
        for (tag, group), num_per in passes:
            if not group:
                continue
            mw = min(REFOLD_MAX_CONTAINERS, len(group))
            batches = [b for b in (group[i::mw] for i in range(mw)) if b]
            print(f"Stage 1b[{tag}]: soluble-MPNN ({num_per}/design over "
                  f"{len(batches)} batches, {len(group)} designs) ...")
            mres = list(mpnn_variants.map(
                batches, kwargs={"num_per": num_per, "omit_aas": omit_aas},
                return_exceptions=True))
            for bi, r in enumerate(mres):
                if not isinstance(r, dict) or "variants" not in r:
                    print(f"  (mpnn {tag} batch {bi} failed: {str(r)[:150]})")
                    continue
                _write_files(run_dir / "mpnn" / f"{tag}_batch_{bi}", r.get("files", {}))
                raw += r["variants"]
        for v in raw:
            if v["binder_seq"] in seen:
                continue  # dedup vs BoltzGen designs + each other
            seen.add(v["binder_seq"])
            parent = next((d for d in bg_designs if d["design_id"] == v["parent_id"]), None)
            variant_designs.append({
                "design_id": f"mpnn_{len(variant_designs):03d}",
                "binder_seq": v["binder_seq"], "binder_chain": v.get("binder_chain", "B"),
                "generator": "soluble_mpnn", "parent": v["parent_id"],
                "generation_mode": gen_mode,  # variants inherit the backbone's mode
                # Whether the BACKBONE was a memorised natural protein, and
                # whether this redesigned sequence still resembles it. The
                # variant is the point of the exercise: parent_natural set with
                # natural_match empty is a successful rescue.
                "parent_natural": (parent or {}).get("natural_match", ""),
                "parent_natural_containment": (parent or {}).get("natural_containment", 0.0),
                "natural_match": _natural_scan(v["binder_seq"])[0],
                "natural_containment": _natural_scan(v["binder_seq"])[1],
                "mpnn_score": v.get("mpnn_score")})
        print(f"  {len(variant_designs)} unique MPNN variants")

    all_designs = bg_designs + variant_designs
    print(f"Total: {len(all_designs)} designs ({len(bg_designs)} BoltzGen + {len(variant_designs)} MPNN)")

    # Persist all binder sequences (FASTA + designs.json) before refold.
    (run_dir / "binders.fasta").write_text("".join(
        f">{d['design_id']} gen={d['generator']} len={len(d['binder_seq'])}\n{d['binder_seq']}\n"
        for d in all_designs))
    (run_dir / "designs.json").write_text(json.dumps(
        [{k: v for k, v in d.items() if k != "complex_cif"} for d in all_designs], indent=2))
    print(f"  wrote binders.fasta + designs.json ({len(all_designs)} sequences)")
    if stage == "generate":
        return

    # Stage 2: refold+score ALL (fanned). Slim inputs (drop complex_cif).
    refold_inputs = [{"design_id": d["design_id"], "binder_seq": d["binder_seq"],
                      "binder_chain": d.get("binder_chain", "B")} for d in all_designs]
    meta = {d["design_id"]: {"generator": d["generator"], "parent": d.get("parent"),
                             "generation_mode": d.get("generation_mode", gen_mode),
                             "mpnn_score": d.get("mpnn_score")} for d in all_designs}
    print(f"Stage 2: refolding + scoring {len(refold_inputs)} designs (fan-out {REFOLD_MAX_CONTAINERS}) ...")
    # Persist EACH design's structure + results.json + DB row as soon as it
    # completes (not after the whole .map() finishes) -- a killed/interrupted
    # client only loses in-flight containers, not everything computed so far.
    # `--stage refold_resume --run-dir <this dir>` picks up any still-missing
    # designs on a subsequent invocation.
    db_config = {"pipeline": "boltzgen+soluble_mpnn->boltz2", "generation_mode": gen_mode,
                "natural_mpnn_per": natural_mpnn_per,
                "n_boltzgen": n_boltzgen, "mpnn_per": mpnn_per, "omit_aas": omit_aas,
                "hotspots": hotspots, "length_min": length_min, "length_max": length_max,
                "run_dir": str(run_dir)}
    clean = []
    for i, r in enumerate(refold_score.map(refold_inputs,
                                           kwargs={"target_seq": target_seq, "use_msa": use_msa},
                                           return_exceptions=True), start=1):
        if not isinstance(r, dict):
            print(f"  ({i}/{len(refold_inputs)}) refold raised: {str(r)[:150]}")
            continue
        did = r.get("design_id")
        m = meta.get(did, {})
        r["generator"] = m.get("generator")
        r["parent"] = m.get("parent")
        r["generation_mode"] = m.get("generation_mode")
        if m.get("mpnn_score") is not None:
            r["mpnn_score"] = m["mpnn_score"]
        ddir = run_dir / "refold" / did
        _write_files(ddir, r.pop("files", {}))
        rel = r.get("structure_relpath")
        r["structure_path"] = str(ddir / rel) if rel else None
        ddir.mkdir(parents=True, exist_ok=True)
        (ddir / "binder.fasta").write_text(f">{did}\n{r.get('binder_seq', '')}\n")
        clean.append(r)
        _persist_to_db(db_path, run_label, target_seq, clean, config=db_config)
        _report(clean, all_designs, target_seq, len(all_designs),
                out_name=str(run_dir / "results.json"), quiet=True)
        print(f"  checkpoint: {i}/{len(refold_inputs)} refolded, {len(clean)} OK so far")

    _report(clean, all_designs, target_seq, len(all_designs), out_name=str(run_dir / "results.json"))
    print(f"\nAll structures + run files under: {run_dir}")
    print(f"DesignDB: {db_path} (run_id={run_label})")


# Run directories are named "<label>_<ts>". They used to carry a "binder_smoke_"
# prefix inherited from scripts/modal_binder_smoke.py, which this module was forked
# from -- accurate for that smoke test, wrong for production campaigns. Runs created
# before 2026-08-28 still carry it, so READERS must accept both (see LEGACY_RUN_PREFIX).
RUN_PREFIX = ""
LEGACY_RUN_PREFIX = "binder_smoke_"


def _local_run_dir(seed, run_label: str = ""):
    """Create runs/<label>_<ts>/ and return (dir, run_label, db_path).

    `seed` is a label seed, historically the raw target sequence. That works for
    a 6-residue smoke target but a real one blows past the 255-byte filename
    limit (a 238-residue eGFP target raised OSError 63 before any GPU work
    started), so it is truncated and disambiguated with a hash of the full seed.
    Pass `run_label` to name the run yourself -- worth doing when several runs
    differ only in their hotspots, since the directory name is otherwise the
    only thing distinguishing them.
    """
    import hashlib
    import os
    import pathlib
    import re
    from datetime import datetime
    repo = pathlib.Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    if run_label:
        stem = re.sub(r"[^A-Za-z0-9_.-]", "_", run_label)[:60]
    else:
        digest = hashlib.sha1(str(seed).encode()).hexdigest()[:8]
        stem = f"{str(seed)[:24]}_{digest}"
    label = f"{RUN_PREFIX}{stem}_{ts}"
    run_dir = repo / "runs" / label
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir, label, repo / "runs" / "designs.db"


def _load_binder_campaign_config(config_path) -> dict:
    """Read the `binder_campaign:` section of config_path; {} if missing/unparseable."""
    import os
    from pathlib import Path
    p = Path(config_path)
    if not p.is_absolute():
        p = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) / config_path
    if not p.exists():
        return {}
    try:
        import yaml
        cfg = yaml.safe_load(p.read_text()) or {}
        return cfg.get("binder_campaign", {}) or {}
    except Exception:  # noqa: BLE001 -- yaml missing/parse error: fall back to {}
        return {}


def _resolve_omit_aas(cli_val, config_path):
    """Resolve the forbidden-residues string. '__cfg__' sentinel => read
    binder_campaign.omit_aas from config_path (default 'C'); any other value
    (incl. '') is an explicit CLI override."""
    if cli_val != "__cfg__":
        return cli_val
    return _load_binder_campaign_config(config_path).get("omit_aas", "C")


def _resolve_hotspots(cli_val, config_path):
    """Resolve the hotspot residue string. '__cfg__' sentinel => read
    binder_campaign.hotspots from config_path (default ""); any other value
    (incl. '') is an explicit CLI override."""
    if cli_val != "__cfg__":
        return cli_val
    return _load_binder_campaign_config(config_path).get("hotspots", "")


def _write_files(dest, files: dict):
    """Write {relpath: bytes} under dest."""
    import pathlib
    dest = pathlib.Path(dest)
    for rel, data in files.items():
        p = dest / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(data if isinstance(data, (bytes, bytearray)) else str(data).encode())


def _composite_scores(results, repo):
    """-> ({design_id: composite_score}, {design_ids passing the phase filters}).

    The scoring loss (`scoring.refold` in configs/default.yaml, implemented in
    modules/scoring.py) was previously only ever invoked programmatically at
    TRIAGE, so no refolded design in the DB carried a composite. This computes it
    here so it lands with the rest of the refold metrics.

    EVERY design is scored, not just the ones the filters keep. `rank()` gates
    before ranking, which on these campaigns scores ~9% of designs and leaves the
    rest unsortable; the filters are therefore evaluated separately and reported
    as a flag, so the whole run can be ordered while the gate information
    survives.

    The composite is min-max normalised across the cohort handed to rank(), so it
    orders designs within a run and not across runs. Because this is called at
    every checkpoint with the designs finished so far, each value is rewritten
    (insert_metrics upserts) as the cohort grows, settling on the full-run
    normalisation when the run ends.

    Never fatal: a scoring failure must not cost a run its GPU-earned metrics.
    """
    import os
    import sys
    try:
        import yaml
        sys.path.insert(0, repo)
        from modules.scoring import ScoringStrategy
        cfg_path = os.path.join(repo, "configs", "default.yaml")
        with open(cfg_path) as fh:
            config = yaml.safe_load(fh)
        phase = ScoringStrategy.from_config(config).for_phase("refold")
        cands = [dict(r) for r in results if r.get("design_id")]
        passed = {c["design_id"] for c in phase._apply_filters(cands)}
        phase.filters = {}
        ranked = phase.rank(cands)
        return ({c["design_id"]: c.get("composite_score") for c in ranked
                 if c.get("composite_score") is not None}, passed)
    except Exception as exc:  # noqa: BLE001
        print(f"  composite scoring skipped: {type(exc).__name__}: {exc}")
        return {}, set()


def _persist_to_db(db_path, run_label, target_seq, results, config):
    """Insert the run + every design (sequence + structure_path) + metrics into
    the DesignDB (src/design_db.py), including the composite ranking score."""
    import os
    import sys
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, os.path.join(repo, "src"))
    from design_db import DesignDB
    ok = [r for r in results if r.get("status") == "ok"]
    composites, passed_filters = _composite_scores(ok, repo)
    with DesignDB(db_path) as db:
        db.register_run(run_label, config=config, target=target_seq, status="completed")
        for r in ok:
            did = r["design_id"]
            db_id = db.insert_design(
                run_id=run_label, lot=did, design_id=did,
                generator=r.get("generator") or "boltzgen",
                stage="refold", binder_seq=r.get("binder_seq"),
                parent_lot=r.get("parent"),  # MPNN variant -> its BoltzGen parent
                structure_path=r.get("structure_path"))
            metrics = {k: r.get(k) for k in
                       ("iptm", "ptm", "plddt", "ipae", "ipsae", "ipsae_AB",
                        "ipsae_BA", "pdockq", "pdockq2", "confidence_score",
                        "binder_len", "mpnn_score")}
            if composites:
                metrics["composite_score"] = composites.get(did)
                metrics["scoring_passed_filters"] = 1.0 if did in passed_filters else 0.0
                metrics["composite_cohort_n"] = float(len(composites))
            db.insert_metrics(db_id, "refold", metrics)
    print(f"DesignDB: inserted {len(ok)} designs + metrics into run '{run_label}'"
          f"{f'; {len(composites)} scored, {len(passed_filters)} pass filters' if composites else ''}")


def _refold_resume(target_seq, run_dir_arg, use_msa: bool = False):
    """Re-refold only the designs in an existing run's designs.json that don't
    yet have a structure on disk (recovers a run whose tail refolds were
    cancelled by the concurrency limit), merging into its results.json + DB."""
    import os
    from pathlib import Path
    repo = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    rd = Path(run_dir_arg) if run_dir_arg else None
    if rd and not rd.is_absolute():
        rd = repo / run_dir_arg
    if rd is None:
        cands = sorted((repo / "runs").glob(f"{RUN_PREFIX}{target_seq}_cofold_*"))
        cands += sorted((repo / "runs").glob(f"{LEGACY_RUN_PREFIX}{target_seq}_cofold_*"))
        rd = cands[-1] if cands else None
    if not rd or not rd.is_dir():
        print(f"resume: run dir not found ({run_dir_arg!r})")
        return
    run_label, db_path = rd.name, repo / "runs" / "designs.db"
    designs = json.loads((rd / "designs.json").read_text())
    existing = {}
    if (rd / "results.json").exists():
        for r in json.loads((rd / "results.json").read_text()).get("results", []):
            existing[r["design_id"]] = r

    def has_struct(did):
        ddir = rd / "refold" / did
        return ddir.is_dir() and any(f.endswith(".cif") for _, _, fs in os.walk(ddir) for f in fs)

    missing = [d for d in designs if not has_struct(d["design_id"])]
    print(f"resume {run_label}: {len(designs)} total, {len(designs) - len(missing)} done, "
          f"{len(missing)} to refold (fan-out {REFOLD_MAX_CONTAINERS})")
    if not missing:
        print("nothing to resume")
        return
    meta = {d["design_id"]: d for d in designs}
    refold_inputs = [{"design_id": d["design_id"], "binder_seq": d["binder_seq"],
                      "binder_chain": d.get("binder_chain", "B")} for d in missing]
    db_config = {
        "pipeline": "boltzgen+soluble_mpnn->boltz2",
        "generation_mode": designs[0].get("generation_mode", "cofold"),
        "n_boltzgen": sum(1 for d in designs if d.get("generator") == "boltzgen"),
        "resumed": True, "run_dir": str(rd)}
    # Persist EACH design as soon as it completes (not after the whole .map()
    # finishes) -- see the matching comment in main()'s Stage 2 for why.
    new_ok, exc = [], 0
    for i, r in enumerate(refold_score.map(refold_inputs,
                                           kwargs={"target_seq": target_seq, "use_msa": use_msa},
                                           return_exceptions=True), start=1):
        if not isinstance(r, dict):
            exc += 1
            if exc <= 3:
                print(f"  still-exception: {str(r)[:200]}")
            continue
        did = r["design_id"]
        m = meta.get(did, {})
        r["generator"] = m.get("generator")
        r["parent"] = m.get("parent")
        r["generation_mode"] = m.get("generation_mode")
        if m.get("mpnn_score") is not None:
            r["mpnn_score"] = m["mpnn_score"]
        ddir = rd / "refold" / did
        _write_files(ddir, r.pop("files", {}))
        rel = r.get("structure_relpath")
        r["structure_path"] = str(ddir / rel) if rel else None
        ddir.mkdir(parents=True, exist_ok=True)
        (ddir / "binder.fasta").write_text(f">{did}\n{r.get('binder_seq', '')}\n")
        if r.get("status") == "ok":
            new_ok.append(r)
            existing[did] = r
            _persist_to_db(db_path, run_label, target_seq, list(existing.values()), config=db_config)
            _report(list(existing.values()), designs, target_seq, len(designs),
                    out_name=str(rd / "results.json"), quiet=True)
        print(f"  ({i}/{len(refold_inputs)}) resumed {did}: status={r.get('status')}")
    print(f"resume: {len(new_ok)} newly refolded OK, {exc} still exceptioned")

    merged_list = list(existing.values())
    _report(merged_list, designs, target_seq, len(designs), out_name=str(rd / "results.json"))
    _persist_to_db(db_path, run_label, target_seq, merged_list, config=db_config)
    print(f"resume complete: {len(merged_list)} total OK in results.json + DesignDB ({run_label})")


def _report(clean, designs, target_seq, n_designs, out_name="binder_smoke_results.json", quiet=False):
    import pathlib
    # Drop the raw file blobs from the JSON summary (they're written to disk).
    for r in clean:
        r.pop("files", None)
    ok = [r for r in clean if r.get("status") == "ok"]
    ok.sort(key=lambda r: (r.get("ipsae") or 0, r.get("iptm") or 0), reverse=True)
    out = {"target": target_seq, "n_designs_requested": n_designs,
           "n_generated": len(designs), "n_refolded_ok": len(ok), "results": ok}
    pathlib.Path(out_name).write_text(json.dumps(out, indent=2, default=str))
    if quiet:
        return
    print(f"\n=== RESULTS: {len(ok)}/{len(designs)} refolded OK ===")
    print(f"{'design':>10} {'len':>4} {'iptm':>6} {'plddt':>6} {'ipae':>6} {'ipsae':>6} {'pdockq':>7}")
    for r in ok[:25]:
        print(f"{r['design_id']:>10} {r.get('binder_len', 0):>4} "
              f"{_f(r.get('iptm')):>6} {_f(r.get('plddt')):>6} {_f(r.get('ipae')):>6} "
              f"{_f(r.get('ipsae')):>6} {_f(r.get('pdockq')):>7}")
    print(f"wrote {out_name}")


def _f(x):
    return f"{x:.3f}" if isinstance(x, (int, float)) else "-"


def _extract_ipsae_metrics(r: dict) -> dict:
    """Extract ipsae_AB/ipsae_BA/ipsae/pdockq/pdockq2 from calculate_ipsae()'s
    nested result, chains A(target)/B(binder) -- same convention/extraction
    path as modules/score_ipsae.py so composite_score (modules/scoring.py,
    configs/default.yaml scoring.refold) is computed consistently between this
    cloud path and the local pipeline. `ipsae` is min(ipsae_AB, ipsae_BA) --
    NOT a max-over-the-whole-blob -- matching src/boltz_io.py's convention.
    Falls back to the old permissive max-over-nested-blob extraction only if
    the expected keys are absent (e.g. an ipsae library schema change)."""
    ipsae_scores = r.get("ipsae_scores", {}) or {}
    pdockq_scores = r.get("pdockq_scores", {}) or {}
    ipsae_d0chn_asym = ipsae_scores.get("ipsae_d0chn_asym", {}) or {}

    ipsae_AB = ipsae_d0chn_asym.get("A", {}).get("B") if "A" in ipsae_d0chn_asym else None
    ipsae_BA = ipsae_d0chn_asym.get("B", {}).get("A") if "B" in ipsae_d0chn_asym else None
    if ipsae_AB is not None and ipsae_BA is not None:
        ipsae = min(ipsae_AB, ipsae_BA)
    elif ipsae_AB is not None or ipsae_BA is not None:
        ipsae = ipsae_AB if ipsae_AB is not None else ipsae_BA
    else:
        ipsae = None

    pdockq_dict = pdockq_scores.get("pDockQ", {})
    pdockq2_dict = pdockq_scores.get("pDockQ2", {})
    pdockq = pdockq_dict.get("A", {}).get("B") if isinstance(pdockq_dict, dict) else pdockq_dict
    pdockq2 = pdockq2_dict.get("A", {}).get("B") if isinstance(pdockq2_dict, dict) else pdockq2_dict

    if ipsae is None and pdockq is None:
        # Schema drift fallback -- keep something rather than all-None.
        return {"ipsae": _max_num(ipsae_scores, cap=1.0), "pdockq": _max_num(pdockq_scores, cap=1.0)}

    out = {"ipsae_AB": round(float(ipsae_AB), 4) if ipsae_AB is not None else None,
           "ipsae_BA": round(float(ipsae_BA), 4) if ipsae_BA is not None else None,
           "ipsae": round(float(ipsae), 4) if ipsae is not None else None,
           "pdockq": round(float(pdockq), 4) if isinstance(pdockq, (int, float)) else None,
           "pdockq2": round(float(pdockq2), 4) if isinstance(pdockq2, (int, float)) else None}
    return out


def _max_num(obj, cap=None):
    """Max numeric leaf in a nested dict/list (for ipsae/pdockq per-pair scores).
    `cap` filters out-of-range values (ipSAE/pDockQ are 0-1)."""
    best = None
    stack = [obj]
    while stack:
        x = stack.pop()
        if isinstance(x, bool):
            continue
        if isinstance(x, (int, float)):
            if cap is not None and x > cap:
                continue
            if best is None or x > best:
                best = round(float(x), 4)
        elif isinstance(x, dict):
            stack.extend(x.values())
        elif isinstance(x, (list, tuple)):
            stack.extend(x)
    return best
