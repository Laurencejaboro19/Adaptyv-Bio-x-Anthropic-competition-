"""PDFA KRAS G12D pMHC (Adaptyv Challenge 5): BindCraft2 -> Boltz-2 refold on 3 targets -> ipSAE + contact check, on Modal.
BindCraft2 (H100) binds the G12D complex and AVOIDS the wild-type complex and the peptide-free HLA+b2m (objective "detarget").
Every harvested sequence is refolded with Boltz-2 (L40S) against all three targets and scored with ipSAE.
Only the G12D target should score high. The controller runs INSIDE Modal with its own spend cap (PC can be off).
Usage (from the kras_modal folder):
  python -m modal run kras_app.py::prep --indir <folder with the 3 PDBs>
  python -m modal run kras_app.py::warm_cache
  python -m modal run --detach kras_app.py::smoke --camp set1
  python -m modal run --detach kras_app.py::start --camp set5 --budget 50 --goal 3 --max-shards 2 --per-round 2 --max-traj 30
  python -m modal run kras_app.py::status
  python -m modal run kras_app.py::logs --tag kras_set1
  python -m modal run kras_app.py::fetch --outdir <folder> --tags "kras_set5"
"""
import csv, glob, io, json, os, pathlib, re, time
import modal

# ============================ SETTINGS ======================================
HLA = ("GSHSMRYFYTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRF"
       "LRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHAAEQQRAYLEGRCVEWLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLTWQRDG"
       "EDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWELSSQPT")
B2M = "IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"
PEP = {"G12D": "VVVGADGVGK", "WT": "VVVGAGGVGK"}          # the peptide-free target ("PF") has no peptide chain
KINDS = ("G12D", "WT", "PF")
assert len(HLA) == 281 and len(B2M) == 99

IN = "/inputs"
FILES = {"G12D": "target_G12D_ABC.pdb", "WT": "avoid_wildtype_ABC.pdb", "PF": "avoid_peptide_free_AB.pdb"}
# hotspot sets (chain letter + residue number in the crystal files; peptide = chain C, residue 6 = Asp6)
SETS = {"set1": "A70,A155,A156,C5,C6,C7", "set2": "C6,C7,C8,A70,A73", "set3": "C5,C6,A147,A152,A155,A156",
        "set4": "C4,C5,C6,C7,C8", "set5": "C4,C5,C6,C7,C8"}
SEEDS = {"set1": 7001, "set2": 7101, "set3": 7201, "set4": 7301, "set5": 7401}
COLD = "A183-276,B1-99"
LO, HI = 60, 100                       # binder length (aa)
DETARGET_IPTM = 0.75                   # loose on purpose: wild-type differs by ONE residue, so a tight ceiling would reject everything
DETARGET_INTERFACE_RESIDUES = 15       # same reason (BC2's default of 3 would reject almost every design)
W_WT, W_PF = -0.3, -1.0                # push away from wild-type (weak) and from peptide-free (strong)

def parse_res(s):
    out = set()
    for t in s.split(","):
        c, body = t[0], t[1:]
        if "-" in body:
            a, b = body.split("-"); out.update((c, i) for i in range(int(a), int(b) + 1))
        else:
            out.add((c, int(body)))
    return out
HOT = {k: parse_res(v) for k, v in SETS.items()}
COLDSET = parse_res(COLD)

FORCED = False                         # forced_targeting restored the wild-type surface at harden and collapsed G12D ipTM
# every BC2 stage floor set to 0.3 (the working "low iptm / plddt cut-off" run); merged over campaign() in design_shard
LOOSE = {"min_plddt_screen": 0.3, "min_plddt_refine": 0.3, "min_plddt_anneal": 0.3, "min_plddt_harden": 0.3, "min_plddt_mutate": 0.3, "min_plddt_final": 0.3, "min_iptm_anneal": 0.3, "min_iptm_harden": 0.3, "min_iptm_mutate": 0.3, "min_iptm_final": 0.3, "min_monomer_plddt_final": 0.3, "min_ptm_final": 0.3, "max_ipae_final": 0.7, "min_hotspot_contact_final": 0.3, "max_detarget_interface_residues_final": 30}
WORKERS_PER_GPU = 2
SHARD_WALL_MIN = 50
PER_BACKBONE = 2
MAX_REFOLD_PER_ROUND = 30
N_SAMPLES = 3                          # Boltz diffusion samples per refold (one sample is noisy)
H100_USD_H, L40S_USD_H = 3.95, 1.95
BC2_COMMIT = "301efdd1937fb963cc40a5b0ecc1bc2f9b2b2d15"
# gates (edit here). G12D must be good; WT and peptide-free must be low.
GATE = {"iptm": 0.65, "ipsae": 0.55, "hot_frac": 0.30, "max_cold": 3, "min_binder_contacts": 10,
        "min_pep_res": 3, "need_asp6": True, "asp6_dist": 4.5,
        "strict_wt": 0.30, "strict_pf": 0.20, "relaxed_delta": 0.15}

app = modal.App("pdfa-kras-bc2")
weights = modal.Volume.from_name("pdfa-bc2-weights", create_if_missing=True)
results = modal.Volume.from_name("pdfa-bc2-results", create_if_missing=True)
bcache = modal.Volume.from_name("pdfa-boltz-cache", create_if_missing=True)
inputs = modal.Volume.from_name("pdfa-kras-inputs", create_if_missing=True)

bc2_image = (
    modal.Image.from_registry("ubuntu:24.04", add_python=None)
    .apt_install("python3", "python3-pip", "python3-dev", "build-essential", "ca-certificates", "git")
    .env({"PIP_BREAK_SYSTEM_PACKAGES": "1", "PIP_NO_CACHE_DIR": "1", "PYTHONUNBUFFERED": "1"})
    .run_commands("ln -sf /usr/bin/python3 /usr/local/bin/python")
    .run_commands("git clone https://github.com/PacesaLab/BindCraft2.git /opt/bindcraft && "
                  "git -C /opt/bindcraft checkout " + BC2_COMMIT)
    .run_commands("python3 -m pip install -e '/opt/bindcraft[cuda13]'")
    .run_commands("python3 -c \"import nvidia, pathlib; print('\\n'.join(sorted(str(p) "
                  "for r in nvidia.__path__ for p in pathlib.Path(r).glob('*/lib'))))\" "
                  "> /etc/ld.so.conf.d/bindcraft-cuda.conf && ldconfig && ldconfig -p | grep -q libcupti")
    .run_commands("cd / && python3 -m bindcraft.selfcheck cuda13 --shipped-only")
    .env({"NVIDIA_VISIBLE_DEVICES": "all", "NVIDIA_DRIVER_CAPABILITIES": "compute,utility"}))

boltz_image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "build-essential", "curl")
    .pip_install("boltz==2.2.1", "numpy", "gemmi==0.6.5", "ipsae==1.0.1")
    .env({"PYTHONUNBUFFERED": "1", "BOLTZ_CACHE": "/cache"}))

cpu_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("numpy", "gemmi==0.6.5")
    .env({"PYTHONUNBUFFERED": "1"}))


# ============================ BC2 CONFIG ====================================
def campaign(camp, seed, project, max_traj):
    return {
        "campaign_name": "KRAS_" + camp, "modality": "binder", "project_folder": project, "campaign_seed": seed,
        "targets": [
            {"name": "G12D", "target_path": IN + "/" + FILES["G12D"], "chains": "A,B,C",
             "hotspots": SETS[camp], "coldspots": COLD, "objective": "target", "weight": 1.0},
            {"name": "WT", "target_path": IN + "/" + FILES["WT"], "chains": "A,B,C",
             "hotspots": SETS[camp], "objective": "detarget", "weight": W_WT},
            {"name": "PEPFREE", "target_path": IN + "/" + FILES["PF"], "chains": "A,B",
             "hotspots": ",".join(h for h in SETS[camp].split(",") if h.startswith("A")),
             "objective": "detarget", "weight": W_PF}],
        "binder_lengths": [LO, HI], "number_of_final_designs": 999,
        "max_trajectories": max_traj, "desperation_trajectories": 40,
        "sequence_candidates": 30, "forced_targeting": FORCED, "termini_accessible": False,
        "aa_bias": {"C": 0},
        "min_plddt_screen": 0.4, "min_plddt_refine": 0.5,
        "min_iptm_final": 0.5, "min_monomer_plddt_final": 0.55,
        "max_coldspot_contact_final": 0.1,
        "max_detarget_iptm_screen": DETARGET_IPTM, "max_detarget_iptm_refine": DETARGET_IPTM,
        "max_detarget_iptm_anneal": DETARGET_IPTM, "max_detarget_iptm_harden": DETARGET_IPTM,
        "max_detarget_iptm_mutate": DETARGET_IPTM, "max_detarget_iptm_final": DETARGET_IPTM,
        "max_detarget_interface_residues_final": DETARGET_INTERFACE_RESIDUES,
        "filters": {"Interface_Residues": {"threshold": 10, "higher": True}},
        "workers_per_gpu": WORKERS_PER_GPU}


# ============================ SCORING =======================================
def own_ipsae(pae, ia, ib, cutoff=15.0):
    import numpy as np
    sub = pae[np.ix_(ia, ib)]
    valid = sub < cutoff
    n0 = valid.sum(1)
    L = np.maximum(n0, 27).astype(float)
    d0 = np.maximum(1.24 * np.cbrt(L - 15) - 1.8, 1.0)
    ptm = 1.0 / (1.0 + (sub / d0[:, None]) ** 2)
    per = np.where(n0 > 0, (ptm * valid).sum(1) / np.maximum(n0, 1), 0.0)
    return float(per.max())


def ipsae_by_chain(pae, cif_path, kind, blen):
    """ipSAE (max over both directions) between the binder (chain D) and each target chain."""
    import tempfile, numpy as np
    lens = {"A": len(HLA), "B": len(B2M)}
    if kind != "PF": lens["C"] = len(PEP[kind])
    out, note = {}, ""
    try:
        from ipsae import core
        core.IpsaeCalculator._save_outputs = lambda *a, **k: None
        tmp = os.path.join(tempfile.mkdtemp(), "x_pae.npz"); np.savez(tmp, pae=pae)
        m = core.IpsaeCalculator(15, 10).calculate(tmp, cif_path)["ipsae_scores"]["ipsae_d0chn_asym"]
        for c in lens:
            out[c] = max(float(m[c]["D"]), float(m["D"][c]))
    except Exception as e:
        note = "package failed, used own: " + repr(e)[:150]
        off, idx = 0, {}
        for c in lens: idx[c] = list(range(off, off + lens[c])); off += lens[c]
        ib = list(range(off, off + blen))
        out = {c: max(own_ipsae(pae, idx[c], ib), own_ipsae(pae, ib, idx[c])) for c in lens}
    return out, note


def contact_info(cif_path, camp, kind):
    import gemmi, numpy as np
    m = gemmi.read_structure(cif_path)[0]
    ch = {c.name: c for c in m}
    if "D" not in ch: raise ValueError("no binder chain D, chains=" + ",".join(ch))
    bx = np.vstack([np.array([a.pos.tolist() for a in r]) for r in ch["D"]])
    br = np.concatenate([np.full(len(r), i + 1) for i, r in enumerate(ch["D"])])
    contacted, bres, pairs, pep_pairs = set(), set(), set(), 0
    asp6 = None
    for c in "ABC":
        if c not in ch: continue
        tx = np.vstack([np.array([a.pos.tolist() for a in r]) for r in ch[c]])
        tr = np.concatenate([np.full(len(r), i + 1) for i, r in enumerate(ch[c])])
        near = np.linalg.norm(tx[:, None, :] - bx[None, :, :], axis=-1) <= 4.5
        ti, bj = np.nonzero(near)
        for t, b in zip(tr[ti], br[bj]):
            pairs.add((c, int(t), int(b)))
        contacted |= {(c, int(x)) for x in tr[near.any(1)]}; bres |= {int(x) for x in br[near.any(0)]}
        if c == "C" and len(ch[c]) >= 6:
            r6 = ch[c][5]
            sc = [a for a in r6 if a.name not in ("N", "CA", "C", "O")] or [a for a in r6 if a.name == "CA"]
            sx = np.array([a.pos.tolist() for a in sc])
            asp6 = float(np.linalg.norm(sx[:, None, :] - bx[None, :, :], axis=-1).min())
    pep_pairs = sum(1 for p in pairs if p[0] == "C")
    hot = HOT[camp]
    hit, cold = len(contacted & hot), len(contacted & COLDSET)
    pep_res = len({r for r in contacted if r[0] == "C"})
    frac = hit / len(hot)
    ok = (frac >= GATE["hot_frac"] and cold <= GATE["max_cold"] and len(bres) >= GATE["min_binder_contacts"]
          and pep_res >= GATE["min_pep_res"]
          and (not GATE["need_asp6"] or (asp6 is not None and asp6 <= GATE["asp6_dist"])))
    return {"binder_len": len(ch["D"]), "n_binder_contacts": len(bres), "hotspots_hit": hit, "hotspots_total": len(hot),
            "hotspot_frac": round(frac, 3), "coldspots_hit": cold, "pep_res_contacted": pep_res,
            "pep_pair_frac": round(pep_pairs / max(len(pairs), 1), 3),
            "asp6_min_dist": None if asp6 is None else round(asp6, 2), "contact_PASS": ok,
            "contacted_residues": " ".join(c + str(i) for c, i in sorted(contacted))}


# ============================ BOLTZ =========================================
MSA = {"hla": "/cache/msa_hla.csv", "b2m": "/cache/msa_b2m.csv"}

def binder_iptm(c, n_chains):
    """Binder (last chain) vs target ipTM from pair_chains_iptm; whole-complex iptm is polluted by the HLA-b2m pair."""
    pc = c.get("pair_chains_iptm") or {}
    last = str(n_chains - 1); vals = []
    for a, row in pc.items():
        for b, v in (row or {}).items():
            if (a == last) != (b == last): vals.append(float(v))
    return max(vals) if vals else float(c.get("iptm", 0.0))


def yaml_for(kind, binder):
    s = "version: 1\nsequences:\n"
    s += "  - protein: {id: A, sequence: " + HLA + ", msa: " + MSA["hla"] + "}\n"
    s += "  - protein: {id: B, sequence: " + B2M + ", msa: " + MSA["b2m"] + "}\n"
    if kind != "PF": s += "  - protein: {id: C, sequence: " + PEP[kind] + ", msa: empty}\n"
    s += "  - protein: {id: D, sequence: " + binder + ", msa: empty}\n"
    return s


def _boltz(yaml_text, extra=(), samples=1):
    import subprocess, tempfile
    work = pathlib.Path(tempfile.mkdtemp()); y = work / "c.yaml"; y.write_text(yaml_text)
    out = work / "out"
    cmd = ["boltz", "predict", str(y), "--out_dir", str(out), "--model", "boltz2", "--accelerator", "gpu",
           "--devices", "1", "--cache", "/cache", "--output_format", "mmcif", "--no_kernels",
           "--write_full_pae", "--diffusion_samples", str(samples)] + list(extra)
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r, out


@app.function(image=boltz_image, gpu="L40S", timeout=3600, volumes={"/cache": bcache})
def warm() -> str:
    """Once: download Boltz-2 weights and store the HLA and b2m MSAs (so refolds never hit the MSA server)."""
    import shutil
    msgs = []
    for key, seq in (("hla", HLA), ("b2m", B2M)):
        if os.path.exists(MSA[key]): msgs.append(key + ": already done"); continue
        y = "version: 1\nsequences:\n  - protein: {id: A, sequence: " + seq + "}\n"
        r, out = _boltz(y, ["--use_msa_server", "--recycling_steps", "1", "--sampling_steps", "20"])
        csvs = sorted(glob.glob(str(out / "**" / "msa" / "*.csv"), recursive=True))
        ok = bool(glob.glob(str(out / "**" / "*model_0.cif"), recursive=True)) and bool(csvs)
        if ok:
            shutil.copy(csvs[0], MSA[key]); bcache.commit()
            msgs.append(key + ": MSA saved (%d bytes)" % os.path.getsize(MSA[key]))
        else:
            msgs.append(key + ": FAILED " + (r.stdout + r.stderr)[-600:])
    return "warm " + ("OK" if all("FAILED" not in m for m in msgs) else "FAILED") + " | " + " | ".join(msgs)


@app.function(image=boltz_image, gpu="L40S", timeout=3600, volumes={"/cache": bcache},
              retries=modal.Retries(max_retries=1, backoff_coefficient=1.0, initial_delay=10.0))
def fold_score(binder: str, design_id: str, kind: str, camp: str) -> dict:
    """Refold binder + one target (G12D, WT or PF) with N_SAMPLES samples and score ipSAE per sample."""
    import numpy as np, tempfile
    t0 = time.time()
    for p in MSA.values():
        if not os.path.exists(p): raise RuntimeError("MSA missing - run warm_cache first")
    r, out = _boltz(yaml_for(kind, binder), ["--recycling_steps", "3", "--sampling_steps", "100"], N_SAMPLES)
    rec = {"design_id": design_id, "kind": kind, "sequence": binder}
    samples, note = [], ""
    for i in range(N_SAMPLES):
        cif = sorted(glob.glob(str(out / "**" / ("*model_%d.cif" % i)), recursive=True))
        cnf = sorted(glob.glob(str(out / "**" / ("confidence*model_%d.json" % i)), recursive=True))
        pae = sorted(glob.glob(str(out / "**" / ("pae*model_%d.npz" % i)), recursive=True))
        if not (cif and cnf and pae): continue
        c = json.load(open(cnf[0])); p = np.load(pae[0])["pae"]
        per, n = ipsae_by_chain(p, cif[0], kind, len(binder)); note = note or n
        samples.append({"iptm": binder_iptm(c, 3 if kind == "PF" else 4), "ipsae": max(per.values()),
                        "ipsae_hla": per.get("A", 0.0), "ipsae_b2m": per.get("B", 0.0), "ipsae_pep": per.get("C", 0.0)})
        if i == 0:
            try: rec.update(contact_info(cif[0], camp, kind))
            except Exception as e: rec["contact_error"] = repr(e)[:250]
            if kind == "G12D": rec["cif_text"] = open(cif[0]).read()
    rec["samples"] = samples; rec["note"] = note; rec["seconds"] = round(time.time() - t0, 1)
    if not samples: rec["error"] = (r.stdout + r.stderr)[-800:]
    return rec


def combine(did, seq, recs):
    """recs: {kind: fold_score result}. Builds one row: G12D must be good, WT and PF must be low."""
    row = {"design_id": did, "sequence": seq, "good": False, "tier": ""}
    bad = [k for k in KINDS if k not in recs or not recs[k].get("samples")]
    if bad:
        row["error"] = "no result for " + ",".join(bad) + ": " + str(recs.get(bad[0], {}).get("error", ""))[:200]
        return row
    g = recs["G12D"]; S = g["samples"]; mean = lambda xs: sum(xs) / len(xs)
    row.update({k: v for k, v in g.items() if k in ("binder_len", "n_binder_contacts", "hotspots_hit", "hotspots_total",
        "hotspot_frac", "coldspots_hit", "pep_res_contacted", "pep_pair_frac", "asp6_min_dist", "contact_PASS",
        "contacted_residues", "contact_error")})
    row["g_iptm"] = round(mean([s["iptm"] for s in S]), 4)
    row["g_ipsae"] = round(mean([s["ipsae"] for s in S]), 4)
    row["g_ipsae_max"] = round(max(s["ipsae"] for s in S), 4)
    row["g_ipsae_pep"] = round(mean([s["ipsae_pep"] for s in S]), 4)
    row["g_ipsae_hla"] = round(mean([s["ipsae_hla"] for s in S]), 4)
    row["g_ipsae_b2m"] = round(mean([s["ipsae_b2m"] for s in S]), 4)
    for k, tag in (("WT", "wt"), ("PF", "pf")):
        row[tag + "_iptm_max"] = round(max(s["iptm"] for s in recs[k]["samples"]), 4)
        row[tag + "_ipsae_max"] = round(max(s["ipsae"] for s in recs[k]["samples"]), 4)
    row["delta"] = round(row["g_ipsae"] - max(row["wt_ipsae_max"], row["pf_ipsae_max"]), 4)
    gc, wc = set(g.get("contacted_residues", "").split()), set(recs["WT"].get("contacted_residues", "").split())
    row["wt_site_overlap"] = round(len(gc & wc) / max(len(gc), 1), 3)
    base = bool(row["contact_PASS"] and row["g_iptm"] >= GATE["iptm"] and row["g_ipsae"] >= GATE["ipsae"])
    if base and row["wt_ipsae_max"] <= GATE["strict_wt"] and row["pf_ipsae_max"] <= GATE["strict_pf"]: row["tier"] = "A"
    elif base and row["delta"] >= GATE["relaxed_delta"]: row["tier"] = "B"
    row["good"] = row["tier"] in ("A", "B")
    return row


# ============================ MODAL FUNCTIONS ===============================
@app.function(image=cpu_image, volumes={"/inputs": inputs}, timeout=600)
def check_inputs() -> str:
    import gemmi
    inputs.reload(); msgs = []; ok = True
    want = {"G12D": {"A": 276, "B": 99, "C": 10}, "WT": {"A": 275, "B": 99, "C": 10}, "PF": {"A": 276, "B": 99}}
    for k, f in FILES.items():
        p = IN + "/" + f
        if not os.path.exists(p): msgs.append(k + ": MISSING " + f); ok = False; continue
        st = gemmi.read_structure(p); got = {c.name: len(c) for c in st[0]}
        peptide = gemmi.one_letter_code([r.name for r in st[0]["C"]]) if "C" in {c.name for c in st[0]} else "-"
        good = got == want[k] and (k == "PF" or peptide == PEP[k])
        ok = ok and good
        msgs.append("%s: chains %s peptide %s %s" % (k, got, peptide, "ok" if good else "<-- UNEXPECTED"))
    return ("PREP OK" if ok else "PREP PROBLEM") + "\n" + "\n".join(msgs)


@app.function(image=bc2_image, gpu="H100", timeout=60 * 150,
              volumes={"/weights": weights, "/results": results, "/inputs": inputs})
def design_shard(camp: str, i: int, tag: str, max_traj: int) -> dict:
    import subprocess, shutil, threading
    seed = SEEDS[camp] + i
    work = pathlib.Path("/tmp/bc2_" + tag + "_" + str(seed)); proj = work / "project"
    proj.mkdir(parents=True, exist_ok=True)
    cache = pathlib.Path("/weights/compile_cache"); cache.mkdir(parents=True, exist_ok=True)
    if not (proj / "compile_cache").exists():
        (proj / "compile_cache").symlink_to(cache)
    cpath = work / "campaign.json"; cpath.write_text(json.dumps({**campaign(camp, seed, str(proj), max_traj), **LOOSE}, indent=2))
    out = pathlib.Path("/results") / tag / ("shard" + format(i, "02d") + "_seed" + str(seed))
    out.mkdir(parents=True, exist_ok=True)
    (out / "campaign.json").write_text(cpath.read_text()); results.commit()
    env = dict(os.environ, BINDCRAFT_WORKERS_PER_GPU=str(WORKERS_PER_GPU), BINDCRAFT_WEIGHTS="/weights")
    t0 = time.time()
    with (out / "campaign.log").open("w") as log:
        p = subprocess.Popen(["bindcraft", "design", str(cpath)], cwd=str(work),
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, env=env)
        timer = threading.Timer(SHARD_WALL_MIN * 60, p.terminate); timer.start()
        for line in p.stdout:
            log.write(line); print("[" + format(i, "02d") + "] " + line, end="")
        rc = p.wait(); timer.cancel()
    mins = (time.time() - t0) / 60
    for root, _, files in os.walk(proj):
        if "compile_cache" in root: continue
        for f in files:
            if f.endswith((".csv", ".fasta", ".json", ".pdb")):
                dst = out / "project" / os.path.relpath(root, proj) / f
                dst.parent.mkdir(parents=True, exist_ok=True)
                try: shutil.copy(pathlib.Path(root) / f, dst)
                except Exception: pass
    (out / "status.json").write_text(json.dumps({"shard": i, "seed": seed, "rc": rc, "minutes": round(mins, 1)}))
    results.commit(); weights.commit()
    return {"shard": i, "rc": rc, "minutes": round(mins, 1)}


@app.function(image=cpu_image, volumes={"/results": results}, timeout=900)
def harvest(tag: str, known: list) -> list:
    results.reload()
    seen, found = set(known), []
    pat = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{25,110}$")
    files = glob.glob("/results/" + tag + "/shard*/project/**/*", recursive=True)
    pref = [f for f in files if re.search(r"final|accepted", f, re.I)]
    use = pref if any(f.endswith((".csv", ".fasta")) for f in pref) else files
    def add(did, s):
        if pat.match(s) and s not in seen:
            seen.add(s); found.append({"design_id": did, "sequence": s})
    for f in use:
        parts = pathlib.Path(f).parts
        sh = parts[3] if len(parts) > 3 else "shard"
        if f.endswith(".fasta"):
            name, buf = None, []
            for line in open(f):
                if line.startswith(">"):
                    if name: add(name, "".join(buf))
                    name, buf = sh + "|" + line[1:].split()[0], []
                else: buf.append(line.strip())
            if name: add(name, "".join(buf))
        elif f.endswith(".csv"):
            try: rows = list(csv.DictReader(open(f)))
            except Exception: continue
            for i, r in enumerate(rows):
                s = next((str(r[k]).strip() for k in r if k and k.lower() in
                          ("sequence", "binder_seq", "seq", "binder_sequence", "binder")), "")
                did = next((str(r[k]) for k in r if k and "design" in k.lower()), pathlib.Path(f).stem + "_" + str(i))
                add(sh + "|" + did, s)
    per, keep = {}, []
    for d in found:
        bb = d["design_id"].rsplit("_candidate", 1)[0]
        if per.get(bb, 0) < PER_BACKBONE: per[bb] = per.get(bb, 0) + 1; keep.append(d)
    print("harvested %d new sequences (%d raw) from %d files" % (len(keep), len(found), len(use)))
    return keep


COLS = ["rank", "good", "tier", "design_id", "binder_len", "sequence", "g_iptm", "g_ipsae", "g_ipsae_max",
        "g_ipsae_pep", "g_ipsae_hla", "g_ipsae_b2m", "wt_iptm_max", "wt_ipsae_max", "pf_iptm_max", "pf_ipsae_max",
        "delta", "hotspots_hit", "hotspots_total", "hotspot_frac", "coldspots_hit", "n_binder_contacts",
        "pep_res_contacted", "pep_pair_frac", "asp6_min_dist", "wt_site_overlap", "contact_PASS",
        "contacted_residues", "contact_error", "error"]

def write_outputs(base, scored):
    rows = [dict(v, design_id=d) for d, v in scored.items()]
    rows.sort(key=lambda r: (not r.get("good"), -(r.get("delta") if r.get("delta") is not None else -9)))
    for i, r in enumerate(rows, 1): r["rank"] = i
    b = io.StringIO(); w = csv.DictWriter(b, fieldnames=COLS, extrasaction="ignore"); w.writeheader()
    w.writerows(rows); (base / "designs.csv").write_text(b.getvalue())
    return sum(1 for r in rows if r.get("good"))


@app.function(image=cpu_image, cpu=2.0, timeout=12 * 3600, volumes={"/results": results})
def controller(camp: str, tag: str, goal: int, max_shards: int, per_round: int, budget: float,
               deadline_h: float, max_traj: int) -> str:
    t0 = time.time(); base = pathlib.Path("/results") / tag; base.mkdir(parents=True, exist_ok=True)
    sp = base / "state.json"; results.reload()
    st = json.loads(sp.read_text()) if sp.exists() else {
        "launched": 0, "spent": 0.0, "seqs": {}, "pending": [], "scored": {}, "shard_cost": []}
    rnd = 0

    def n_good(): return sum(1 for s in st["scored"].values() if s.get("good"))

    def save(state_name, why=None):
        sp.write_text(json.dumps(st)); g = write_outputs(base, st["scored"])
        (base / "progress.json").write_text(json.dumps({
            "state": state_name, "camp": camp, "round": rnd, "shards_launched": st["launched"],
            "sequences_seen": len(st["seqs"]), "pending": len(st["pending"]),
            "refolded": len(st["scored"]), "good": g, "goal": goal,
            "est_spend_usd": round(st["spent"], 2), "budget_cap_usd": budget, "stop_reason": why,
            "elapsed_h": round((time.time() - t0) / 3600, 2), "updated": time.strftime("%Y-%m-%d %H:%M:%S")}, indent=1))
        results.commit()

    def process_pending():
        take = st["pending"][:MAX_REFOLD_PER_ROUND]; st["pending"] = st["pending"][MAX_REFOLD_PER_ROUND:]
        if not take: return
        save("refolding")
        jobs = [(st["seqs"][d], d, k, camp) for d in take for k in KINDS]
        outs = list(fold_score.starmap(jobs, return_exceptions=True))
        by = {d: {} for d in take}
        for (seq, d, k, _), o in zip(jobs, outs):
            if isinstance(o, Exception): o = {"design_id": d, "kind": k, "sequence": seq, "error": repr(o)[:300]}
            else: st["spent"] += o.get("seconds", 300) / 3600 * L40S_USD_H
            by[d][k] = o
        for d in take:
            st["scored"][d] = combine(d, st["seqs"][d], by[d])
            if st["scored"][d].get("good") and "cif_text" in by[d]["G12D"]:
                p = base / "refold" / d.replace("|", "_"); p.mkdir(parents=True, exist_ok=True)
                (p / "G12D_complex.cif").write_text(by[d]["G12D"]["cif_text"])

    save("starting")
    if st["pending"]: process_pending(); save("resumed")
    stop = None
    while True:
        rnd += 1
        avg = (sum(st["shard_cost"]) / len(st["shard_cost"])) if st["shard_cost"] else 7.0
        if n_good() >= goal: stop = "goal reached"; break
        if st["launched"] >= max_shards: stop = "shard ceiling reached"; break
        if deadline_h - (time.time() - t0) / 3600 <= 1.0: stop = "deadline reached"; break
        k = min(per_round, max_shards - st["launched"])
        while k > 0 and st["spent"] + k * (avg * 1.15 + 3.0) > budget: k -= 1
        if k <= 0: stop = "budget cap reached (est $%.2f)" % st["spent"]; break
        idx = list(range(st["launched"], st["launched"] + k))
        print("=== round %d: shards %s | good %d/%d | est $%.2f" % (rnd, idx, n_good(), goal, st["spent"]))
        save("designing")
        rows = list(design_shard.starmap([(camp, i, tag, max_traj) for i in idx], return_exceptions=True))
        st["launched"] += k
        for r in rows:
            mins = SHARD_WALL_MIN + 10.0 if isinstance(r, Exception) else r["minutes"]
            c = mins / 60 * H100_USD_H; st["spent"] += c; st["shard_cost"].append(c)
        if all(isinstance(r, Exception) or r.get("rc", 1) != 0 for r in rows):
            print("WARNING: every shard in this round ended with an error/non-zero exit; still harvesting what exists")
        save("harvesting")
        for d in harvest.remote(tag, list(st["seqs"].values())):
            st["seqs"][d["design_id"]] = d["sequence"]; st["pending"].append(d["design_id"])
        if not st["pending"] and all(isinstance(r, Exception) or r.get("rc", 1) != 0 for r in rows):
            stop = "all shards failed and nothing to refold - run: logs"; break
        while st["pending"]:
            process_pending(); save("scored")
            if n_good() >= goal: break
    save("done", stop)
    return "STOPPING: %s | shards %d | refolded %d | good %d | est $%.2f" % (
        stop, st["launched"], len(st["scored"]), n_good(), st["spent"])


@app.function(image=cpu_image, volumes={"/results": results}, timeout=600)
def read_progress(tag: str) -> str:
    results.reload(); p = pathlib.Path("/results") / tag / "progress.json"
    return p.read_text() if p.exists() else "no progress file yet for " + tag

@app.function(image=cpu_image, volumes={"/results": results}, timeout=600)
def read_file(tag: str, name: str) -> bytes:
    results.reload(); p = pathlib.Path("/results") / tag / name
    return p.read_bytes() if p.exists() else b""

@app.function(image=cpu_image, volumes={"/results": results}, timeout=900)
def tail_logs(tag: str, n: int = 25) -> str:
    results.reload(); out = []
    for f in sorted(glob.glob("/results/" + tag + "/shard*/campaign.log")):
        out.append("--- " + f.split("/")[-2]); out.extend(open(f, errors="ignore").read().splitlines()[-n:])
    return "\n".join(out) or "no shard logs yet for " + tag

@app.function(image=cpu_image, volumes={"/results": results}, timeout=600)
def inspect_tag(tag: str) -> str:
    results.reload(); out = []
    for f in sorted(glob.glob("/results/" + tag + "/shard*/**/*", recursive=True)):
        if os.path.isfile(f):
            out.append("%8d  %s" % (os.path.getsize(f), f[len("/results/"):]))
            if f.endswith(".csv"):
                out.append("          HEADER: " + open(f, errors="ignore").readline().strip()[:300])
    return "\n".join(out[:120]) or "nothing under " + tag


# ============================ ENTRYPOINTS ===================================
@app.local_entrypoint()
def prep(indir: str):
    """Upload the 3 target PDBs to a Modal volume and check chains and peptides."""
    with inputs.batch_upload(force=True) as b:
        for f in FILES.values():
            p = os.path.join(indir, f)
            if not os.path.exists(p): raise SystemExit("missing file: " + p)
            b.put_file(p, "/" + f)
    print(check_inputs.remote())

@app.local_entrypoint()
def warm_cache():
    print(warm.remote())

@app.local_entrypoint()
def smoke(camp: str = "set1"):
    """1 shard x 30 trajectories, then refold on all 3 targets and score. Runs in the cloud (spawn).
    Same tag resumes saved state and does nothing: change "kras_smoke9999_" below to run it again."""
    print(controller.spawn(camp, "kras_smoke9999_" + camp, 1, 1, 1, 20.0, 3.0, 30))

@app.local_entrypoint()
def start(camp: str = "set1", budget: float = 40.0, goal: int = 8, max_shards: int = 12,
          per_round: int = 4, max_traj: int = 25, deadline_h: float = 20.0):
    tag = "kras_" + camp
    call = controller.spawn(camp, tag, goal, max_shards, per_round, budget, deadline_h, max_traj)
    print("CONTROLLER STARTED IN THE CLOUD  tag:", tag, " goal:", goal, " budget cap (est): $", budget)
    print("  call id:", call.object_id)
    print("You can close this window and switch the PC off now.")

@app.local_entrypoint()
def status():
    for t in ("kras_set1", "kras_set2", "kras_set4", "kras_set5"):
        print(t, read_progress.remote(t))

@app.local_entrypoint()
def inspect_out(tag: str = "kras_set5"):
    print(inspect_tag.remote(tag))

@app.local_entrypoint()
def logs(tag: str = "kras_set1"):
    print(tail_logs.remote(tag))

@app.local_entrypoint()
def fetch(outdir: str, tags: str = "kras_set5"):
    out = pathlib.Path(outdir); out.mkdir(parents=True, exist_ok=True)
    allrows = []
    for t in tags.split(","):
        b = read_file.remote(t, "designs.csv")
        if not b: print("not yet:", t); continue
        (out / (t + "_designs.csv")).write_bytes(b)
        rows = list(csv.DictReader(io.StringIO(b.decode())))
        for r in rows: r["campaign"] = t.replace("kras_", "")
        allrows += rows; print("wrote", t + "_designs.csv", len(rows), "rows")
    if allrows:
        keys = ["campaign"] + COLS
        with open(out / "kras_all_designs.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore"); w.writeheader(); w.writerows(allrows)
        print("wrote", out / "kras_all_designs.csv", len(allrows), "rows")


@app.local_entrypoint()
def refold_test():
    """Refold a dummy 78-aa sequence on G12D, WT and PF and print the combined row (checks the Boltz-2 + ipSAE path)."""
    dummy = "SPEELLKKAEELIKRAAELAKKGDLEGALKLLEEALRLAKESGNPELIARVLLLLAEIYLKLGDKEKAKEYAEEALKL"
    jobs = [(dummy, "dummy", k, "set5") for k in KINDS]
    outs = list(fold_score.starmap(jobs))
    print(combine("dummy", dummy, dict(zip(KINDS, outs))))


# ============================ RESCORE (appended) ============================
# Refold chosen designs on G12D, WT and PF with more samples; keep per-sample ipTM, ipSAE, peptide-chain ipSAE and Asp6 distance.
@app.function(image=boltz_image, gpu="L40S", timeout=3600, volumes={"/cache": bcache},
              retries=modal.Retries(max_retries=1, backoff_coefficient=1.0, initial_delay=10.0))
def fold_n(binder: str, design_id: str, kind: str, camp: str, n: int) -> dict:
    import numpy as np
    t0 = time.time()
    for p in MSA.values():
        if not os.path.exists(p): raise RuntimeError("MSA missing - run warm_cache first")
    r, out = _boltz(yaml_for(kind, binder), ["--recycling_steps", "3", "--sampling_steps", "100"], n)
    rec = {"design_id": design_id, "kind": kind, "sequence": binder}
    samples = []
    for i in range(n):
        cif = sorted(glob.glob(str(out / "**" / ("*model_%d.cif" % i)), recursive=True))
        cnf = sorted(glob.glob(str(out / "**" / ("confidence*model_%d.json" % i)), recursive=True))
        pae = sorted(glob.glob(str(out / "**" / ("pae*model_%d.npz" % i)), recursive=True))
        if not (cif and cnf and pae): continue
        c = json.load(open(cnf[0])); p = np.load(pae[0])["pae"]
        per, _ = ipsae_by_chain(p, cif[0], kind, len(binder))
        s = {"iptm": binder_iptm(c, 3 if kind == "PF" else 4), "ipsae": max(per.values()),
             "ipsae_pep": per.get("C", 0.0), "ipsae_hla": per.get("A", 0.0), "asp6": None}
        try: s["asp6"] = contact_info(cif[0], camp, kind).get("asp6_min_dist")
        except Exception: pass
        samples.append(s)
    rec["samples"] = samples; rec["seconds"] = round(time.time() - t0, 1)
    if not samples: rec["error"] = (r.stdout + r.stderr)[-600:]
    return rec


def _agg(rec):
    S = (rec or {}).get("samples") or []
    if not S: return None
    mean = lambda k: sum(s[k] for s in S) / len(S)
    d6 = sorted(s["asp6"] for s in S if s.get("asp6") is not None)
    return {"n": len(S), "iptm": mean("iptm"), "ipsae": mean("ipsae"), "ipsae_max": max(s["ipsae"] for s in S),
            "ipsae_min": min(s["ipsae"] for s in S), "pep": mean("ipsae_pep"),
            "asp6_min": d6[0] if d6 else None, "asp6_med": d6[len(d6) // 2] if d6 else None}


def rescore_row(did, seq, recs):
    row = {"design_id": did, "sequence": seq, "binder_len": len(seq)}
    A = {k: _agg(recs.get(k)) for k in KINDS}
    if any(v is None for v in A.values()):
        row["error"] = "missing: " + ",".join(k for k, v in A.items() if v is None); return row
    for k, tag in (("G12D", "g"), ("WT", "wt"), ("PF", "pf")):
        for f, v in A[k].items(): row[tag + "_" + f] = None if v is None else round(v, 4)
    row["delta_mean"] = round(A["G12D"]["ipsae"] - max(A["WT"]["ipsae"], A["PF"]["ipsae"]), 4)
    row["delta_worst"] = round(A["G12D"]["ipsae"] - max(A["WT"]["ipsae_max"], A["PF"]["ipsae_max"]), 4)
    row["pep_delta"] = round(A["G12D"]["pep"] - A["WT"]["pep"], 4)
    return row


# ---- detached rescore: runs inside Modal, writes to the results volume ----
RESROOT = "/results"
RESCORE_KEYS = ["rank", "campaign", "design_id", "binder_len", "sequence",
                "g_iptm", "g_ipsae", "g_ipsae_min", "g_ipsae_max", "g_pep", "g_asp6_min", "g_asp6_med",
                "wt_iptm", "wt_ipsae", "wt_ipsae_max", "wt_pep", "wt_asp6_min",
                "pf_iptm", "pf_ipsae", "pf_ipsae_max", "delta_mean", "delta_worst", "pep_delta", "g_n", "error"]


@app.function(image=cpu_image, cpu=2.0, timeout=4 * 3600, volumes={"/results": results})
def rescore_ctl(rows: list, tag: str, samples: int, camp: str) -> str:
    base = pathlib.Path(RESROOT) / tag; base.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    def prog(state, note=""):
        (base / "progress.json").write_text(json.dumps({"state": state, "designs": len(rows), "samples_each": samples,
            "note": note, "elapsed_min": round((time.time() - t0) / 60, 1), "updated": time.strftime("%Y-%m-%d %H:%M:%S")}, indent=1))
        results.commit()
    prog("folding")
    jobs = [(r["sequence"], r["design_id"], k, camp, samples) for r in rows for k in KINDS]
    outs = list(fold_n.starmap(jobs, return_exceptions=True))
    by, spent = {}, 0.0
    for (seq, d, k, _, _), o in zip(jobs, outs):
        if isinstance(o, Exception): o = {"error": repr(o)[:200]}
        else: spent += o.get("seconds", 0) / 3600 * L40S_USD_H
        by.setdefault(d, {})[k] = o
    res = [rescore_row(r["design_id"], r["sequence"], by.get(r["design_id"], {})) for r in rows]
    for x, r in zip(res, rows): x["campaign"] = r.get("campaign", "")
    res.sort(key=lambda x: -(x.get("delta_mean") if x.get("delta_mean") is not None else -9))
    for i, x in enumerate(res, 1): x["rank"] = i
    b = io.StringIO(); w = csv.DictWriter(b, fieldnames=RESCORE_KEYS, extrasaction="ignore"); w.writeheader(); w.writerows(res)
    (base / "rescored.csv").write_text(b.getvalue())
    prog("done", "compute est (L40S) $%.2f" % spent)
    return "done"


@app.local_entrypoint()
def rescore_start(csv_path: str, tag: str = "kras_rescore1", samples: int = 8, top: int = 20, max_usd: float = 15.0, camp: str = "set5"):
    rows = [dict(design_id=r["design_id"], sequence=r["sequence"], campaign=r.get("campaign", ""))
            for r in csv.DictReader(open(csv_path, encoding="utf-8-sig")) if r.get("sequence")][:top]
    bound = len(rows) * 3 * 0.25
    print("designs:", len(rows), " folds:", len(rows) * 3, " samples each:", samples, " worst-case est $%.2f (cap $%.2f)" % (bound, max_usd))
    if bound > max_usd: raise SystemExit("estimated upper bound above max_usd; lower --top/--samples or raise --max-usd")
    call = rescore_ctl.spawn(rows, tag, samples, camp)
    print("RESCORE STARTED IN THE CLOUD  tag:", tag, " call id:", call.object_id)
    print("You can close this window and switch the PC off now. Later run rescore_get.")


@app.local_entrypoint()
def rescore_get(outdir: str, tag: str = "kras_rescore1"):
    print(read_progress.remote(tag))
    b = read_file.remote(tag, "rescored.csv")
    if not b: print("not finished yet"); return
    out = pathlib.Path(outdir); out.mkdir(parents=True, exist_ok=True)
    (out / (tag + ".csv")).write_bytes(b)
    rows = list(csv.DictReader(io.StringIO(b.decode())))
    print("%-4s %-4s %-6s %-6s %-6s %-6s %-6s %-7s %-7s %-6s" % ("rank", "len", "g_ipt", "g_ips", "wt_ips", "pf_ips", "d_mean", "d_worst", "pep_d", "asp6"))
    for x in rows:
        if x.get("error"): print(x["rank"], x["design_id"], x["error"]); continue
        print("%-4s %-4s %-6s %-6s %-6s %-6s %-6s %-7s %-7s %-6s" % (x["rank"], x["binder_len"], x["g_iptm"], x["g_ipsae"],
              x["wt_ipsae"], x["pf_ipsae"], x["delta_mean"], x["delta_worst"], x["pep_delta"], x["g_asp6_min"]))
    print("saved:", out / (tag + ".csv"))


# ============================ PROTEINMPNN ROUND (appended) ============================
# Seeds (chosen designs) -> Boltz-2 refold on G12D to get a backbone -> ProteinMPNN redesigns binder chain D against the
# fixed target chains -> every variant refolded on G12D, WT and PF (stage 1) -> best ones rescored with more samples (stage 2).
MPNN_DIR = "/opt/ProteinMPNN"
mpnn_image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git")
    .pip_install("torch", index_url="https://download.pytorch.org/whl/cpu")
    .pip_install("numpy", "gemmi==0.6.5")
    .run_commands("git clone https://github.com/dauparas/ProteinMPNN.git /opt/ProteinMPNN && "
                  "git -C /opt/ProteinMPNN checkout 8907e66")
    .env({"PYTHONUNBUFFERED": "1"}))


@app.function(image=boltz_image, gpu="L40S", timeout=3600, volumes={"/cache": bcache})
def backbone_fold(binder: str, design_id: str, camp: str, n: int) -> dict:
    """Refold one binder on G12D with n samples and return the CIF of the sample with the best ipSAE (the backbone to redesign)."""
    import numpy as np
    for p in MSA.values():
        if not os.path.exists(p): raise RuntimeError("MSA missing - run warm_cache first")
    r, out = _boltz(yaml_for("G12D", binder), ["--recycling_steps", "3", "--sampling_steps", "100"], n)
    best = None
    for i in range(n):
        cif = sorted(glob.glob(str(out / "**" / ("*model_%d.cif" % i)), recursive=True))
        pae = sorted(glob.glob(str(out / "**" / ("pae*model_%d.npz" % i)), recursive=True))
        if not (cif and pae): continue
        per, _ = ipsae_by_chain(np.load(pae[0])["pae"], cif[0], "G12D", len(binder))
        sc = max(per.values())
        if best is None or sc > best["ipsae"]:
            best = {"design_id": design_id, "sample": i, "ipsae": round(sc, 4), "cif_text": open(cif[0]).read()}
    if best is None: raise RuntimeError("no Boltz output: " + (r.stdout + r.stderr)[-500:])
    return best


@app.function(image=mpnn_image, cpu=4.0, timeout=1800)
def mpnn_design(cif_text: str, native: str, per_temp: int, temps: str, seed: int) -> list:
    """ProteinMPNN on binder chain D of a predicted complex, target chains fixed. Returns [{seq, T, score, recovery}]."""
    import gemmi, subprocess, tempfile
    work = pathlib.Path(tempfile.mkdtemp()); cif = work / "c.cif"; cif.write_text(cif_text)
    st = gemmi.read_structure(str(cif)); st.remove_hydrogens()
    ch = {c.name: c for c in st[0]}
    if "D" not in ch: raise RuntimeError("no chain D, chains=" + ",".join(ch))
    got = gemmi.one_letter_code([r.name for r in ch["D"]])
    if got != native: raise RuntimeError("chain D sequence differs from the binder: " + got[:25])
    pdb = work / "c.pdb"; st.write_pdb(str(pdb)); outd = work / "out"
    cmd = ["python", "protein_mpnn_run.py", "--pdb_path", str(pdb), "--pdb_path_chains", "D", "--out_folder", str(outd),
           "--num_seq_per_target", str(per_temp), "--sampling_temp", temps, "--seed", str(seed), "--batch_size", "1",
           "--omit_AAs", "CX", "--suppress_print", "1"]
    r = subprocess.run(cmd, cwd=MPNN_DIR, capture_output=True, text=True)
    fa = outd / "seqs" / "c.fa"
    if not fa.exists(): raise RuntimeError("ProteinMPNN failed: " + (r.stdout + r.stderr)[-600:])
    recs, hdr = [], None
    for line in fa.read_text().splitlines():
        if line.startswith(">"): hdr = line[1:]
        elif hdr is not None:
            m = re.search(r"T=([0-9.]+), sample=(\d+), score=([0-9.]+).*seq_recovery=([0-9.]+)", hdr)
            if m: recs.append({"seq": line.strip(), "T": float(m.group(1)), "score": float(m.group(3)), "recovery": float(m.group(4))})
            hdr = None
    return recs


def _ident(a, b):
    return round(sum(1 for x, y in zip(a, b) if x == y) / len(a), 3) if len(a) == len(b) else 0.0


def _score(row):
    return (row["delta_mean"] + row["delta_worst"]) / 2


MPNN_KEYS = ["rank", "parent", "is_parent", "design_id", "mpnn_T", "identity_to_parent", "binder_len", "sequence", "stage1_score",
             "score", "g_iptm", "g_ipsae", "g_ipsae_min", "g_ipsae_max", "g_pep", "g_asp6_min", "g_asp6_med",
             "wt_iptm", "wt_ipsae", "wt_ipsae_max", "wt_pep", "wt_asp6_min", "pf_iptm", "pf_ipsae", "pf_ipsae_max",
             "delta_mean", "delta_worst", "pep_delta", "g_n", "error"]


def _write_csv(path, rows):
    b = io.StringIO(); w = csv.DictWriter(b, fieldnames=MPNN_KEYS, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    pathlib.Path(path).write_text(b.getvalue())


@app.function(image=cpu_image, cpu=2.0, timeout=6 * 3600, volumes={"/results": results})
def mpnn_ctl(seeds: list, tag: str, per_temp: int, temps: str, keep_per_seed: int, s1_samples: int,
             s2_top: int, s2_samples: int, max_usd: float, camp: str) -> str:
    base = pathlib.Path(RESROOT) / tag; base.mkdir(parents=True, exist_ok=True)
    t0 = time.time(); note = {"spent_est": 0.0}
    def prog(state, msg=""):
        (base / "progress.json").write_text(json.dumps({"state": state, "msg": msg, "seeds": [s["name"] for s in seeds],
            "est_spend_usd": round(note["spent_est"], 2), "cap_usd": max_usd,
            "elapsed_min": round((time.time() - t0) / 60, 1), "updated": time.strftime("%Y-%m-%d %H:%M:%S")}, indent=1))
        results.commit()
    def run_folds(items, n):
        jobs = [(s, d, k, camp, n) for (d, s) in items for k in KINDS]
        outs = list(fold_n.starmap(jobs, return_exceptions=True))
        by = {}
        for (s, d, k, _, _), o in zip(jobs, outs):
            if isinstance(o, Exception): o = {"error": repr(o)[:200]}
            else: note["spent_est"] += o.get("seconds", 0) / 3600 * L40S_USD_H
            by.setdefault(d, {})[k] = o
        return {d: rescore_row(d, s, by.get(d, {})) for (d, s) in items}

    prog("backbones")
    bbs = list(backbone_fold.starmap([(s["sequence"], s["name"], camp, 5) for s in seeds], return_exceptions=True))
    meta, items = {}, []
    for s, bb in zip(seeds, bbs):
        if isinstance(bb, Exception): prog("failed", "backbone fold failed for " + s["name"] + ": " + repr(bb)[:300]); return "failed"
        note["spent_est"] += 0.5
        meta[s["name"] + "_parent"] = {"parent": s["name"], "is_parent": "yes", "mpnn_T": "", "identity_to_parent": 1.0, "seq": s["sequence"]}
        items.append((s["name"] + "_parent", s["sequence"]))
    prog("mpnn")
    des = list(mpnn_design.starmap([(bb["cif_text"], s["sequence"], per_temp, temps, 100 + i)
                                    for i, (s, bb) in enumerate(zip(seeds, bbs))], return_exceptions=True))
    for s, rec in zip(seeds, des):
        if isinstance(rec, Exception): prog("failed", "ProteinMPNN failed for " + s["name"] + ": " + repr(rec)[:400]); return "failed"
        seen, j = {s["sequence"]}, 0
        for x in rec:
            if x["seq"] in seen or len(x["seq"]) != len(s["sequence"]): continue
            seen.add(x["seq"]); j += 1
            if j > keep_per_seed: break
            did = "%s_v%02d" % (s["name"], j)
            meta[did] = {"parent": s["name"], "is_parent": "", "mpnn_T": x["T"], "identity_to_parent": _ident(x["seq"], s["sequence"]), "seq": x["seq"]}
            items.append((did, x["seq"]))
    (base / "variants.csv").write_text("design_id,parent,mpnn_T,identity_to_parent,sequence\n" + "\n".join(
        "%s,%s,%s,%s,%s" % (d, meta[d]["parent"], meta[d]["mpnn_T"], meta[d]["identity_to_parent"], s) for d, s in items))
    prog("stage1", "%d sequences (incl. %d parents), %d samples each" % (len(items), len(seeds), s1_samples))
    r1 = run_folds(items, s1_samples)
    rows1 = []
    for d, s in items:
        x = r1[d]; x.update({k: v for k, v in meta[d].items() if k != "seq"}); x["stage1_score"] = None if "error" in x else round(_score(x), 4)
        rows1.append(x)
    rows1.sort(key=lambda x: -(x["stage1_score"] if x["stage1_score"] is not None else -9))
    for i, x in enumerate(rows1, 1): x["rank"] = i
    _write_csv(base / "stage1_all.csv", rows1)
    need = (s2_top + len(seeds)) * 3 * 0.05
    if note["spent_est"] + need > max_usd:
        prog("stopped", "budget: stage 1 done (est $%.2f), stage 2 skipped; see stage1_all.csv" % note["spent_est"]); return "stopped"
    ok = [x for x in rows1 if x["stage1_score"] is not None and not x["is_parent"]]
    ok.sort(key=lambda x: (-(float(x["g_ipsae"]) >= 0.5), -x["stage1_score"]))
    pick = [x["design_id"] for x in ok[:s2_top]] + [d for d in meta if meta[d]["is_parent"]]
    seqmap = dict(items)
    prog("stage2", "%d sequences, %d samples each" % (len(pick), s2_samples))
    r2 = run_folds([(d, seqmap[d]) for d in pick], s2_samples)
    s1 = {x["design_id"]: x["stage1_score"] for x in rows1}
    rows2 = []
    for d in pick:
        x = r2[d]; x.update({k: v for k, v in meta[d].items() if k != "seq"}); x["stage1_score"] = s1.get(d)
        x["score"] = None if "error" in x else round(_score(x), 4); rows2.append(x)
    rows2.sort(key=lambda x: -(x["score"] if x["score"] is not None else -9))
    for i, x in enumerate(rows2, 1): x["rank"] = i
    _write_csv(base / "stage2_all.csv", rows2)
    top = [x for x in rows2 if not x["is_parent"] and x["score"] is not None and float(x["g_ipsae"]) >= 0.5][:20]
    for i, x in enumerate(top, 1): x["rank"] = i
    _write_csv(base / "top20_new.csv", top)
    prog("done", "est spend $%.2f" % note["spent_est"]); return "done"


@app.local_entrypoint()
def mpnn_start(csv_path: str, tag: str = "kras_mpnn1", n_seeds: int = 2, per_temp: int = 17, temps: str = "0.1 0.2 0.3",
               keep_per_seed: int = 50, s1_samples: int = 3, s2_top: int = 30, s2_samples: int = 8,
               max_usd: float = 25.0, camp: str = "set5"):
    rows = [r for r in csv.DictReader(open(csv_path, encoding="utf-8-sig")) if r.get("sequence")][:n_seeds]
    seeds = [{"name": "S%d" % (i + 1), "sequence": r["sequence"]} for i, r in enumerate(rows)]
    n1 = (len(seeds) * keep_per_seed + len(seeds)) * 3; n2 = (s2_top + len(seeds)) * 3
    bound = (n1 + n2) * 0.05 + 0.5 * len(seeds)
    print("seeds:", len(seeds), " variants per seed:", keep_per_seed, " stage1 folds:", n1, " stage2 folds:", n2,
          " rough est up to $%.2f (cap $%.2f)" % (bound, max_usd))
    if bound > max_usd: raise SystemExit("rough estimate above max_usd; lower the numbers or raise --max-usd")
    call = mpnn_ctl.spawn(seeds, tag, per_temp, temps, keep_per_seed, s1_samples, s2_top, s2_samples, max_usd, camp)
    print("PROTEINMPNN ROUND STARTED IN THE CLOUD  tag:", tag, " call id:", call.object_id)
    print("You can close this window and switch the PC off now. Later run mpnn_get.")


@app.local_entrypoint()
def mpnn_get(outdir: str, tag: str = "kras_mpnn1"):
    print(read_progress.remote(tag))
    out = pathlib.Path(outdir); out.mkdir(parents=True, exist_ok=True)
    for name in ("variants.csv", "stage1_all.csv", "stage2_all.csv", "top20_new.csv"):
        b = read_file.remote(tag, name)
        if b: (out / (tag + "_" + name)).write_bytes(b); print("saved", tag + "_" + name)
    b = read_file.remote(tag, "top20_new.csv")
    if not b: print("final table not ready yet"); return
    print("%-4s %-7s %-4s %-5s %-6s %-6s %-6s %-6s %-6s %-7s %-6s %-6s" % ("rank", "parent", "len", "ident", "g_ipt", "g_ips", "wt_ips", "pf_ips", "d_mean", "d_worst", "pep_d", "asp6"))
    for x in csv.DictReader(io.StringIO(b.decode())):
        print("%-4s %-7s %-4s %-5s %-6s %-6s %-6s %-6s %-6s %-7s %-6s %-6s" % (x["rank"], x["parent"], x["binder_len"], x["identity_to_parent"],
              x["g_iptm"], x["g_ipsae"], x["wt_ipsae"], x["pf_ipsae"], x["delta_mean"], x["delta_worst"], x["pep_delta"], x["g_asp6_min"]))


MPNN_KEYS_FX = MPNN_KEYS[:5] + ["n_fixed"] + MPNN_KEYS[5:]


def _write_csv_fx(path, rows):
    b = io.StringIO(); w = csv.DictWriter(b, fieldnames=MPNN_KEYS_FX, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    pathlib.Path(path).write_text(b.getvalue())


# ============================ PROTEINMPNN, ASP6 + HOTSPOT RESIDUES KEPT FIXED (appended) ============================
@app.function(image=mpnn_image, cpu=4.0, timeout=1800)
def mpnn_design_hot(cif_text: str, native: str, per_temp: int, temps: str, seed: int, cutoff: float, asp_cutoff: float, hot: str) -> dict:
    import gemmi, subprocess, tempfile
    import numpy as np
    work = pathlib.Path(tempfile.mkdtemp()); cif = work / "c.cif"; cif.write_text(cif_text)
    st = gemmi.read_structure(str(cif)); st.remove_hydrogens()
    ch = {c.name: c for c in st[0]}
    if "D" not in ch or "C" not in ch: raise RuntimeError("need chains C and D, chains=" + ",".join(ch))
    got = gemmi.one_letter_code([r.name for r in ch["D"]])
    if got != native: raise RuntimeError("chain D sequence differs from the binder: " + got[:25])
    hot_idx = [int(x) for x in hot.split(",")]
    if len(ch["C"]) < max(hot_idx): raise RuntimeError("peptide chain too short")
    hot_xyz = np.array([a.pos.tolist() for j, r in enumerate(ch["C"], 1) if j in hot_idx for a in r])
    asp_xyz = np.array([a.pos.tolist() for a in ch["C"][5]])
    fixed = []
    for i, r in enumerate(ch["D"], 1):
        x = np.array([a.pos.tolist() for a in r])
        near_hot = (np.linalg.norm(x[:, None, :] - hot_xyz[None, :, :], axis=-1) <= cutoff).any()
        near_asp = (np.linalg.norm(x[:, None, :] - asp_xyz[None, :, :], axis=-1) <= asp_cutoff).any()
        if near_hot or near_asp: fixed.append(i)
    pdb = work / "c.pdb"; st.write_pdb(str(pdb)); outd = work / "out"
    fx = work / "fixed.jsonl"; fx.write_text(json.dumps({"c": {"A": [], "B": [], "C": [], "D": fixed}}) + "\n")
    cmd = ["python", "protein_mpnn_run.py", "--pdb_path", str(pdb), "--pdb_path_chains", "D", "--out_folder", str(outd),
           "--num_seq_per_target", str(per_temp), "--sampling_temp", temps, "--seed", str(seed), "--batch_size", "1",
           "--omit_AAs", "CX", "--fixed_positions_jsonl", str(fx), "--suppress_print", "1"]
    r = subprocess.run(cmd, cwd=MPNN_DIR, capture_output=True, text=True)
    fa = outd / "seqs" / "c.fa"
    if not fa.exists(): raise RuntimeError("ProteinMPNN failed: " + (r.stdout + r.stderr)[-600:])
    recs, hdr = [], None
    for line in fa.read_text().splitlines():
        if line.startswith(">"): hdr = line[1:]
        elif hdr is not None:
            m = re.search(r"T=([0-9.]+), sample=(\d+), score=([0-9.]+).*seq_recovery=([0-9.]+)", hdr)
            if m: recs.append({"seq": line.strip(), "T": float(m.group(1)), "score": float(m.group(3)), "recovery": float(m.group(4))})
            hdr = None
    return {"recs": recs, "fixed": fixed}


@app.function(image=cpu_image, cpu=2.0, timeout=6 * 3600, volumes={"/results": results})
def mpnn_ctl_hot(seeds: list, tag: str, per_temp: int, temps: str, keep_per_seed: int, s1_samples: int,
                s2_top: int, s2_samples: int, max_usd: float, camp: str, cutoff: float, asp_cutoff: float, hot: str) -> str:
    base = pathlib.Path(RESROOT) / tag; base.mkdir(parents=True, exist_ok=True)
    t0 = time.time(); note = {"spent_est": 0.0}
    def prog(state, msg=""):
        (base / "progress.json").write_text(json.dumps({"state": state, "msg": msg, "seeds": [s["name"] for s in seeds],
            "est_spend_usd": round(note["spent_est"], 2), "cap_usd": max_usd,
            "elapsed_min": round((time.time() - t0) / 60, 1), "updated": time.strftime("%Y-%m-%d %H:%M:%S")}, indent=1))
        results.commit()
    def run_folds(items, n):
        jobs = [(s, d, k, camp, n) for (d, s) in items for k in KINDS]
        outs = list(fold_n.starmap(jobs, return_exceptions=True))
        by = {}
        for (s, d, k, _, _), o in zip(jobs, outs):
            if isinstance(o, Exception): o = {"error": repr(o)[:200]}
            else: note["spent_est"] += o.get("seconds", 0) / 3600 * L40S_USD_H
            by.setdefault(d, {})[k] = o
        return {d: rescore_row(d, s, by.get(d, {})) for (d, s) in items}

    prog("backbones")
    bbs = list(backbone_fold.starmap([(s["sequence"], s["name"], camp, 5) for s in seeds], return_exceptions=True))
    meta, items = {}, []
    for s, bb in zip(seeds, bbs):
        if isinstance(bb, Exception): prog("failed", "backbone fold failed for " + s["name"] + ": " + repr(bb)[:300]); return "failed"
        note["spent_est"] += 0.5
        meta[s["name"] + "_parent"] = {"parent": s["name"], "is_parent": "yes", "mpnn_T": "", "identity_to_parent": 1.0, "n_fixed": "", "seq": s["sequence"]}
        items.append((s["name"] + "_parent", s["sequence"]))
    prog("mpnn")
    des = list(mpnn_design_hot.starmap([(bb["cif_text"], s["sequence"], per_temp, temps, 300 + i, cutoff, asp_cutoff, hot)
                                          for i, (s, bb) in enumerate(zip(seeds, bbs))], return_exceptions=True))
    for s, res in zip(seeds, des):
        if isinstance(res, Exception): prog("failed", "ProteinMPNN failed for " + s["name"] + ": " + repr(res)[:400]); return "failed"
        seen, j = {s["sequence"]}, 0
        for x in res["recs"]:
            if x["seq"] in seen or len(x["seq"]) != len(s["sequence"]): continue
            seen.add(x["seq"]); j += 1
            if j > keep_per_seed: break
            did = "%s_h%02d" % (s["name"], j)
            meta[did] = {"parent": s["name"], "is_parent": "", "mpnn_T": x["T"], "identity_to_parent": _ident(x["seq"], s["sequence"]),
                         "n_fixed": len(res["fixed"]), "seq": x["seq"]}
            items.append((did, x["seq"]))
    (base / "variants.csv").write_text("design_id,parent,mpnn_T,identity_to_parent,n_fixed,sequence\n" + "\n".join(
        "%s,%s,%s,%s,%s,%s" % (d, meta[d]["parent"], meta[d]["mpnn_T"], meta[d]["identity_to_parent"], meta[d]["n_fixed"], s) for d, s in items))
    prog("stage1", "%d sequences (incl. %d parents), %d samples each" % (len(items), len(seeds), s1_samples))
    r1 = run_folds(items, s1_samples)
    rows1 = []
    for d, s in items:
        x = r1[d]; x.update({k: v for k, v in meta[d].items() if k != "seq"}); x["stage1_score"] = None if "error" in x else round(_score(x), 4)
        rows1.append(x)
    rows1.sort(key=lambda x: -(x["stage1_score"] if x["stage1_score"] is not None else -9))
    for i, x in enumerate(rows1, 1): x["rank"] = i
    _write_csv_fx(base / "stage1_all.csv", rows1)
    need = (s2_top + len(seeds)) * 3 * 0.05
    if note["spent_est"] + need > max_usd:
        prog("stopped", "budget: stage 1 done (est $%.2f), stage 2 skipped; see stage1_all.csv" % note["spent_est"]); return "stopped"
    ok = [x for x in rows1 if x["stage1_score"] is not None and not x["is_parent"]]
    ok.sort(key=lambda x: (-(float(x["g_ipsae"]) >= 0.5), -x["stage1_score"]))
    pick = [x["design_id"] for x in ok[:s2_top]] + [d for d in meta if meta[d]["is_parent"]]
    seqmap = dict(items)
    prog("stage2", "%d sequences, %d samples each" % (len(pick), s2_samples))
    r2 = run_folds([(d, seqmap[d]) for d in pick], s2_samples)
    s1 = {x["design_id"]: x["stage1_score"] for x in rows1}
    rows2 = []
    for d in pick:
        x = r2[d]; x.update({k: v for k, v in meta[d].items() if k != "seq"}); x["stage1_score"] = s1.get(d)
        x["score"] = None if "error" in x else round(_score(x), 4); rows2.append(x)
    rows2.sort(key=lambda x: -(x["score"] if x["score"] is not None else -9))
    for i, x in enumerate(rows2, 1): x["rank"] = i
    _write_csv_fx(base / "stage2_all.csv", rows2)
    top = [x for x in rows2 if not x["is_parent"] and x["score"] is not None and float(x["g_ipsae"]) >= 0.5][:20]
    for i, x in enumerate(top, 1): x["rank"] = i
    _write_csv_fx(base / "top20_new.csv", top)
    prog("done", "est spend $%.2f" % note["spent_est"]); return "done"


@app.local_entrypoint()
def mpnn_hot_start(csv_path: str, rows: str = "1,2", tag: str = "kras_mpnn3", per_temp: int = 17, temps: str = "0.1 0.15 0.2",
                  keep_per_seed: int = 50, s1_samples: int = 3, s2_top: int = 30, s2_samples: int = 8,
                  max_usd: float = 25.0, camp: str = "set5", cutoff: float = 5.0, asp_cutoff: float = 7.0, hot: str = "4,5,6,7,8"):
    allrows = [r for r in csv.DictReader(open(csv_path, encoding="utf-8-sig")) if r.get("sequence")]
    want = [int(x) for x in rows.split(",")]
    seeds = [{"name": "S%d" % n, "sequence": allrows[n - 1]["sequence"]} for n in want]
    n1 = (len(seeds) * keep_per_seed + len(seeds)) * 3; n2 = (s2_top + len(seeds)) * 3
    bound = (n1 + n2) * 0.05 + 0.5 * len(seeds)
    print("seeds (CSV rows):", rows, " variants per seed:", keep_per_seed, " fixed = binder residues within", cutoff, "A of peptide residues", hot, "or within", asp_cutoff, "A of residue 6 (Asp6)",
          " stage1 folds:", n1, " stage2 folds:", n2, " rough est up to $%.2f (cap $%.2f)" % (bound, max_usd))
    if bound > max_usd: raise SystemExit("rough estimate above max_usd; lower the numbers or raise --max-usd")
    call = mpnn_ctl_hot.spawn(seeds, tag, per_temp, temps, keep_per_seed, s1_samples, s2_top, s2_samples, max_usd, camp, cutoff, asp_cutoff, hot)
    print("PROTEINMPNN (HOTSPOT + ASP6 FIXED) STARTED IN THE CLOUD  tag:", tag, " call id:", call.object_id)
    print("You can close this window and switch the PC off now. Later run mpnn_get --tag", tag)
