"""PDFA GLP1R (Adaptyv Challenge 3): BindCraft2 -> Boltz-2 refold -> ipSAE + contact check, on Modal.
Design: BC2 on H100 (bridge 30-50 aa, deep 35-60 aa). Refold: Boltz-2 (open-source boltz 2.2.1) on L40S (460 tokens; EGFR/TNF used L4).
Controller runs INSIDE Modal with its own spend cap, so the PC can be switched off.
Usage (all from the glp1r_modal folder; see README.md for what each one does):
  python -m modal run glp1r_app.py::warm_cache
  python -m modal run --detach glp1r_app.py::smoke --camp bridge
  python -m modal run --detach glp1r_app.py::start --camp bridge --budget 45
  python -m modal run --detach glp1r_app.py::start --camp bridge --goal 14 --max-shards 40 --per-round 6 --budget 70
  python -m modal run glp1r_app.py::status
  python -m modal run glp1r_app.py::logs --tag glp1r_bridge
  python -m modal run glp1r_app.py::inspect_out --tag glp1r_bridge
  python -m modal run glp1r_app.py::fetch --outdir <folder> --tags glp1r_bridge
Only the bridge campaign has been run. The deep campaign is in the code but has not been tested.
"""
import csv, glob, io, json, os, pathlib, re, shlex, time
import modal

# ============================ SETTINGS ======================================
RECEPTOR = ("RPQGATVSLWETVQKWREYRRQCQRSLTEDPPPATDLFCNRTFDEYACWPDGEPGSFVNVSCPWYLPWASSVPQGHVYRFCTAEGLWLQKDNSSLPWRDL"
            "SECEESKRGERSSPEEQLLFLYIIYTVGYALSFSALVIASAILLGFRHLHCTRNYIHLNLFASFILRALSVFIKDAALKWMYSTAAQQHQWDGLLSYQDSLS"
            "CRLVFLLMQYCVAANYYWLLVEGVYLYTLLAFSVLSEQWIFRLYVSIGWGVPLLFVVPWGIVKYLYEDEGCWTRNSNMNYWLIIRLPILFAIGVNFLIFVRVI"
            "CIVVSKLKANLMCKTDIKCRLAKSTLTLIPLLGTHEVIFAFVMDEHARGTLRFIKLFTELSFTSFQGLMVAILYCFVNNEVQLEFRKSWERWRLEHLHIQRDSS")
UNIPROT_FIRST = 24                      # receptor residue i (0-based) = UniProt i + 24
GROOVE = "A32,A39,A69,A88,A91,A121"
RIM = "A197,A201,A205,A210,A211,A214,A215,A221,A298,A299,A300,A306"
LOWRIM = "A298,A299,A300,A306"
DEEPH = "A230,A233,A234,A309,A310"
COLD = ("A63,A82,A115,A146,A207,A208,A279,A282,A283,A286,A290,A303,A304,"
        "A307,A311,A371,A373,A374,A376,A378,A381,A382")
# camp: (min len, max len, base seed, BC2 hotspots, BC2 min_hotspot_contact_final)
CAMPS = {"bridge": (30, 50, 9001, GROOVE + "," + RIM, 0.55),
         "deep":   (35, 60, 9501, LOWRIM + "," + DEEPH, 0.70)}
nums = lambda s: {int(t[1:]) for t in s.split(",")}
HOT = {"bridge": nums(GROOVE + "," + RIM), "deep": nums(LOWRIM + "," + DEEPH)}
COLDSET = nums(COLD)

WORKERS_PER_GPU = 3
SHARD_WALL_MIN = 100            # stop BC2 after this many minutes and keep what it has
PER_BACKBONE = 2
MAX_REFOLD_PER_ROUND = 40
H100_USD_H, L40S_USD_H = 3.95, 1.95
# "good" = contact check PASS and these confidence gates (the contact rules are the ones we used on the 100)
GATE = {"iptm": 0.60, "ipsae": 0.50, "hot_frac": 0.30, "max_cold": 2, "min_binder_contacts": 10}
BC2_COMMIT = "301efdd1937fb963cc40a5b0ecc1bc2f9b2b2d15"

# 6X18 chain R, renamed A, only these residues (identical to the crop used in the plan: 203 residues)
KEEP = [(29, 129), (136, 146), (193, 235), (278, 311), (370, 382), (384, 384)]
CROP_SCRIPT = (
    "import sys, os, urllib.request\n"
    "KEEP = " + repr(KEEP) + "\n"
    "src, dst = sys.argv[1], sys.argv[2]\n"
    "txt = urllib.request.urlopen(src).read().decode() if src.startswith('http') else open(src).read()\n"
    "keep = set()\n"
    "for a, b in KEEP: keep.update(range(a, b + 1))\n"
    "out = [l[:21] + 'A' + l[22:] for l in txt.splitlines()\n"
    "       if l.startswith('ATOM') and l[21] == 'R' and l[16] in ' A' and l[26] == ' ' and int(l[22:26]) in keep]\n"
    "os.makedirs(os.path.dirname(dst), exist_ok=True)\n"
    "open(dst, 'w').write('\\n'.join(out) + '\\nEND\\n')\n"
    "n = len({l[22:26] for l in out}); print('crop residues', n, 'atoms', len(out))\n"
    "assert n == 203, n\n")

app = modal.App("pdfa-glp1r-bc2")
weights = modal.Volume.from_name("pdfa-bc2-weights", create_if_missing=True)
results = modal.Volume.from_name("pdfa-bc2-results", create_if_missing=True)
bcache = modal.Volume.from_name("pdfa-boltz-cache", create_if_missing=True)

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
    .run_commands("echo " + __import__("base64").b64encode(CROP_SCRIPT.encode()).decode() + " | base64 -d > /tmp/crop.py && python3 /tmp/crop.py https://files.rcsb.org/download/6X18.pdb /adapter/GLP1R_crop.pdb")
    .env({"NVIDIA_VISIBLE_DEVICES": "all", "NVIDIA_DRIVER_CAPABILITIES": "compute,utility"}))

boltz_image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "build-essential", "curl")
    .pip_install("boltz==2.2.1", "numpy")
    .env({"PYTHONUNBUFFERED": "1", "BOLTZ_CACHE": "/cache"}))

cpu_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("numpy", "gemmi==0.6.5", "ipsae==1.0.1")
    .env({"PYTHONUNBUFFERED": "1"}))


# ============================ BC2 CONFIG ====================================
def campaign(camp, seed, project, max_traj):
    lo, hi, _, hot, floor = CAMPS[camp]
    return {
        "campaign_name": "GLP1R_" + camp, "modality": "binder", "project_folder": project,
        "campaign_seed": seed,
        "targets": [{"name": "GLP1R", "target_path": "/adapter/GLP1R_crop.pdb", "chains": "A",
                     "hotspots": hot, "coldspots": COLD, "objective": "target", "weight": 1.0}],
        "binder_lengths": [lo, hi], "number_of_final_designs": 999,
        "max_trajectories": max_traj, "desperation_trajectories": 40,
        "sequence_candidates": 30,
        "forced_targeting": False, "termini_accessible": False,
        "min_hotspot_contact_final": floor, "max_coldspot_contact_final": 0.1,
        "aa_bias": {"C": 0},
        "min_iptm_final": 0.55, "min_monomer_plddt_final": 0.65,
        "filters": {"Interface_Residues": {"threshold": 10, "higher": True}, "C_Terminus_Away_Cosine": {"threshold": 0.0, "higher": True}},
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


def pkg_ipsae(pae, cif_path):
    """Official-style ipSAE (d0chn, PAE cutoff 15, dist cutoff 10), max over A->B and B->A."""
    import tempfile, numpy as np
    from ipsae import core
    core.IpsaeCalculator._save_outputs = lambda *a, **k: None
    tmp = os.path.join(tempfile.mkdtemp(), "x_pae.npz")
    np.savez(tmp, pae=pae)
    m = core.IpsaeCalculator(15, 10).calculate(tmp, cif_path)["ipsae_scores"]["ipsae_d0chn_asym"]
    ch = sorted({str(a) for a in m} | {str(b) for a in m for b in m[a]})
    rec = "A" if "A" in ch else ch[0]
    bnd = [c for c in ch if c != rec][0]
    return float(m[rec][bnd]), float(m[bnd][rec])


def contact_check(cif_path, camp):
    import gemmi, numpy as np
    st = gemmi.read_structure(cif_path); st.setup_entities(); m = st[0]
    rec, bnd = m["A"], m["B"]
    names = {i + UNIPROT_FIRST: r.name for i, r in enumerate(rec)}
    if names.get(29) != "THR" or len(rec) != len(RECEPTOR):
        raise ValueError("receptor numbering check failed: res29=%s len=%d" % (names.get(29), len(rec)))
    rx = np.vstack([np.array([a.pos.tolist() for a in r]) for r in rec])
    rr = np.concatenate([np.full(len(r), i + UNIPROT_FIRST) for i, r in enumerate(rec)])
    bx = np.vstack([np.array([a.pos.tolist() for a in r]) for r in bnd])
    br = np.concatenate([np.full(len(r), i + 1) for i, r in enumerate(bnd)])
    near = np.linalg.norm(rx[:, None, :] - bx[None, :, :], axis=-1) <= 4.5
    crec = set(rr[near.any(1)].tolist()); cbnd = set(br[near.any(0)].tolist())
    L = len(bnd); hot = HOT[camp]
    hit = len(crec & hot); cold = len(crec & COLDSET)
    cterm = any(i in cbnd for i in range(L - 2, L + 1))
    frac = hit / len(hot)
    ok = (frac >= GATE["hot_frac"] and cold <= GATE["max_cold"]
          and len(cbnd) >= GATE["min_binder_contacts"] and not cterm)
    return {"binder_len": L, "n_binder_contacts": len(cbnd), "hotspots_hit": hit, "hotspots_total": len(hot),
            "hotspot_frac": round(frac, 3), "coldspots_hit": cold, "cterm_in_interface": cterm,
            "contact_PASS": ok, "contacted_residues": " ".join(str(x) for x in sorted(crec))}


# ============================ MODAL FUNCTIONS ===============================
@app.function(image=bc2_image, gpu="H100", timeout=60 * 150,
              volumes={"/weights": weights, "/results": results})
def design_shard(camp: str, i: int, tag: str, max_traj: int) -> dict:
    import subprocess, shutil, threading
    seed = CAMPS[camp][2] + i
    work = pathlib.Path("/tmp/bc2_" + tag + "_" + str(seed)); proj = work / "project"
    proj.mkdir(parents=True, exist_ok=True)
    cache = pathlib.Path("/weights/compile_cache"); cache.mkdir(parents=True, exist_ok=True)
    if not (proj / "compile_cache").exists():
        (proj / "compile_cache").symlink_to(cache)
    cpath = work / "campaign.json"; cpath.write_text(json.dumps(campaign(camp, seed, str(proj), max_traj), indent=2))
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


def _boltz(yaml_text, extra=()):
    import subprocess, tempfile
    work = pathlib.Path(tempfile.mkdtemp()); y = work / "c.yaml"; y.write_text(yaml_text)
    out = work / "out"
    cmd = ["boltz", "predict", str(y), "--out_dir", str(out), "--model", "boltz2", "--accelerator", "gpu",
           "--devices", "1", "--cache", "/cache", "--output_format", "mmcif", "--no_kernels",
           "--write_full_pae", "--diffusion_samples", "1"] + list(extra)
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r, out


MSA_PATH = "/cache/msa_glp1r.csv"      # receptor MSA, computed once by warm() (TNF pattern)

def _yaml(binder):
    return ("version: 1\nsequences:\n"
            "  - protein: {id: A, sequence: " + RECEPTOR + ", msa: " + MSA_PATH + "}\n"
            "  - protein: {id: B, sequence: " + binder + ", msa: empty}\n")


@app.function(image=boltz_image, gpu="L40S", timeout=3600, volumes={"/cache": bcache})
def warm() -> str:
    """Once: download Boltz-2 weights and compute + store the receptor MSA (so refolds never hit the MSA server)."""
    import shutil
    if os.path.exists(MSA_PATH): return "warm OK (already done)"
    y = "version: 1\nsequences:\n  - protein: {id: A, sequence: " + RECEPTOR + "}\n"
    r, out = _boltz(y, ["--use_msa_server", "--recycling_steps", "1", "--sampling_steps", "20"])
    csvs = sorted(glob.glob(str(out / "**" / "msa" / "*.csv"), recursive=True))
    ok = bool(glob.glob(str(out / "**" / "*model_0.cif"), recursive=True)) and bool(csvs)
    if ok:
        shutil.copy(csvs[0], MSA_PATH); bcache.commit()
        return "warm OK, receptor MSA saved (%d bytes)" % os.path.getsize(MSA_PATH)
    return "warm FAILED " + (r.stdout + r.stderr)[-800:]


@app.function(image=boltz_image, gpu="L40S", timeout=3600, volumes={"/cache": bcache},
              retries=modal.Retries(max_retries=1, backoff_coefficient=1.0, initial_delay=10.0))
def refold(binder_seq: str, design_id: str) -> dict:
    t0 = time.time()
    if not os.path.exists(MSA_PATH): raise RuntimeError("receptor MSA missing - run warm_cache first")
    r, out = _boltz(_yaml(binder_seq), ["--recycling_steps", "3", "--sampling_steps", "100"])
    cif = sorted(glob.glob(str(out / "**" / "*model_0.cif"), recursive=True))
    cnf = sorted(glob.glob(str(out / "**" / "confidence*model_0.json"), recursive=True))
    pae = sorted(glob.glob(str(out / "**" / "pae*model_0.npz"), recursive=True))
    rec = {"design_id": design_id, "sequence": binder_seq, "seconds": round(time.time() - t0, 1)}
    if cif and cnf and pae:
        rec["conf"] = json.load(open(cnf[0])); rec["cif_text"] = open(cif[0]).read()
        rec["pae_bytes"] = open(pae[0], "rb").read()
    else:
        rec["error"] = (r.stdout + r.stderr)[-800:]
    return rec


@app.function(image=cpu_image, cpu=2.0, timeout=1800)
def score_one(rec: dict, camp: str) -> dict:
    import tempfile, numpy as np
    out = {k: rec[k] for k in ("design_id", "sequence", "seconds") if k in rec}
    if "cif_text" not in rec:
        out["error"] = rec.get("error", "no structure"); out["good"] = False; return out
    cp = tempfile.mkdtemp() + "/c.cif"; open(cp, "w").write(rec["cif_text"])
    pae = np.load(io.BytesIO(rec["pae_bytes"]))["pae"]
    c = rec["conf"]
    for k in ("iptm", "ptm", "complex_plddt", "complex_iplddt", "complex_ipde", "confidence_score"):
        if k in c: out[k] = round(float(c[k]), 4)
    n = len(RECEPTOR); ia = list(range(n)); ib = list(range(n, pae.shape[0]))
    out["ipsae_own"] = round(max(own_ipsae(pae, ia, ib), own_ipsae(pae, ib, ia)), 4)
    try:
        ab, ba = pkg_ipsae(pae, cp)
        out["ipsae_A_to_B"], out["ipsae_B_to_A"] = round(ab, 4), round(ba, 4)
        out["ipsae"] = round(max(ab, ba), 4)
    except Exception as e:
        out["ipsae"] = out["ipsae_own"]; out["ipsae_note"] = "package failed, used own: " + repr(e)[:150]
    try:
        out.update(contact_check(cp, camp))
    except Exception as e:
        out["contact_error"] = repr(e)[:250]; out["contact_PASS"] = False
    out["good"] = bool(out.get("contact_PASS") and out.get("iptm", 0) >= GATE["iptm"]
                       and out["ipsae"] >= GATE["ipsae"])
    return out


@app.function(image=cpu_image, volumes={"/results": results}, timeout=900)
def harvest(tag: str, known: list) -> list:
    results.reload()
    seen, found = set(known), []
    pat = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{25,80}$")
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


COLS = ["rank", "good", "design_id", "binder_len", "sequence", "iptm", "ptm", "ipsae", "ipsae_A_to_B",
        "ipsae_B_to_A", "ipsae_own", "complex_plddt", "complex_iplddt", "complex_ipde", "confidence_score",
        "hotspots_hit", "hotspots_total", "hotspot_frac", "coldspots_hit", "n_binder_contacts",
        "cterm_in_interface", "contact_PASS", "contacted_residues", "contact_error", "ipsae_note", "error"]

def write_outputs(base, camp, scored):
    rows = [dict(v, design_id=d) for d, v in scored.items()]
    rows.sort(key=lambda r: (not r.get("good"), -(r.get("ipsae") or -1)))
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
        sp.write_text(json.dumps(st)); g = write_outputs(base, camp, st["scored"])
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
        folds = list(refold.starmap([(st["seqs"][d], d) for d in take], return_exceptions=True))
        recs = []
        for d, f in zip(take, folds):
            if isinstance(f, Exception): f = {"design_id": d, "sequence": st["seqs"][d], "error": repr(f)[:300]}
            else: st["spent"] += f.get("seconds", 300) / 3600 * L40S_USD_H
            recs.append(f)
        save("scoring")
        scores = list(score_one.map(recs, [camp] * len(recs), return_exceptions=True))
        for f, s in zip(recs, scores):
            d = f["design_id"]
            if isinstance(s, Exception): s = {"design_id": d, "error": repr(s)[:300], "good": False}
            st["spent"] += 0.01
            st["scored"][d] = s
            if s.get("good") and "cif_text" in f:
                p = base / "refold" / d.replace("|", "_"); p.mkdir(parents=True, exist_ok=True)
                (p / "complex.cif").write_text(f["cif_text"])

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
    """List what BindCraft2 actually wrote, plus the header line of each csv (to check the harvest step)."""
    results.reload(); out = []
    for f in sorted(glob.glob("/results/" + tag + "/shard*/**/*", recursive=True)):
        if os.path.isfile(f):
            out.append("%8d  %s" % (os.path.getsize(f), f[len("/results/"):]))
            if f.endswith(".csv"):
                out.append("          HEADER: " + open(f, errors="ignore").readline().strip()[:300])
    return "\n".join(out[:120]) or "nothing under " + tag


# ============================ ENTRYPOINTS ===================================
@app.local_entrypoint()
def warm_cache():
    print(warm.remote())

@app.local_entrypoint()
def smoke(camp: str = "bridge"):
    """1 shard x 20 trajectories, then refold + score whatever it produced (~$1-3). Keep the PC on until it prints.
    A second run with the same tag resumes the saved state and does nothing: change "glp1r_smoke_" to run it again."""
    print(controller.remote(camp, "glp1r_smoke_" + camp, 1, 1, 1, 20.0, 3.0, 20))
    print(read_progress.remote("glp1r_smoke_" + camp))

@app.local_entrypoint()
def start(camp: str = "bridge", budget: float = 40.0, goal: int = 0, max_shards: int = 12,
          per_round: int = 4, max_traj: int = 25, deadline_h: float = 20.0):
    goal = goal or (14 if camp == "bridge" else 6)
    tag = "glp1r_" + camp
    call = controller.spawn(camp, tag, goal, max_shards, per_round, budget, deadline_h, max_traj)
    print("CONTROLLER STARTED IN THE CLOUD  tag:", tag, " goal:", goal, " budget cap (est): $", budget)
    print("  call id:", call.object_id)
    print("You can close this window and switch the PC off now.")

@app.local_entrypoint()
def status():
    for t in ("glp1r_smoke_bridge", "glp1r_bridge", "glp1r_deep"):
        print(t, read_progress.remote(t))

@app.local_entrypoint()
def inspect_out(tag: str = "glp1r_smoke_bridge"):
    print(inspect_tag.remote(tag))

@app.local_entrypoint()
def logs(tag: str = "glp1r_bridge"):
    print(tail_logs.remote(tag))

@app.local_entrypoint()
def fetch(outdir: str, tags: str = "glp1r_bridge,glp1r_deep"):
    out = pathlib.Path(outdir); out.mkdir(parents=True, exist_ok=True)
    allrows = []
    for t in tags.split(","):
        b = read_file.remote(t, "designs.csv")
        if not b: print("not yet:", t); continue
        (out / (t + "_designs.csv")).write_bytes(b)
        rows = list(csv.DictReader(io.StringIO(b.decode())))
        for r in rows: r["campaign"] = t.replace("glp1r_", "")
        allrows += rows; print("wrote", t + "_designs.csv", len(rows), "rows")
    if allrows:
        keys = ["campaign"] + COLS
        with open(out / "glp1r_bc2_refolded_ipsae.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore"); w.writeheader(); w.writerows(allrows)
        print("wrote", out / "glp1r_bc2_refolded_ipsae.csv", len(allrows), "rows")
