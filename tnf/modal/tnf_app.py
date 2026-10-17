"""PDFA TNF-alpha dual-species (human + mouse) pH-switch binder campaign on Modal.
Design: BindCraft2 on A/B dimers (H100).  Check: Boltz-2 on FULL trimer + binder (L4).
Score: inline ipSAE / interface PAE / groove hotspots and PROPKA-Wyman His linkage (pH switch) on BOTH species.
Controller runs INSIDE Modal (spawned) so the PC can be switched off."""
import json, os, csv, io, pathlib, time, re, random, modal

# ============================ SETTINGS ======================================
HUMAN_AB, MOUSE_AB = "tnf_human_AB.pdb", "tnf_mouse_AB.pdb"
HUMAN_TRI, MOUSE_TRI = "tnf_human_trimer_12on.pdb", "tnf_mouse_trimer_12on.pdb"

HOT = {"human": "A29-29,A32-33,A143-147,B86-87,B90-92",
       "mouse": "A29-29,A32-33,A142-146,B85-86,B89-91"}
COLD = {
 "human": "A55-55,A57-57,A71-71,A73-73,A75-75,A82-82,A91-99,A102-103,A118-125,A151-157,B13-13,B15-15,B34-34,B36-36,B57-57,B59-59,B61-61,B63-63,B69-69,B98-98,B100-106,B112-117,B119-119,B143-143,B146-149,B151-157",
 "mouse": "A55-55,A57-57,A72-72,A74-74,A81-81,A90-98,A100-103,A117-124,A151-156,B13-13,B15-15,B34-34,B36-36,B57-57,B59-59,B61-61,B63-63,B97-97,B102-105,B111-116,B118-118,B145-148,B150-156"}

BINDER_LENGTHS = [70, 120]
AA_BIAS = {"C": 0, "H": 2.5, "D": 0.8, "E": 0.8, "K": 1.0, "R": 1.0}
MAX_TRAJ_PER_SHARD = 10
WORKERS_PER_GPU = 2
BASE_SEED = 9000
BC2_COMMIT = "301efdd1937fb963cc40a5b0ecc1bc2f9b2b2d15"

# ---- controller defaults (all overridable from the command line)
GOAL = 24                  # stop when this many designs pass every gate (you submit the best 20)
MAX_SHARDS = 14
SHARDS_PER_ROUND = 3
BUDGET_USD = 85.0          # controller's own spend cap (estimate); Modal credit is the real cap
DEADLINE_HOURS = 14.0
MAX_REFOLD_SEQS_PER_ROUND = 60
PER_BACKBONE = 2
SCAN_N = 5                 # His-scan: point mutants to His tried per near-miss design (learned from the colleague's run)
H100_USD_H, L4_USD_H, CPU_USD_H = 3.95, 0.80, 0.30

# ---- score gates (defaults; adjustable later via --gate-json)
GATE = {"iptm": 0.75, "ipsae": 0.60, "cpae": 7.0, "bplddt": 0.80,
        "hot": 0.40, "cold": 0.30, "chains": 2, "his_fold": 0.30}
PH_LO, PH_HI = 6.0, 7.4

HERE = pathlib.Path(__file__).resolve().parent
PDBS = [HUMAN_AB, MOUSE_AB, HUMAN_TRI, MOUSE_TRI]
app = modal.App("pdfa-tnf-dual")
weights = modal.Volume.from_name("pdfa-bc2-weights", create_if_missing=True)
results = modal.Volume.from_name("pdfa-bc2-results", create_if_missing=True)
bcache = modal.Volume.from_name("pdfa-boltz-cache", create_if_missing=True)

def _add(img):
    for p in PDBS:
        img = img.add_local_file(HERE / p, "/adapter/" + p, copy=True)
    return img

bc2_image = _add(
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

boltz_image = _add(
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "build-essential", "curl")
    .pip_install("boltz==2.2.1", "biopython", "numpy")
    .env({"PYTHONUNBUFFERED": "1", "BOLTZ_CACHE": "/cache"}))

cpu_image = _add(
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("propka==3.5.1", "biopython", "numpy", "pandas")
    .env({"PYTHONUNBUFFERED": "1"}))


# ============================ BC2 CONFIG ====================================
def campaign(seed, project):
    def tgt(sp, path):
        return {"name": "TNF_" + sp, "target_path": "/adapter/" + path, "chains": "A,B",
                "hotspots": HOT[sp], "coldspots": COLD[sp], "objective": "target", "weight": 1.0}
    return {
        "campaign_name": "TNF_C2_dual", "modality": "binder", "project_folder": project,
        "campaign_seed": seed, "number_of_final_designs": 999,
        "max_trajectories": MAX_TRAJ_PER_SHARD, "desperation_trajectories": 40,
        "sequence_candidates": 30,
        "targets": [tgt("human", HUMAN_AB), tgt("mouse", MOUSE_AB)],
        "binder_lengths": BINDER_LENGTHS, "forced_targeting": False,
        "termini_accessible": False, "aa_bias": AA_BIAS,
        "min_hotspot_contact_final": 0.4, "max_coldspot_contact_final": 0.1,
        "min_receptor_chains_contacted_final": 2,
        "min_iptm_final": 0.55, "max_ipae_final": 0.50, "min_ptm_final": 0.55,
        "min_monomer_plddt_final": 0.70,
        "filters": {"Interface_Residues": {"threshold": 12, "higher": True},
                    "Binder_RMSD": {"threshold": 3.5, "higher": False}},
        "workers_per_gpu": WORKERS_PER_GPU}


# ============================ SCORING (pure python) =========================
T31 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V'}

def pdb_chain_info(path, chain="A"):
    seq, nums = [], []
    for l in open(path).read().splitlines():
        if l.startswith("ATOM") and l[12:16].strip() == "CA" and l[21] == chain:
            seq.append(T31[l[17:20]]); nums.append(int(l[22:26]))
    return "".join(seq), nums

def parse_ranges(s):
    out = {}
    for tok in s.split(","):
        tok = tok.strip()
        if not tok: continue
        c, a, b = tok[0], tok[1:].split("-")[0], tok[1:].split("-")[-1]
        out.setdefault(c, set()).update(range(int(a), int(b) + 1))
    return out

def ipsae_asym(pae, ia, ib, cutoff=15.0):
    import numpy as np
    sub = pae[np.ix_(ia, ib)]
    valid = sub < cutoff
    n0 = valid.sum(1)
    L = np.maximum(n0, 27).astype(float)
    d0 = np.maximum(1.24 * np.cbrt(L - 15) - 1.8, 1.0)
    ptm = 1.0 / (1.0 + (sub / d0[:, None]) ** 2)
    per = np.where(n0 > 0, (ptm * valid).sum(1) / np.maximum(n0, 1), 0.0)
    return float(per.max())

def run_propka(pdb_path):
    import subprocess, os
    d = os.path.dirname(pdb_path)
    subprocess.run(["propka3", os.path.basename(pdb_path)], cwd=d, capture_output=True, timeout=900)
    return parse_propka(pdb_path[:-4] + ".pka")

def parse_propka(pka_path):
    out = {}
    if not os.path.exists(pka_path): return out
    on = False
    for line in open(pka_path):
        if "SUMMARY OF THIS PREDICTION" in line: on = True; continue
        if not on: continue
        f = line.split()
        if len(f) >= 4 and len(f[0]) == 3 and f[1].lstrip("-").isdigit() and f[2].isalpha():
            try: out[(f[2], int(f[1]), f[0])] = float(f[3])
            except ValueError: pass
    return out

def wyman(pc, pf, only_his=False, lo=PH_LO, hi=PH_HI):
    """log10 [Ka(hi)/Ka(lo)] = -integral(lo..hi) dn dpH ; dn = sum(f_complex - f_free)."""
    import numpy as np
    keys = [k for k in set(pc) & set(pf) if (k[2] == "HIS" or not only_his)]
    if not keys: return 0.0
    fp = lambda pk, ph: 1.0 / (1.0 + 10.0 ** (ph - pk))
    g = np.linspace(lo, hi, 29)
    dn = np.array([sum(fp(pc[k], p) - fp(pf[k], p) for k in keys) for p in g])
    tz = getattr(np, "trapezoid", None) or np.trapz
    return float(-tz(dn, g))

def structural_score(cif_path, pae, conf, species):
    import numpy as np
    from Bio.PDB import MMCIFParser
    tri = HUMAN_TRI if species == "human" else MOUSE_TRI
    _, nums = pdb_chain_info("/adapter/" + tri, "A")
    idx = {n: i for i, n in enumerate(nums)}
    n = len(nums)
    m = MMCIFParser(QUIET=True).get_structure("x", cif_path)[0]
    ch = {c.id: [r for r in c if r.id[0] == " "] for c in m}
    if sorted(ch) != ["A", "B", "C", "D"] or any(len(ch[c]) != n for c in "ABC"):
        return {"error": "unexpected chains " + str({c: len(v) for c, v in ch.items()})}
    nb = len(ch["D"])
    ib = list(range(3 * n, 3 * n + nb)); it = list(range(3 * n))
    o = {"binder_length": nb}
    o["iptm"] = float(conf.get("iptm", 0)); o["ptm"] = float(conf.get("ptm", 0))
    o["ipsae"] = round(max(ipsae_asym(pae, it, ib), ipsae_asym(pae, ib, it)), 4)
    bp = [r["CA"].get_bfactor() for r in ch["D"] if "CA" in r]
    o["binder_plddt"] = round(float(np.mean(bp)) / (100.0 if np.mean(bp) > 1.5 else 1.0), 4)
    # contacts (heavy atoms < 5 A), vectorised
    tx = np.vstack([np.array([at.coord for at in r]) for c in "ABC" for r in ch[c]])
    tr = np.concatenate([np.full(len(r), ci * n + i) for ci, c in enumerate("ABC") for i, r in enumerate(ch[c])])
    bx = np.vstack([np.array([at.coord for at in r]) for r in ch["D"]])
    br = np.concatenate([np.full(len(r), j) for j, r in enumerate(ch["D"])])
    near = np.linalg.norm(tx[:, None, :] - bx[None, :, :], axis=-1) < 5.0
    ti_set = set(tr[np.where(near.any(1))[0]].tolist()); bi_set = set(br[np.where(near.any(0))[0]].tolist())
    contact = {c: set(t - ci * n for t in ti_set if t // n == ci) for ci, c in enumerate("ABC")}
    o["chains_contacted"] = sum(1 for c in "ABC" if len(contact[c]) >= 3)
    o["binder_contact"] = sorted(bi_set)
    ti = sorted(ti_set); bi = [3 * n + j for j in sorted(bi_set)]
    if ti and bi:
        o["contact_pae"] = round(float((pae[np.ix_(ti, bi)].mean() + pae[np.ix_(bi, ti)].mean()) / 2), 3)
    else:
        o["contact_pae"] = 99.0
    # grooves: (X,Y) = (A,B),(B,C),(C,A); hotspots 'A' residues sit on X, 'B' residues on Y
    hot, cold = parse_ranges(HOT[species]), parse_ranges(COLD[species])
    best = None
    for X, Y in (("A", "B"), ("B", "C"), ("C", "A")):
        def frac(spec):
            tot = hit = 0
            for role, chain in (("A", X), ("B", Y)):
                for rn in spec.get(role, ()):
                    if rn in idx:
                        tot += 1; hit += idx[rn] in contact[chain]
            return hit / tot if tot else 0.0
        hf, cf = frac(hot), frac(cold)
        if best is None or (hf, -cf) > (best[1], -best[2]): best = (X + Y, hf, cf)
    o["groove"], o["hot_frac"], o["cold_frac"] = best[0], round(best[1], 3), round(best[2], 3)
    return o

def structural_pass(o):
    return ("error" not in o and o["iptm"] >= GATE["iptm"] and o["ipsae"] >= GATE["ipsae"]
            and o["contact_pae"] <= GATE["cpae"] and o["binder_plddt"] >= GATE["bplddt"]
            and o["hot_frac"] >= GATE["hot"] and o["cold_frac"] <= GATE["cold"]
            and o["chains_contacted"] >= GATE["chains"])

def ph_score(cif_path):
    import tempfile, numpy as np
    from concurrent.futures import ThreadPoolExecutor
    from Bio.PDB import MMCIFParser, PDBIO, Select
    m = MMCIFParser(QUIET=True).get_structure("x", cif_path)[0]
    td = tempfile.mkdtemp()
    def write(name, chains):
        class S(Select):
            def accept_chain(s, c): return c.id in chains
            def accept_residue(s, r): return r.id[0] == " "
        io = PDBIO(); io.set_structure(m); p = td + "/" + name + ".pdb"; io.save(p, S()); return p
    paths = [write("cx", "ABCD"), write("tg", "ABC"), write("bd", "D")]
    with ThreadPoolExecutor(3) as ex:
        pc, pt, pb = list(ex.map(run_propka, paths))
    pf = {}; pf.update(pt); pf.update(pb)
    o = {"his_log10_fold": round(wyman(pc, pf, True), 3), "all_log10_fold": round(wyman(pc, pf, False), 3)}
    # His shifts + carboxylate proximity
    carb = np.array([a.coord for r in m.get_residues() if r.get_resname() in ("ASP", "GLU")
                     for a in r if a.get_id() in ("OD1", "OD2", "OE1", "OE2")])
    shifts, nsw, nclean = [], 0, 0
    for k in sorted(set(pc) & set(pf)):
        if k[2] != "HIS": continue
        d = pc[k] - pf[k]
        if abs(d) < 0.5: continue
        r = m[k[0]][(" ", k[1], " ")] if (" ", k[1], " ") in m[k[0]] else None
        dc = 99.0
        if r is not None and len(carb) and "NE2" in r and "ND1" in r:
            ring = np.array([r["ND1"].coord, r["NE2"].coord])
            dc = float(np.linalg.norm(carb[:, None, :] - ring[None, :, :], axis=-1).min())
        good = d <= -0.5 and 6.5 <= pf[k] <= 8.0
        clean = good and dc >= 5.0
        nsw += good; nclean += clean
        shifts.append("%s%d:%.1f>%.1f%s" % (k[0], k[1], pf[k], pc[k], "*" if clean else ("+" if good else "")))
    o["n_switch_his"], o["n_clean_switch_his"] = nsw, nclean
    o["his_shifts"] = " ".join(shifts)
    return o

def ph_pass(o):
    return o.get("his_log10_fold", -9) >= GATE["his_fold"]


# ============================ MODAL FUNCTIONS ===============================
@app.function(image=bc2_image, gpu="H100", timeout=60 * 180,
              volumes={"/weights": weights, "/results": results})
def design_shard(i: int, tag: str) -> dict:
    import subprocess, shutil
    seed = BASE_SEED + i
    work = pathlib.Path("/tmp/bc2_" + tag + "_" + str(seed)); proj = work / "project"
    proj.mkdir(parents=True, exist_ok=True)
    cache = pathlib.Path("/weights/compile_cache"); cache.mkdir(parents=True, exist_ok=True)
    if not (proj / "compile_cache").exists():
        (proj / "compile_cache").symlink_to(cache)
    cpath = work / "campaign.json"; cpath.write_text(json.dumps(campaign(seed, str(proj)), indent=2))
    out = pathlib.Path("/results") / tag / ("shard" + format(i, "02d") + "_seed" + str(seed))
    out.mkdir(parents=True, exist_ok=True)
    (out / "campaign.json").write_text(cpath.read_text()); results.commit()
    env = dict(os.environ, BINDCRAFT_WORKERS_PER_GPU=str(WORKERS_PER_GPU), BINDCRAFT_WEIGHTS="/weights")
    t0 = time.time()
    with (out / "campaign.log").open("w") as log:
        p = subprocess.Popen(["bindcraft", "design", str(cpath)], cwd=str(work),
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, env=env)
        for line in p.stdout:
            log.write(line); print("[" + format(i, "02d") + "] " + line, end="")
        rc = p.wait()
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


@app.function(image=boltz_image, gpu="L4", timeout=3600, volumes={"/cache": bcache})
def warm() -> str:
    """Once: download Boltz weights, compute & store the two target MSAs (so refolds never hit the MSA server)."""
    import glob, shutil
    if all(os.path.exists("/cache/msa_" + sp + ".csv") for sp in ("human", "mouse")): return "already warm"
    msgs = []
    for sp, tri in (("human", HUMAN_TRI), ("mouse", MOUSE_TRI)):
        seq, _ = pdb_chain_info("/adapter/" + tri, "A")
        y = "version: 1\nsequences:\n  - protein: {id: [A, B, C], sequence: " + seq + "}\n"
        r, out = _boltz(y, ["--use_msa_server", "--recycling_steps", "1", "--sampling_steps", "20"])
        csvs = sorted(glob.glob(str(out / "**" / "msa" / "*.csv"), recursive=True))
        ok = bool(glob.glob(str(out / "**" / "*model_0.cif"), recursive=True)) and csvs
        if ok: shutil.copy(csvs[0], "/cache/msa_" + sp + ".csv")
        msgs.append(sp + ": " + ("OK msa saved (%d bytes)" % os.path.getsize("/cache/msa_" + sp + ".csv")
                                 if ok else "FAILED " + (r.stdout + r.stderr)[-600:]))
    bcache.commit()
    return "\n".join(msgs)


@app.function(image=boltz_image, gpu="L4", timeout=3600, volumes={"/cache": bcache},
              retries=modal.Retries(max_retries=1, backoff_coefficient=1.0, initial_delay=10.0))
def refold(binder_seq: str, species: str, design_id: str) -> dict:
    import glob
    t0 = time.time()
    seq, _ = pdb_chain_info("/adapter/" + (HUMAN_TRI if species == "human" else MOUSE_TRI), "A")
    msa = "/cache/msa_" + species + ".csv"
    if not os.path.exists(msa): raise RuntimeError("target MSA missing - run the warm step first")
    y = ("version: 1\nsequences:\n"
         "  - protein: {id: [A, B, C], sequence: " + seq + ", msa: " + msa + "}\n"
         "  - protein: {id: D, sequence: " + binder_seq + ", msa: empty}\n")
    r, out = _boltz(y, ["--recycling_steps", "3", "--sampling_steps", "100"])
    cif = sorted(glob.glob(str(out / "**" / "*model_0.cif"), recursive=True))
    cnf = sorted(glob.glob(str(out / "**" / "confidence*model_0.json"), recursive=True))
    pae = sorted(glob.glob(str(out / "**" / "pae*model_0.npz"), recursive=True))
    rec = {"design_id": design_id, "species": species, "sequence": binder_seq,
           "seconds": round(time.time() - t0, 1)}
    if cif and cnf and pae:
        rec["conf"] = json.load(open(cnf[0])); rec["cif_text"] = open(cif[0]).read()
        rec["pae_bytes"] = open(pae[0], "rb").read()
    else:
        rec["error"] = (r.stdout + r.stderr)[-800:]
    return rec


@app.function(image=cpu_image, cpu=4.0, timeout=3600)
def score_one(rec: dict, gate: dict = None) -> dict:
    """Structural gate first; PROPKA/Wyman only if the complex passes (saves time)."""
    import tempfile, io as _io, numpy as np
    if gate: GATE.update(gate)
    out = {k: rec[k] for k in ("design_id", "species", "sequence", "seconds") if k in rec}
    if "cif_text" not in rec:
        out["error"] = rec.get("error", "no structure"); out["pass"] = False; return out
    d = tempfile.mkdtemp(); cp = d + "/c.cif"; open(cp, "w").write(rec["cif_text"])
    pae = np.load(_io.BytesIO(rec["pae_bytes"]))["pae"]
    s = structural_score(cp, pae, rec["conf"], rec["species"]); out.update(s)
    out["struct_pass"] = structural_pass(s)
    # pH switch (His-linkage 6.0 -> 7.4) required on BOTH human and mouse complexes
    if out["struct_pass"]:
        try: out.update(ph_score(cp))
        except Exception as e: out["ph_error"] = repr(e)[:200]
    out["pass"] = bool(out["struct_pass"] and ph_pass(out))
    return out


@app.function(image=cpu_image, volumes={"/results": results}, timeout=900)
def harvest(tag: str, known: list) -> list:
    import glob, csv as _csv
    results.reload()
    seen, out = set(known), []
    pat = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{50,250}$")
    def add(did, s):
        if pat.match(s) and s not in seen:
            seen.add(s); out.append({"design_id": did, "sequence": s})
    for f in glob.glob("/results/" + tag + "/shard*/**/*.fasta", recursive=True):
        name, buf = None, []
        for line in open(f):
            if line.startswith(">"):
                if name: add(name, "".join(buf))
                name, buf = pathlib.Path(f).parts[3] + "|" + line[1:].split()[0], []
            else: buf.append(line.strip())
        if name: add(name, "".join(buf))
    for f in glob.glob("/results/" + tag + "/shard*/**/*.csv", recursive=True):
        try: rows = list(_csv.DictReader(open(f)))
        except Exception: continue
        for i, r in enumerate(rows):
            s = next((str(r[k]).strip() for k in r if k and k.lower() in
                      ("sequence", "binder_seq", "seq", "binder_sequence")), "")
            did = next((str(r[k]) for k in r if k and "design" in k.lower()), pathlib.Path(f).stem + "_" + str(i))
            add(pathlib.Path(f).parts[3] + "|" + did, s)
    per, keep = {}, []
    for d in out:
        bb = d["design_id"].rsplit("_candidate", 1)[0]
        if per.get(bb, 0) < PER_BACKBONE: per[bb] = per.get(bb, 0) + 1; keep.append(d)
    print("harvested %d new sequences (%d raw)" % (len(keep), len(out)))
    return keep


def backbone_of(did): return did.rsplit("_candidate", 1)[0]

def write_outputs(base, state, final=False):
    both = []
    for did, sp in state["scored"].items():
        h, m = sp.get("human"), sp.get("mouse")
        if h and m and h.get("pass") and m.get("pass"):
            r = {"design_id": did, "sequence": state["seqs"].get(did, ""), "length": h.get("binder_length"),
                 "min_ipsae": min(h["ipsae"], m["ipsae"]),
                 "min_his_log10_fold": min(h["his_log10_fold"], m["his_log10_fold"])}
            r["score"] = round(r["min_ipsae"] + 0.2 * min(1.5, r["min_his_log10_fold"]), 4)
            for sn, v in (("human", h), ("mouse", m)):
                for k, x in v.items():
                    if k not in ("design_id", "species", "sequence"): r[sn + "_" + k] = x
            both.append(r)
    both.sort(key=lambda r: -r["score"])
    cnt, rank = {}, 0
    for r in both:
        bb = backbone_of(r["design_id"]); cnt[bb] = cnt.get(bb, 0) + 1
        r["pick"] = "YES" if (cnt[bb] <= 2 and rank < 20) else ""
        rank += r["pick"] == "YES"
    def dump(name, rows):
        if not rows: (base / name).write_text(""); return
        keys = ["design_id", "pick", "score", "min_ipsae", "min_his_log10_fold", "length", "sequence"]
        keys += sorted({k for r in rows for k in r} - set(keys))
        b = io.StringIO(); w = csv.DictWriter(b, fieldnames=keys, extrasaction="ignore"); w.writeheader()
        w.writerows(rows); (base / name).write_text(b.getvalue())
    dump("designs.csv", both)
    allrows = [dict(v, design_id=d, species=sp) for d, s in state["scored"].items() for sp, v in s.items()]
    dump("all_scored.csv", allrows)
    (base / "top20.fasta").write_text("".join(">%s\n%s\n" % (r["design_id"], r["sequence"])
                                              for r in both if r["pick"] == "YES"))
    return len(both)


@app.function(image=cpu_image, cpu=2.0, timeout=12 * 3600, volumes={"/results": results})
def controller(tag: str, goal: int, max_shards: int, per_round: int, budget: float,
               deadline_h: float, gate_json: str = "") -> str:
    if gate_json: GATE.update(json.loads(gate_json))
    t0 = time.time(); base = pathlib.Path("/results") / tag; base.mkdir(parents=True, exist_ok=True)
    sp_path = base / "state.json"; results.reload()
    st = json.loads(sp_path.read_text()) if sp_path.exists() else {
        "launched": 0, "spent": 0.0, "seqs": {}, "pending": [], "scored": {}, "shard_cost": []}
    rnd = 0; stop = None

    def save(state_name, stop_reason=None):
        sp_path.write_text(json.dumps(st))
        n_good = write_outputs(base, st)
        (base / "progress.json").write_text(json.dumps({
            "state": state_name, "round": rnd, "shards_launched": st["launched"],
            "sequences_seen": len(st["seqs"]), "pending": len(st["pending"]),
            "scored_designs": len(st["scored"]), "good_all_gates": n_good, "goal": goal,
            "est_spend_usd": round(st["spent"], 2), "budget_cap_usd": budget,
            "elapsed_h": round((time.time() - t0) / 3600, 2), "stop_reason": stop_reason,
            "gate": GATE, "updated": time.strftime("%Y-%m-%d %H:%M:%S")}, indent=1))
        results.commit()

    def n_good():
        return sum(1 for s in st["scored"].values()
                   if s.get("human", {}).get("pass") and s.get("mouse", {}).get("pass"))

    def process_pending():
        take = st["pending"][:MAX_REFOLD_SEQS_PER_ROUND]; st["pending"] = st["pending"][MAX_REFOLD_SEQS_PER_ROUND:]
        if not take: return
        for stage in ("human", "mouse"):
            ids = [d for d in take if stage == "human" or st["scored"].get(d, {}).get("human", {}).get("pass")]
            if not ids: break
            save("refolding_" + stage)
            folds = list(refold.starmap([(st["seqs"][d], stage, d) for d in ids], return_exceptions=True))
            recs = []
            for d, f in zip(ids, folds):
                if isinstance(f, Exception): f = {"design_id": d, "species": stage, "error": repr(f)[:300]}
                else: st["spent"] += f.get("seconds", 200) / 3600 * L4_USD_H
                recs.append(f)
            keep_dir = base / "refold"
            save("scoring_" + stage)
            scores = list(score_one.map(recs, [GATE] * len(recs), return_exceptions=True))
            for f, s in zip(recs, scores):
                d = f["design_id"]
                if isinstance(s, Exception): s = {"design_id": d, "species": stage, "error": repr(s)[:300], "pass": False}
                st["spent"] += 0.01
                if s.get("struct_pass") and not s.get("pass") and "__" not in d and d in st["seqs"]:
                    # near-miss: good binder but weak pH switch -> try His point mutants at contact residues
                    seq = st["seqs"][d]; have = set(st["seqs"].values())
                    pos = [p for p in s.get("binder_contact", []) if p < len(seq) and seq[p] not in "HGPC"]
                    random.Random(d).shuffle(pos)
                    for p in pos[:SCAN_N]:
                        ms = seq[:p] + "H" + seq[p + 1:]
                        if ms not in have:
                            nm = "%s__%s%dH" % (d, seq[p], p + 1); st["seqs"][nm] = ms; st["pending"].append(nm)
                st["scored"].setdefault(d, {})[stage] = {k: v for k, v in s.items() if k != "cif_text"}
                if s.get("struct_pass") and "cif_text" in f:
                    p = keep_dir / d.replace("|", "_") / stage; p.mkdir(parents=True, exist_ok=True)
                    (p / "complex.cif").write_text(f["cif_text"])

    save("starting")
    w = warm.remote()           # one-time Boltz weights + target MSAs (skipped if already done)
    if "FAILED" in w:
        save("done", "warm failed"); return "STOPPING: warm failed: " + w
    # finish anything left from a previous run first
    if st["pending"]: process_pending(); save("resumed")
    while True:
        rnd += 1
        avg = (sum(st["shard_cost"]) / len(st["shard_cost"])) if st["shard_cost"] else 5.0
        if n_good() >= goal: stop = "goal reached"; break
        if st["launched"] >= max_shards: stop = "shard ceiling reached"; break
        if deadline_h - (time.time() - t0) / 3600 <= 1.6: stop = "deadline reached"; break
        k = min(per_round, max_shards - st["launched"])
        while k > 0 and st["spent"] + k * (avg * 1.15 + 3.0) > budget: k -= 1
        if k <= 0: stop = "budget cap reached (est $%.2f)" % st["spent"]; break
        idx = list(range(st["launched"], st["launched"] + k))
        print("=== round %d: shards %s | good %d/%d | est $%.2f" % (rnd, idx, n_good(), goal, st["spent"]))
        save("designing")
        rows = list(design_shard.map(idx, [tag] * k, return_exceptions=True))
        st["launched"] += k
        if all(isinstance(r, Exception) or r.get("rc", 1) != 0 for r in rows):
            for r in rows:
                st["spent"] += (60.0 if isinstance(r, Exception) else r["minutes"]) / 60 * H100_USD_H
            stop = "ALL shards in this round failed - run: logs"; break
        for r in rows:
            mins = 60.0 if isinstance(r, Exception) else r["minutes"]
            c = mins / 60 * H100_USD_H; st["spent"] += c; st["shard_cost"].append(c)
        save("harvesting")
        new = harvest.remote(tag, list(st["seqs"].values()))
        for d in new:
            st["seqs"][d["design_id"]] = d["sequence"]; st["pending"].append(d["design_id"])
        while st["pending"]:
            process_pending(); save("scored")
            if n_good() >= goal: break
    save("done", stop)
    return "STOPPING: %s | shards %d | good %d | est $%.2f" % (stop, st["launched"], n_good(), st["spent"])


@app.function(image=cpu_image, volumes={"/results": results}, timeout=600)
def read_progress(tag: str) -> str:
    results.reload(); p = pathlib.Path("/results") / tag / "progress.json"
    return p.read_text() if p.exists() else "no progress file yet"

@app.function(image=cpu_image, volumes={"/results": results}, timeout=600)
def read_file(tag: str, name: str) -> bytes:
    results.reload(); p = pathlib.Path("/results") / tag / name
    return p.read_bytes() if p.exists() else b""

@app.function(image=cpu_image, volumes={"/results": results}, timeout=900)
def tail_logs(tag: str, n: int = 25) -> str:
    import glob
    results.reload(); out = []
    for f in sorted(glob.glob("/results/" + tag + "/shard*/campaign.log")):
        out.append("--- " + f.split("/")[-2]); out.extend(open(f, errors="ignore").read().splitlines()[-n:])
    return "\n".join(out) or "no shard logs yet"
