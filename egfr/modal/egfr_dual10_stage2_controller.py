import json, os, csv, io, pathlib, time, modal

# ============== BINDCRAFT2 SETTINGS ========================================
HOTSPOTS  = "406,409,410,411"
COLDSPOTS = "328,337,353,388,389,420,443,468,471,473,504"
BINDER_LENGTHS       = [80, 120]
AA_BIAS              = {"C": 0, "H": 2.0, "D": 1.3, "E": 1.3, "K": 0.8, "R": 0.8}
FORCED_TARGETING     = False
MIN_HOTSPOT_CONTACT  = 0.5
MAX_COLDSPOT_CONTACT = 0.1
TERMINI_ACCESSIBLE   = False
SEQUENCE_CANDIDATES  = 30
DESPERATION_TRAJ     = 40
MAX_TRAJ_PER_SHARD   = 30
N_FINAL_DESIGNS      = 999
INTERFACE_RESIDUES   = 12
WORKERS_PER_GPU      = 3
BASE_SEED            = 7000
DUAL_TARGET          = True

# ============== CONTROLLER =================================================
GOAL = 10
SHARDS_PER_ROUND = 4
MAX_SHARDS = 16
BUDGET_USD = 150.0
DEADLINE_HOURS = 10.0
REFOLD_PER_BACKBONE = 3
GATE = {"iptm": 0.80, "ipsae": 0.70, "ipae": 0.35}
H100_USD_H, L4_USD_H = 3.95, 0.80

HUMAN, MOUSE = "egfr_d3ext_human.pdb", "egfr_d3ext_mouse.pdb"
HERE = pathlib.Path(__file__).resolve().parent

app = modal.App("pdfa-egfr-dual")
weights = modal.Volume.from_name("pdfa-bc2-weights", create_if_missing=True)
results = modal.Volume.from_name("pdfa-bc2-results", create_if_missing=True)
bcache  = modal.Volume.from_name("pdfa-boltz-cache", create_if_missing=True)
BC2_COMMIT = "301efdd1937fb963cc40a5b0ecc1bc2f9b2b2d15"

bc2_image = (
    modal.Image.from_registry("ubuntu:24.04", add_python=None)
    .apt_install("python3","python3-pip","python3-dev","build-essential","ca-certificates","git")
    .env({"PIP_BREAK_SYSTEM_PACKAGES":"1","PIP_NO_CACHE_DIR":"1","PYTHONUNBUFFERED":"1"})
    .run_commands("ln -sf /usr/bin/python3 /usr/local/bin/python")
    .run_commands("git clone https://github.com/PacesaLab/BindCraft2.git /opt/bindcraft && "
                  "git -C /opt/bindcraft checkout " + BC2_COMMIT)
    .run_commands("python3 -m pip install -e '/opt/bindcraft[cuda13]'")
    .run_commands("python3 -m pip install 'gemmi==0.6.5' 'ipsae==1.0.1' numpy")
    .run_commands("python3 -c \"import nvidia, pathlib; print('\\n'.join(sorted(str(p) "
                  "for r in nvidia.__path__ for p in pathlib.Path(r).glob('*/lib'))))\" "
                  "> /etc/ld.so.conf.d/bindcraft-cuda.conf && ldconfig && ldconfig -p | grep -q libcupti")
    .run_commands("cd / && python3 -m bindcraft.selfcheck cuda13 --shipped-only")
    .env({"NVIDIA_VISIBLE_DEVICES":"all","NVIDIA_DRIVER_CAPABILITIES":"compute,utility"})
    .add_local_file(HERE / HUMAN, "/adapter/" + HUMAN, copy=True)
    .add_local_file(HERE / MOUSE, "/adapter/" + MOUSE, copy=True)
)

boltz_image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git","build-essential","curl")
    .pip_install("boltz==2.2.1","ipsae==1.0.1","biopython","numpy")
    .env({"PYTHONUNBUFFERED":"1","BOLTZ_CACHE":"/cache"})
    .add_local_file(HERE / HUMAN, "/adapter/" + HUMAN, copy=True)
    .add_local_file(HERE / MOUSE, "/adapter/" + MOUSE, copy=True)
)

cpu_image = (modal.Image.debian_slim(python_version="3.11")
             .pip_install("propka","biopython","numpy","pandas")
             .env({"PYTHONUNBUFFERED":"1"}))


def campaign(seed, project):
    def tgt(n, p):
        return {"name": n, "target_path": "/adapter/" + p, "chains": "A",
                "hotspots": HOTSPOTS, "coldspots": COLDSPOTS,
                "objective": "target", "weight": 1}
    targets = [tgt("EGFR_d3_human", HUMAN)]
    if DUAL_TARGET:
        targets.append(tgt("EGFR_d3_mouse", MOUSE))
    return {"campaign_name":"EGFR_d3_H409_dual","modality":"binder","project_folder":project,
            "campaign_seed":seed,"number_of_final_designs":N_FINAL_DESIGNS,
            "max_trajectories":MAX_TRAJ_PER_SHARD,"desperation_trajectories":DESPERATION_TRAJ,
            "sequence_candidates":SEQUENCE_CANDIDATES,"targets":targets,
            "binder_lengths":BINDER_LENGTHS,"forced_targeting":FORCED_TARGETING,
            "min_hotspot_contact_final":MIN_HOTSPOT_CONTACT,
            "max_coldspot_contact_final":MAX_COLDSPOT_CONTACT,
            "termini_accessible":TERMINI_ACCESSIBLE,"aa_bias":AA_BIAS,
            "filters":{"Interface_H_Count":{"threshold":0,"higher":True,"mandatory":False},
                       "Interface_Residues":{"threshold":INTERFACE_RESIDUES,"higher":True}},
            "workers_per_gpu":WORKERS_PER_GPU}

T31 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

def seq_of(p):
    return "".join(T31[l[17:20]] for l in open(p).read().splitlines()
                   if l.startswith("ATOM") and l[12:16].strip() == "CA")


@app.function(image=bc2_image, gpu="H100", timeout=60*180,
              volumes={"/weights": weights, "/results": results},
              retries=modal.Retries(max_retries=1, backoff_coefficient=1.0))
def design_shard(i: int, tag: str) -> dict:
    import subprocess, shutil
    seed = BASE_SEED + i
    work = pathlib.Path("/tmp/bc2_" + tag + "_" + str(seed)); proj = work / "project"
    proj.mkdir(parents=True, exist_ok=True)
    cache = pathlib.Path("/weights/compile_cache"); cache.mkdir(parents=True, exist_ok=True)
    if not (proj / "compile_cache").exists():
        (proj / "compile_cache").symlink_to(cache)
    cfg = campaign(seed, str(proj))
    cpath = work / "campaign.json"; cpath.write_text(json.dumps(cfg, indent=2))
    out = pathlib.Path("/results") / tag / ("shard" + format(i,"02d") + "_seed" + str(seed))
    out.mkdir(parents=True, exist_ok=True)
    (out / "campaign.json").write_text(cpath.read_text()); results.commit()
    env = dict(os.environ, BINDCRAFT_WORKERS_PER_GPU=str(WORKERS_PER_GPU),
               BINDCRAFT_WEIGHTS="/weights")
    t0 = time.time()
    with (out / "campaign.log").open("w") as log:
        p = subprocess.Popen(["bindcraft","design",str(cpath)], cwd=str(work),
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             text=True, bufsize=1, env=env)
        for line in p.stdout:
            log.write(line); print("[" + format(i,"02d") + "] " + line, end="")
        rc = p.wait()
    mins = (time.time() - t0) / 60
    for root, _, files in os.walk(proj):
        if "compile_cache" in root: continue
        for f in files:
            if f.endswith((".csv",".fasta",".json",".pdb")):
                dst = out / "project" / os.path.relpath(root, proj) / f
                dst.parent.mkdir(parents=True, exist_ok=True)
                try: shutil.copy(pathlib.Path(root)/f, dst)
                except Exception: pass
    (out/"status.json").write_text(json.dumps({"shard":i,"seed":seed,"rc":rc,
                                               "minutes":round(mins,1)}))
    results.commit(); weights.commit()
    return {"shard": i, "rc": rc, "minutes": round(mins, 1)}


@app.function(image=boltz_image, gpu="L4", timeout=3600, volumes={"/cache": bcache})
def refold(binder_seq: str, species: str, design_id: str) -> dict:
    import subprocess, tempfile, glob
    tgt = seq_of("/adapter/" + (HUMAN if species == "human" else MOUSE))
    work = pathlib.Path(tempfile.mkdtemp()); y = work / "c.yaml"
    y.write_text("version: 1\nsequences:\n"
                 + "  - protein: {id: A, sequence: " + tgt + "}\n"
                 + "  - protein: {id: B, sequence: " + binder_seq + ", msa: empty}\n")
    out = work / "out"
    r = subprocess.run(["boltz","predict",str(y),"--out_dir",str(out),"--model","boltz2",
                        "--accelerator","gpu","--devices","1","--recycling_steps","1",
                        "--sampling_steps","50","--diffusion_samples","1",
                        "--output_format","mmcif","--cache","/cache","--use_msa_server"],
                       capture_output=True, text=True)
    bcache.commit()
    cif = sorted(glob.glob(str(out/"**"/"*model_0.cif"), recursive=True))
    cnf = sorted(glob.glob(str(out/"**"/"confidence*model_0.json"), recursive=True))
    rec = {"design_id": design_id, "species": species, "rc": r.returncode,
           "sequence": binder_seq}
    if cnf: rec.update(json.load(open(cnf[0])))
    if cif: rec["cif_text"] = open(cif[0]).read()
    else:   rec["log"] = (r.stdout + r.stderr)[-1000:]
    return rec


SCORER = r"""
import sys, json, subprocess, tempfile, os, copy
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.SASA import ShrakeRupley
ANCHOR = 105
COO = {("ASP","OD1"),("ASP","OD2"),("GLU","OE1"),("GLU","OE2")}
CAT = {("LYS","NZ"),("ARG","NH1"),("ARG","NH2"),("ARG","NE")}

def pkas(model, chains):
    class S(Select):
        def accept_chain(s, c): return c.id in chains
        def accept_residue(s, r): return r.id[0] == " "
    with tempfile.TemporaryDirectory() as td:
        io = PDBIO(); io.set_structure(model); io.save(td + "/x.pdb", S())
        subprocess.run(["propka3","x.pdb"], cwd=td, capture_output=True, timeout=900)
        out, p = {}, td + "/x.pka"
        if not os.path.exists(p): return out
        for line in open(p):
            f = line.split()
            if len(f) >= 4 and len(f[0]) == 3 and f[1].lstrip("-").isdigit():
                try: out[(f[2], int(f[1]), f[0])] = float(f[3])
                except ValueError: pass
        return out

def wyman(pc, pf, lo=6.4, hi=7.4):
    k = set(pc) & set(pf)
    if not k: return None, 0
    fp = lambda pk, ph: 1.0/(1.0 + 10.0**(ph - pk))
    g = np.linspace(lo, hi, 41)
    dn = np.array([sum(fp(pc[x], p) - fp(pf[x], p) for x in k) for p in g])
    tz = getattr(np, "trapezoid", None) or np.trapz
    return float(tz(dn, g)), len(k)

def score(cif, T="A", B="B"):
    m = MMCIFParser(QUIET=True).get_structure("x", cif)[0]
    o = {}
    o["binder_length"] = len([r for r in m[B] if r.id[0] == " "])
    his = m[T].child_dict.get((" ", ANCHOR, " "))
    if his is not None and his.get_resname() == "HIS":
        ring = np.array([his[a].coord for a in ("ND1","NE2") if a in his])
        ba, bc, part = 99.0, 99.0, None
        for r in m[B]:
            for at in r:
                k = (r.get_resname(), at.get_id())
                d = float(np.linalg.norm(at.coord - ring, axis=-1).min())
                if k in COO and d < ba: ba, part = d, r.get_resname() + str(r.id[1])
                if k in CAT: bc = min(bc, d)
        sr = ShrakeRupley()
        free = copy.deepcopy(m)
        for c in [c.id for c in free]:
            if c != T: free.detach_child(c)
        sr.compute(free, level="R"); sr.compute(m, level="R")
        a = free[T][(" ", ANCHOR, " ")].sasa; b = m[T][(" ", ANCHOR, " ")].sasa
        o.update({"H409_carboxylate_partner": part,
                  "H409_d_carboxylate_A": round(ba,2) if part else None,
                  "H409_d_lysarg_A": round(bc,2) if bc < 99 else None,
                  "H409_burial_on_binding": round((a-b)/a,3) if a > 1 else None,
                  "H409_contacted": bool(a - b > 1)})
        tot = sum(max(0.0, free[T][r.id].sasa - m[T][r.id].sasa)
                  for r in m[T] if r.id[0] == " " and r.id in free[T])
        o["interface_buried_area_A2"] = round(tot, 1)
    pc = pkas(m, {T, B}); pf = {}
    pf.update(pkas(m, {T})); pf.update(pkas(m, {B}))
    lg, n = wyman(pc, pf)
    o["log10_Ka_6.4_over_7.4"] = round(lg,3) if lg is not None else None
    o["fold_Ka_6.4_over_7.4"] = round(10**lg,3) if lg is not None else None
    o["n_titratable_groups"] = n
    k = (T, ANCHOR, "HIS")
    if k in pc and k in pf:
        o["H409_pka_free"] = round(pf[k],2); o["H409_pka_complex"] = round(pc[k],2)
        o["H409_dpka"] = round(pc[k]-pf[k],2)
    return o

if __name__ == "__main__":
    print(json.dumps(score(sys.argv[1])))
"""

@app.function(image=cpu_image, cpu=2.0, timeout=1800)
def score_one(cif_text: str, design_id: str, species: str) -> dict:
    import subprocess, tempfile
    d = tempfile.mkdtemp()
    pathlib.Path(d + "/s.py").write_text(SCORER)
    pathlib.Path(d + "/c.cif").write_text(cif_text)
    r = subprocess.run(["python", d + "/s.py", d + "/c.cif"],
                       capture_output=True, text=True, timeout=1500)
    rec = {"design_id": design_id, "species": species}
    for line in r.stdout.splitlines():
        if line.startswith("{"):
            try: rec.update(json.loads(line))
            except Exception: pass
    if r.returncode: rec["error"] = (r.stdout + r.stderr)[-300:]
    return rec


@app.function(image=cpu_image, volumes={"/results": results}, timeout=900)
def harvest(tag: str, already: list) -> list:
    import glob, csv as _csv, re
    results.reload()
    seen, out = set(already), []
    pat = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{50,250}$")
    def add(did, s):
        if pat.match(s) and s not in seen:
            seen.add(s); out.append({"design_id": did, "sequence": s})
    for f in glob.glob("/results/" + tag + "/**/*.fasta", recursive=True):
        name, buf = None, []
        for line in open(f):
            if line.startswith(">"):
                if name: add(name, "".join(buf))
                name, buf = line[1:].split()[0], []
            else: buf.append(line.strip())
        if name: add(name, "".join(buf))
    for f in glob.glob("/results/" + tag + "/**/*.csv", recursive=True):
        try: rows = list(_csv.DictReader(open(f)))
        except Exception: continue
        for i, r in enumerate(rows):
            s = next((str(r[k]).strip() for k in r if k and k.lower() in
                      ("sequence","binder_seq","seq","binder_sequence")), "")
            did = next((str(r[k]) for k in r if k and "design" in k.lower()),
                       pathlib.Path(f).stem + "_" + str(i))
            add(did, s)
    by_bb, keep = {}, []
    for d in out:
        bb = d["design_id"].rsplit("_candidate", 1)[0]
        by_bb.setdefault(bb, []).append(d)
    for bb, lst in by_bb.items():
        keep.extend(lst[:REFOLD_PER_BACKBONE])
    print("harvested " + str(len(keep)) + " new sequences")
    return keep


def passes(rec):
    try:
        return (float(rec.get("iptm", 0)) >= GATE["iptm"]
                and float(rec.get("ipsae", rec.get("iptm", 0))) >= GATE["ipsae"]
                and float(rec.get("ipae", 0)) <= GATE["ipae"]
                and bool(rec.get("H409_contacted", True)))
    except Exception:
        return False


# =================== THE CONTROLLER — runs IN THE CLOUD =====================
@app.function(image=cpu_image, cpu=2.0, timeout=12*3600,
              volumes={"/results": results})
def controller(tag: str, goal: int, max_shards: int, per_round: int,
               budget: float, deadline_h: float) -> str:
    t0 = time.time(); spent = 0.0
    done_seqs, cross_ok, seqs, launched, rnd = [], {}, {}, 0, 0
    base = pathlib.Path("/results") / tag
    base.mkdir(parents=True, exist_ok=True)

    def progress(state, stop=None):
        (base / "progress.json").write_text(json.dumps({
            "state": state, "round": rnd, "shards_launched": launched,
            "designs_refolded": len(done_seqs), "good": len(cross_ok),
            "goal": goal, "est_spend_usd": round(spent, 2),
            "elapsed_h": round((time.time()-t0)/3600, 2),
            "stop_reason": stop,
            "updated": time.strftime("%Y-%m-%d %H:%M:%S")}, indent=1))
        results.commit()

    progress("starting")
    while True:
        rnd += 1
        left_t = deadline_h - (time.time() - t0) / 3600
        if len(cross_ok) >= goal: stop = "goal reached: " + str(len(cross_ok)); break
        if launched >= max_shards: stop = "shard ceiling reached"; break
        if spent >= budget:        stop = "budget spent"; break
        if left_t <= 0.5:          stop = "deadline reached"; break

        k = min(per_round, max_shards - launched)
        idx = list(range(launched, launched + k))
        print("=== round " + str(rnd) + ": " + str(k) + " shards | have "
              + str(len(cross_ok)) + "/" + str(goal) + " | $" + format(spent,".2f"))
        progress("designing")
        rows = list(design_shard.map(idx, [tag]*k)); launched += k
        mins = sum(r["minutes"] for r in rows); spent += mins/60 * H100_USD_H
        progress("harvesting")

        new = harvest.remote(tag, done_seqs)
        if not new:
            print("no new sequences"); progress("idle"); continue
        for d in new:
            done_seqs.append(d["sequence"]); seqs[d["design_id"]] = d["sequence"]
        jobs = [(d["sequence"], sp, d["design_id"]) for d in new for sp in ("human","mouse")]
        progress("refolding")
        folds = list(refold.starmap(jobs))
        spent += len(folds) * (3.5/60) * L4_USD_H
        good = [f for f in folds if "cif_text" in f]

        rbase = base / "refold"
        for f in folds:
            d = rbase / f["design_id"] / f["species"]; d.mkdir(parents=True, exist_ok=True)
            if "cif_text" in f: (d / "complex.cif").write_text(f["cif_text"])
            (d / "confidence.json").write_text(json.dumps(
                {k2: v2 for k2, v2 in f.items() if k2 != "cif_text"}, indent=1))
        results.commit()

        progress("scoring")
        scores = list(score_one.starmap(
            [(f["cif_text"], f["design_id"], f["species"]) for f in good]))
        by = {}
        for f in good:
            by[(f["design_id"], f["species"])] = {k2: v2 for k2, v2 in f.items()
                                                  if k2 != "cif_text"}
        for s in scores:
            by.setdefault((s["design_id"], s["species"]), {}).update(s)
        for d in {k2[0] for k2 in by}:
            h, m = by.get((d,"human")), by.get((d,"mouse"))
            if h and m and passes(h) and passes(m):
                cross_ok[d] = {"human": h, "mouse": m}
        print("good designs now " + str(len(cross_ok)) + "/" + str(goal))
        progress("scored")

    L = "log10_Ka_6.4_over_7.4"
    rowsout = []
    for d, v in cross_ok.items():
        r = {"design_id": d, "sequence": seqs.get(d, ""),
             "length": v["human"].get("binder_length")}
        for sp in ("human","mouse"):
            for k2, v2 in v[sp].items():
                if k2 not in ("design_id","species","sequence","cif_text","log"):
                    r[sp + "_" + k2] = v2
        hL = r.get("human_" + L) or 0; mL = r.get("mouse_" + L) or 0
        r["tier"] = "A" if (hL > 0 and mL > 0) else ("B" if hL > 0 else "C")
        r["submit"] = "YES" if r["tier"] == "A" else "hold"
        rowsout.append(r)
    rowsout.sort(key=lambda r: (r["tier"], -(r.get("human_" + L) or -99)))
    for i, r in enumerate(rowsout, 1): r["rank"] = i
    fixed = ["rank","tier","submit","design_id","length","sequence"]
    keys = fixed + sorted({k2 for r in rowsout for k2 in r} - set(fixed))
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=keys); w.writeheader(); w.writerows(rowsout)
    (base / "designs.csv").write_text(buf.getvalue())
    progress("done", stop)
    return ("STOPPING: " + stop + " | shards " + str(launched) + " | good "
            + str(len(cross_ok)) + " | $" + format(spent, ".2f"))


@app.function(image=cpu_image, volumes={"/results": results}, timeout=600)
def read_progress(tag: str) -> str:
    results.reload()
    p = pathlib.Path("/results") / tag / "progress.json"
    return p.read_text() if p.exists() else "no progress file yet"


@app.function(image=cpu_image, volumes={"/results": results}, timeout=600)
def read_csv_bytes(tag: str) -> bytes:
    results.reload()
    p = pathlib.Path("/results") / tag / "designs.csv"
    return p.read_bytes() if p.exists() else b""


# ----------------------------- entrypoints ---------------------------------
@app.local_entrypoint()
def start(tag: str = "dual10", goal: int = GOAL, max_shards: int = MAX_SHARDS,
          per_round: int = SHARDS_PER_ROUND, budget: float = BUDGET_USD,
          deadline_h: float = DEADLINE_HOURS):
    call = controller.spawn(tag, goal, max_shards, per_round, budget, deadline_h)
    print("CONTROLLER STARTED IN THE CLOUD")
    print("  call id : " + str(call.object_id))
    print("  tag     : " + tag)
    print("\nYou can close this window and shut the laptop now.")
    print("\nTomorrow:")
    print("  python -m modal run egfr_dual10.py::status --tag " + tag)
    print("  python -m modal run egfr_dual10.py::fetch  --tag " + tag)


@app.local_entrypoint()
def status(tag: str = "dual10"):
    print(read_progress.remote(tag))


@app.local_entrypoint()
def fetch(tag: str = "dual10"):
    data = read_csv_bytes.remote(tag)
    if not data:
        print("no CSV yet - run status to see where it is")
        return
    pathlib.Path(tag + "_designs.csv").write_bytes(data)
    print("wrote " + tag + "_designs.csv (" + str(len(data)) + " bytes)")
