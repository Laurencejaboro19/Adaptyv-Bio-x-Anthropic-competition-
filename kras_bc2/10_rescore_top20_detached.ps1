# Block 10: detached rescore of the top 20. Appends fold_n / rescore_ctl / rescore_start / rescore_get to kras_app.py (skipped if already there;
# the kras_app.py in this folder already contains them), then refolds each design alone on G12D, WT and peptide-free with 8 samples inside Modal.
# Ran 2026-11-01: 20 designs, 60 folds, 8.8 min, est compute $2.72. Safe to close PowerShell after 'RESCORE STARTED'.
$f = "C:\Users\HP\Downloads\kras_modal\kras_app.py"
if (-not (Select-String -Path $f -Pattern "def fold_n" -Quiet)) {
$addA = @'


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
'@
Add-Content -Path $f -Value $addA -Encoding utf8
Write-Host "part A added"
} else { Write-Host "part A already in file" }
if (-not (Select-String -Path $f -Pattern "def rescore_ctl" -Quiet)) {
$addB = @'


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
'@
Add-Content -Path $f -Value $addB -Encoding utf8
Write-Host "part B added"
} else { Write-Host "part B already in file" }
Set-Location "C:\Users\HP\Downloads\kras_modal"
python -m modal run --detach kras_app.py::rescore_start --csv-path "C:\Users\HP\Downloads\kras_final\kras_submit_top20_v3_no_diversity.csv" --tag kras_rescore1 --samples 8 --top 20 --max-usd 15
