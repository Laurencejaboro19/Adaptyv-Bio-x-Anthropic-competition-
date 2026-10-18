$work = "$env:USERPROFILE\Downloads\glp1r_refold"
$ext  = "$work\extracted"
New-Item -ItemType Directory -Force -Path $ext | Out-Null

# 1) which designs passed contact
$pass = Import-Csv "$work\contact_check.csv" | Where-Object { $_.PASS -match '^(True|true|1|yes)$' } | ForEach-Object { $_.design }
"PASS designs: " + ($pass -join ", ")

# 2) extract only those archives
foreach ($n in $pass) {
    $a = "$work\archives\$n.tar.gz"
    $d = "$ext\$n"
    if (-not (Test-Path "$d\prediction\sample_0_pae.npz")) {
        New-Item -ItemType Directory -Force -Path $d | Out-Null
        tar -xzf $a -C $d
    }
}

# 3) deps
python -m pip install --quiet ipsae numpy

# 4) scorer
@'
import sys, os, csv, tempfile, traceback
import numpy as np
from ipsae import core
core.IpsaeCalculator._save_outputs = lambda *a, **k: None
work = sys.argv[1]
names = sys.argv[2:]
cc = {r["design"]: r for r in csv.DictReader(open(os.path.join(work, "contact_check.csv")))}
out = []
for name in names:
    d = os.path.join(work, "extracted", name, "prediction")
    cif = os.path.join(d, "sample_0_predicted_structure.cif")
    npz = os.path.join(d, "sample_0_pae.npz")
    try:
        z = np.load(npz)
        if "pae" in z.files:
            pae = z["pae"]
        else:
            c = [k for k in z.files if z[k].ndim == 2 and z[k].shape[0] == z[k].shape[1]]
            if not c:
                raise ValueError("no square 2D array in npz; keys=%s" % z.files)
            pae = z[c[0]]
        tmp = os.path.join(tempfile.mkdtemp(), "x_pae.npz")
        np.savez(tmp, pae=pae)
        r = core.IpsaeCalculator(15, 10).calculate(tmp, cif)["ipsae_scores"]
        m = r["ipsae_d0chn_asym"]
        chains = sorted({str(a) for a in m} | {str(b) for a in m for b in m[a]})
        rec = "A" if "A" in chains else chains[0]
        bnd = [x for x in chains if x != rec][0]
        ab = float(m[rec][bnd]); ba = float(m[bnd][rec])
        row = {"design": name, "ipsae_A_to_B": round(ab, 4), "ipsae_B_to_A": round(ba, 4), "ipsae_max": round(max(ab, ba), 4)}
    except Exception as e:
        print("FAILED", name, type(e).__name__, e); traceback.print_exc()
        row = {"design": name, "ipsae_A_to_B": "", "ipsae_B_to_A": "", "ipsae_max": ""}
    c = cc.get(name, {})
    for k in ("hotspot_frac", "hotspots_hit", "coldspots_hit", "n_binder_contacts", "cterm_in_interface"):
        row[k] = c.get(k, "")
    out.append(row)
out.sort(key=lambda r: -(r["ipsae_max"] if r["ipsae_max"] != "" else -1))
fn = os.path.join(work, "ipsae_pass4.csv")
with open(fn, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
print("\nwrote", fn)
for r in out: print(r)
'@ | Set-Content -Encoding UTF8 "$work\score_pass.py"

# 5) run
python "$work\score_pass.py" $work @pass
