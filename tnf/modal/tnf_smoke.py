"""Cheap pipeline check (~10 min, well under $1): Boltz setup + refold + score of 1 known binder on BOTH species.
The binder is only a positive control - it is NOT part of your submission."""
import modal
APP = "pdfa-tnf-dual"
CTL = "EEEEILERLKEIEKRLWAFLEKFWERHPDANSVAYEFLWDQLVKLHHMLKNPEEYSIEEILEILRKILEALKHILEKYGLEKEPEYQELIKEIEEVLE"
F = lambda n: modal.Function.from_name(APP, n)
print("1/2 Boltz setup (weights + target MSAs) ...")
w = F("warm").remote(); print(w)
if "FAILED" in w: raise SystemExit("SETUP FAILED - send me the text above")
print("2/2 refold + score on human and mouse ...")
recs = [F("refold").remote(CTL, sp, "ctl") for sp in ("human", "mouse")]
keep = ("species", "iptm", "ipsae", "contact_pae", "binder_plddt", "chains_contacted", "groove",
        "hot_frac", "cold_frac", "struct_pass", "his_log10_fold", "his_shifts", "pass", "error", "ph_error")
for r in F("score_one").map(recs, [{}, {}]):
    print({k: r[k] for k in keep if k in r})
print("SMOKE DONE")
