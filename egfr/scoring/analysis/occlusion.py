"""Which domain-III surface residues are occluded in the full receptor / by EGF / by known antibodies?
Superpose our domain III onto each EGFR chain in 1NQL (tethered), 1IVO (EGF-bound extended), 1YY9
(cetuximab Fab), 5SX4 (panitumumab Fab) via CA of residues 311-480 and measure min distance from each
of our exposed residues to anything that is NOT domain III of that chain."""
import json, numpy as np
from Bio.PDB import MMCIFParser, PDBParser, Superimposer, NeighborSearch
P = "/private/tmp/claude-501/-Users-mike-Code-Raptorbind/cf683a21-06b2-4caa-9b02-fd5de9766761/scratchpad/pdb/"
ref = PDBParser(QUIET=True).get_structure("h", "/Users/mike/Code/Raptorbind/inputs/targets/Adaptyv EGFR/egfr_d3ext_human.pdb")[0]["A"]
d = json.load(open("patches.json")); rs = {int(k): v for k, v in d["rsasa_h"].items()}
exposed = sorted(i for i in rs if rs[i] >= 0.2)
GLY = {"NAG", "BMA", "MAN", "FUC", "GAL", "NDG"}
def domain(n):
    return "I" if n <= 165 else "II" if n <= 310 else "III" if n <= 480 else "IV"
results = {}
for pid in ("1NQL", "1IVO", "1YY9", "5SX4"):
    st = MMCIFParser(QUIET=True).get_structure(pid, P + pid + ".cif"); model = st[0]
    for chain in model:
        common = [n for n in range(311, 481) if n in chain and n in ref and "CA" in chain[n]]
        if len(common) < 120: continue
        ident = np.mean([chain[n].get_resname() == ref[n].get_resname() for n in common])
        if ident < 0.9: continue
        sup = Superimposer(); sup.set_atoms([chain[n]["CA"] for n in common], [ref[n]["CA"] for n in common])
        rot, tran = sup.rotran
        others = []
        for ch in model:
            for res in ch:
                het = res.id[0].strip()
                if res.get_resname() in GLY:
                    others += [(a, "glycan") for a in res]; continue
                if het: continue
                if ch.id == chain.id:
                    if 311 <= res.id[1] <= 480: continue
                    others += [(a, f"dom{domain(res.id[1])}") for a in res]
                else:
                    others += [(a, f"chain{ch.id}") for a in res]
        ns = NeighborSearch([a for a, _ in others]); label = {a: l for a, l in others}
        occ = {}
        for i in exposed:
            best = None
            for a in ref[i]:
                c = np.dot(a.coord, rot) + tran
                for b in ns.search(c, 5.0):
                    dist = np.linalg.norm(c - b.coord)
                    if best is None or dist < best[0]: best = (dist, label[b])
            if best: occ[i] = (round(float(best[0]), 1), best[1])
        key = f"{pid}:{chain.id}"
        results[key] = {"rms": round(sup.rms, 2), "ident": round(float(ident), 3), "n_common": len(common), "occ": occ}
        print(f"{key}: CA rmsd {sup.rms:.2f} A over {len(common)} res, identity {ident:.3f}; occluded exposed residues (<5 A): {len(occ)}/{len(exposed)}")
        bylab = {}
        for i, (dd, l) in occ.items(): bylab.setdefault(l, []).append(i)
        for l, ids in sorted(bylab.items()): print(f"    {l:>10}: {ids}")
json.dump(results, open("occlusion.json", "w"))
