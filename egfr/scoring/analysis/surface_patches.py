"""Surface-patch analysis of EGFR domain III (human vs mouse threaded model).

For every solvent-exposed residue as a patch center, collect exposed residues within
PATCH_R of its CB and score: conserved His, charged residues, substitutions, glycan
proximity, hydrophobics, convexity.
"""
import sys, json, math
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.SASA import ShrakeRupley

HUMAN = "/Users/mike/Code/Raptorbind/inputs/targets/Adaptyv EGFR/egfr_d3ext_human.pdb"
MOUSE = "/Users/mike/Code/Raptorbind/inputs/targets/Adaptyv EGFR/egfr_d3ext_mouse.pdb"
PATCH_R = float(sys.argv[1]) if len(sys.argv) > 1 else 11.0

# Tien et al. 2013 theoretical max SASA
MAXSASA = {'ALA':129,'ARG':274,'ASN':195,'ASP':193,'CYS':167,'GLN':225,'GLU':223,'GLY':104,
           'HIS':224,'ILE':197,'LEU':201,'LYS':236,'MET':224,'PHE':240,'PRO':159,'SER':155,
           'THR':172,'TRP':285,'TYR':263,'VAL':174}
ONE = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V'}

p = PDBParser(QUIET=True)
hs = p.get_structure("h", HUMAN)[0]["A"]
ms = p.get_structure("m", MOUSE)[0]["A"]
sr = ShrakeRupley()
sr.compute(hs, level="R"); sr.compute(ms, level="R")

hseq = {r.id[1]: ONE[r.get_resname()] for r in hs}
mseq = {r.id[1]: ONE[r.get_resname()] for r in ms}
subs = {i for i in hseq if hseq[i] != mseq[i]}
# N-glycosylation sequons N-X-S/T (X != P)
def sequons(seq):
    ids = sorted(seq)
    out = set()
    for k, i in enumerate(ids[:-2]):
        if seq[i] == 'N' and seq[ids[k+1]] != 'P' and seq[ids[k+2]] in 'ST':
            out.add(i)
    return out
hgly, mgly = sequons(hseq), sequons(mseq)

rsasa_h = {r.id[1]: r.sasa / MAXSASA[r.get_resname()] for r in hs}
rsasa_m = {r.id[1]: r.sasa / MAXSASA[r.get_resname()] for r in ms}
EXPOSED = 0.20
exposed = {i for i in hseq if rsasa_h[i] >= EXPOSED or rsasa_m[i] >= EXPOSED}

def cb(res):
    return res["CB"].coord if "CB" in res else res["CA"].coord
coords = {r.id[1]: cb(r) for r in hs}
ca = {r.id[1]: r["CA"].coord for r in hs}
allca = np.array(list(ca.values()))
atoms = [a for a in hs.get_atoms()]
ns = NeighborSearch(atoms)

# glycan-shadow: residues with any atom within 8 A of a sequon Asn ND2
gly_shadow = set()
for g in hgly | mgly:
    nd2 = hs[g]["ND2"] if "ND2" in hs[g] else hs[g]["CA"]
    for a in ns.search(nd2.coord, 8.0):
        gly_shadow.add(a.get_parent().id[1])

print(f"human sequons (N-glycan): {sorted(hgly)}  mouse: {sorted(mgly)}")
print(f"substitutions ({len(subs)}): {sorted(subs)}")
print(f"His human: {[i for i in hseq if hseq[i]=='H']}  His mouse: {[i for i in mseq if mseq[i]=='H']}")
print(f"exposed residues (rSASA>={EXPOSED}): {len(exposed)} of {len(hseq)}")
print("per-residue His exposure:")
for i in sorted(i for i in hseq if hseq[i]=='H' or mseq[i]=='H'):
    print(f"  {hseq[i]}{i}{'' if i not in subs else '->'+mseq[i]}  rSASA_h={rsasa_h[i]:.2f} rSASA_m={rsasa_m[i]:.2f} B={hs[i]['CA'].bfactor:.0f}")

rows = []
for c in sorted(exposed):
    members = [j for j in exposed if np.linalg.norm(coords[j] - coords[c]) <= PATCH_R]
    if len(members) < 8:
        continue
    cons = [j for j in members if j not in subs]
    his_cons = [j for j in members if hseq[j] == 'H' and mseq[j] == 'H']
    his_h = [j for j in members if hseq[j] == 'H']
    his_m = [j for j in members if mseq[j] == 'H']
    charged = [j for j in members if hseq[j] in 'KRDE']
    basic = [j for j in members if hseq[j] in 'KR']
    acidic = [j for j in members if hseq[j] in 'DE']
    hyd = [j for j in members if hseq[j] in 'AILMFVWY']
    sub_in = [j for j in members if j in subs]
    gly_in = [j for j in members if j in gly_shadow]
    conv = int(np.sum(np.linalg.norm(allca - ca[c], axis=1) <= 10.0))
    bmean = float(np.mean([hs[j]['CA'].bfactor for j in members]))
    sasa = float(sum(hs[j].sasa for j in members))
    rows.append(dict(center=c, center_aa=hseq[c], n=len(members), his_cons=len(his_cons),
                     his_cons_ids=his_cons, his_h=len(his_h), his_m=len(his_m),
                     charged=len(charged), basic=len(basic), acidic=len(acidic), hyd=len(hyd),
                     subs=len(sub_in), subs_ids=sub_in, gly=len(gly_in), conv=conv, bmean=bmean,
                     sasa=sasa, members=members))

# Score: want conserved His (2-4), few substitutions, few glycan-shadowed, some charge, some hydrophobics
for r in rows:
    r['score'] = (min(r['his_cons'], 3) * 2.0 - r['subs'] * 1.5 - r['gly'] * 0.5
                  + 0.3 * min(r['charged'], 6) + 0.3 * min(r['hyd'], 6) - 0.05 * max(r['conv'] - 25, 0))
rows.sort(key=lambda r: -r['score'])
print(f"\nTop patches (radius {PATCH_R} A around CB, exposed residues only):")
print(f"{'center':>7} {'n':>3} {'HisC':>4} {'HisH':>4} {'HisM':>4} {'chg':>3} {'bas':>3} {'aci':>3} {'hyd':>3} {'sub':>3} {'gly':>3} {'conv':>4} {'B':>5} {'score':>6}  his_ids / sub_ids")
for r in rows[:25]:
    print(f"{r['center_aa']}{r['center']:>6} {r['n']:>3} {r['his_cons']:>4} {r['his_h']:>4} {r['his_m']:>4} {r['charged']:>3} {r['basic']:>3} {r['acidic']:>3} {r['hyd']:>3} {r['subs']:>3} {r['gly']:>3} {r['conv']:>4} {r['bmean']:>5.0f} {r['score']:>6.2f}  {r['his_cons_ids']} / {r['subs_ids']}")
json.dump({"rows": rows, "rsasa_h": rsasa_h, "rsasa_m": rsasa_m, "subs": sorted(subs),
           "hgly": sorted(hgly), "mgly": sorted(mgly), "gly_shadow": sorted(gly_shadow),
           "hseq": hseq, "mseq": mseq},
          open("/private/tmp/claude-501/-Users-mike-Code-Raptorbind/cf683a21-06b2-4caa-9b02-fd5de9766761/scratchpad/analysis/patches.json", "w"), default=lambda o: int(o) if isinstance(o, (np.integer,)) else float(o))
