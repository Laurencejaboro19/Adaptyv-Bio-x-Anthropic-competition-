"""
MODULE: scripts/score_ph_switch.py
Score a target+binder complex for a bind-at-low-pH / release-at-neutral-pH histidine switch,
and compare its target footprint with a reference footprint (e.g. BindCraft2's design-time
contacts). Pure post-processing: Biopython + numpy + PROPKA3; no GPU.

DIRECTION. The switch wanted here is tighter binding at pH 6.4 than at 7.4, so every
titratable group should bind a proton MORE tightly in the complex than free. In Wyman
linkage terms, d log10 Ka / d pH = -dnH, where dnH is the number of protons taken up on
binding; integrating dnH over [6.4, 7.4] gives the predicted log10 Ka(6.4)/Ka(7.4). PROPKA3
pKa values of the complex and of each separated chain supply dnH(pH). This is an estimate
(PROPKA is empirical, the structure is a prediction), used to RANK designs, not to predict
a measured ratio.

Per interface His the geometry is also read directly, because the pKa shift is driven by
what sits across the interface: a carboxylate O within 4 A of ND1/NE2 raises the pKa (wanted);
a Lys NZ / Arg NH within 6 A lowers it (forbidden); full burial without a partner lowers it.

INDEX:
  load_structure -- mmCIF or PDB -> Bio.PDB model
  chain_residues -- ordered residues of a chain (amino acids only)
  interface_residues -- residue ids on both chains within `cutoff` A (heavy atoms)
  plddt_by_residue -- mean B-factor (pLDDT in Boltz-2 output) per residue
  his_geometry -- per interface His: partner distances, burial, verdict
  propka_pkas -- {(chain, resnum, restype): pKa} for complex and separated chains
  proton_linkage -- dnH(pH) and the integrated log10 Ka ratio over [6.4, 7.4]
  footprint_compare -- Jaccard between two sets of target residue numbers
  parse_bc2_footprint -- "N1,V3,H105" -> {305, 307, 409} with the file's first residue
  score_complex -- everything above for one complex, as a flat dict
  main -- CLI: one or many complexes -> JSON lines

USAGE:
  python scripts/score_ph_switch.py complex.cif --target A --binder B --target-offset 304 \
      --anchor 409 --bc2-footprint "T102,H105,G106,Q107" [--pae pae.npz]
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import math
import sys
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch, PDBIO, PDBParser, Select
from Bio.PDB.SASA import ShrakeRupley

# Tien et al. 2013 theoretical maximum SASA, for relative burial.
MAX_SASA = {'ALA': 129, 'ARG': 274, 'ASN': 195, 'ASP': 193, 'CYS': 167, 'GLN': 225, 'GLU': 223,
            'GLY': 104, 'HIS': 224, 'ILE': 197, 'LEU': 201, 'LYS': 236, 'MET': 224, 'PHE': 240,
            'PRO': 159, 'SER': 155, 'THR': 172, 'TRP': 285, 'TYR': 263, 'VAL': 174}
ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E',
       'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
       'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
CARBOXYLATE_O = {('ASP', 'OD1'), ('ASP', 'OD2'), ('GLU', 'OE1'), ('GLU', 'OE2')}
CATION_N = {('LYS', 'NZ'), ('ARG', 'NH1'), ('ARG', 'NH2'), ('ARG', 'NE')}
AROMATIC = {'PHE', 'TYR', 'TRP'}
# PROPKA group types that titrate; sign of the charge they carry when protonated.
BASIC_GROUPS = {'HIS', 'LYS', 'ARG', 'N+'}
ACIDIC_GROUPS = {'ASP', 'GLU', 'TYR', 'CYS', 'C-'}


def load_structure(path: str):
    p = Path(path)
    parser = MMCIFParser(QUIET=True) if p.suffix.lower() in ('.cif', '.mmcif') else PDBParser(QUIET=True)
    return parser.get_structure(p.stem, str(p))[0]


def chain_residues(model, chain_id: str) -> list:
    return [r for r in model[chain_id] if r.id[0] == ' ' and r.get_resname() in MAX_SASA]


def interface_residues(model, target: str, binder: str, cutoff: float = 4.5) -> Tuple[Set[int], Set[int]]:
    """Residue numbers (as in the file) on each chain with a heavy atom within `cutoff` A."""
    t_atoms = [a for r in chain_residues(model, target) for a in r if a.element != 'H']
    b_atoms = [a for r in chain_residues(model, binder) for a in r if a.element != 'H']
    ns = NeighborSearch(t_atoms)
    t_set, b_set = set(), set()
    for a in b_atoms:
        hits = ns.search(a.coord, cutoff)
        if hits:
            b_set.add(a.get_parent().id[1])
            t_set.update(h.get_parent().id[1] for h in hits)
    return t_set, b_set


def plddt_by_residue(model, chain_id: str) -> Dict[int, float]:
    return {r.id[1]: float(np.mean([a.bfactor for a in r])) for r in chain_residues(model, chain_id)}


def _sasa_per_residue(model, chain_ids: Iterable[str]) -> Dict[Tuple[str, int], float]:
    """Relative SASA per residue computed on a copy containing only `chain_ids`."""
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model
    s = Structure('tmp'); m = Model(0); s.add(m)
    for cid in chain_ids:
        m.add(model[cid].copy())
    ShrakeRupley().compute(m, level='R')
    out = {}
    for cid in chain_ids:
        for r in m[cid]:
            if r.get_resname() in MAX_SASA:
                out[(cid, r.id[1])] = r.sasa / MAX_SASA[r.get_resname()]
    return out


def his_geometry(model, target: str, binder: str, t_if: Set[int], b_if: Set[int],
                 cutoff_salt: float = 4.0, cutoff_cation: float = 6.0) -> List[dict]:
    """For every His in either chain's interface: what sits across the interface."""
    bound = _sasa_per_residue(model, (target, binder))
    free = {}
    free.update(_sasa_per_residue(model, (target,)))
    free.update(_sasa_per_residue(model, (binder,)))
    rows = []
    for own, other, own_if in ((target, binder, t_if), (binder, target, b_if)):
        other_atoms = [a for r in chain_residues(model, other) for a in r if a.element != 'H']
        ns = NeighborSearch(other_atoms)
        for r in chain_residues(model, own):
            if r.get_resname() != 'HIS' or r.id[1] not in own_if:
                continue
            ring_n = [r[n] for n in ('ND1', 'NE2') if n in r]
            if not ring_n:
                continue
            d_carb = d_cat = d_arom = d_polar = math.inf
            carb_partner = cat_partner = None
            for n in ring_n:
                for a in ns.search(n.coord, 8.0):
                    res = a.get_parent(); key = (res.get_resname(), a.get_id())
                    d = float(np.linalg.norm(n.coord - a.coord))
                    if key in CARBOXYLATE_O and d < d_carb:
                        d_carb, carb_partner = d, f"{ONE[res.get_resname()]}{res.id[1]}"
                    elif key in CATION_N and d < d_cat:
                        d_cat, cat_partner = d, f"{ONE[res.get_resname()]}{res.id[1]}"
                    elif res.get_resname() in AROMATIC and a.get_id() not in ('N', 'C', 'O', 'CA', 'CB') and d < d_arom:
                        d_arom = d
                    elif a.element in ('O', 'N') and key not in CARBOXYLATE_O and key not in CATION_N and d < d_polar:
                        d_polar = d
            r_free = free.get((own, r.id[1]), float('nan'))
            r_bound = bound.get((own, r.id[1]), float('nan'))
            burial = 1.0 - (r_bound / r_free) if r_free and r_free > 0 else float('nan')
            if d_cat <= cutoff_cation:
                verdict = 'bad: Lys/Arg within 6 A'
            elif d_carb <= cutoff_salt:
                verdict = 'good: carboxylate salt bridge'
            elif d_polar <= 3.5 or d_arom <= 4.5:
                verdict = 'ok: polar/aromatic partner'
            elif burial > 0.9:
                verdict = 'bad: buried without partner'
            else:
                verdict = 'neutral: solvent-facing'
            rows.append({
                'chain': own, 'resnum': r.id[1],
                'd_carboxylate': None if math.isinf(d_carb) else round(d_carb, 2), 'carboxylate': carb_partner,
                'd_cation': None if math.isinf(d_cat) else round(d_cat, 2), 'cation': cat_partner,
                'd_aromatic': None if math.isinf(d_arom) else round(d_arom, 2),
                'd_polar': None if math.isinf(d_polar) else round(d_polar, 2),
                'rsasa_free': round(r_free, 2), 'rsasa_bound': round(r_bound, 2),
                'burial_on_binding': round(burial, 2) if not math.isnan(burial) else None,
                'verdict': verdict,
            })
    return rows


class _ChainSelect(Select):
    def __init__(self, chains):
        self.chains = set(chains)

    def accept_chain(self, chain):
        return chain.id in self.chains

    def accept_residue(self, residue):
        return residue.id[0] == ' '

    def accept_atom(self, atom):
        return atom.element != 'H' and (not atom.is_disordered() or atom.get_altloc() in (' ', 'A'))


def propka_pkas(model, chains: Iterable[str]) -> Dict[Tuple[str, int, str], float]:
    """PROPKA3 pKa per titratable group for the chains given, written out as a PDB."""
    import propka.run  # noqa: WPS433 -- optional dependency, imported here on purpose
    io_ = PDBIO(); io_.set_structure(model)
    with tempfile.TemporaryDirectory() as td:
        pdb = Path(td) / 'sel.pdb'
        io_.save(str(pdb), _ChainSelect(chains))
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            mol = propka.run.single(str(pdb), optargs=['--quiet'], write_pka=False)
    out = {}
    for g in mol.conformations['AVR'].groups:
        if g.residue_type in BASIC_GROUPS | ACIDIC_GROUPS:
            out[(g.atom.chain_id, int(g.atom.res_num), g.residue_type)] = float(g.pka_value)
    return out


def proton_linkage(pk_complex: Dict, pk_free: Dict, ph_low: float = 6.4, ph_high: float = 7.4,
                   restrict: Optional[Set[Tuple[str, int]]] = None) -> dict:
    """dnH(pH) = protons bound by the complex minus by the free chains, summed over groups;
    log10 Ka(ph_low)/Ka(ph_high) = integral of dnH over [ph_low, ph_high]."""
    def frac_prot(pka, ph):
        return 1.0 / (1.0 + 10.0 ** (ph - pka))

    keys = set(pk_complex) & set(pk_free)
    if restrict is not None:
        keys = {k for k in keys if (k[0], k[1]) in restrict}

    def dnh(ph):
        return sum(frac_prot(pk_complex[k], ph) - frac_prot(pk_free[k], ph) for k in keys)

    grid = np.linspace(ph_low, ph_high, 41)
    vals = np.array([dnh(p) for p in grid])
    _trapz = getattr(np, 'trapezoid', None) or np.trapz    # numpy 2 renamed it
    log_ratio = float(_trapz(vals, grid))   # positive = tighter at ph_low
    shifts = {f"{k[0]}:{k[2]}{k[1]}": round(pk_complex[k] - pk_free[k], 2) for k in keys
              if abs(pk_complex[k] - pk_free[k]) >= 0.3}
    return {'dnH_6.4': round(float(dnh(ph_low)), 3), 'dnH_7.4': round(float(dnh(ph_high)), 3),
            'log10_Ka_ratio_6.4_over_7.4': round(log_ratio, 3), 'n_groups': len(keys),
            'pka_shifts_ge_0.3': shifts}


def footprint_compare(a: Set[int], b: Set[int]) -> dict:
    inter, union = a & b, a | b
    return {'jaccard': round(len(inter) / len(union), 3) if union else None,
            'shared': sorted(inter), 'only_a': sorted(a - b), 'only_b': sorted(b - a)}


def parse_bc2_footprint(text: str, first_resnum: int = 1) -> Set[int]:
    """BindCraft2 writes 'L325,Q384,H409' in the INPUT FILE's numbering (the PLA2G2A file
    started at 1, which is why an earlier note said it renumbers; it does not). With the
    EGFR 305-515 file the numbers are already mature numbering, so first_resnum stays 1."""
    out = set()
    for tok in (text or '').replace('/', ',').split(','):
        tok = tok.strip()
        digits = ''.join(ch for ch in tok if ch.isdigit())
        if digits:
            out.add(int(digits) + first_resnum - 1)
    return out


def _mean_contact_pae(pae_path: str, model, target: str, binder: str, t_if: Set[int], b_if: Set[int]) -> Optional[float]:
    """Mean PAE over (target interface residue, binder interface residue) pairs, both directions."""
    try:
        pae = np.load(pae_path)['pae']
    except Exception:  # noqa: BLE001
        return None
    nt = len(chain_residues(model, target)); nb = len(chain_residues(model, binder))
    if pae.shape[0] != nt + nb:
        return None
    first_t = chain_residues(model, target)[0].id[1]; first_b = chain_residues(model, binder)[0].id[1]
    ti = [i - first_t for i in sorted(t_if)]; bi = [nt + j - first_b for j in sorted(b_if)]
    if not ti or not bi:
        return None
    block = np.concatenate([pae[np.ix_(ti, bi)].ravel(), pae[np.ix_(bi, ti)].ravel()])
    return round(float(block.mean()), 2)


def score_complex(path: str, target: str = 'A', binder: str = 'B', target_offset: int = 0,
                  anchor: Optional[int] = None, bc2_footprint: str = '', pae: Optional[str] = None,
                  plddt_min: float = 70.0, cutoff: float = 4.5) -> dict:
    """Everything for one complex. Residue numbers are reported in file numbering +
    `target_offset` (304 for the EGFR 305-515 construct written from 1)."""
    model = load_structure(path)
    t_if, b_if = interface_residues(model, target, binder, cutoff)
    pl_t = plddt_by_residue(model, target); pl_b = plddt_by_residue(model, binder)
    conf_b = {i for i in b_if if pl_b.get(i, 0) >= plddt_min}
    # High-confidence footprint: target residues contacted by confident binder residues.
    hc_t = set()
    if conf_b:
        t_atoms = [a for r in chain_residues(model, target) for a in r if a.element != 'H']
        ns = NeighborSearch(t_atoms)
        for r in chain_residues(model, binder):
            if r.id[1] in conf_b:
                for a in r:
                    for h in ns.search(a.coord, cutoff):
                        if pl_t.get(h.get_parent().id[1], 0) >= plddt_min:
                            hc_t.add(h.get_parent().id[1])
    t_num = {i + target_offset for i in t_if}; hc_num = {i + target_offset for i in hc_t}
    out = {
        'complex': str(path), 'n_target_if': len(t_if), 'n_binder_if': len(b_if),
        'target_footprint': sorted(t_num), 'target_footprint_hc': sorted(hc_num),
        'binder_if_mean_plddt': round(float(np.mean([pl_b[i] for i in b_if])), 1) if b_if else None,
        'binder_mean_plddt': round(float(np.mean(list(pl_b.values()))), 1),
    }
    if pae:
        out['mean_contact_pae'] = _mean_contact_pae(pae, model, target, binder, t_if, b_if)
    if anchor is not None:
        out['anchor'] = anchor
        out['anchor_contacted'] = anchor in t_num
        out['anchor_contacted_hc'] = anchor in hc_num
    if bc2_footprint:
        ref = parse_bc2_footprint(bc2_footprint)      # already in the target file's numbering
        out['footprint_vs_bc2'] = footprint_compare(t_num, ref)
        out['footprint_hc_vs_bc2'] = footprint_compare(hc_num, ref)
    his = his_geometry(model, target, binder, t_if, b_if)
    for h in his:
        h['label'] = f"{h['chain']}:H{h['resnum'] + (target_offset if h['chain'] == target else 0)}"
    out['interface_his'] = his
    out['n_interface_his'] = len(his)
    out['n_interface_his_good'] = sum(1 for h in his if h['verdict'].startswith('good'))
    out['n_interface_his_bad'] = sum(1 for h in his if h['verdict'].startswith('bad'))
    # Anchor-specific geometry: the His409 salt bridge the switch is built on.
    if anchor is not None:
        a_rows = [h for h in his if h['chain'] == target and h['resnum'] + target_offset == anchor]
        out['anchor_geometry'] = a_rows[0] if a_rows else None
    # His-tag decoy guard: share of buried target area carried by the anchor residue.
    if anchor is not None:
        free_t = _sasa_per_residue(model, (target,)); bound_t = _sasa_per_residue(model, (target, binder))
        buried = {k: max(0.0, (free_t[k] - bound_t.get(k, free_t[k])) * MAX_SASA[model[target][k[1]].get_resname()])
                  for k in free_t}
        total = sum(buried.values())
        a_key = (target, anchor - target_offset)
        out['target_buried_area'] = round(total, 1)
        out['anchor_area_share'] = round(buried.get(a_key, 0.0) / total, 3) if total else None
    # PROPKA linkage: all groups, and interface His only.
    try:
        pk_c = propka_pkas(model, (target, binder))
        pk_f = {}; pk_f.update(propka_pkas(model, (target,))); pk_f.update(propka_pkas(model, (binder,)))
        out['linkage_all_groups'] = proton_linkage(pk_c, pk_f)
        his_keys = {(h['chain'], h['resnum']) for h in his}
        out['linkage_interface_his'] = proton_linkage(pk_c, pk_f, restrict=his_keys)
        if anchor is not None:
            k = (target, anchor - target_offset, 'HIS')
            out['anchor_pka_complex'] = pk_c.get(k); out['anchor_pka_free'] = pk_f.get(k)
    except Exception as e:  # noqa: BLE001
        out['propka_error'] = str(e)[:200]
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('INDEX:')[0])
    ap.add_argument('complexes', nargs='+')
    ap.add_argument('--target', default='A'); ap.add_argument('--binder', default='B')
    ap.add_argument('--target-offset', type=int, default=0)
    ap.add_argument('--anchor', type=int, default=None)
    ap.add_argument('--bc2-footprint', default='')
    ap.add_argument('--pae', default=None)
    ap.add_argument('--plddt-min', type=float, default=70.0)
    args = ap.parse_args(argv)
    for c in args.complexes:
        row = score_complex(c, args.target, args.binder, args.target_offset, args.anchor,
                            args.bc2_footprint, args.pae, args.plddt_min)
        print(json.dumps(row))


if __name__ == '__main__':
    sys.exit(main())
