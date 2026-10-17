# PDFA method 2 - design and validation notes

**Target and epitope.** Human TNF-alpha soluble trimer (P01375 77-233; PDB 1TNF chains A-C with the
native Asp143 side chain restored from 1A8M) and mouse TNF-alpha (P06804 80-235; PDB 2TNF). We mapped
the TNFR2 footprint from 3ALQ onto the human trimer: the receptor-binding groove straddles two
protomers (side 1: residues 20-33, 63, 67, 113-115, 143-149; side 2: 71-97, 128). Human/mouse
conservation (79 % identity) was computed per residue; the 71-73 loop, 85 and 89 are divergent
(the reason adalimumab and infliximab do not bind mouse TNF), so the design hotspots were restricted
to conserved, solvent-exposed, receptor-contacting residues within 13.5 A of the groove centre:
K90/V91/N92 on one protomer and R32/A33/F144/A145 on its neighbour.

**Design.** BindCraft2 (AF2-multimer hallucination + ProteinMPNN redesign) against the two protomers
that form the groove (1TNF chains A+B; the third protomer is >= 18 A from every hotspot and untouched
by the receptor), chain-prefixed hotspots, `min_receptor_chains_contacted_final 2` so every candidate
straddles the groove, 80-120 aa (plus 60-90 aa shards), 30 MPNN sequences per backbone, histidine
favoured in the composition prior (`aa_bias H 2.0`, Cys excluded), no species coldspots. A server-side
controller on Modal kept 10-14 H100 shards running (~1,100 trajectories over two runs), streamed every
MPNN candidate, re-gated it, and refolded the best 4 per backbone. The design-time screen was i_pTM >= 0.55 /
i_pAE <= 0.50 (pTM >= 0.55, unbound pLDDT >= 0.70, RMSD <= 3.5 A, >= 12 interface residues, hotspot
contact >= 0.5, two receptor chains) and the Boltz-2 trimer refold made the call.

**Validation (source of truth).** Every candidate was refolded with stock Boltz-2 2.2.1 (MSA for the
target, single sequence for the binder) against THREE copies of the mature sequence plus the binder,
for human AND mouse. Gates per species, calibrated on controls (three folded EGFR binders: merged
ipSAE 0.00, interface PAE 20-21 A; a generic helix: 0.50 / 7 A): binder-target pair ipTM >= 0.80,
merged-target ipSAE >= 0.70, mean interface PAE <= 8 A, and >= 2 groove hotspots contacted on 2
protomers in any of the three symmetric grooves. Boltz's complex ipAE is not used (protomer-protomer
contacts dominate it). Because single-sample refolds of the same sequence moved by up to 0.6 ipSAE in
mouse, every shortlisted design was refolded twice more per species; a species "passes" when the
majority of its refolds pass, and the ranking uses the mean ipSAE minus the replicate spread.

**pH selectivity (bind at 7.4, release at 6.0).** The mechanism we design for is an interface
histidine whose pKa DROPS on binding: protonation at pH 6.0 is then penalised in the complex
(electrostatic repulsion from a target Lys/Arg -- R32 and K90 are hotspots, K128/K65/R82 line the
groove -- or desolvation of a buried neutral His), so the complex dissociates in the endosome while
the neutral His at pH 7.4 packs normally. The opposite geometry (His across a target carboxylate,
E146/D143 at the groove centre) would raise the pKa and tighten binding at low pH, so it is
penalised. Each refolded complex was scored with PROPKA 3.5 on the complex and on the separated
partners (the free trimer kept intact); the Wyman linkage integral of the proton-uptake difference
over pH 6.0-7.4 gives a predicted log10 Ka(6.0)/Ka(7.4), reported sign-corrected as the "pH switch"
(positive = weaker at 6.0). Per-His geometry (nearest carboxylate O / Lys-Arg N, burial on binding)
was read directly from the model. A design-time prior from the same scores on the AF2 complex
ordered the refold queue; it never gated.

**Histidine scan.** For each cross-reactive backbone a geometric scan proposed His substitutions, in
order: binder Asp/Glu salt-bridged to a target Lys/Arg (the charge swap: neutral His keeps the contact
at 7.4, protonated His is repelled at 6.0), positions within 7 A of a target Lys/Arg and >= 6 A from any
target carboxylate, buried positions without a carboxylate partner; target cations that differ in mouse
(R6, R31, R131, R138) were demoted. Singles and pairwise doubles were refolded in both species and kept
when they stayed cross-reactive with a clean His mechanism.

**Ranking and selection.** Tiers follow the challenge priority: (1) cross-reactive with a clean His
switch in both species (>= 1 interface His next to a target Lys/Arg or buried, none paired with a target
carboxylate, positive His-only linkage); (2) cross-reactive, clean switch on human; (3) cross-reactive;
(4)/(5) human only. Within a tier: mean species ipSAE minus half the species gap minus a quarter of the
replicate spread plus 0.3 x clip(His switch, 0, 1). Diversity: a cap per BindCraft2 backbone and a
pairwise identity ceiling between picks; no Cys; 60-120 aa. The all-groups PROPKA linkage is reported
but not ranked on: it is negative for every groove binder, dominated by the target's Asp/Glu shifts.

**Code.** Raptorbind (this repository): `scripts/bc2_campaigns.py` (campaign `tnf_groove`),
`scripts/modal_bc2_egfr.py` (design shards), `scripts/modal_bc2_egfr_overnight.py` (controller),
`scripts/refold_multichain.py` (trimer refold), `scripts/score_ph_switch.py` (His geometry + PROPKA
linkage), `scripts/bc2_submission.py` (ranking). Analysis: `agent/projects/tnf_alpha_conditional/`.
