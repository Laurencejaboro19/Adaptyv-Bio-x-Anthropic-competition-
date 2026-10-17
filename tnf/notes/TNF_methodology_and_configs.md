# PDFA TNF-α — methodology

PDFA - Conditional TNF-alpha binders (Anthropic x Adaptyv, Challenge 2). 15 miniprotein designs (70-118 aa, no cysteines, all single_chain) built to bind human and mouse TNF-alpha at pH 7.4 and release at pH 6.0 through an engineered histidine switch.

DESIGN RULE
Every design must be predicted to bind BOTH human and mouse TNF-alpha (soluble trimer) at pH 7.4 AND show the histidine pH switch on BOTH species. This is stricter than the challenge scoring, which measures pH selectivity on human only. Residue numbers are mature-chain numbers (human P01375 77-233; mouse P06804 80-235; the equivalent residue is Asp143 in human and Asp142 in mouse; both targets carry the native Asp, matching the wet-lab construct).

TWO PDFA METHODS (both: BindCraft2 -> Boltz-2 -> PROPKA/Wyman)
Both methods use BindCraft2 (AF2-multimer hallucination + ProteinMPNN) against the two protomers that form the receptor groove, Boltz-2 2.2.1 refolding of every sequence with the full trimer (the refolded complex is the single source of truth for all scoring and selection), PROPKA 3.5 + Wyman linkage for the histidine pH switch, and a histidine scan to place switch histidines at the interface.

PDFA METHOD 1 - dual-species design (9 of the 15: 7 direct BindCraft2 designs + 2 His-scan variants)
- Targets: human and mouse TNF-alpha A/B protomer dimers designed against SIMULTANEOUSLY (two targets in one BindCraft2 campaign, weight 1.0 each); scoring on the full trimers.
- Hotspots (inter-protomer receptor groove): human A29, A32-33, A143-147, B86-87, B90-92; mouse A29, A32-33, A142-146, B85-86, B89-91.
- Coldspots (surface to avoid): human A55, A57, A71, A73, A75, A82, A91-99, A102-103, A118-125, A151-157, B13, B15, B34, B36, B57, B59, B61, B63, B69, B98, B100-106, B112-117, B119, B143, B146-149, B151-157; mouse A55, A57, A72, A74, A81, A90-98, A100-103, A117-124, A151-156, B13, B15, B34, B36, B57, B59, B61, B63, B97, B102-105, B111-116, B118, B145-148, B150-156.
- BindCraft2 (PacesaLab commit 301efdd): binder 70-120 aa; aa_bias C 0, H +2.5, D/E +0.8, K/R +1.0; 10 trajectories per worker batch, 40 desperation trajectories, 30 sequence candidates, 2 workers/GPU; min hotspot contact 0.4, max coldspot contact 0.1, >= 2 receptor chains contacted; i_pTM >= 0.55, pTM >= 0.55, i_pAE <= 0.50, monomer pLDDT >= 0.70; interface residues >= 12, binder RMSD <= 3.5 A; forced_targeting off, termini_accessible off. Seeds 9000 + shard index.
- Boltz-2 2.2.1: full trimer (A,B,C) + binder, recycling 3, sampling steps 100, 1 diffusion sample, cached MSAs; human refold first, mouse refold for human passers.
- Gates per species: ipTM >= 0.75, ipSAE (PAE cutoff 15) >= 0.60, contact PAE <= 7 A, binder pLDDT >= 0.80, hotspot fraction >= 0.40, coldspot fraction <= 0.30, >= 2 chains contacted, His fold >= 0.30 on BOTH species.
- His-scan: near-misses received up to 5 single-His mutants at contact residues, refolded and scored the same way.
- Compute: resumable cloud controller on Modal (H100 design shards, L4 refolds, budget governor); 12 shards, 106 sequences scored, 33 passed all gates.

PDFA METHOD 2 - human-first design, mouse check, replicate refolds (6 of the 15: 1 direct BindCraft2 design + 5 His-scan variants)
- Target: human TNF-alpha, PDB 1TNF chains A+B (the third protomer is >= 18 A from every hotspot), native Asp143 restored from 1A8M; mouse TNF-alpha (P06804, PDB 2TNF) for the cross-species refold.
- Hotspots: K90, V91, N92 on one protomer and R32, A33, F144, A145 on its neighbour (conserved, solvent-exposed, receptor-contacting residues within 13.5 A of the groove centre, mapped from the TNFR2 footprint in 3ALQ; the divergent 71-73 loop and residues 85/89 avoided). Coldspots: none.
- BindCraft2: binder 80-120 aa (plus 60-90 aa shards), 30 ProteinMPNN sequences per backbone, aa_bias H +2.0, Cys excluded, >= 2 receptor chains contacted; design filters i_pTM >= 0.55, i_pAE <= 0.50, pTM >= 0.55, unbound pLDDT >= 0.70, RMSD <= 3.5 A, >= 12 interface residues, hotspot contact >= 0.5.
- Compute: Modal controller keeping 10-14 H100 shards running, ~1,100 trajectories over two runs; best 4 sequences per backbone refolded; 154 candidates from 28 backbones.
- Boltz-2 2.2.1: THREE copies of the target + binder, human AND mouse; gates per species ipTM >= 0.80, merged-target ipSAE >= 0.70, mean interface PAE <= 8 A, >= 2 groove hotspots contacted on 2 protomers. Every shortlisted design refolded 3 times per species; a species passes when the majority of refolds pass; ranking uses mean ipSAE minus replicate spread.
- pH switch: PROPKA 3.5 + Wyman linkage on the refolded complex and separated partners; target mechanism is an interface His whose pKa drops on binding (repulsion from a target Lys/Arg such as R32/K90, or desolvation of a buried neutral His); a His facing a target carboxylate (E146/D143) is penalised.
- His-scan: substitutions proposed in order of (1) binder Asp/Glu salt-bridged to a target Lys/Arg (charge swap), (2) positions within 7 A of a target Lys/Arg and >= 6 A from any carboxylate, (3) buried positions without a carboxylate partner; mouse-divergent target cations demoted; singles and doubles refolded on both species and kept if cross-reactive with a clean His mechanism.
- Internal tiers: 1 cross-reactive with clean His switch on both species, 2 clean on human, 3 cross-reactive, 4/5 human only.

PH-SWITCH SCORE (both methods)
PROPKA 3.5 pKa values of every histidine in the bound complex, the free target and the free binder; Wyman linkage over pH 6.0-7.4: log10[Ka(7.4)/Ka(6.0)] = -integral(dn dpH), histidines only ("His fold"). Positive = stronger binding at pH 7.4. A switch histidine has free pKa 6.5-8.0 and a binding shift <= -0.5 (clean if no carboxylate within 5 A).

COMBINING THE METHODS, RANKING AND SELECTION
Candidates from both methods were rescored on one scale and gated identically (ipTM >= 0.75, ipSAE >= 0.60, contact PAE <= 7 A, pLDDT >= 0.80, His fold >= 0.30 on both species): 34 designs passed. They were ranked by pH switch first: the LOWER of the human and mouse His fold (highest first), so the best-ranked designs show a strong switch on both species. Binding metrics were gates, not ranking terms. Designs were screened with the Adaptyv binder-prescreen and passed the competition novelty check (>= 3/4). The 15 submitted designs are 9 from method 1 and 6 from method 2, listed in rank order (His fold, lower of the two species, 0.95 down to 0.36).

RESULTS
All 15 pass both-species binding and pH gates. Human His fold 0.36-1.02 (median 0.52); mouse 0.37-1.33 (median 0.72). Lowest ipSAE across species 0.72-0.87; ipTM >= 0.88 on both species. The designs come from 8 parent backbone families.

USE OF CLAUDE
Claude Sonnet 5.5 (Anthropic) was used to write the configurations, the Modal run code, the scoring, ranking and design selection, and this report.

SUPPORTING FILES
Full configs for both methods; method 2 design and validation notes; Modal pipeline code and input structures (method 1); the 15 designs with their metrics.


---

# Appendix: full configurations

# PDFA TNF-α — full configurations for PDFA method 1 and PDFA method 2
Configs, run code, scoring, ranking, selection and this report were produced with Claude Sonnet 5.5 (Anthropic) and run by the PDFA team.
Residue numbers are mature-chain numbering (human P01375 77–233, mouse P06804 80–235). In the designed dimers chain A and chain B are two neighbouring protomers of the trimer.

## A. PDFA method 1 — dual-species BindCraft2 on Modal (app `pdfa-tnf-dual`)
### Targets
- Design inputs: `tnf_human_AB.pdb` (human, chains A,B; residues 12–157, 146 aa/chain, Asp143) and `tnf_mouse_AB.pdb` (mouse, chains A,B; residues 12–156, 145 aa/chain, equivalent Asp142). AlphaFold-predicted structures.
- Scoring inputs: `tnf_human_trimer_12on.pdb`, `tnf_mouse_trimer_12on.pdb` (full trimers A,B,C) so binders are checked against all three grooves (A,B), (B,C), (C,A).
### Hotspots (inter-protomer receptor groove)
- Human: `A29-29,A32-33,A143-147,B86-87,B90-92`
- Mouse: `A29-29,A32-33,A142-146,B85-86,B89-91`
### Coldspots (residues the binder must avoid; outer surface and divergent loops)
- Human: `A55-55,A57-57,A71-71,A73-73,A75-75,A82-82,A91-99,A102-103,A118-125,A151-157,B13-13,B15-15,B34-34,B36-36,B57-57,B59-59,B61-61,B63-63,B69-69,B98-98,B100-106,B112-117,B119-119,B143-143,B146-149,B151-157`
- Mouse: `A55-55,A57-57,A72-72,A74-74,A81-81,A90-98,A100-103,A117-124,A151-156,B13-13,B15-15,B34-34,B36-36,B57-57,B59-59,B61-61,B63-63,B97-97,B102-105,B111-116,B118-118,B145-148,B150-156`
### BindCraft2 campaign (PacesaLab, commit 301efdd1937fb963cc40a5b0ecc1bc2f9b2b2d15, `[cuda13]`)
```
campaign_name TNF_C2_dual; modality binder; number_of_final_designs 999
targets: [TNF_human, TNF_mouse], chains "A,B", objective "target", weight 1.0 each,
         hotspots/coldspots as above
binder_lengths [70,120]; aa_bias {C:0, H:2.5, D:0.8, E:0.8, K:1.0, R:1.0}
max_trajectories 10 (per worker batch); desperation_trajectories 40; sequence_candidates 30
forced_targeting false; termini_accessible false; workers_per_gpu 2
min_hotspot_contact_final 0.4; max_coldspot_contact_final 0.1
min_receptor_chains_contacted_final 2
min_iptm_final 0.55; min_ptm_final 0.55; max_ipae_final 0.50; min_monomer_plddt_final 0.70
filters: Interface_Residues >= 12; Binder_RMSD <= 3.5 A
```
Compute: H100 design shards, seed = 9000 + shard index, 12 shards run (5 rounds of up to 3).
### Boltz-2 refold (source of truth)
Boltz-2 2.2.1 on L4: full trimer (A,B,C) + binder (D); `--recycling_steps 3 --sampling_steps 100 --diffusion_samples 1 --no_kernels --output_format mmcif --write_full_pae`; precomputed MSAs cached; human first, mouse only for designs passing the human structural gate; 1 refold per species. At most 60 sequences refolded per round, at most 2 per parent backbone.
### Scoring and pH switch
- ipSAE (PAE cutoff 15), contact PAE, 5 Å contacts, hotspot / coldspot fractions over the three grooves, chains contacted, binder pLDDT, ipTM.
- PROPKA 3.5.1 on complex, target alone and binder alone; Wyman linkage over pH 6.0–7.4 (trapezoid); `his_log10_fold` (histidines only) and `all_log10_fold`; switch His = free pKa 6.5–8.0 and shift ≤ −0.5; clean = no carboxylate within 5 Å.
### Gates (both species)
ipTM ≥ 0.75; ipSAE ≥ 0.60; contact PAE ≤ 7 Å; binder pLDDT ≥ 0.80; hotspot fraction ≥ 0.40; coldspot fraction ≤ 0.30; chains contacted ≥ 2; His fold ≥ 0.30.
### His-scan
Near-misses (structure gate passed, pH gate not) with no existing substitution get up to 5 single-His mutants at binder contact residues, refolded and scored the same way.
### Controller, cost and resilience
Controller runs in Modal (12 h timeout, PC can be off); resumable state; budget cap USD 85 (estimated); H100 USD 3.95/h, L4 USD 0.80/h; goal 24 passing designs; stops on goal, shard ceiling (14), deadline, or budget; circuit breaker if every shard in a round fails.
### Outcome
12 shards, 106 sequences scored, 33 passed all gates, est. spend USD 72.26, 7.4 h. Entered the final 15: 9 designs (7 direct BindCraft2 + 2 His-scan variants).

## B. PDFA method 2 — human-first BindCraft2, mouse check, replicate refolds (campaign `tnf_groove_ab`)
(Full notes in `pdfa_method2_notes.md`.)
### Targets and hotspots
- Design target: human TNF-α, PDB 1TNF chains A+B (the two protomers forming the groove; third protomer ≥ 18 Å from every hotspot), native Asp143 restored from 1A8M. Mouse TNF-α (P06804 80–235; PDB 2TNF) used for the cross-species refold only.
- Hotspots (chain-prefixed): **K90, V91, N92** on one protomer and **R32, A33, F144, A145** on its neighbour. Chosen as conserved, solvent-exposed, receptor-contacting residues within 13.5 Å of the groove centre (mapped from TNFR2 footprint in 3ALQ; divergent loop 71–73, residues 85 and 89 avoided).
- Coldspots: **none** (no species coldspots).
### BindCraft2 settings
```
binder length 80-120 aa (plus 60-90 aa shards); 30 ProteinMPNN sequences per backbone
aa_bias H 2.0, Cys excluded
min_receptor_chains_contacted_final 2
i_pTM >= 0.55 ; i_pAE <= 0.50 ; pTM >= 0.55 ; unbound pLDDT >= 0.70
binder RMSD <= 3.5 A ; interface residues >= 12 ; hotspot contact >= 0.5
```
Design-time filters were set to the values above and the Boltz-2 trimer refold made the final call.
### Compute
Server-side Modal controller keeping 10–14 H100 shards running, ~1,100 trajectories over two runs; every MPNN candidate re-gated; best 4 per backbone refolded.
### Boltz-2 validation
Stock Boltz-2 2.2.1 (MSA for the target, single sequence for the binder) against **three copies** of the mature target + binder, for human and mouse. Gates per species: binder–target ipTM ≥ 0.80; merged-target ipSAE ≥ 0.70; mean interface PAE ≤ 8 Å; ≥ 2 groove hotspots contacted on 2 protomers in any of the three grooves. Gates calibrated on controls (three folded EGFR binders: ipSAE 0.00; generic helix: 0.50). Every shortlisted design refolded two more times per species (3 total); a species passes when the majority of refolds pass; ranking uses mean ipSAE minus replicate spread.
### pH switch
PROPKA 3.5 on complex and separated partners; Wyman linkage of proton uptake over pH 6.0–7.4 gives predicted log10 Ka(6.0)/Ka(7.4), sign-corrected as the "pH switch" (positive = weaker at pH 6.0). Mechanism: interface His whose pKa drops on binding (repulsion from target Lys/Arg, e.g. R32, K90; or desolvation of buried neutral His); His across a target carboxylate (E146/D143) is penalised.
### His-scan
Per cross-reactive backbone: (1) binder Asp/Glu salt-bridged to target Lys/Arg (charge swap), (2) positions within 7 Å of Lys/Arg and ≥ 6 Å from any carboxylate, (3) buried positions without carboxylate partner; mouse-divergent target cations (R6, R31, R131, R138) demoted. Singles and doubles refolded on both species, kept when cross-reactive with clean His mechanism.
### Tiers and ranking inside method 2
Tier 1 cross-reactive with clean His switch on both species; tier 2 cross-reactive, clean on human; tier 3 cross-reactive; 4/5 human only. Score = mean ipSAE − ½·species gap − ¼·replicate spread + 0.3·clip(His switch, 0, 1); per-backbone cap and identity ceiling.
### Outcome
154 candidates from 28 backbones (tags tnf_ab_prod1, tnf_ab_prod2, his_scan2, his_scan3); all human ipSAE ≥ 0.70. Entered the final 15: 6 designs (1 direct BindCraft2 + 5 His-scan variants).

## C. Combining PDFA methods 1 and 2
Both candidate sets were rescored on one scale (ipSAE, ipTM, contact PAE, PROPKA His fold for human and mouse) and gated identically; 34 designs passed. Ranking: **lower of the human and mouse His fold, highest first**. Binding metrics were gates, not ranking terms. Then the Adaptyv binder-prescreen (known-TNF arm only) and the competition novelty check (≥ 3/4): 15 designs submitted.

| | Method 1 | Method 2 | Total |
|---|---|---|---|
| Designs in the 15 | 9 (7 direct + 2 His-scan) | 6 (1 direct + 5 His-scan) | 15 |
| Refolds per species | 1 | 3 | |
| Species at design time | human + mouse | human only (mouse checked) | |
