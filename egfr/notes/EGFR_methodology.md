# Protein Design for Africa (PDFA)
## Challenge 1 — conditional EGFR binder: methodology, validation and per-design metrics

**Team:** Protein Design for Africa (PDFA), a pan-African non-profit training African scientists in de novo protein design.
**Track:** 2 · **Submitted:** 15 de novo single-chain miniproteins · **Molecule class:** `single_chain`
**All target preparation, both design campaigns, all validation, scoring, ranking and selection was carried out by the single PDFA team.**
**Compute:** Modal serverless cloud — NVIDIA H100 (design), NVIDIA L4 (co-folding), CPU (pKa). Total recorded spend $139.03.
**AI assistance:** Claude Opus 5 (Anthropic) — Modal orchestration code for both campaigns, pH-switch scoring implementation, results analysis and drafting of this document. All scientific decisions, target selection, hotspot and coldspot definitions, gate thresholds and the final ranking were made by PDFA.

> **All evidence here is computational.** No experimental validation has been performed. The pH preference is a PROPKA3 / Wyman-linkage estimate on a predicted structure, used to **rank** designs — not a measured affinity ratio.

---

## 0. Summary

PDFA ran **two independent design campaigns** against the same epitope, with identical hotspots, identical coldspots and an identical validation pipeline. They differ in exactly one controlled variable: whether the binder was optimised against the human structure alone, or against human and mouse simultaneously.

| | Stage 1 | Stage 2 |
|---|---|---|
| Design | BindCraft2, human EGFR only | BindCraft2 multi-target, human **and** mouse at once |
| Cross-species | tested after design, by refolding | **optimised during design** |
| BC2 trajectories | 456 | 42 |
| MPNN candidates | 4,891 | 165 |
| Backbones | 231 generated / 119 screened / 65 refolded | 7 accepted |
| Designs validated | 21 | 21 |
| Boltz-2 refolds | 284 (0 failed) | 42 (0 failed) |
| GPU | NVIDIA H100 + L4 | NVIDIA H100 + L4 |
| Compute cost | $104.06 | $34.97 |
| In the submitted 15 | 9 | 6 |

---

## 1. Target definition — identical in both stages

**Construct.** Human EGFR domain III, mature numbering **305–515**, chain A, **211 residues, no gaps**. PDFA trimmed this from the full ectodomain rather than loading an accession, so the design target matches the domain actually presented in the assay. The mouse ortholog was prepared to the identical span and numbering.

**Anchor residue: His409** (mature numbering; UniProt 433). In the cropped coordinate file this is residue index **105**, i.e. mature = 304 + index. PDFA's scoring code uses `ANCHOR = 105` with `--target-offset 304`.

### 1.1 Target sequences as supplied to BindCraft2 and Boltz-2

Verified directly from the coordinate files: 211 residues each, residues 305–515 continuous with no missing density, **87.2% sequence identity** between species. All four hotspot residues are **Thr406, His409, Gly410, Gln411 in both species** — identical residue types at identical numbering.

**Human EGFR domain III, chain A, residues 305–515 (211 aa)**
```
CEGPCRKVCNGIGIGEFKDSLSINATNIKHFKNCTSISGDLHILPVAFRGDSFTHTPPLD
PQELDILKTVKEITGFLLIQAWPENRTDLHAFENLEIIRGRTKQHGQFSLAVVSLNITSL
GLRSLKEISDGDVIISGNKNLCYANTINWKKLFGTSGQKTKIISNRGENSCKATGQVCHA
LCSPEGCWGPEPRDCVSCRNVSRGRECVDKC
```

**Mouse EGFR domain III, chain A, residues 305–515 (211 aa)**
```
CDGPCRKVCNGIGIGEFKDTLSINATNIKHFKYCTAISGDLHILPVAFKGDSFTRTPPLD
PRELEILKTVKEITGFLLIQAWPDNWTDLHAFENLEIIRGRTKQHGQFSLAVVGLNITSL
GLRSLKEISDGDVIISGNRNLCYANTINWKKLFGTPNQKTKIMNNRAEKDCKAVNHVCNP
LCSSEGCWGPEPRDCVSCQNVSRGRECVEKC
```

### 1.2 Why His409 — three requirements met by one residue

| Requirement | Evidence established by PDFA |
|---|---|
| **Functional epitope** | Direct EGF contact: **2.9 Å to Leu47 of EGF**, measured in 1IVO. A binder here competes with the natural ligand. |
| **Cross-species** | All four hotspots identical in human and mouse (confirmed above). Independently folded models of the two species agree to **0.03 Å Cα RMSD** across the hotspot set. |
| **pH-switchable** | Histidine, free pKa measured at **6.21–6.37** across all PDFA designs — between the two assay pH values, so its protonation state genuinely differs at 6.5 and 7.4. No other titratable side chain has a pKa in this window. |
| **Accessibility** | **65% solvent exposed** (SASA 145.1 Å²), **18.4 Å** from the nearest N-glycosylation sequon (Asn328) — available on the glycosylated, tethered assay target. |

### 1.3 Hotspots and coldspots (both stages, both species)

```
hotspots   406, 409, 410, 411
coldspots  328, 337, 353, 388, 389, 420, 443, 468, 471, 473, 504
```

**Hotspots.** Four contiguous residues centred on the anchor — enough that BindCraft2 builds a real interface patch rather than a point contact, few enough that the switch residue dominates the interface energetics.

**Coldspots, two classes, both of which would invalidate a design on the real target:**

- **Species-divergent surface residues** — 337, 353, 388, 443, 468, 471, 473. Contacting these makes the interface species-specific and breaks the K<sub>D</sub>-ratio requirement.
- **N-glycosylation sequons (N-X-S/T)** — 328, 337, 389, 420, 504. Asn337 is both. Glycans are present on the expressed target but absent from the design model, so a binder contacting them would be predicted to bind and then fail experimentally.

**A calibration PDFA performed and acted on.** An earlier coldspot set included residues 324, 418 and 467. PDFA measured their distances to the nearest hotspot at **13.0, 13.3 and 10.6 Å** and ran a dedicated ablation campaign (`coldfar_seed5003`) which showed that coldspots inside roughly 15 Å of the hotspot collapse interface i_pTM to as low as **0.08** — the exclusion zone swallows the epitope itself. Those three were removed. The final set is **11 coldspots, not 14**, and that decision rests on PDFA's own ablation data rather than on assumption.

---

## 2. Stage 1 — BindCraft2 design against human EGFR, cross-species tested afterwards

### 2.1 Design configuration (verbatim from PDFA's campaign file)

```json
{
  "campaign_name": "EGFR_d3_H409",
  "modality": "binder",
  "targets": [{
    "name": "EGFR_d3_human",
    "target_path": "egfr_d3ext_human.pdb",
    "chains": "A",
    "hotspots":  "406,409,410,411",
    "coldspots": "328,337,353,388,389,420,443,468,471,473,504"
  }],
  "binder_lengths": [80, 120],
  "forced_targeting": false,
  "min_hotspot_contact_final": 0.5,
  "max_coldspot_contact_final": 0.1,
  "termini_accessible": false,
  "aa_bias": { "C": 0, "H": 2.0, "D": 1.3, "E": 1.3, "K": 0.8, "R": 0.8 },
  "filters": {
    "Interface_H_Count":  { "threshold": 0,  "higher": true, "mandatory": false },
    "Interface_Residues": { "threshold": 12, "higher": true }
  },
  "sequence_candidates": 30,
  "desperation_trajectories": 40,
  "workers_per_gpu": 3,
  "max_trajectories": 30,
  "number_of_final_designs": 999,
  "campaign_seed": 7000
}
```

### 2.2 Why each setting — PDFA's reasoning

| Setting | Value | Reasoning |
|---|---|---|
| `aa_bias H` | 2.0 | The single most consequential setting. Doubling histidine frequency lets the binder contribute its *own* titratable groups to the interface, not merely respond to the target's His409. |
| `aa_bias D, E` | 1.3 | Raises the chance an Asp or Glu lands within salt-bridge range of the protonated His409 — the physical basis of the switch. |
| `aa_bias K, R` | 0.8 | A Lys or Arg near His409 **suppresses** its pKa, driving the switch backwards. The mirror image of the carboxylate bias. |
| `aa_bias C` | 0 | Free cysteines scramble disulfides and aggregate during expression. All 15 submitted designs contain **zero cysteines**. |
| `min_hotspot_contact_final` | 0.5 | At least half the hotspots must be contacted, so the epitope argument binds the output. |
| `max_coldspot_contact_final` | 0.1 | Hard ceiling on glycan and divergent-residue contact. |
| `Interface_Residues` | ≥ 12 | Doubles as a **His-tag decoy guard**: the assay construct carries a His-tag, and a binder merely gripping a run of histidines cannot reach twelve genuine interface residues. |
| `Interface_H_Count` | ≥0, non-mandatory | Records interface histidine count for ranking without rejecting designs on it. |
| `binder_lengths` | 80–120 | Below ~60 aa it is hard to bury enough surface for nanomolar affinity; above ~130 folding reliability and expression both fall. |
| `forced_targeting` | false | Hotspot satisfaction enforced through filters rather than a hard geometric constraint; PDFA found this gives higher-quality interfaces at the same epitope. |
| `termini_accessible` | false | Relaxed after PDFA observed it over-constrained topology at these lengths; terminal accessibility checked post hoc instead. |
| `sequence_candidates` | 30 | Thirty ProteinMPNN redesigns per accepted backbone, re-predicted and filtered. |
| `desperation_trajectories` | 40 | Fallback exploration budget if acceptance stalls. |

### 2.3 Stage 1 funnel — PDFA's measured numbers

```
456   BindCraft2 trajectories      (32 shards, seeds 7000-7015 and 8000-8015)
4891  ProteinMPNN candidates
231   backbones generated
119   backbones passed the design screen
65    backbones refolded against BOTH species with Boltz-2
284   Boltz-2 refold jobs, 0 failed
65    passed the human gate
21    passed BOTH species   <-- cross-species is the real bottleneck (32%)
15    of those 21 also pH-positive on human

Wall clock 0.89 h across parallel shards | recorded compute $104.06
Stop condition: 'goal reached: 21 backbones pass gate cross'
```

**The finding that shaped Stage 2.** 65 designs passed against human and only 21 survived mouse. Two-thirds of human-validated designs failed cross-reactivity *despite identical hotspot residues in both species*. Cross-species transfer was not a formality — which is precisely what motivated PDFA to design against both species simultaneously.

### 2.4 Stage 1 validation

Every accepted sequence was refolded **independently against human and mouse** with Boltz-2. The target chain receives a real MSA; the binder chain is given `msa: empty` so its fold is predicted de novo rather than recalled from the design. Boltz-2 did not design these sequences — AlphaFold2 inside BindCraft2 did — so this is a genuinely independent structural check, not self-scoring.

Stage 1 additionally computed **ipSAE** (Dunbrack) from the Boltz PAE matrices, which outperforms raw i_pTM for small binders.

```
Stage 1 gate:
  i_pTM >= 0.80  AND  ipSAE >= 0.70  AND  i_pAE <= 0.35  AND  His409 in footprint
  ... in EACH species independently
```

---

## 3. Stage 2 — BindCraft2 multi-target design against human AND mouse simultaneously

### 3.1 What changed, and what deliberately did not

Stage 2 is a controlled comparison. Hotspots, coldspots, amino-acid bias, contact thresholds, binder lengths, MPNN depth and the entire validation pipeline are **unchanged from Stage 1**. The single variable is the target specification:

```python
"targets": [
  { "name": "EGFR_d3_human", "target_path": "egfr_d3ext_human.pdb", "chains": "A",
    "hotspots": "406,409,410,411",
    "coldspots": "328,337,353,388,389,420,443,468,471,473,504",
    "objective": "target", "weight": 1 },
  { "name": "EGFR_d3_mouse", "target_path": "egfr_d3ext_mouse.pdb", "chains": "A",
    "hotspots": "406,409,410,411",
    "coldspots": "328,337,353,388,389,420,443,468,471,473,504",
    "objective": "target", "weight": 1 }
]
```

**This is multi-state design, not a bridging binder.** BindCraft2 evaluates the binder against each target structure separately and combines the losses at equal weight. The binder is never asked to contact both proteins at once; it is asked to be simultaneously good against two folds. Equal weights (1:1) because the objective is a K<sub>D</sub> **ratio near 1** — allowing either species to dominate the loss would defeat the purpose.

Every BindCraft2 stage reports both species. A representative accepted trajectory from PDFA's log:

```
passed screen  i_pTM.human=0.87  i_pTM.mouse=0.87  pLDDT.human=0.91  pLDDT.mouse=0.90
passed refine  i_pTM.human=0.87  i_pTM.mouse=0.87  pLDDT.human=0.92  pLDDT.mouse=0.91
passed anneal  i_pTM.human=0.89  i_pTM.mouse=0.88  pLDDT.human=0.93  pLDDT.mouse=0.93
passed harden  i_pTM.human=0.88  i_pTM.mouse=0.87  pLDDT.human=0.92  pLDDT.mouse=0.91
passed mutate  i_pTM.human=0.86  i_pTM.mouse=0.86  pLDDT.human=0.89  pLDDT.mouse=0.89
binder hallucination successful -> 30 MPNN redesigns
3 of 30 redesigns passed, keeping the best 1 by i_pDAE: candidate 5 (0.29)
```

### 3.2 Stage 2 funnel

```
42    BindCraft2 dual-target trajectories
165   MPNN candidates evaluated
7     backbones accepted  (13% trajectory acceptance)
21    sequences carried forward (3 MPNN variants per backbone)
42    Boltz-2 refold jobs (21 designs x 2 species) - 42 structures, 0 failed
42    PROPKA3 + Wyman scores

Recorded compute $34.97
```

Dual-target design costs roughly **2× the GPU time per trajectory** because every trajectory is folded against two structures. PDFA measured ~9 trajectories per container-hour on NVIDIA H100 against ~21 for single-target.

### 3.3 Boltz-2 co-folding — exact invocation (both stages)

```bash
boltz predict complex.yaml --model boltz2 --accelerator gpu --devices 1 \
      --recycling_steps 1 --sampling_steps 50 --diffusion_samples 1 \
      --output_format mmcif --use_msa_server --no_kernels
```

```yaml
version: 1
sequences:
  - protein: { id: A, sequence: <EGFR domain III, human or mouse> }
  - protein: { id: B, sequence: <binder>, msa: empty }
```

`msa: empty` on chain B is the critical detail: the binder gets no evolutionary information, so Boltz-2 must fold it from sequence alone.

### 3.4 Compute infrastructure

Both campaigns ran entirely on **Modal** serverless cloud, orchestrated by Python written by PDFA with Claude Opus 5. Design ran on **NVIDIA H100** GPUs with three BindCraft2 workers per GPU; Boltz-2 co-folding ran on **NVIDIA L4** GPUs; PROPKA3 and the Wyman integration ran on CPU containers. Model weights (AlphaFold2 5.3 GB, Boltz-2 checkpoint and CCD library) were cached in persistent Modal volumes. The orchestrator is a long-lived detached controller that survives client disconnection, with per-shard trajectory caps, budget ceilings and wall-clock deadlines.

---

## 4. pH-switch scoring — identical in both stages

PDFA wrote the scoring module (`score_ph_switch.py`). It runs on CPU, post-prediction, on every complex and on the separated chains.

### 4.1 The thermodynamics

The switch wanted is **tighter binding at pH 6.5 than at 7.4**, so every titratable group should bind a proton more tightly in the complex than when free. In Wyman linkage terms:

```
d log10 Ka / d pH  =  -dnH

  where  dnH(pH) = protons bound by the complex - protons bound by the free chains

  log10 Ka(6.4)/Ka(7.4) = integral of dnH over pH 6.4 -> 7.4   (41-point trapezoid)
```

PROPKA3 supplies the pKa of every titratable group in the complex and in each separated chain. **A positive integral means the design is predicted to bind more tightly at acidic pH.**

### 4.2 The geometry, read directly from each predicted structure

| Observation | Effect on His pKa | PDFA verdict label |
|---|---|---|
| Carboxylate O (Asp OD1/OD2, Glu OE1/OE2) within range of ND1/NE2 | **raises** — wanted | `good: carboxylate salt bridge` |
| Polar or aromatic partner, no carboxylate | mild | `ok: polar/aromatic partner` |
| Lys NZ or Arg NH/NE within 6 Å | **lowers** — forbidden | `bad: Lys/Arg within 6 Å` |
| Buried on binding with no partner | lowers | `neutral: solvent-facing` |

Burial is computed with Shrake–Rupley SASA against Tien et al. (2013) maximum values; interfaces are defined at a **4.5 Å heavy-atom cutoff**.

### 4.3 A result PDFA generated that corrects the obvious heuristic

The intuitive rule — *put a carboxylate close to the histidine* — **is not sufficient, and PDFA has the data to show it**:

- One PDFA design places a glutamate **2.44 Å** from His409 — closer than any submitted design — and its ΔpKa is **+0.05**: effectively no switch.
- Another buries His409 at 0.74 occupancy and drives its pKa **down by 4.2 units** — the opposite of the mechanism.
- Across PDFA's 68-complex rate-validation set, **30 of 68** interface histidines were flagged *bad: Lys/Arg within 6 Å* despite many having a carboxylate present, and only **42 of 68** achieved a positive log₁₀ Ka ratio.
- In an earlier PDFA campaign, **24 of 34** designs drove His409's pKa *down* on binding.

The switch depends on the **whole electrostatic environment** of the site, not one contact distance. PDFA therefore selects on the **measured ΔpKa and the Wyman integral**, using geometry only as a diagnostic. This is why the submitted ranking is not sorted by carboxylate proximity.

Observed fold-changes of 1.1–15× per pH unit are consistent with one-to-few protons under Henderson–Hasselbalch limits (≤10× per proton per pH unit).

---

## 5. Selection of the submitted designs

Both pools (21 + 21) were merged and passed through **three hard gates**. All three must be satisfied for Tier A:

```
BINDING           human i_pTM >= 0.85
CROSS-REACTIVITY  mouse i_pTM >= 0.78  AND  |human - mouse| <= 0.15
pH SWITCH         delta-pKa(His409) > 0  AND  fold Ka(6.4)/Ka(7.4) > 1   in BOTH species
```

Designs passing all three form **Tier A**, ordered by the fold-change on the **weaker** species. A design switching 30× on human and 1.1× on mouse has not satisfied a ratio objective, and averaging would hide that. Designs passing binding and cross-reactivity but with no verified switch form **Tier B**, ordered by the weaker species' i_pTM.

PDFA submits **15 designs — 7 Tier A and 8 Tier B — spanning 13 distinct backbones**, ranked 1–15 by the scheme above. Every design also clears the competition's automated novelty check; the per-design novelty score is reported alongside its metrics in section 7.

### Three decisions PDFA made deliberately

**1. pH outranks affinity.** PDFA's strongest binder — rank 9, i_pTM **0.949 human / 0.939 mouse**, the best pair in the whole candidate pool — is ranked *below eight weaker designs* because its ΔpKa is **−0.63**. It binds beautifully and switches backwards. Ranking it first would have been metric-chasing.

**2. ipSAE was excluded from the ranking.** Stage 1 computed it; Stage 2 did not. Including it would have systematically favoured Stage 1 on *instrumentation* rather than on quality. Both pools were ranked only on quantities measured identically for both. ipSAE is reported for Stage 1 designs but carries zero weight.

**3. Backbone diversity enforced.** The 15 submitted designs span **13 distinct backbones**. Only one backbone appears more than once (ranks 2, 3, 5), and only because all three of its MPNN variants pass the pH gate independently — and all three also clear the novelty filter, one of them at 4/4.

---

## 6. The central result of the two-stage comparison

Co-optimising against both species produced **markedly tighter cross-species parity**:

| | Best human−mouse i_pTM gap |
|---|---|
| **Stage 1** — human-only design | 0.070 |
| **Stage 2** — dual-target design | **0.005** |

Stage 2's three submitted Tier-A designs have gaps of **0.005, 0.014 and 0.034**. Its best binder sits at 0.949 / 0.939 — a gap of **0.010** at the top of the affinity range.

Stage 1 required 456 trajectories to find 21 cross-reactive designs (4.6%). Stage 2 found 7 cross-reactive backbones from 42 trajectories (17%), and **every one of its 21 sequences passed the cross-species gate**, because cross-reactivity was an optimisation objective rather than a filter.

There is a second signal worth noting. The only design scoring **4/4 on the competition's novelty check is a Stage 2 dual-target design**, and all three variants of that backbone clear the check. PDFA notes this as suggestive, not established — the sample is small — but it is consistent with multi-target optimisation pushing into less-explored sequence space.

> **PDFA's methodological claim, measurable in the submitted data:** for a cross-reactive binder, design against both orthologs rather than designing against one and screening the other.

---

## 7. All 15 submitted designs with every score

Novelty score from the competition's automated check (3/4 required). i_pTM and ipSAE from Boltz-2 co-folding, per species. ΔpKa = pKa(His409 in complex) − pKa(His409 free). Fold Ka = Ka(6.4)/Ka(7.4) from the Wyman integral. d(COO⁻) = His409 to nearest binder carboxylate oxygen.

### 7.1 Binding, cross-reactivity and novelty

| # | Name | aa | Stage | Novelty | i_pTM h | i_pTM m | gap | ipSAE h | ipSAE m | Tier |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | PDFA_EGFR_01_l86_d33b2b22_c4 | 86 | 1 human | **3/4** | 0.945 | 0.875 | 0.070 | 0.9242 | 0.7911 | A |
| 2 | PDFA_EGFR_02_l117_f0b94c6e_c2 | 117 | 2 dual | **3/4** | 0.925 | 0.920 | 0.005 | — | — | A |
| 3 | PDFA_EGFR_03_l117_f0b94c6e_c1 | 117 | 2 dual | **3/4** | 0.930 | 0.916 | 0.014 | — | — | A |
| 4 | PDFA_EGFR_04_l93_482b9c69_c4 | 93 | 1 human | **3/4** | 0.919 | 0.817 | 0.102 | 0.9142 | 0.7345 | A |
| 5 | PDFA_EGFR_05_l117_f0b94c6e | 117 | 2 dual | **4/4** | 0.921 | 0.887 | 0.034 | — | — | A |
| 6 | PDFA_EGFR_06_l107_0e0c2791_c24 | 107 | 1 human | **3/4** | 0.919 | 0.804 | 0.115 | 0.9015 | 0.7168 | A |
| 7 | PDFA_EGFR_07_l114_c95cacd1_c2 | 114 | 1 human | **3/4** | 0.921 | 0.809 | 0.112 | 0.8974 | 0.7092 | A |
| 8 | PDFA_EGFR_08_l111_dcbc7fb3_c2 | 111 | 1 human | **3/4** | 0.938 | 0.895 | 0.043 | 0.9395 | 0.8784 | B |
| 9 | PDFA_EGFR_09_l97_9cb614f2 | 97 | 2 dual | **3/4** | 0.949 | 0.939 | 0.010 | — | — | B |
| 10 | PDFA_EGFR_10_l99_edc75802_c4 | 99 | 1 human | **3/4** | 0.947 | 0.884 | 0.063 | 0.9358 | 0.8119 | B |
| 11 | PDFA_EGFR_11_l117_7f11b0bb_c3 | 117 | 1 human | **3/4** | 0.934 | 0.890 | 0.044 | 0.9277 | 0.8471 | B |
| 12 | PDFA_EGFR_12_l94_33634cd0_c22 | 94 | 1 human | **3/4** | 0.928 | 0.891 | 0.037 | 0.9044 | 0.8309 | B |
| 13 | PDFA_EGFR_13_l108_fd369d4f_c2 | 108 | 2 dual | **3/4** | 0.914 | 0.895 | 0.019 | — | — | B |
| 14 | PDFA_EGFR_14_l109_90db780a_c3 | 109 | 1 human | **3/4** | 0.916 | 0.880 | 0.036 | 0.9022 | 0.8497 | B |
| 15 | PDFA_EGFR_15_l99_616d5218_c2 | 99 | 2 dual | **3/4** | 0.909 | 0.883 | 0.026 | — | — | B |

### 7.2 pH switch and interface

| # | Name | ΔpKa h | ΔpKa m | fold Ka h | fold Ka m | pKa free | pKa cplx | partner | d(COO⁻) Å | buried Å² | charge | helix |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | PDFA_EGFR_01_l86_d33b2b22_c4 | 0.92 | 0.92 | 9.27 | 9.42 | 6.36 | 7.27 | D36 | 3.06 | 1034.7 | -9.0 | 0.83 |
| 2 | PDFA_EGFR_02_l117_f0b94c6e_c2 | 0.96 | 0.93 | 4.47 | 4.03 | 6.23 | 7.19 | GLU43 | 2.6 | 888.7 | -12.0 | 0.91 |
| 3 | PDFA_EGFR_03_l117_f0b94c6e_c1 | 0.9 | 0.82 | 3.82 | 3.05 | 6.25 | 7.15 | GLU43 | 3.11 | 965.2 | -15.9 | 0.91 |
| 4 | PDFA_EGFR_04_l93_482b9c69_c4 | 1.89 | 1.89 | 4.79 | 2.84 | 6.27 | 8.16 | E25 | 2.77 | 888.1 | -4.0 | 0.86 |
| 5 | PDFA_EGFR_05_l117_f0b94c6e | 0.83 | 0.78 | 3.07 | 2.41 | 6.26 | 7.09 | GLU43 | 3.16 | 979.2 | -12.0 | 0.91 |
| 6 | PDFA_EGFR_06_l107_0e0c2791_c24 | 0.42 | 0.42 | 2.42 | 3.62 | 6.26 | 6.69 | E95 | 2.93 | 1182.4 | -4.9 | 0.85 |
| 7 | PDFA_EGFR_07_l114_c95cacd1_c2 | 0.96 | 0.96 | 15.38 | 2.10 | 6.25 | 7.21 | E44 | 3.28 | 1435.7 | -10.9 | 0.82 |
| 8 | PDFA_EGFR_08_l111_dcbc7fb3_c2 | 0.66 | 0.66 | 1.17 | 0.57 | 6.21 | 6.87 | D100 | 3.12 | 970.5 | -9.9 | 0.86 |
| 9 | PDFA_EGFR_09_l97_9cb614f2 | -0.63 | -0.57 | 0.71 | 0.96 | 6.37 | 5.74 | ASP10 | 4.72 | 1205.8 | -1.0 | 0.91 |
| 10 | PDFA_EGFR_10_l99_edc75802_c4 | -3.65 | -3.65 | 1.01 | 0.38 | 6.26 | 2.61 | E68 | 4.22 | 1138.0 | -2.0 | 0.82 |
| 11 | PDFA_EGFR_11_l117_7f11b0bb_c3 | -2.39 | -2.39 | 0.98 | 0.23 | 6.27 | 3.88 | D39 | 4.27 | 1137.4 | -9.9 | 0.87 |
| 12 | PDFA_EGFR_12_l94_33634cd0_c22 | -1.59 | -1.59 | 0.29 | 0.77 | 6.3 | 4.71 | E55 | 7.26 | 983.4 | -5.9 | 0.91 |
| 13 | PDFA_EGFR_13_l108_fd369d4f_c2 | -1.6 | -1.56 | 0.52 | 0.47 | 6.24 | 4.64 | ASP88 | 3.91 | 948.9 | 0.07 | 0.82 |
| 14 | PDFA_EGFR_14_l109_90db780a_c3 | -2.18 | -2.18 | 1.85 | 2.40 | 6.28 | 4.11 | D12 | 4.24 | 861.7 | -7.0 | 0.84 |
| 15 | PDFA_EGFR_15_l99_616d5218_c2 | -0.08 | 0.0 | 0.39 | 0.63 | 6.26 | 6.18 | ASP63 | 2.74 | 1015.5 | -9.9 | 0.8 |

### 7.3 Sequences

```
>PDFA_EGFR_01_l86_d33b2b22_c4 | rank 1 | 86 aa | Tier A | novelty 3/4 | Stage 1 human-only
ALELAEQLVEVSQEIIEEYKDTMDPEDLVYAYTFDDFLRWLVEYMKRHKEWSMEELYARLKEFASLYGERMQYWLDLVLSRYKEKL
>PDFA_EGFR_02_l117_f0b94c6e_c2 | rank 2 | 117 aa | Tier A | novelty 3/4 | Stage 2 dual-target
LEEKLEEVKFKVRYLVKRRTYDFEFLEQIGLPEEERRRLFEAELERFARELKEVLSKFAEELGNEKGLKEVEELVEELKERLREAWEESFEKVLEVLEEVPKRLEELAEKIYEESLK
>PDFA_EGFR_03_l117_f0b94c6e_c1 | rank 3 | 117 aa | Tier A | novelty 3/4 | Stage 2 dual-target
LEKKLEEVKFKIRYLVKRFRYDFEFLHQIGLPEEERRRLFEEELERFAKELKEVLLEFAEYLGNEEGKEKVEELVEELEEELRRAWEESHEKVLEVLEEAEERLTKLAEEIYEESLK
>PDFA_EGFR_04_l93_482b9c69_c4 | rank 4 | 93 aa | Tier A | novelty 3/4 | Stage 1 human-only
MKEVEKEALKLLYELKREVYRAKYELDFMEMTQVYDPAYWEKIDPKIEPLLEKLEKFLKENKIENEELKKLIEEVLKLYEEVKERKKEMEKSQ
>PDFA_EGFR_05_l117_f0b94c6e | rank 5 | 117 aa | Tier A | novelty 4/4 | Stage 2 dual-target
MEKQYKRVEFKVRYYVKRFRYDAMFLEQIGLPKKELNRLWRGEMRKFQQRMRHWLMKFAVEMGDKKGREQVREATGRTYTQLEEAWNKDFKEVLRVLDRVPVEWVGLAHRIYKRFQK
>PDFA_EGFR_06_l107_0e0c2791_c24 | rank 6 | 107 aa | Tier A | novelty 3/4 | Stage 1 human-only
SEEEEIKKKLEELFEKLEKAVNKAVEAHKNHTWSEKKPEVDKEWNEVYEEMKKFAKELKGKYTVEEIYDILYEKFKEKYKDQPELVYFFEWFRYEWRLRELGEQMKP
>PDFA_EGFR_07_l114_c95cacd1_c2 | rank 7 | 114 aa | Tier A | novelty 3/4 | Stage 1 human-only
SELLEELVELARELYKVLWDIIHNGKPKERTDDWEKMLKLMVLETELYLLMDDIVRMIREGKETTPEFRKTVEEALKTLKEYKKLLEELNFEIPEEHKEILERTIEVLEKVLEG
>PDFA_EGFR_08_l111_dcbc7fb3_c2 | rank 8 | 111 aa | Tier B | novelty 3/4 | Stage 1 human-only
SMAELQEKLRRMAELLEELVKLLLSLVVEEPSGLSKEEIIKKIEELLKEIEELMEEIRRLHRESDSPFSDADMEYLEELLEEVIEAFKKKEWVQAHFRADYLLQLVKRMAR
>PDFA_EGFR_09_l97_9cb614f2 | rank 9 | 97 aa | Tier B | novelty 3/4 | Stage 2 dual-target
MYWAWRKLSDMLFMMKVDTIEYSARAKKAMKSKSPRRQKVVEKVREQIESQFHEFLRHLKKYGPMMRPKLAKSMRYYGRWYYQEAMEEVKAIVRKYQ
>PDFA_EGFR_10_l99_edc75802_c4 | rank 10 | 99 aa | Tier B | novelty 3/4 | Stage 1 human-only
KLEELLEKLEKLIEEAEKKSKLNKRTLFFYIVRSKQFIEHMKKAKTKEEQIEALREYVEFVEEQLRMEKLMQVKWTEEEKEFKKKLEEIIEELKKLIEE
>PDFA_EGFR_11_l117_7f11b0bb_c3 | rank 11 | 117 aa | Tier B | novelty 3/4 | Stage 1 human-only
MKELAERIYERYVAMAQKMTEVLSDEIFHKLTFLSLYADVVLRWALRHQTPEALERAIRELERLRKEMEEIVKELEEKLPSLDEEEAEKAKKEIEFLREQIETIEEIIEELRKLKEL
>PDFA_EGFR_12_l94_33634cd0_c22 | rank 12 | 94 aa | Tier B | novelty 3/4 | Stage 1 human-only
DLEERKKEIVEEYERKLKEFEEKIKELSKVADKESLEKLEIEFLHFWLTGGVLVERLKRARTEPRFEHLVGMLEEEIAKLEKKIDEVKKEAEEK
>PDFA_EGFR_13_l108_fd369d4f_c2 | rank 13 | 108 aa | Tier B | novelty 3/4 | Stage 2 dual-target
SELKKLEEEKKKLLEEIKKAAEELLRLLEEASRETGKPIDRKEVEELLKNLSDPSLSVWEKAVNMFKLVPYILRFRHHPRFPELAAVDQKISDKLRELEKVEDKEKKL
>PDFA_EGFR_14_l109_90db780a_c3 | rank 14 | 109 aa | Tier B | novelty 3/4 | Stage 1 human-only
SLPYQVMKLKYDLQFSLEIEDFSKEVLELHLKQIDEILKLFEETPEYKENPEIKKLMEELKKIKEELEKKLKEGTLSEEEAKKYYEEVLKLLEELIPLLRKIEEEKQKK
>PDFA_EGFR_15_l99_616d5218_c2 | rank 15 | 99 aa | Tier B | novelty 3/4 | Stage 2 dual-target
MATEHLRFFVELAEEMGDNIPEEVKKPLKELLEYVEKEGDEVVTEEMYELLRKVYHLLWWIWDKVRLMTAWGEPVVLDEEVVKKMVEYLEEIKELLKKK
```

---

## 8. Shared properties of all 15

| Property | Value across all 15 |
|---|---|
| Length | 86–117 aa (all within the 10–250 limit) |
| Topology | all-α, **β-sheet fraction 0.00** |
| Cysteines | **0 in every design** — no disulfide scrambling, no aggregation risk |
| Net charge at pH 7.4 | −15.9 to +0.07 |
| Helix fraction | 0.80–0.91 |
| Interface residues | 15–36 |
| Interface buried area | 861.7–1435.7 Å² |
| His409 contacted | **True in both species, all 15** |
| Hotspot contact fraction | ≥ 0.5 |
| Coldspot contact fraction | ≤ 0.09 |

**All-α by design.** PDFA noted Anthropic's own measurement that β-containing EGFR binders validated at **12.8%** against **28.4%** for all-α, and switched the β-sheet filter off rather than requiring it. Every submitted design is pure helix.

**Assay-condition check PDFA performed.** Ionic strength is matched at ~170 mM between the two pH conditions, so pH is the only variable. The Debye length of ~7.4 Å screens long-range electrostatics but **not** the 2.6–3.3 Å contact ion pair the switch depends on. The epitope patch is net +2 against binders at median −10, so the electrostatic complementarity is favourable.

---

## 9. Limitations, stated plainly by PDFA

1. **All evidence here is computational.** No experimental validation has been performed on these 15 designs. The pH preference is a PROPKA3 / Wyman-linkage estimate on a predicted structure, used to **rank** designs — it is not a measured affinity ratio.
2. **Stage 1 reports a single ΔpKa** rather than one per species. Its per-species fold-Ka values are independent, but its cross-species pKa check is less rigorous than Stage 2's, where PROPKA3 ran separately on each species' complex.
3. **ipSAE exists only for Stage 1**, and was therefore excluded from ranking entirely rather than applied asymmetrically.
4. **Ranks 2, 3 and 5 share one backbone.** They are included because all three pass the pH gate independently and all three clear the novelty filter, but they are not three independent tests of that fold.
5. **Boltz-2 ran with one diffusion sample and one recycling step.** Deeper sampling would tighten the confidence estimates; this was an explicit compute trade-off.
6. **PROPKA3 is an empirical pKa predictor** applied to a predicted structure. Two layers of approximation sit between these numbers and a measured K<sub>D</sub> ratio.
7. **The 0.03 Å Cα RMSD and 2.9 Å EGF contact** are single measurements on single predicted or crystallographic structures, not ensemble averages.

---

## 10. Reproducibility and attribution

| Component | Version / identifier |
|---|---|
| BindCraft2 | v1.0.1, commit `301efdd1937fb963cc40a5b0ecc1bc2f9b2b2d15` (PacesaLab) |
| ProteinMPNN | internal to BindCraft2 — 30 sequence candidates per backbone |
| AlphaFold2 | internal to BindCraft2 — design-time structure prediction |
| Boltz-2 | v2.2.1 |
| PROPKA | 3 |
| ipSAE | v1.0.1 (Stage 1 only) |
| Compute platform | Modal serverless cloud |
| GPUs | **NVIDIA H100** (BindCraft2 design, 3 workers/GPU) · **NVIDIA L4** (Boltz-2 co-folding) · CPU containers (PROPKA3, Wyman integration) |
| AI assistance | **Claude Opus 5 (Anthropic)** — Modal orchestration code for both campaigns, pH-switch scoring implementation, results analysis, drafting of this document |
| Stage 1 seeds | 7000–7015 and 8000–8015 (32 shards) |
| Stage 2 seed | 7000 |
| Total recorded compute | $139.03 |
| Team | Protein Design for Africa (PDFA) — single team, all stages |

All campaign configurations, controller logs, cost records, predicted structures, PAE matrices and scoring code are retained by PDFA and available on request.
