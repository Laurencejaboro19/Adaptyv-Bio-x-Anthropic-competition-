# Design strategy

## Shared approach

1. **Target preparation.** Trim the target to the domain or groove that matters, keep numbering consistent between species, and define hotspots (residues we want touched) and coldspots (surface we want avoided).
2. **Generate with BindCraft2.** AF2-multimer hallucination plus ProteinMPNN, with the hotspots and coldspots above, histidine favoured in the amino-acid prior (cysteine excluded), and a design-time screen on i_pTM, i_pAE, pTM, unbound pLDDT, interface size and hotspot contact.
3. **Refold with Boltz-2 (source of truth).** Every candidate sequence is refolded as a full complex, and all scoring and selection use the refolded structure, not the generator output. Both species are refolded where cross-reactivity is required.
4. **Score the pH switch.** PROPKA3 on the complex and on the separated partners, then the Wyman linkage integral of the proton-uptake difference over the pH window gives a predicted log10 Ka ratio between the two pH values. Interface histidine geometry (nearest target carboxylate or Lys/Arg, burial on binding) is read directly from the model.
5. **Histidine scan.** Near-misses get geometry-guided histidine substitutions that are refolded and re-scored.
6. **Rank.** Tiered by the challenge priorities (cross-reactive with a clean histidine switch first), then by interface confidence with a penalty for replicate spread, with a diversity cap per backbone.

## Challenge 1: EGFR (conditional, human and mouse)

- **Target:** EGFR domain III, mature residues 305-515 (211 aa), anchored on His409. Mouse prepared to the same span and numbering (87.2% identity). Hotspots are Thr406, His409, Gly410 and Gln411, identical in both species.
- **Two campaigns, same epitope and validation:** Stage 1 optimised against human EGFR only (456 trajectories, 21 validated designs). Stage 2 optimised against human and mouse at once (42 trajectories, 21 validated designs). 15 designs were submitted: 9 from Stage 1 and 6 from Stage 2.
- **pH direction:** tighter binding at pH 6.4 than at 7.4 (the challenge page specifies pH 6.5 vs 7.4). The scoring code integrates over pH 6.4 to 7.4.
- **Finding:** coldspots placed within about 15 A of the hotspot collapsed i_pTM, shown by an ablation, and carboxylate proximity alone was not enough for a good histidine switch (in a 68-complex validation, 30 had an interface histidine flagged bad because of a nearby Lys/Arg, even where carboxylates were present).

## Challenge 2: TNF-alpha (conditional, human and mouse)

- **Target:** the receptor-binding groove, which straddles two protomers of the soluble trimer. Hotspots are restricted to conserved, solvent-exposed, receptor-contacting residues, avoiding the human/mouse-divergent 71-73 loop.
- **Mechanism designed for:** binding at pH 7.4 and release at pH 6.0 through an interface histidine whose pKa drops on binding (repulsion from a target Lys/Arg, or desolvation of a buried neutral histidine). A histidine facing a target carboxylate is penalised.
- **Two methods, 15 designs submitted:** Method 1 designed against human and mouse simultaneously (9 of the 15). Method 2 was human-first with a mouse check and replicate refolds (6 of the 15).

Details and exact settings are in the notes under `egfr/notes/` and `tnf/notes/`.
