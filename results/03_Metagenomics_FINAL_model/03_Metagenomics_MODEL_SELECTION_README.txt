METAGENOMICS — MODEL SELECTION
================================

Three models were compared to handle the sex imbalance
(IBS-C was 21 female / 1 male; sex was imbalanced, chi-sq p=0.035).

MODEL 1 (archived: model1_Group_only)
  Design: ~ Group        | 75 subjects
  Result: 79 significant species (C:7, D:72)
  Note: No covariate adjustment. Exploratory only.

MODEL 2 (archived: model2_Group_Age_Gender)
  Design: ~ Group + Age + Gender | 75 subjects
  Result: 13 significant species (C:2, D:11)
  Note: Adding Gender over-constrained the model because sex
        is collinear with Group (IBS-C ~all female). Rejected.

MODEL 3 = FINAL (folder: 03_Metagenomics_FINAL_model)
  Design: ~ Group + Age  | 74 subjects (1 male IBS-C excluded)
  Result: 49 significant species (C:6, D:43)
  Rationale: Excluding the single IBS-C male makes that group
        all-female, matching the transcriptomic layer's structure,
        so sex no longer confounds and age-adjustment suffices.
  Main finding: Streptococcus (18 species, all up in IBS-D) - SIBO pattern.

Organized on: 2026-07-27


ADDENDUM (2026-08-10) - WHY THE MODEL WITH FEWER FINDINGS WAS CHOSEN
=====================================================================
Model 1 (Group only) reported 79 significant species; the final model
(Group + Age, one IBS-C male excluded) reports 49. Selecting the model
with fewer discoveries requires explicit justification.

1. Model 1 is unadjusted. Age is an established correlate of gut
   microbiome composition, so an unadjusted comparison cannot separate
   group effects from age structure. Its larger yield is expected of an
   unadjusted model and is not evidence of greater sensitivity.

2. Model 2 (Group + Age + Gender) collapsed to 13 species because sex was
   collinear with group: IBS-C comprised 21 females and 1 male
   (chi-square p = 0.035). Adjusting for a variable that is nearly
   determined by group removes the group effect itself. The collapse
   reflects collinearity, not absence of signal.

3. The final model resolves this structurally rather than statistically.
   Excluding the single IBS-C male makes that group entirely female,
   matching the transcriptomic layer, so sex is no longer a confounder
   and age adjustment alone is sufficient. This is a documented
   analytical decision, not missing data.

4. Direction of bias: the final model is the CONSERVATIVE choice. Of the
   three, it is neither the most permissive (Model 1) nor the
   over-constrained one (Model 2). Findings reported here are therefore
   unlikely to be artefacts of insufficient adjustment.

5. Robustness of the main finding: verified quantitatively. Model 1
   (unadjusted) recovers 28 significant Streptococcus species, ALL assigned
   to IBS-D and ALL elevated, with none in IBS-C and none in the opposite
   direction. The final model reduces this to 18 species - the expected
   effect of adjustment - without altering the pattern. The headline
   finding is therefore independent of model choice. All rejected models
   are archived in _sensitivity/ for audit.

COST OF THE DECISION: excluding one participant reduces the metagenomic
  sample from 75 to 74 and means no male IBS-C patient is represented.
  Findings specific to IBS-C therefore apply to female patients only,
  and this restriction is stated as a limitation.


ADDENDUM (2026-08-10) - IBS-C SIGNATURE: SCATTERED BUT SUBSTANTIAL
===================================================================

Prompted by the ordered coefficient plot (Fig3), which showed that IBS-C species
are NOT confined to the bottom of the ranking as informally assumed.

THE SIX IBS-C SPECIES (verified from sig_final):
  Fimadaptatus faecigallinarum      +2.67  q=0.044  Christensenellales / CAG
  Ruthenibacterium sp013316265      +2.41  q=0.016  Oscillospirales / Ruminococcaceae
  Aphodomonas merdavium             +1.98  q=0.041  Christensenellales / CAG
  Roslinia caecavium                +1.88  q=0.041  Monoglobales / Monoglobaceae
  Fournierella excrementigallinarum +1.30  q=0.047  Oscillospirales / Ruminococcaceae
  Prevotella sp900769055            -2.08  q=0.016  Bacteroidales / Bacteroidaceae

1. IBS-C IS NOT A WEAK LAYER. Five of six species are ELEVATED, and the single
   largest coefficient in the entire analysis (+2.67) belongs to IBS-C, not IBS-D.
   The difference between subtypes is one of COHERENCE, not effect size.

2. TAXONOMIC COHERENCE - QUANTIFIED:
     IBS-D (43 spp.): Lactobacillales 23, Lachnospirales 13, Oscillospirales 3,
                      Coriobacteriales 2, Bacteroidales 1, Erysipelotrichales 1
                      -> 84% in just two orders.
     IBS-C (6 spp.) : Christensenellales 2, Oscillospirales 2, Monoglobales 1,
                      Bacteroidales 1 -> four orders for six species.
   IBS-D is a clustered signature; IBS-C is dispersed.

3. CORRECTION OF AN EARLIER READING: it was initially suggested that five of the
   six IBS-C species belong to one family, inferred from their poultry-derived
   epithets (faecigallinarum, merdavium, caecavium, excrementigallinarum). This is
   WRONG - they span three distinct orders. The naming similarity reflects the
   isolation source of the GTDB reference genomes, not phylogenetic proximity.

4. WHAT IS ACTUALLY SHARED - functional, not taxonomic: all five elevated species
   are strictly anaerobic fermentative Clostridia from families associated with
   short-chain fatty acid production (Ruminococcaceae, Christensenellaceae,
   Monoglobaceae). This is notable because the metabolomic layer shows SCFAs
   REDUCED by 27-32% in the same subtype.
   INTERPRETATION IS NOT CLAIMED. A higher abundance of potential producers
   alongside lower faecal concentrations is only apparently contradictory: faecal
   SCFA levels reflect the balance of production and absorption, and slower transit
   in IBS-C allows more absorption. This is a hypothesis consistent with both
   layers, NOT a demonstrated mechanism, and no test of it is reported here.

CAVEAT: re-reading of existing results; no new analysis, no addition to the test count.
