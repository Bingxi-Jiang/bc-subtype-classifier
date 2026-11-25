# Interpretation of consistently selected genes

Working notes for the report's Discussion. Evidence levels: **data** = computed from our TCGA-BRCA cohort; **lit** = published literature. Interpretations are hypotheses about why a gene is informative for reproducing PAM50 labels, not causal or clinical claims. Every citation below was checked against its Europe PMC / NCBI record (Nov 2026); "abstract" means only the abstract was read.

## Which genes

- **L1 minimal signature:** genes in the L1 top 25 in at least 20 of 25 training splits (`docs/validation/selection-frequency.tsv`).
- **ANOVA comparison:** genes in the ANOVA top 10 in all 25 splits.

Per-subtype values are mean z-scores of log2(RSEM+1) expression over all 945 patients (>0 = above the cohort average). ★ = PAM50 gene.

## Four functional groups chosen by L1

| Group | Genes (L1 splits/25) | Data: expression pattern | Role in classification |
| --- | --- | --- | --- |
| A. Luminal / ER | FOXA1★ (25), MLPH★ (25), XBP1 (24) | Luminal A/B ≈ +0.4, Basal ≈ −1.9 | Luminal vs Basal |
| B. Basal | FOXC1★ (25), SFRP1★ (25), TRIM29 (24) | Basal +1.0 to +1.8 | Confirms Basal-like |
| C. Proliferation | CENPA (25), CENPN (25), CEP55★ (25), AURKA (25), TIMELESS (25) | Luminal A lowest (≈ −0.7); Luminal B, HER2-E, Basal higher | Luminal A vs Luminal B |
| D. HER2-enriched | ERBB2★ (25), STARD3 (21), MPHOSPH6 (25), ZP2 (25), SILC1/LOC150622 (24), G6PD (24), KRT20 (20) | HER2-E +1.5 to +1.8 | Identifies HER2-enriched |

ANOVA's stable top 10 (FOXA1, ESR1, MLPH, AGR3, FOXC1, SPDEF, XBP1, …) are almost entirely group A plus FOXC1: all variations of the Luminal-vs-Basal axis. This matches ANOVA's weak Luminal B F1 at small k (0.60 at k=10 vs 0.70 for L1; LR).

## Group D: two mechanisms

**Data, correlation.** Spearman correlation with ERBB2 within HER2-enriched patients (n=78): STARD3 0.89, GRB7 0.91; MPHOSPH6 0.14, ZP2 −0.12, SILC1 0.31, G6PD −0.13, KRT20 −0.01. The non-amplicon genes are also weakly correlated with each other (|r| ≤ 0.5). ZP2 and KRT20 are near zero in many HER2-E tumors and high in a subset (90th percentile ≈ 8 log2 units vs a median of 0 in other subtypes).

**Data, ERBB2 alone is not enough.** Using the 95th percentile of ERBB2 in Luminal A as a cutoff, 72% of HER2-E tumors are ERBB2-high (Luminal B 12%, Luminal A 5%, Basal 2%), and only 51% of ERBB2-high tumors are HER2-E (n=109). This is a descriptive expression cutoff, not clinical HER2 testing.

1. **ERBB2 amplicon (STARD3, GRB7).** In ERBB2-amplified tumors, ERBB2–C17orf37–GRB7 forms the minimal common 17q12–q21 amplicon. STARD3 lies in the broader amplicon and was amplified with upregulated expression in 82% of amplified samples (lit: Sircoulomb et al. 2010, *BMC Cancer* 10:539, PMC2958950, full text). GRB7 modulates ERBB2 signaling as an ERBB2-dependent oncogene (Saito et al. 2012, *FEBS Lett* 586:1708, PMID 22584052, abstract). L1 keeps ERBB2 every time but STARD3 only in 21/25 splits and GRB7 rarely, consistent with L1 dropping near-duplicate information.
2. **ERBB2-independent HER2-E markers (MPHOSPH6, ZP2, SILC1, G6PD, KRT20).** High in HER2-E but uncorrelated with ERBB2. The PAM50 HER2-enriched class is an expression phenotype, not the same as clinical HER2 positivity: in TCGA and METABRIC, 64.6% of HER2-E tumors were clinically HER2-positive, and all intrinsic subtypes occurred regardless of clinical HER2 status (lit: Prat et al. 2014, *JNCI* 106(8):dju152, PMC4151853). Our 72% ERBB2-high figure is consistent with this. A classifier therefore needs HER2-E signals beyond ERBB2. Individual evidence:
   - **G6PD:** by IHC in 348 breast cancers, pentose-phosphate-pathway proteins differed by molecular subtype, with G6PDH higher in the HER-2 type (lit: Choi et al. 2018, *Dis Markers* 2018:9369358, PMC5845514, full text). Separately, high G6PD/TKT expression was associated with worse survival in public datasets (Benito et al. 2017, *Oncotarget* 8:106693, PMC5739767). This is the best-supported gene in this group.
   - **KRT20:** CK20 was found more often than previously reported in breast cancer and "may be associated with the molecular apocrine subtype"; the paper also suggests it "might" indicate poor prognosis (lit: Shao et al. 2012, *Virchows Arch* 461:313, PMID 22851038, abstract). Molecular apocrine tumors are often HER2-positive, so this is indirect support only.
   - **ZP2:** zona pellucida glycoprotein (oocyte). The Human Protein Atlas places it in an uncharacterized, breast-cancer-associated expression cluster with low confidence ([HPA ZP2](https://www.proteinatlas.org/ENSG00000103310-ZP2/cancer)). **No breast cancer study found.**
   - **MPHOSPH6:** **no breast cancer literature found.**
   - **SILC1** (Entrez 150622, formerly LOC150622/LINC01105; long intergenic non-coding RNA at 2p25.2): the dataset uses the outdated LOC symbol. **No breast cancer literature found.**

   For ZP2, MPHOSPH6 and SILC1 we can only report a reproducible statistical association with the HER2-E label in this cohort. Possible explanations include HER2-E subgroup biology, technical/batch effects, or chance, and the data here cannot distinguish them. They are candidates for validation in an independent cohort (e.g. METABRIC).

## Groups A–C: well established

- **FOXA1 (A):** a luminal marker; necessary for estrogen-receptor binding and transcriptional activation of luminal genes, with mRNA levels proportional to ER (lit: Metovic et al. 2022, *Cancers* 14:4699, PMC9564251, full text).
- **XBP1 (A):** GATA3, HNF3A (FOXA1) and XBP1 were identified by microarray as transcription-factor genes co-expressed with ESR1 in breast cancer (lit: Lacroix & Leclercq 2004, *Mol Cell Endocrinol* 219:1–7, PMID 15149721, abstract).
- **FOXC1 (B):** consistently overexpressed in basal-like vs other molecular subtypes (lit: Ray et al. 2010, *Cancer Res* 70:3870, PMID 20406990, abstract); a FOXC1-based assay identified basal-like cancer with AUC 0.88 (Jensen et al. 2015, *JNCI* 107(8), PMC4554196, abstract).
- **SFRP1 (B):** significantly overexpressed in triple-negative vs other subtypes (lit: Bernemann et al. 2014, *Mol Cancer* 13:174, PMID 25033833, abstract).
- **TRIM29 (B):** listed among basal-epithelial genes in the NCI definition of basal-like breast carcinoma ([MedGen C3642347](https://www.ncbi.nlm.nih.gov/medgen/C3642347)); its regulation clusters by PAM50 subtype, with robust expression in the ER-negative basal-like group (Avraham et al. 2014, *PLoS One* 9:e91805, PMC3961270, full text).
- **Proliferation (C):** in 3,093 tumors, AURKA outperformed Ki67 and other proliferation markers as an independent prognostic marker in ER-positive breast cancer (lit: Ali et al. 2012, *Br J Cancer* 106:1798, PMC3365239). No source was found specifically for CENPA, CENPN or TIMELESS in the Luminal A/B split; we report their expression pattern (data) only. That Luminal B shows higher proliferation-gene expression than Luminal A is standard background (cite Parker et al. 2009 from `docs/REFERENCES.md`).

## Points for the Discussion

1. L1 recovers one or more markers for each subtype-defining axis without biological input; ANOVA concentrates on the single strongest axis. This explains L1 reaching full-gene macro-F1 with 25 genes versus 100 for ANOVA.
2. L1 picks one representative among redundant genes (ERBB2 over GRB7, CENPA/AURKA over MKI67/BIRC5, FOXA1 over ESR1). Individual L1 genes are less stable across splits (Jaccard 0.64 vs 0.88 at k=25), but the *groups* are stable. Report stability at the group level.
3. Selected genes that are not PAM50 members (e.g. STARD3, AURKA, XBP1) are often known biology from the same pathways, so non-overlap with PAM50 does not mean novel biology.
4. Because HER2-E is not equivalent to ERBB2 status (about one-third of HER2-E tumors are not clinically HER2-positive), L1 adds ERBB2-independent HER2-E markers. G6PD has supporting literature; ZP2, MPHOSPH6 and SILC1 do not, so present them as observations needing external validation, not as discoveries.

## Citation check record (Nov 2026)

Removed or corrected after checking: PMC6042281 (Desai et al. 2018) reports the share of HER2-positive cases that are HER2-E, not the reverse, so it was replaced by Prat et al. 2014. PMC7850037 (Prat et al. 2020) covers dual HER2 blockade in HER2-positive patients only and was dropped. PMC2716509 (Albergaria et al. 2009) does not mention XBP1 and was replaced by Lacroix & Leclercq 2004. An unverified Kao & Pollack 2006 citation was replaced by Sircoulomb et al. 2010, which was read in full text.
