# GSE68086 — Review of `00b_series_df_dictionary.ipynb`, gene-filter decision, and comparison with the source paper

*Written 2026-09-14. Reviewer notes on the state of the pipeline after commit `5058598`. Every number below was computed on the actual data; the read-only script [`functions/00c_filter_comparison.R`](functions/00c_filter_comparison.R) reproduces all of them (run from the project root: `Rscript "functions/00c_filter_comparison.R"`). The notebook and the pipeline scripts 01–04 were **not** modified by this review.*

---

## Summary

The notebook is sound and its conclusions hold up against the data. Three things need fixing before any DE model, and the filter question has a clear answer with measured evidence:

- keep `filterByExpr(group = ...)` in 04,
- do **not** add a manual gene filter on top,
- move the one manual step that *is* justified (dropping the author-excluded sample `VU398-GBM`) into 02, with a comment quoting the GEO note as the reason.

---

## 1. Review of 00b

### What is right (all verified)

- Parsing the six `!Sample_characteristics_ch1*` columns **by key** instead of by position. 245 samples use one key order, the 40 later-submitted (MGH) samples use a rotated order, so positional indexing mixes `batch` with `tissue`.
- The join to the counts matrix is `!Sample_source_name_ch1`, not the GSM accession (285/285 match by source name, 0 by GSM).
- The batch × group crosstab and the confounding it shows (Batch05/06 contain no healthy controls).
- The QC flags hidden in `!Sample_description`: 2 samples excluded by the authors for too few intron-spanning reads (`VU258-CRC`, `VU398-GBM`) and 5 HD samples merged *in silico*.
- Deriving depth (`lib_size`, `genes_detected`) from the counts, since the series matrix carries no depth field.

### What needs attention

1. **The gotcha "drop the 2 author-excluded samples" is not implemented anywhere.**
   `VU398-GBM` (Batch02) passes the subset rule in 02 and is inside the 133-sample object.

   | | VU398-GBM | GBM median in subset |
   |---|---|---|
   | library size | 396,494 reads (rank 1 of 133, the shallowest) | 1,634,691 |
   | genes detected | 4,643 | 8,874 |
   | median Spearman correlation to other samples | 0.636 | 0.680 (5th percentile 0.598) |

   Its correlation is within the normal range, so the reason to drop it is **provenance**, not outlier behaviour: the authors excluded it, and keeping it makes our population differ from the paper's 283. Dropping it moves `filterByExpr` from 6,994 to 7,082 genes; all 6,994 survive, so nothing is lost.

2. **Two sample sheets are in play.**
   `01_sample sheet creation.R` writes `GSE68086_sample_sheet.csv` (columns `sample, group, batch, subclass, patient`). `02_make subset.R` reads `GSE68086_sample_sheet_claude.csv`, a different file with different column names (`source`, `cancer type`, …). 01 is also not runnable as committed: the GEOquery lines are commented out, so `raw` is never defined.
   → Pick one source of truth. 01 should write `sample, gsm, group, batch, subclass, patient, submission_date, author_excluded`; 02 should read that file.

3. **Stale outputs in the notebook.** The cell that reads the counts matrix (`counts = pd.read_csv(COUNTS, index_col=0)`) displays the depth summary that belongs to the following cell. Re-run top to bottom before committing. The `crosstab_batchslice` cell assigns a new column into a slice of `crosstab` (`crosstab.iloc[:, 1:4]`); use `.copy()` to avoid pandas' chained-assignment behaviour.

4. **Small text mismatches.** The dictionary says 5 samples lack a tumour tag; the check finds 7 (the extra two are `VU383Platelet-hiseq` and `VU394Platelet-hiseq`). Worth adding to the dictionary: GEO's 285 samples = the paper's 283 + the 2 excluded, which is why the paper reports CRC 41 and GBM 39 while GEO gives 42 and 40.

---

## 2. The filter decision, with numbers

All numbers on the 133-sample Lung / GBM / HC, Batch02–04 subset produced by 02.

| Rule | Genes kept |
|---|---|
| `filterByExpr(y, group = group)` | **6,994** |
| `filterByExpr(y, design = ~ batch + group)` | 8,079 |
| manual: ≥10 reads in ≥38 samples (04 before the switch) | 6,109 |
| CPM > 1 in ≥38 samples (`subset/make_subset.py`) | 9,666 |
| ≥5 reads in ≥10 % of samples ("paper-like") | 9,903 |
| any read at all | 26,177 (of 57,736) |

### Why `filterByExpr(group=)` is the right choice

**It is fair to the shallow control libraries.**

| group | median library size |
|---|---|
| HC | 1,037,596 |
| GBM | 1,634,691 |
| Lung | 1,971,280 |

HC depth is about half of Lung. A raw threshold of 10 reads is therefore twice as hard for an HC sample to pass. `filterByExpr` sets its cutoff in CPM: 10 / median-library-size × 1e6 = **5.86 CPM** here, which is ≈ 6 reads in a median HC library and ≈ 12 reads in a median Lung library. It asks the same *relative* question of every sample.

Evidence: the two rules agree on 6,108 genes. The **886 genes** that `filterByExpr` keeps and the raw rule drops are HC-expressed genes:

| | HC | GBM | Lung |
|---|---|---|---|
| mean logCPM of those 886 genes | 2.23 | 1.62 | 1.00 |
| median number of samples with ≥10 reads | 15 | 9 | 5 |

The raw rule silently biases the analysis **against genes that are down in cancer** — exactly the class the paper reports 793 of. Only 1 gene goes the other way.

**`group=` rather than `design=`.** With `design`, edgeR derives the minimum group size from the hat matrix, i.e. from the batch × group cells: 1/max(hat) ≈ **19 samples**. A gene expressed in a single batch-by-group cell then gets through, which is precisely the kind of batch-driven gene we want to avoid. With `group=`, the minimum is the smallest group (38) after edgeR's tolerance (`large.n = 10`, `min.prop = 0.7` → 29.6 samples). Our contrasts are group-level, so the group rule matches the question.

**Not the looser rules.** CPM > 1 and "≥5 in ≥10 %" keep about 3,000 more genes that reach a few reads in a minority of samples. In a matrix that is 85 % zeros, those genes destabilise dispersion estimation and add multiple-testing burden without being testable. If a sensitivity analysis is wanted later, re-run the DE under "≥5 in ≥10 %" and check whether the top list changes.

**No manual gene filter on top.** Nothing in the data argues for one. The one manual step that is justified is at the *sample* level (`VU398-GBM`), and it belongs in 02 with the GEO note quoted as the reason.

**Consequence for 04.** After dropping `VU398-GBM` the subset is 132 samples (GBM 37) and `filterByExpr(group=)` keeps 7,082 genes; update the comment in `04_DESeq2.R` accordingly. `filterByExpr(counts(dds), group = dds$group)` on a plain matrix is fine (library sizes are taken as column sums).

### Composition facts to carry forward

These shape how HC-vs-cancer results should be read.

| measure (median by group) | HC | GBM | Lung |
|---|---|---|---|
| TMM normalisation factor (after filterByExpr) | 1.22 | 0.89 | 0.82 |
| ribosomal-protein (RPS/RPL) share of reads | 6.3 % | 2.0 % | 1.3 % |
| fraction of reads retained by filterByExpr | 98.9 % | 99.4 % | 99.6 % |

- Normalisation is absorbing a **group-wide composition shift** (HC factors well above 1, Lung below 1). Cancer-vs-HC lists will be long and will partly reflect library character, not only tumour biology. The Lung-vs-GBM contrast is where disease specificity lives.
- Top-10 genes by total reads carry 28 % of all reads: B2M 7.9 %, PPBP 6.0 %, TMSB4X 4.5 %, ACTB 2.8 %, CLU 1.8 %, F13A1 1.4 %, FTL 1.2 %, SPARC 1.2 %, HBB 1.0 %, SH3BGRL3 0.9 %. All platelet/blood transcripts, which confirms the matrix is what it claims.
- **Mitochondrial genes are absent** from the matrix (0 of 5 canonical MT gene IDs present): intron-spanning counting removes intronless genes. An MT-fraction QC is therefore impossible; use the ribosomal share as the complexity/composition covariate instead.
- Annotation: 23,875 of 57,736 Ensembl-75 IDs have no symbol in `org.Hs.eg.db`, but only 233 of the 6,994 kept genes lack one. Ensembl v75 (hg19) is the vintage to use for biotype if needed.
- **Sex is inferable but noisy.** Not recorded in GEO. Y-linked reads (RPS4Y1 + DDX3Y + UTY) are > 0 in 90 of 133 samples (HC 38/55, GBM 29/38, Lung 23/40). XIST is detected in only 36 samples and cannot be used. A zero-Y call in a very shallow library is uncertain.

---

## 3. What the paper did, and what to bring back

Source paper: Best MG et al., *Cancer Cell* 2015; 28:666–676 (PMID 26525104, [PMC4644263](https://pmc.ncbi.nlm.nih.gov/articles/PMC4644263/)).

| aspect | Best et al. 2015 | this pipeline |
|---|---|---|
| population | 283 samples (55 HC, 228 cancer); 2 excluded for too few spliced reads | 285 in GEO; one excluded sample still in the subset |
| depth | mean ≈ 22 M raw reads/sample; intron-spanning reads only | median 2.0 M intron-spanning reads/sample in the matrix |
| gene filter | "transcripts with low expression (<5 reads in all samples) were excluded" → 5,003 genes | `filterByExpr(group=)`, ≈ 7,000 on the subset |
| normalisation | "as described in the Supplemental Procedures" (one summary attributes it to edgeR; the supplement could not be opened to confirm) | TMM (edgeR) / median-of-ratios (DESeq2) |
| cancer vs HC | 1,453 up and 793 down of 5,003 at FDR < 0.001 (45 % of genes) | to be measured in 05 |
| batch | counts voom-transformed "using sequencing batch and sample group as variables" (for the pathway/CAGE analysis) | `~ batch + group` |
| sex | male HC excluded from the breast-cancer clustering "to avoid sample bias due to gender-specific platelet mRNA profiles" | infer sex from Y genes, test as covariate |
| age | HC aged 21–64, not age-matched to patients | not in GEO metadata, cannot be adjusted |
| classifier | SVM with LOOCV (e1071); ANOVA feature selection with logCPM > 3 and FDR thresholds; 1,072 genes; training n = 175, validation n = 108; accuracy 95 % / 96 %; multiclass 71 % | later, with glmnet |
| biology reported | up: vesicle-mediated transport, cytoskeletal protein binding; down: RNA processing and splicing | to check by enrichment |

### Ideas to take from this

- **Expect huge HC-vs-cancer lists.** The paper found 45 % of its genes differentially expressed. That is a global shift, consistent with the TMM and ribosomal differences above. Treat Lung-vs-HC and GBM-vs-HC as "library character + biology"; treat Lung-vs-GBM as the disease-specific list.
- **The paper's 5,003 cannot be reproduced from its text.** On the same 283 samples: `filterByExpr(group=, 7 groups)` gives 9,744; "≥5 reads in ≥10 % of samples" gives 9,893; "≥5 reads in ≥50 % of samples" gives 5,817. Their filter was stricter than anything considered here. In the authors' 2017 follow-up ([Best et al., Cancer Cell 2017, PMC6381325](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6381325/)) genes with "<30 intron-spanning reads in >90 % of the cohort" were excluded (4,722 genes) and confounders (patient age, blood-storage time) were corrected with RUVg. That tells you which confounders the authors themselves considered real once they revisited the design.
- **Batch as a covariate is the authors' own choice too**, which supports the `~ batch + group` design.
- **Use their enrichment result as an external sanity check.** If the Lung-vs-HC list shows vesicle-transport genes up and RNA-processing/splicing genes down, the pipeline is at least reproducing the paper.
- **A published critique** (Chakraborty S., bioRxiv 2017, single author, not peer reviewed, [doi 10.1101/146134](https://www.biorxiv.org/content/10.1101/146134.full.pdf)) argues that the paper's discriminator sets contained very low-count genes (TRAT1, MET; standard deviation far above the mean) alongside bona fide platelet genes (F13A1). Whatever one thinks of its conclusion, it points to two checks worth building in: rank by **shrunken** fold changes, and verify that DE lists are not simply the top platelet genes and their complement.

---

## 4. What comes next, in order

1. **01** — make it runnable and the single sample-sheet writer, with an `author_excluded` flag (from `!Sample_description`) and `submission_date`. **02** — read that file, drop the 2 flagged samples (subset becomes 132, GBM 37), and store `lib_size`, `genes_detected`, `rp_fraction` and a Y-based `sex` call in `ss` so 03/04 stop recomputing them.
2. **03** — add an MDS plot on logCPM (colour = group, shape = batch), the ribosomal share by group, and a sample-correlation heatmap. Lowest-correlation samples at present: `HD-45-1`, `VU280-GBM`, `HD-10`, `HD-38-2`, `HD-1`, `Lung-008` (all within range, none an outlier).
3. **04** — keep `filterByExpr(group=)`, write the depth-fairness rationale as comments, update the gene-count comment (7,082), and build a `DGEList` with TMM (`normLibSizes`) next to the `dds` so TMM factors and DESeq2 size factors can be compared by group.
4. **05** — edgeR quasi-likelihood: `estimateDisp(robust = TRUE)`, `glmQLFit(robust = TRUE)`, contrasts Lung−HC, GBM−HC, Lung−GBM plus the global group test. Look at p-value histograms first, then the fraction of genes DE (compare with the paper's 45 %), then the enrichment sanity check. DESeq2 follows as 06 and the tool comparison as 07.
5. **If a classifier is the end goal** — hold out ≈ 30 % of patients *before* step 04 (the paper split 175 / 108) and keep any DE-based gene selection inside the cross-validation folds.

---

## Appendix — files touched by this review

- Added: [`functions/00c_filter_comparison.R`](functions/00c_filter_comparison.R) — read-only, prints every number in sections 1–3; run from the project root.
- Added: this document.
- Not modified: `functions/00b_series_df_dictionary.ipynb`, `functions/01–04`, any data file.

Sources: [Best et al. 2015, Cancer Cell (PMC4644263)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4644263/) · [Best et al. 2017, Cancer Cell (PMC6381325)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6381325/) · [Chakraborty 2017, bioRxiv critique](https://www.biorxiv.org/content/10.1101/146134.full.pdf).
