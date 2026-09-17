# What the raw data looks like — GSE68086, 133-sample subset

Findings behind the four plots produced by `functions/03_plot raw data.R`.

**What was looked at.** The 133 samples chosen in step 02: healthy controls (HC, 55
samples), glioblastoma (GBM, 38) and lung cancer (Lung, 40), all from Batch02, Batch03
and Batch04. Nothing has been normalised and no genes have been removed yet, so these
are raw read counts. Everything below is about the *quality and shape of the data*, not
about biology.

Figures are in `data/derived/figures/`. Re-create everything with:

```bash
Rscript "functions/03_plot raw data.R"
```

---

## 1. Samples were sequenced very unevenly, and it follows the batch

![depth by batch and group](data/derived/figures/01_depth_by_batch_group.png)

The deepest sample has **18 times** more reads than the shallowest (7.27 million against
0.40 million; median 1.71 million).

Almost all of that spread comes from which batch a sample was processed in:

| Batch | Median reads |
|---|---|
| Batch02 | 2.56 M |
| Batch03 | 2.13 M |
| Batch04 | 0.99 M |

That difference is very unlikely to be chance (p = 4 × 10⁻¹⁰).

**The part that matters most:** the healthy controls are the shallowest group. Their
median is 1.04 million reads against 1.97 million for lung cancer. Two things cause it —
44% of the HC samples sit in Batch04, the shallow batch (against 29% of GBM), and HC
inside Batch04 is shallower still:

| Median reads | Batch02 | Batch03 | Batch04 |
|---|---|---|---|
| HC | 2.52 M | 2.69 M | **0.59 M** |
| GBM | 2.06 M | 1.38 M |1.36 M  |
| Lung | 2.64 M | 3.12 M | 1.56 M |

Seven of the eight shallowest samples in the whole subset are healthy controls from
Batch04.

**"Shallowest group" is not the same as "smallest group".** HC is in fact the *largest*
group here, and still the shallowest, because the two words count different things.
Group size counts people; depth counts reads inside one person's sample:

| | Samples in group | Median reads **per sample** | All reads in group |
|---|---|---|---|
| HC | **55** (most) | **1,037,596** (fewest) | 89,157,880 |
| GBM | 38 | 1,634,691 | 67,123,400 |
| Lung | 40 | 1,971,280 | 92,952,890 |

HC has 15 more samples than Lung and still fewer reads in total, because each HC sample
is so much thinner.

in practice: so the depth of the columns decides which rows look real, and the 24 healthy controls sitting at 0.59M in Batch04 are the columns most likely to turn a real gene into a row of zeros.

**How depth is worked out.** The counts table has one column per sample. A sample's depth
is simply that column added up — all 57,736 gene values summed into one number. In the
script that is one line, `ss$lib_size <- colSums(counts)`. "Depth", "library size",
"total reads" and "the sum of that sample's raw counts" are four names for the same
number. Group size never touches the counts table at all; it is just how many rows of
the sample sheet carry that group name.

**Why this matters.** If you compared cancer against healthy without accounting for the
batch, part of what you found would simply be that the healthy samples were sequenced
less. It would look like a biological difference and it would not be one. This is the
reason step 04 fits `~ batch + group` instead of group alone, and the reason step 02
insisted that all nine group-and-batch combinations contain samples — if a group were
missing from a batch, the model could not tell the two effects apart.

---

## 2. More reads did not find more genes

![genes detected against depth](data/derived/figures/02_genes_vs_depth.png)

You would expect a sample with more reads to detect more genes. It does not happen here.
The correlation between reads and genes found is **−0.01** (p = 0.92) — in other words,
no relationship at all. Samples find between 3,817 and 12,889 genes, median 9,074, and
which end of that range a sample lands on has nothing to do with how deeply it was
sequenced.

The reason is that extra reads keep landing on genes that were already found. Platelet
samples contain a small number of very abundant transcripts (see finding 3), so
sequencing deeper mostly counts those same transcripts again.

**Why this matters.** Sequencing depth is not a measure of sample quality in this
dataset. A shallow sample is not a worse sample. So removing samples because they are in
the bottom slice by depth would throw away usable data — and, given finding 1, it would
remove mostly healthy controls, damaging the group we compare everything against.

---

## 3. A small number of genes hold most of the reads

![read concentration](data/derived/figures/03_read_concentration.png)

Ranking genes by how many reads they collected. These figures treat all 133 samples as
one pool, which is the **grey line** on the plot:

| | Share of all reads |
|---|---|
| Top 1 gene | 7.9% |
| Top 10 genes | 28.7% |
| Top 100 genes | 57.6% |

Only **57 genes** are needed to account for half of every read. 1,434 genes cover 90%,
and 6,268 cover 99%. The remaining ~51,000 genes share the last 1%.

The three coloured lines are the same calculation done inside each group, and they reach
the halfway mark at different points:

| Curve | Genes needed for half of its reads |
|---|---|
| HC (blue) | 76 |
| GBM (orange) | 56 |
| Lung (aqua) | 42 |
| All 133 pooled (grey) | 57 |

So read a number off the curve you mean. The pooled figure is not the middle of the
three: pooling counts reads, not samples, and the deeper groups contribute more reads,
which pulls the grey line up towards Lung. Lung's reads are the most concentrated of the
three, HC's the least.

**Why this matters.** It explains finding 2 — there is not much new to find, so extra
reads add little. It also means normalisation deserves attention: when a handful of
genes carries a third of the signal, a shift in just those genes moves a sample's whole
total.

---

## 4. Over half the genes are empty, and concentration varies a lot per sample

![top 10 share per sample](data/derived/figures/04_top10_share_per_sample.png)

**31,559 of the 57,736 genes (54.7%) have zero reads in all 133 samples.** They are rows
of nothing. Only 8,446 genes appear in at least half the samples.

Concentration also differs strongly between individual samples. The share of a sample's
reads sitting in the ten most abundant genes runs from **11.9% to 50.3%** (median 28.5%).
By group the medians are close — HC 26.7%, GBM 29.5%, Lung 29.3% — so this is
sample-to-sample variation, not a group difference.

**Why this matters.** Removing low-information genes is necessary, not optional: more
than half the table is empty before you start. Step 04 does this with `filterByExpr`,
which keeps 6,994 genes.

---

## What follows from all this

1. **Keep batch in the model.** Depth is driven by batch, and the healthy controls are
   concentrated in the shallowest batch. Dropping batch would let a processing
   difference be read as a disease difference.
2. **Do not remove samples based on depth.** Depth does not track quality here
   (finding 2), and a depth cut would mostly remove healthy controls (finding 1).
3. **Do remove low-information genes.** Over half of them are empty (finding 4).
4. **Normalise for library size.** An 18-fold range cannot be compared raw.
5. **Keep an eye on the most abundant genes.** With 10 genes holding 28.7% of reads,
   they have real influence over the size factors used to normalise (finding 3).

## Open question

Samples differ a lot in how concentrated their reads are (11.9% to 50.3%, finding 4) and
this is not explained by group or by depth. It may be a real difference in sample
quality that library-size normalisation alone does not correct. Worth checking whether
those samples stand out once the data is normalised and plotted in step 04.
