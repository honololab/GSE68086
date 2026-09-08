"""Subset GSE68086 TEP counts to two diseases + healthy controls."""
import pandas as pd, numpy as np, re, os

DISEASES = ["GBM", "Lung"]        # <-- swap here, e.g. ["Pancreatic", "CRC"]
CONTROL  = "Healthy"
OUTDIR   = "subset"

RULES = [("Healthy",r"\bhd\b|^hd|control|healthy"),("Breast",r"breast|brca"),
         ("CRC",r"crc|coad|colon|rectal"),("GBM",r"gbm|glio"),("NSCLC",r"nsclc|lung"),
         ("Pancreatic",r"panc|paad"),("Hepatobiliary",r"liver|chol|hepat"),
         ("Unknown",r"type-unknown")]

def label(s):
    sl = s.lower()
    for n, p in RULES:
        if re.search(p, sl): return n
    return "Unlabeled"

def site(s):
    sl = s.lower()
    if sl.startswith("mgh"): return "MGH"
    if re.match(r"^vu\d|^vu-|^vumc", sl): return "VUMC"
    return "unprefixed"

def subtype(s):
    m = re.findall(r"(?i)\b(her2|kras|egfr|alk|braf|pik3ca|wt)\b", s)
    return m[0].upper() if m else "NA"

os.makedirs(OUTDIR, exist_ok=True)
df = pd.read_csv("GSE68086_TEP_data_matrix.csv", index_col=0)

keep_classes = DISEASES + [CONTROL]
lab = pd.Series({c: label(c) for c in df.columns})
sel = lab[lab.isin(keep_classes)].index
sub = df[sel]

meta = pd.DataFrame({
    "sample_id": sel,
    "group":     [label(c) for c in sel],
    "site":      [site(c) for c in sel],
    "subtype":   [subtype(c) for c in sel],
    "lib_size":  sub.sum(axis=0).values,
}).set_index("sample_id")
meta["condition"] = np.where(meta.group == CONTROL, "control", "cancer")

# --- edgeR-style CPM filter: CPM > 1 in >= smallest-group-many samples ---
cpm = sub / sub.sum(axis=0) * 1e6
n_min = meta.group.value_counts().min()
keep = (cpm > 1).sum(axis=1) >= n_min
filt = sub[keep]

sub.to_csv(f"{OUTDIR}/counts_raw.csv")
filt.to_csv(f"{OUTDIR}/counts_filtered.csv")
meta.to_csv(f"{OUTDIR}/metadata.csv")

print(f"samples kept : {sub.shape[1]} of {df.shape[1]}")
print(f"genes  raw   : {sub.shape[0]:,}")
print(f"genes  filter: {filt.shape[0]:,}  (CPM>1 in >={n_min} samples)")
print("\ngroup x site:")
print(pd.crosstab(meta.group, meta.site).to_string())
print("\nlib size by group:")
print(meta.groupby("group").lib_size.describe()[["min","50%","max"]].round(0).to_string())
print(f"\nwrote -> {OUTDIR}/counts_raw.csv, counts_filtered.csv, metadata.csv")
