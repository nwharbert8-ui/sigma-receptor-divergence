"""
01_genome_wide_correlations.py
==============================
Computes genome-wide Pearson correlations for all 7 target genes against
16,212 expressed genes in GTEx v8 Brain – Frontal Cortex (BA9; n=209).

Methodology (Manuscript Section 2.1–2.2):
  - Gene filtering: median TPM >= 1.0 (retains 16,212 genes)
  - Transformation: log2(TPM + 1)
  - Metric: Pearson correlation coefficient
  - Network definition: top 5% (811–812 genes per network)

Target genes:
  SIGMAR1 (sigma-1 receptor, ER chaperone) — primary target
  TMEM97  (sigma-2 receptor, cholesterol) — primary target
  PELO    (ribosome rescue)
  LTN1    (RQC ubiquitin ligase)
  NEMF    (CAT-tail addition)
  EIF2S1  (ISR master regulator)
  HSPA5   (BiP/GRP78, ER chaperone)

Input:  data/gtex_v8_tpm.parquet + data/gtex_v8_sample_attrs.parquet
Output: results/{gene}_genome_wide_correlations.csv
        results/{gene}_top5pct.csv
"""

import os
import numpy as np
import pandas as pd
from scipy import stats

# ── Configuration (relative paths for portability) ─────────────────────────
DATA_DIR = "data"
GTEX_TPM = os.path.join(DATA_DIR, "gtex_v8_tpm.parquet")
GTEX_ATTR = os.path.join(DATA_DIR, "gtex_v8_sample_attrs.parquet")
OUTPUT_DIR = "results"

TARGET_GENES = ['SIGMAR1', 'TMEM97', 'PELO', 'LTN1', 'NEMF', 'EIF2S1', 'HSPA5']
BRAIN_REGION = "Brain - Frontal Cortex (BA9)"
MIN_TPM_MEDIAN = 1.0     # Manuscript: "Genes with median TPM < 1.0 were excluded"
TOP_PERCENT = 5

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Load Data ──────────────────────────────────────────────────────────────
print("=" * 60)
print("Script 01: Genome-Wide Co-expression Analysis")
print("=" * 60)
print(f"\nLoading GTEx v8 data...")

tpm = pd.read_parquet(GTEX_TPM) if GTEX_TPM.endswith('.parquet') else pd.read_csv(GTEX_TPM, sep='\t', index_col=0)
attr = pd.read_parquet(GTEX_ATTR) if GTEX_ATTR.endswith('.parquet') else pd.read_csv(GTEX_ATTR, sep='\t')

smtsd_col = [c for c in attr.columns if 'SMTSD' in c.upper() or 'tissue' in c.lower()][0]
sampid_col = [c for c in attr.columns if 'SAMPID' in c.upper() or 'sample' in c.lower()][0]

# ── Filter to BA9 Samples ─────────────────────────────────────────────────
brain_samples = attr[attr[smtsd_col] == BRAIN_REGION][sampid_col].tolist()
brain_cols = [c for c in tpm.columns if c in brain_samples]
brain_tpm = tpm[brain_cols]

# ── Gene Filtering: Median TPM >= 1.0 ─────────────────────────────────────
gene_medians = brain_tpm.median(axis=1)
expressed = gene_medians[gene_medians >= MIN_TPM_MEDIAN].index
brain_log = np.log2(brain_tpm.loc[expressed] + 1)

print(f"  Region: {BRAIN_REGION}")
print(f"  Samples: {len(brain_cols)}")
print(f"  Genes passing median TPM >= {MIN_TPM_MEDIAN}: {len(expressed)}")

# ── Compute Genome-Wide Correlations ──────────────────────────────────────
for gene in TARGET_GENES:
    if gene not in brain_log.index:
        print(f"\n  WARNING: {gene} not found in expression matrix — skipping")
        continue

    print(f"\n  Computing correlations for {gene}...")
    target_vals = brain_log.loc[gene].values
    other_genes = [g for g in brain_log.index if g != gene]

    results = []
    for other in other_genes:
        r, pval = stats.pearsonr(target_vals, brain_log.loc[other].values)
        results.append({'gene': other, 'pearson_r': r, 'p_value': pval})

    df = pd.DataFrame(results).sort_values('pearson_r', ascending=False).reset_index(drop=True)
    df['rank'] = range(1, len(df) + 1)

    full_path = os.path.join(OUTPUT_DIR, f"{gene}_genome_wide_correlations.csv")
    df.to_csv(full_path, index=False)

    threshold = np.percentile(df['pearson_r'], 100 - TOP_PERCENT)
    top5 = df[df['pearson_r'] >= threshold].copy()
    top5_path = os.path.join(OUTPUT_DIR, f"{gene}_top5pct.csv")
    top5.to_csv(top5_path, index=False)

    print(f"    Full ranking: {len(df)} genes → {full_path}")
    print(f"    Top {TOP_PERCENT}%: {len(top5)} genes (r >= {threshold:.4f})")
    print(f"    Top partner: {df.iloc[0]['gene']} (r = {df.iloc[0]['pearson_r']:.4f})")

print("\n" + "=" * 60)
print("✓ Script 01 complete")
print("=" * 60)
