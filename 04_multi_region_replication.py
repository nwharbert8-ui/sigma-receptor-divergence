"""
04_multi_region_replication.py
==============================
Validates co-expression network stability across 5 GTEx brain regions.

Methodology (Manuscript Section 2.5, Results 3.6):
  - Independent genome-wide correlations per region
  - Cross-region Spearman rank correlation (10 pairs)
  - SIGMAR1 range: ρ = 0.811–0.938
  - Gene filtering: median TPM >= 1.0 (consistent with Script 01)

Output: {gene}_cross_region_matrix.csv — 5×5 Spearman ρ
"""

import os
import numpy as np
import pandas as pd
from scipy import stats

DATA_DIR = "data"
GTEX_TPM = os.path.join(DATA_DIR, "gtex_v8_tpm.parquet")
GTEX_ATTR = os.path.join(DATA_DIR, "gtex_v8_sample_attrs.parquet")
OUTPUT_DIR = "results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

TARGET_GENES = ['SIGMAR1', 'TMEM97']
MIN_TPM_MEDIAN = 1.0  # Matches manuscript: "median TPM < 1.0 were excluded"

BRAIN_REGIONS = [
    "Brain - Frontal Cortex (BA9)",
    "Brain - Putamen (basal ganglia)",
    "Brain - Hippocampus",
    "Brain - Nucleus accumbens (basal ganglia)",
    "Brain - Anterior cingulate cortex (BA24)",
]
REGION_SHORT = ['BA9', 'Putamen', 'Hippocampus', 'Nuc_Acc', 'BA24']

print("=" * 60)
print("Script 04: Multi-Region Replication")
print("=" * 60)

print("\nLoading GTEx v8 data...")
tpm = pd.read_parquet(GTEX_TPM) if GTEX_TPM.endswith('.parquet') else pd.read_csv(GTEX_TPM, sep='\t', index_col=0)
attr = pd.read_parquet(GTEX_ATTR) if GTEX_ATTR.endswith('.parquet') else pd.read_csv(GTEX_ATTR, sep='\t')

smtsd_col = [c for c in attr.columns if 'SMTSD' in c.upper() or 'tissue' in c.lower()][0]
sampid_col = [c for c in attr.columns if 'SAMPID' in c.upper() or 'sample' in c.lower()][0]


def get_region_correlations(gene, region):
    """Compute genome-wide Pearson correlations for a gene in one brain region."""
    samples = attr[attr[smtsd_col] == region][sampid_col].tolist()
    cols = [c for c in tpm.columns if c in samples]
    if len(cols) < 30:
        print(f"    WARNING: {region} has only {len(cols)} samples")
        return None

    region_tpm = tpm[cols]
    gene_medians = region_tpm.median(axis=1)
    expressed = gene_medians[gene_medians >= MIN_TPM_MEDIAN].index
    region_expr = np.log2(region_tpm.loc[expressed] + 1)

    if gene not in region_expr.index:
        print(f"    WARNING: {gene} not expressed in {region}")
        return None

    target = region_expr.loc[gene].values
    correlations = {}
    for other in region_expr.index:
        if other == gene:
            continue
        r, _ = stats.pearsonr(target, region_expr.loc[other].values)
        correlations[other] = r

    return pd.Series(correlations)


for gene in TARGET_GENES:
    print(f"\n{'─' * 50}")
    print(f"Cross-region analysis: {gene}")
    print(f"{'─' * 50}")

    region_profiles = {}
    for i, region in enumerate(BRAIN_REGIONS):
        print(f"  Computing {REGION_SHORT[i]} (n={len(attr[attr[smtsd_col] == region])})...")
        profile = get_region_correlations(gene, region)
        if profile is not None:
            region_profiles[REGION_SHORT[i]] = profile

    regions_available = list(region_profiles.keys())
    n = len(regions_available)
    cross_matrix = np.ones((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            common = region_profiles[regions_available[i]].index.intersection(
                     region_profiles[regions_available[j]].index)
            rho, _ = stats.spearmanr(
                region_profiles[regions_available[i]].loc[common],
                region_profiles[regions_available[j]].loc[common])
            cross_matrix[i, j] = round(rho, 3)
            cross_matrix[j, i] = round(rho, 3)

    cross_df = pd.DataFrame(cross_matrix, index=regions_available, columns=regions_available)
    outpath = os.path.join(OUTPUT_DIR, f"{gene}_cross_region_matrix.csv")
    cross_df.to_csv(outpath)

    triu = cross_matrix[np.triu_indices(n, k=1)]
    print(f"\n  Range: ρ = {triu.min():.3f}–{triu.max():.3f}")
    print(f"  Saved: {outpath}")
    print(cross_df.to_string())

print("\n" + "=" * 60)
print("✓ Script 04 complete")
print("=" * 60)
