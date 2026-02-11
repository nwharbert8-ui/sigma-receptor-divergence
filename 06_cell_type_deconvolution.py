"""
06_cell_type_deconvolution.py
=============================
Estimates cell-type proportions from 53 marker genes across 6 brain cell
types, then computes partial correlations controlling for composition.

Methodology (Manuscript Section 2.5, Results 3.7):
  - 53 established marker genes across 6 cell types
  - Proportions estimated as mean log2(TPM+1) of markers per sample
  - Partial correlation via residualization (OLS regression)
  - Validation: Spearman ρ > 0.95 between raw and deconvolved rankings
  - Top 5% threshold drops from r ≥ 0.831 to r ≥ 0.402

Output:
  SIGMAR1_deconvolved_correlations.csv
  deconvolution_comparison.csv
"""

import os
import numpy as np
import pandas as pd
from scipy import stats

DATA_DIR = "data"
GTEX_TPM = os.path.join(DATA_DIR, "gtex_v8_tpm.parquet")
GTEX_ATTR = os.path.join(DATA_DIR, "gtex_v8_sample_attrs.parquet")
OUTPUT_DIR = "results"
BRAIN_REGION = "Brain - Frontal Cortex (BA9)"
MIN_TPM_MEDIAN = 1.0

# ── Cell-Type Marker Genes (Darmanis et al. 2015; Lake et al. 2018) ──────
CELL_TYPE_MARKERS = {
    'neurons': [
        'RBFOX3', 'SYT1', 'SNAP25', 'SLC17A7', 'GAD1', 'GAD2',
        'STMN2', 'SYN1', 'NRGN', 'THY1',
    ],
    'astrocytes': [
        'GFAP', 'AQP4', 'SLC1A2', 'SLC1A3', 'ALDH1L1', 'GJA1',
        'SOX9', 'NDRG2', 'S100B',
    ],
    'oligodendrocytes': [
        'MBP', 'MOG', 'PLP1', 'MAG', 'CLDN11', 'CNP', 'OLIG1',
        'OLIG2', 'SOX10',
    ],
    'microglia': [
        'CX3CR1', 'P2RY12', 'TMEM119', 'CSF1R', 'AIF1', 'CD68',
        'HEXB', 'TREM2', 'C1QA',
    ],
    'endothelial': [
        'CLDN5', 'VWF', 'FLT1', 'PECAM1', 'CDH5', 'ERG',
        'ESAM', 'ENG',
    ],
    'OPC': [
        'PDGFRA', 'CSPG4', 'VCAN', 'GPR17', 'OLIG2', 'SOX6',
        'PCDH15', 'NEU4',
    ],
}

all_markers = []
for v in CELL_TYPE_MARKERS.values():
    all_markers.extend(v)
unique_markers = list(set(all_markers))
print(f"Total unique markers: {len(unique_markers)}")
assert len(unique_markers) == 53, f"Expected 53 markers, got {len(unique_markers)}"

print("=" * 60)
print("Script 06: Cell-Type Deconvolution Sensitivity")
print("=" * 60)

# ── Load Data ──────────────────────────────────────────────────────────────
print("\nLoading GTEx v8 data...")
tpm = pd.read_parquet(GTEX_TPM) if GTEX_TPM.endswith('.parquet') else pd.read_csv(GTEX_TPM, sep='\t', index_col=0)
attr = pd.read_parquet(GTEX_ATTR) if GTEX_ATTR.endswith('.parquet') else pd.read_csv(GTEX_ATTR, sep='\t')

smtsd_col = [c for c in attr.columns if 'SMTSD' in c.upper() or 'tissue' in c.lower()][0]
sampid_col = [c for c in attr.columns if 'SAMPID' in c.upper() or 'sample' in c.lower()][0]

brain_samples = attr[attr[smtsd_col] == BRAIN_REGION][sampid_col].tolist()
brain_cols = [c for c in tpm.columns if c in brain_samples]
brain_tpm = tpm[brain_cols]

gene_medians = brain_tpm.median(axis=1)
expressed = gene_medians[gene_medians >= MIN_TPM_MEDIAN].index
brain_log = np.log2(brain_tpm.loc[expressed] + 1)

print(f"  Samples: {len(brain_cols)}, Genes: {len(expressed)}")

# ── Estimate Cell-Type Proportions ─────────────────────────────────────────
print("\nEstimating cell-type proportions...")
proportions = pd.DataFrame(index=brain_cols)
for cell_type, markers in CELL_TYPE_MARKERS.items():
    available = [m for m in markers if m in brain_log.index]
    if available:
        proportions[cell_type] = brain_log.loc[available].mean(axis=0)
        print(f"  {cell_type:20s}: {len(available)}/{len(markers)} markers available")
    else:
        proportions[cell_type] = 0.0
        print(f"  {cell_type:20s}: 0/{len(markers)} markers available — skipped")

# ── Compute Partial Correlations (Residualization) ─────────────────────────
print("\nComputing partial correlations for SIGMAR1...")
gene = 'SIGMAR1'
target = brain_log.loc[gene].values
prop_matrix = proportions.values

# Regress cell-type proportions from target
X = np.column_stack([np.ones(len(target)), prop_matrix])
beta_target = np.linalg.lstsq(X, target, rcond=None)[0]
target_resid = target - X @ beta_target

results = []
for other in brain_log.index:
    if other == gene:
        continue
    other_vals = brain_log.loc[other].values
    beta_other = np.linalg.lstsq(X, other_vals, rcond=None)[0]
    other_resid = other_vals - X @ beta_other
    r, p = stats.pearsonr(target_resid, other_resid)
    results.append({'gene': other, 'partial_r': r, 'p_value': p})

df_deconv = pd.DataFrame(results).sort_values('partial_r', ascending=False).reset_index(drop=True)
df_deconv['rank'] = range(1, len(df_deconv) + 1)

outpath = os.path.join(OUTPUT_DIR, f"{gene}_deconvolved_correlations.csv")
df_deconv.to_csv(outpath, index=False)

threshold_deconv = np.percentile(df_deconv['partial_r'], 95)
print(f"\n  Deconvolved top 5% threshold: r >= {threshold_deconv:.4f}")
print(f"  (vs. raw threshold r >= 0.831)")

# ── Compare Raw vs. Deconvolved Rankings ───────────────────────────────────
raw_path = os.path.join(OUTPUT_DIR, f"{gene}_genome_wide_correlations.csv")
if os.path.exists(raw_path):
    df_raw = pd.read_csv(raw_path).set_index('gene')
    df_deconv_idx = df_deconv.set_index('gene')
    common = df_raw.index.intersection(df_deconv_idx.index)
    rho, _ = stats.spearmanr(df_raw.loc[common, 'pearson_r'],
                              df_deconv_idx.loc[common, 'partial_r'])
    print(f"  Rank preservation (Spearman ρ): {rho:.4f}")

    # Check VCP persistence
    raw_top5 = set(df_raw.nlargest(len(df_raw) // 20, 'pearson_r').index)
    deconv_top5 = set(df_deconv.nlargest(len(df_deconv) // 20, 'partial_r')['gene'])
    vcp_in_raw = 'VCP' in raw_top5
    vcp_in_deconv = 'VCP' in deconv_top5
    print(f"  VCP in raw top 5%: {vcp_in_raw}")
    print(f"  VCP in deconvolved top 5%: {vcp_in_deconv}")

    pd.DataFrame([{
        'analysis': 'cell_type_deconvolution',
        'target': gene,
        'raw_threshold': 0.831,
        'adjusted_threshold': round(threshold_deconv, 4),
        'rank_preservation_rho': round(rho, 4),
        'VCP_persists': vcp_in_deconv,
        'n_markers': len(unique_markers),
        'n_cell_types': 6,
    }]).to_csv(os.path.join(OUTPUT_DIR, "deconvolution_comparison.csv"), index=False)

print("\n" + "=" * 60)
print("✓ Script 06 complete")
print("=" * 60)
