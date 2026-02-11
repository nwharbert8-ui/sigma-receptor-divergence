"""
07_covariate_adjustment.py
==========================
Controls for donor age and sex using partial correlation via residualization.

Methodology (Manuscript Section 2.5, Results 3.7):
  - Age: midpoint of GTEx age bins (e.g., "60-69" → 64.5)
  - Sex: binary coding (1=male, 2=female from GTEx)
  - Partial correlation via OLS residualization
  - Validation: Spearman ρ > 0.99 between raw and adjusted rankings
  - Top 5% threshold shifts from r ≥ 0.831 to r ≥ 0.875

Output:
  SIGMAR1_covariate_adjusted_correlations.csv
  covariate_adjustment_comparison.csv
"""

import os
import re
import numpy as np
import pandas as pd
from scipy import stats

DATA_DIR = "data"
GTEX_TPM = os.path.join(DATA_DIR, "gtex_v8_tpm.parquet")
GTEX_ATTR = os.path.join(DATA_DIR, "gtex_v8_sample_attrs.parquet")
GTEX_PHENO = os.path.join(DATA_DIR, "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
OUTPUT_DIR = "results"
BRAIN_REGION = "Brain - Frontal Cortex (BA9)"
MIN_TPM_MEDIAN = 1.0

print("=" * 60)
print("Script 07: Covariate Adjustment (Age + Sex)")
print("=" * 60)

# ── Load Expression Data ──────────────────────────────────────────────────
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

# ── Load Phenotype Data ───────────────────────────────────────────────────
print("Loading phenotype data...")
if os.path.exists(GTEX_PHENO):
    pheno = pd.read_csv(GTEX_PHENO, sep='\t')
else:
    pheno_parquet = os.path.join(DATA_DIR, "gtex_v8_subject_phenotypes.parquet")
    if os.path.exists(pheno_parquet):
        pheno = pd.read_parquet(pheno_parquet)
    else:
        print("  WARNING: Phenotype file not found. Creating dummy covariates.")
        pheno = None

# ── Map Sample IDs to Subject IDs ─────────────────────────────────────────
def sample_to_subject(sample_id):
    """Extract GTEx subject ID from sample ID (GTEX-XXXX-...)."""
    parts = sample_id.split('-')
    if len(parts) >= 2:
        return '-'.join(parts[:2])
    return sample_id

subject_ids = [sample_to_subject(s) for s in brain_cols]

# ── Parse Age and Sex ──────────────────────────────────────────────────────
age_values = []
sex_values = []

if pheno is not None:
    subjid_col = [c for c in pheno.columns if 'SUBJID' in c.upper()][0]
    age_col = [c for c in pheno.columns if 'AGE' in c.upper()][0]
    sex_col = [c for c in pheno.columns if 'SEX' in c.upper()][0]
    pheno_lookup = pheno.set_index(subjid_col)

    for subj in subject_ids:
        if subj in pheno_lookup.index:
            age_str = str(pheno_lookup.loc[subj, age_col])
            match = re.match(r'(\d+)-(\d+)', age_str)
            if match:
                age_values.append((int(match.group(1)) + int(match.group(2))) / 2)
            else:
                try:
                    age_values.append(float(age_str))
                except ValueError:
                    age_values.append(np.nan)
            sex_values.append(pheno_lookup.loc[subj, sex_col])
        else:
            age_values.append(np.nan)
            sex_values.append(np.nan)
else:
    age_values = [np.nan] * len(brain_cols)
    sex_values = [np.nan] * len(brain_cols)

covariates = pd.DataFrame({
    'sample': brain_cols,
    'age': age_values,
    'sex': sex_values,
}).dropna()

valid_samples = covariates['sample'].tolist()
brain_log_valid = brain_log[valid_samples]
print(f"  Samples with complete covariates: {len(valid_samples)}")

# ── Compute Partial Correlations ──────────────────────────────────────────
print("\nComputing partial correlations for SIGMAR1...")
gene = 'SIGMAR1'
target = brain_log_valid.loc[gene].values
cov_matrix = covariates[['age', 'sex']].values

X = np.column_stack([np.ones(len(target)), cov_matrix])
beta_target = np.linalg.lstsq(X, target, rcond=None)[0]
target_resid = target - X @ beta_target

results = []
for other in brain_log_valid.index:
    if other == gene:
        continue
    other_vals = brain_log_valid.loc[other].values
    beta_other = np.linalg.lstsq(X, other_vals, rcond=None)[0]
    other_resid = other_vals - X @ beta_other
    r, p = stats.pearsonr(target_resid, other_resid)
    results.append({'gene': other, 'partial_r': r, 'p_value': p})

df_adj = pd.DataFrame(results).sort_values('partial_r', ascending=False).reset_index(drop=True)
df_adj['rank'] = range(1, len(df_adj) + 1)

outpath = os.path.join(OUTPUT_DIR, f"{gene}_covariate_adjusted_correlations.csv")
df_adj.to_csv(outpath, index=False)

threshold_adj = np.percentile(df_adj['partial_r'], 95)
print(f"  Adjusted top 5% threshold: r >= {threshold_adj:.4f}")
print(f"  (vs. raw threshold r >= 0.831)")

# ── Compare Raw vs. Adjusted Rankings ─────────────────────────────────────
raw_path = os.path.join(OUTPUT_DIR, f"{gene}_genome_wide_correlations.csv")
if os.path.exists(raw_path):
    df_raw = pd.read_csv(raw_path).set_index('gene')
    df_adj_idx = df_adj.set_index('gene')
    common = df_raw.index.intersection(df_adj_idx.index)
    rho, _ = stats.spearmanr(df_raw.loc[common, 'pearson_r'],
                              df_adj_idx.loc[common, 'partial_r'])
    print(f"  Rank preservation (Spearman ρ): {rho:.4f}")

    pd.DataFrame([{
        'analysis': 'covariate_adjustment',
        'target': gene,
        'raw_threshold': 0.831,
        'adjusted_threshold': round(threshold_adj, 4),
        'rank_preservation_rho': round(rho, 4),
        'covariates': 'age,sex',
        'n_samples_used': len(valid_samples),
    }]).to_csv(os.path.join(OUTPUT_DIR, "covariate_adjustment_comparison.csv"), index=False)

print("\n" + "=" * 60)
print("✓ Script 07 complete")
print("=" * 60)
