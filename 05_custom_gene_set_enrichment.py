"""
05_custom_gene_set_enrichment.py
================================
Tests 5 curated gene sets against SIGMAR1 and TMEM97 networks using
one-sided Fisher's exact test with the expressed gene universe as background.

Gene sets (Manuscript Section 2.4):
  1. MAM-mitochondrial (14 genes)
  2. Sigma receptors (4 genes)
  3. ER stress/UPR (14 genes)
  4. Methylation pathway (18 genes)
  5. Vascular markers [negative control] (18 genes)

Multi-region testing: each set tested per region for replication.

Output:
  SIGMAR1_custom_enrichment.csv
  TMEM97_custom_enrichment.csv
  custom_enrichment_by_region.csv
"""

import os
import numpy as np
import pandas as pd
from scipy import stats

RESULTS_DIR = "results"
DATA_DIR = "data"
GTEX_TPM = os.path.join(DATA_DIR, "gtex_v8_tpm.parquet")
GTEX_ATTR = os.path.join(DATA_DIR, "gtex_v8_sample_attrs.parquet")

GENOME_SIZE = 16212
TOP_PERCENT = 5

# ── Curated Gene Sets (Manuscript Section 2.4) ────────────────────────────
CUSTOM_GENE_SETS = {
    'mam_mitochondrial': [
        'VDAC1', 'VDAC2', 'VDAC3', 'MFN1', 'MFN2', 'RHOT1', 'RHOT2',
        'ITPR1', 'ITPR2', 'ITPR3', 'VAPB', 'RMDN3', 'PACS2', 'FATE1',
    ],
    'sigma_receptors': [
        'SIGMAR1', 'TMEM97', 'PGRMC1', 'NPC1',
    ],
    'er_stress_upr': [
        'HSPA5', 'HSP90B1', 'DDIT3', 'ATF4', 'ATF6', 'ERN1', 'EIF2AK3',
        'XBP1', 'DNAJB9', 'HERPUD1', 'EDEM1', 'CALR', 'CANX', 'P4HB',
    ],
    'methylation_pathway': [
        'MAT1A', 'MAT2A', 'MAT2B', 'AHCY', 'AHCYL1', 'AHCYL2', 'MTR',
        'MTRR', 'MTHFR', 'BHMT', 'BHMT2', 'CBS', 'CTH', 'GNMT', 'PEMT',
        'NNMT', 'INMT', 'DNMT1',
    ],
    'vascular_markers': [
        'PECAM1', 'CDH5', 'VWF', 'FLT1', 'KDR', 'ENG', 'CLDN5', 'ESAM',
        'ERG', 'TIE1', 'TEK', 'ANGPT1', 'ANGPT2', 'NOS3', 'MCAM',
        'PODXL', 'EMCN', 'ROBO4',
    ],
}

# ── Assert guards: sizes must match manuscript ─────────────────────────────
assert len(CUSTOM_GENE_SETS['mam_mitochondrial']) == 14, "MAM set must be 14 genes"
assert len(CUSTOM_GENE_SETS['sigma_receptors']) == 4, "Sigma set must be 4 genes"
assert len(CUSTOM_GENE_SETS['er_stress_upr']) == 14, "UPR set must be 14 genes"
assert len(CUSTOM_GENE_SETS['methylation_pathway']) == 18, "Methylation set must be 18 genes"
assert len(CUSTOM_GENE_SETS['vascular_markers']) == 18, "Vascular set must be 18 genes"

BRAIN_REGIONS = [
    "Brain - Frontal Cortex (BA9)",
    "Brain - Putamen (basal ganglia)",
    "Brain - Hippocampus",
    "Brain - Nucleus accumbens (basal ganglia)",
    "Brain - Anterior cingulate cortex (BA24)",
]
REGION_SHORT = ['BA9', 'Putamen', 'Hippocampus', 'Nuc_Acc', 'BA24']


def fisher_enrichment(network_genes, custom_set, universe_size):
    """One-sided Fisher's exact test for enrichment."""
    network = set(network_genes)
    custom = set(custom_set)
    a = len(network & custom)
    b = len(network - custom)
    c = len(custom - network)
    d = universe_size - len(network | custom)
    odds, p = stats.fisher_exact([[a, b], [c, d]], alternative='greater')

    expected = len(network) * len(custom) / universe_size
    fold = a / expected if expected > 0 else 0

    return {
        'overlap': a,
        'network_size': len(network),
        'set_size': len(custom),
        'fold_enrichment': round(fold, 2),
        'fisher_OR': round(odds, 2),
        'fisher_p': p,
        'overlapping_genes': ','.join(sorted(network & custom)),
    }


print("=" * 60)
print("Script 05: Custom Gene Set Enrichment (5 sets)")
print("=" * 60)

# ── Primary analysis (BA9) for SIGMAR1 and TMEM97 ─────────────────────────
for gene in ['SIGMAR1', 'TMEM97']:
    top5_path = os.path.join(RESULTS_DIR, f"{gene}_top5pct.csv")
    if not os.path.exists(top5_path):
        print(f"\n  SKIP: {top5_path} not found")
        continue

    network_genes = pd.read_csv(top5_path)['gene'].tolist()
    print(f"\n  {gene} (BA9, {len(network_genes)} genes):")

    records = []
    for set_name, gene_list in CUSTOM_GENE_SETS.items():
        result = fisher_enrichment(network_genes, gene_list, GENOME_SIZE)
        result['gene_set'] = set_name
        result['target'] = gene
        records.append(result)
        sig = '*' if result['fisher_p'] < 0.05 else ' '
        print(f"    {set_name:22s}: {result['fold_enrichment']:5.1f}x  "
              f"p={result['fisher_p']:.3f} {sig}  [{result['overlapping_genes'][:50]}]")

    pd.DataFrame(records).to_csv(
        os.path.join(RESULTS_DIR, f"{gene}_custom_enrichment.csv"), index=False)

# ── Multi-region replication (SIGMAR1 only, matching manuscript) ───────────
print(f"\n{'─' * 50}")
print(f"Multi-region custom enrichment: SIGMAR1")
print(f"{'─' * 50}")

tpm = pd.read_parquet(GTEX_TPM) if GTEX_TPM.endswith('.parquet') else pd.read_csv(GTEX_TPM, sep='\t', index_col=0)
attr = pd.read_parquet(GTEX_ATTR) if GTEX_ATTR.endswith('.parquet') else pd.read_csv(GTEX_ATTR, sep='\t')
smtsd_col = [c for c in attr.columns if 'SMTSD' in c.upper() or 'tissue' in c.lower()][0]
sampid_col = [c for c in attr.columns if 'SAMPID' in c.upper() or 'sample' in c.lower()][0]

region_results = []
for i, region in enumerate(BRAIN_REGIONS):
    samples = attr[attr[smtsd_col] == region][sampid_col].tolist()
    cols = [c for c in tpm.columns if c in samples]
    region_tpm = tpm[cols]
    gene_medians = region_tpm.median(axis=1)
    expressed = gene_medians[gene_medians >= 1.0].index
    region_expr = np.log2(region_tpm.loc[expressed] + 1)
    region_universe = len(expressed)

    if 'SIGMAR1' not in region_expr.index:
        continue

    target = region_expr.loc['SIGMAR1'].values
    corrs = {}
    for g in region_expr.index:
        if g == 'SIGMAR1':
            continue
        r, _ = stats.pearsonr(target, region_expr.loc[g].values)
        corrs[g] = r

    corr_series = pd.Series(corrs).sort_values(ascending=False)
    threshold = np.percentile(corr_series.values, 100 - TOP_PERCENT)
    region_net = set(corr_series[corr_series >= threshold].index)

    for set_name, gene_list in CUSTOM_GENE_SETS.items():
        result = fisher_enrichment(region_net, gene_list, region_universe)
        result['region'] = REGION_SHORT[i]
        result['gene_set'] = set_name
        region_results.append(result)

    print(f"  {REGION_SHORT[i]}: {len(region_net)} genes, universe={region_universe}")

pd.DataFrame(region_results).to_csv(
    os.path.join(RESULTS_DIR, "custom_enrichment_by_region.csv"), index=False)

print("\n" + "=" * 60)
print("✓ Script 05 complete")
print("=" * 60)
