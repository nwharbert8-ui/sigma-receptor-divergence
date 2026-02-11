"""
02_network_comparison.py
========================
Computes all 21 pairwise network comparisons across the 7-gene panel and
exports unique/shared gene lists for SIGMAR1 vs TMEM97.

Methodology (Manuscript Section 2.3):
  - Jaccard similarity, Fisher's exact (one-sided), Spearman rank ρ
  - Gene universe: 16,212 expressed genes
  - Exports: SIGMAR1-unique, TMEM97-unique, shared gene lists

Output:
  network_comparisons_all_pairs.csv  — 21-row table (includes Table 1 data)
  SIGMAR1_unique_top5pct.csv         — 665 genes unique to SIGMAR1
  TMEM97_unique_top5pct.csv          — 666 genes unique to TMEM97
  Sigma_shared_top5pct.csv           — 146 shared genes
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from itertools import combinations

RESULTS_DIR = "results"
TARGET_GENES = ['SIGMAR1', 'TMEM97', 'PELO', 'LTN1', 'NEMF', 'EIF2S1', 'HSPA5']
GENOME_SIZE = 16212

print("=" * 60)
print("Script 02: Pairwise Network Comparison (21 pairs)")
print("=" * 60)

networks = {}
full_rankings = {}

for gene in TARGET_GENES:
    top5_path = os.path.join(RESULTS_DIR, f"{gene}_top5pct.csv")
    full_path = os.path.join(RESULTS_DIR, f"{gene}_genome_wide_correlations.csv")
    if os.path.exists(top5_path):
        df = pd.read_csv(top5_path)
        networks[gene] = set(df['gene'].tolist())
        print(f"  Loaded {gene}: {len(networks[gene])} genes")
    if os.path.exists(full_path):
        df = pd.read_csv(full_path)
        full_rankings[gene] = df.set_index('gene')['pearson_r']

# ── All 21 Pairwise Comparisons ───────────────────────────────────────────
results = []
for g1, g2 in combinations(TARGET_GENES, 2):
    if g1 not in networks or g2 not in networks:
        continue
    s1, s2 = networks[g1], networks[g2]
    inter = s1 & s2
    union = s1 | s2
    jaccard = len(inter) / len(union) if union else 0

    a, b, c = len(inter), len(s1 - s2), len(s2 - s1)
    d = GENOME_SIZE - len(union)
    odds, fisher_p = stats.fisher_exact([[a, b], [c, d]], alternative='greater')

    rho = np.nan
    if g1 in full_rankings and g2 in full_rankings:
        common = full_rankings[g1].index.intersection(full_rankings[g2].index)
        if len(common) > 100:
            rho, _ = stats.spearmanr(full_rankings[g1].loc[common], full_rankings[g2].loc[common])

    results.append({
        'gene_a': g1, 'gene_b': g2,
        'set_a_size': len(s1), 'set_b_size': len(s2),
        'shared_genes': a, 'jaccard': round(jaccard, 4),
        'fisher_OR': round(odds, 2), 'fisher_p': fisher_p,
        'rank_spearman': round(rho, 4) if not np.isnan(rho) else np.nan,
    })
    print(f"  {g1:8s}–{g2:8s}: J={jaccard:.3f}  shared={a:4d}  ρ={rho:.3f}")

comp_df = pd.DataFrame(results).sort_values('jaccard', ascending=False)
comp_df.to_csv(os.path.join(RESULTS_DIR, "network_comparisons_all_pairs.csv"), index=False)
print(f"\n  Saved: network_comparisons_all_pairs.csv ({len(comp_df)} pairs)")

# ── Export SIGMAR1–TMEM97 Gene Lists ───────────────────────────────────────
if 'SIGMAR1' in networks and 'TMEM97' in networks:
    s1, s2 = networks['SIGMAR1'], networks['TMEM97']
    for name, genes in [('SIGMAR1_unique_top5pct', s1 - s2),
                         ('TMEM97_unique_top5pct', s2 - s1),
                         ('Sigma_shared_top5pct', s1 & s2)]:
        pd.DataFrame({'gene': sorted(genes)}).to_csv(
            os.path.join(RESULTS_DIR, f"{name}.csv"), index=False)
        print(f"  {name}: {len(genes)} genes")

print("\n" + "=" * 60)
print("✓ Script 02 complete")
print("=" * 60)
