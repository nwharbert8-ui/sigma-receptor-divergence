"""
03_go_enrichment.py
===================
GO/KEGG/Reactome enrichment via gProfiler for 5 gene sets.

Gene sets (Manuscript Section 2.4):
  1. SIGMAR1 full top 5% (811 genes)
  2. TMEM97 full top 5% (812 genes)
  3. SIGMAR1-unique (665 genes)
  4. TMEM97-unique (666 genes)
  5. Shared (146 genes)

Output: {set_name}_GO_enrichment.csv per gene set
"""

import os
import pandas as pd
import requests

RESULTS_DIR = "results"
ORGANISM = 'hsapiens'
SOURCES = ['GO:BP', 'GO:MF', 'GO:CC', 'REAC', 'KEGG']
SIGNIFICANCE_THRESHOLD = 0.05

GENE_SETS = {
    'SIGMAR1_full':    os.path.join(RESULTS_DIR, "SIGMAR1_top5pct.csv"),
    'TMEM97_full':     os.path.join(RESULTS_DIR, "TMEM97_top5pct.csv"),
    'SIGMAR1_unique':  os.path.join(RESULTS_DIR, "SIGMAR1_unique_top5pct.csv"),
    'TMEM97_unique':   os.path.join(RESULTS_DIR, "TMEM97_unique_top5pct.csv"),
    'Sigma_shared':    os.path.join(RESULTS_DIR, "Sigma_shared_top5pct.csv"),
}


def run_gprofiler(gene_list, organism='hsapiens', sources=None):
    """Query gProfiler g:GOSt API for functional enrichment."""
    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
    payload = {
        'organism': organism,
        'query': gene_list,
        'sources': sources or SOURCES,
        'user_threshold': SIGNIFICANCE_THRESHOLD,
        'significance_threshold_method': 'g_SCS',
        'no_evidences': False,
    }
    try:
        response = requests.post(url, json=payload, timeout=120)
    except requests.exceptions.RequestException as e:
        print(f"    API error: {e}")
        return pd.DataFrame()

    if response.status_code != 200:
        return pd.DataFrame()

    data = response.json()
    if 'result' not in data or len(data['result']) == 0:
        return pd.DataFrame()

    records = []
    for term in data['result']:
        records.append({
            'source': term.get('source', ''),
            'term_id': term.get('native', ''),
            'name': term.get('name', ''),
            'p_value': term.get('p_value', 1.0),
            'term_size': term.get('term_size', 0),
            'query_size': term.get('query_size', 0),
            'intersection_size': term.get('intersection_size', 0),
            'precision': term.get('precision', 0),
            'recall': term.get('recall', 0),
        })
    return pd.DataFrame(records).sort_values('p_value')


print("=" * 60)
print("Script 03: GO / Pathway Enrichment (5 gene sets)")
print("=" * 60)

for name, filepath in GENE_SETS.items():
    if not os.path.exists(filepath):
        print(f"\n  SKIP: {filepath} not found")
        continue

    genes = pd.read_csv(filepath).iloc[:, 0].tolist()
    print(f"\n  {name} ({len(genes)} genes)...")
    enrichment = run_gprofiler(genes)

    if len(enrichment) > 0:
        outpath = os.path.join(RESULTS_DIR, f"{name}_GO_enrichment.csv")
        enrichment.to_csv(outpath, index=False)
        print(f"    Saved: {outpath} ({len(enrichment)} terms)")
        for _, row in enrichment[enrichment['source'] == 'GO:BP'].head(3).iterrows():
            print(f"      {row['name'][:55]}: p = {row['p_value']:.2e}")
    else:
        print(f"    No significant terms found")

print("\n" + "=" * 60)
print("✓ Script 03 complete")
print("=" * 60)
