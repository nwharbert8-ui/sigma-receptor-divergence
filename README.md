# Sigma-1 and Sigma-2 Receptor Co-Expression Divergence in Human Brain

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.XXXXXXX-blue)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Google Colab](https://img.shields.io/badge/Run%20in-Google%20Colab-F9AB00?logo=googlecolab)](https://colab.research.google.com/)
[![GTEx v8](https://img.shields.io/badge/Data-GTEx%20v8-green)](https://gtexportal.org/)

**Reproducible analysis pipeline** for the manuscript:

> **Sigma-1 and Sigma-2 Receptors Exhibit Divergent Genome-Wide Co-Expression Networks in Human Brain Despite Shared Subcellular Localization**
>
> Drake H. Harbert  
> *Journal of Neurochemistry* (submitted 2026)

---

## Overview

SIGMAR1 (sigma-1) and TMEM97 (sigma-2) are pharmacologically related endoplasmic reticulum–resident receptors that share subcellular localization at mitochondria-associated membranes (MAMs). Both are implicated in neurodegeneration, neuroprotection, and pain signaling, and they bind overlapping ligand classes. These similarities have led to the longstanding assumption that the two receptors participate in coordinated biological programs.

This study challenges that assumption. Using genome-wide co-expression analysis of GTEx v8 RNA-seq data across five human brain regions, we demonstrate that SIGMAR1 and TMEM97 maintain **fundamentally divergent transcriptional networks** (Jaccard index = 0.112; 88.8% non-overlap) despite their shared subcellular localization. This divergence is robust to cell-type deconvolution, covariate adjustment, and replicates across all five brain regions examined.

### Key Findings

- **Divergent networks**: SIGMAR1 and TMEM97 top 5% co-expression networks share only ~11% of genes (Jaccard = 0.112), with 649 genes unique to each receptor
- **Cross-region robustness**: Network architectures replicate across BA9, putamen, hippocampus, nucleus accumbens, and BA24 (Spearman ρ = 0.857–0.966)
- **SIGMAR1 enrichment for ribosome quality control (RQC)**: VCP, NPLOC4, and TCF25 — extraction-arm components of the RQC pathway — are exclusive to the SIGMAR1 network, while the detection arm (LTN1) shows near-zero network overlap (J = 0.001)
- **Sensitivity-validated**: Divergence persists after cell-type composition adjustment (rank ρ > 0.76) and age/sex covariate regression (rank ρ > 0.99)

---

## Repository Structure

```
sigma-divergence-brain/
│
├── MS3_Sigma_Divergence_Pipeline.py   # Complete analysis pipeline (23 cells)
├── README.md                          # This file
├── LICENSE                            # MIT License
├── requirements.txt                   # Python dependencies
├── CITATION.cff                       # Citation metadata
├── .gitignore                         # Git ignore rules
│
├── docs/
│   └── METHODS.md                     # Detailed statistical methods
│
└── outputs/                           # Generated outputs (git-ignored)
    └── .gitkeep
```

---

## Quick Start

### Option 1: Google Colab (Recommended)

1. Open [Google Colab](https://colab.research.google.com/)
2. Upload `MS3_Sigma_Divergence_Pipeline.py` or copy cells sequentially into a new notebook
3. Run all cells in order — the pipeline will:
   - Mount your Google Drive
   - Download GTEx v8 data (~2.6 GB, cached after first run)
   - Execute all analyses
   - Save results, figures, and supplementary tables to `MyDrive/MS3_JNC_Submission/`

**Runtime**: ~20–30 minutes (12 GB RAM, standard Colab)

### Option 2: Local Execution

```bash
git clone https://github.com/nwharbert8-ui/sigma-divergence-brain.git
cd sigma-divergence-brain
pip install -r requirements.txt
```

> **Note**: The pipeline is designed for Google Colab with 12 GB RAM. Local execution requires modification of the Drive mount cell and sufficient memory for the GTEx TPM matrix.

---

## Pipeline Architecture

The pipeline consists of 23 sequential cells organized into four phases:

### Phase 1: Data Acquisition & Preprocessing (Cells 1–6)

| Cell | Description |
|------|-------------|
| 1 | Environment setup, Google Drive mount, dependency installation |
| 2 | Locked parameter configuration (gene sets, brain regions, thresholds) |
| 3 | GTEx v8 TPM, sample attributes, and subject phenotype download |
| 4 | Sample metadata loading, brain region mapping, covariate extraction |
| 5 | Memory-efficient TPM loading (brain-only columns, 12 GB safe) |
| 6 | Region-specific matrix extraction, gene filtering (median TPM ≥ 1.0), log₂ transform |

### Phase 2: Core Analysis (Cells 7–14)

| Cell | Description |
|------|-------------|
| 7 | Genome-wide Pearson correlations for all 7 target genes in BA9 |
| 8 | Top 5% network extraction, all 21 pairwise comparisons (Jaccard, Fisher's, Spearman) |
| 9 | Top co-expression partner validation |
| 10 | Custom gene set enrichment (6 sets × Fisher's exact test) |
| 11 | gProfiler functional enrichment (GO:BP/MF/CC, KEGG, Reactome) |
| 12 | Cell-type deconvolution sensitivity (53 markers, 6 cell types) |
| 13 | Covariate adjustment (age + sex partial correlations) |
| 14 | Multi-region replication across 5 brain regions |

### Phase 3: Outputs (Cells 15–18)

| Cell | Description |
|------|-------------|
| 15 | Gene list export (top 5%, shared, unique, genome-wide rankings) |
| 16 | Figure 1: Venn diagram, rank scatter, top partner comparison |
| 17 | Figure 2: GO/Reactome/KEGG enrichment bar charts |
| 18 | Figure 3: Cross-region replication heatmaps |

### Phase 4: Supplementary Tables (Cells 19–23)

| Cell | Description |
|------|-------------|
| 19 | **Table S1**: Genome-wide rankings (16,224 genes × both targets) — `.xlsx` |
| 20 | **Table S2**: gProfiler enrichment (5 gene set queries, multi-sheet) — `.xlsx` |
| 21 | **Table S3**: Custom gene set Fisher's exact tests — `.docx` + `.xlsx` |
| 22 | **Table S4**: Shared gene characterization with dual-target ranks — `.docx` + `.xlsx` |
| 23 | Final summary with all manuscript numbers |

---

## Data Sources

| Resource | Version | Access |
|----------|---------|--------|
| GTEx RNA-seq TPM | v8 (2017-06-05) | [GTEx Portal](https://gtexportal.org/) |
| GTEx Sample Attributes | v8 | Downloaded automatically by pipeline |
| GTEx Subject Phenotypes | v8 | Downloaded automatically by pipeline |
| gProfiler | Current release | [g:Profiler API](https://biit.cs.ut.ee/gprofiler/) |

All data are downloaded automatically on first execution. The GTEx TPM file (~2.6 GB compressed) is cached locally for subsequent runs.

### Brain Regions Analyzed

| Abbreviation | GTEx SMTSD Label | Typical n |
|-------------|-------------------|-----------|
| BA9 | Brain - Frontal Cortex (BA9) | ~209 |
| Putamen | Brain - Putamen (basal ganglia) | ~170 |
| Hippocampus | Brain - Hippocampus | ~170 |
| NAcc | Brain - Nucleus accumbens (basal ganglia) | ~202 |
| BA24 | Brain - Anterior cingulate cortex (BA24) | ~157 |

---

## Custom Gene Sets

Six curated gene sets are tested for enrichment in both receptor networks:

| Gene Set | n | Biological Rationale |
|----------|---|---------------------|
| MAM-mitochondrial | 14 | Mitochondria-associated membrane markers (shared localization) |
| Sigma receptor network | 4 | Established sigma receptor interactors (PGRMC1, NPC1) |
| ER stress/UPR | 14 | Unfolded protein response (SIGMAR1 chaperone function) |
| Methylation pathway | 18 | One-carbon metabolism and SAM cycle |
| Vascular markers | 18 | Endothelial-specific genes (**negative control**) |
| Ribosome quality control | 12 | RQC pathway including detection and extraction arms |

---

## Output Files

After pipeline completion, the following directory structure is created on Google Drive:

```
MyDrive/MS3_JNC_Submission/
├── MS3_Results/
│   ├── all_21_pairwise_comparisons.csv
│   ├── SIGMAR1_top5pct.csv
│   ├── TMEM97_top5pct.csv
│   ├── SIGMAR1_TMEM97_shared.csv
│   ├── SIGMAR1_unique.csv
│   ├── TMEM97_unique.csv
│   ├── SIGMAR1_genome_wide_rankings.csv
│   ├── TMEM97_genome_wide_rankings.csv
│   ├── [5 additional target ranking files]
│   ├── SIGMAR1_custom_enrichment.csv
│   ├── TMEM97_custom_enrichment.csv
│   └── gProfiler_*.csv (5 enrichment result files)
│
├── MS3_Figures/
│   ├── Figure1_Divergent_Networks.{png,pdf}
│   ├── Figure2_GO_Enrichment.{png,pdf}
│   └── Figure3_MultiRegion_Replication.{png,pdf}
│
└── MS3_Supplementary_Tables/
    ├── Table_S1_Genome_Wide_Rankings.xlsx
    ├── Table_S2_gProfiler_Enrichment.xlsx
    ├── Table_S3_Custom_Gene_Set_Enrichment.{docx,xlsx}
    └── Table_S4_Shared_Genes.{docx,xlsx}
```

---

## Statistical Methods

- **Co-expression**: Pearson correlation on log₂(TPM + 1) transformed expression values
- **Network definition**: Top 5% of genome-wide correlation coefficients (ceil-rounded)
- **Overlap metrics**: Jaccard index, Fisher's exact test (one-sided, greater), Spearman rank correlation of genome-wide profiles
- **Enrichment**: Fisher's exact test for custom gene sets; g:SCS-corrected gProfiler for GO/KEGG/Reactome (custom background = all expressed genes)
- **Cell-type adjustment**: OLS regression of cell-type marker means (6 types, 53 markers from Darmanis et al. 2015 and Lake et al. 2018), partial correlation on residuals
- **Covariate adjustment**: OLS regression of age (midpoint-coded) and sex, partial correlation on residuals
- **Multi-region replication**: Independent genome-wide correlation computation per region, cross-region Spearman ρ on shared gene rankings

See [`docs/METHODS.md`](docs/METHODS.md) for complete statistical specification.

---

## Reproducibility

This pipeline produces deterministic results from the GTEx v8 public release. All parameters are locked in Cell 2 and cannot drift between runs. The only source of minor variation is the gProfiler API, which uses a live database — functional enrichment term names and p-values may shift slightly with database updates, but biological conclusions are stable.

**Random seed**: `np.random.seed(42)` is set for the Figure 1B scatter plot subsampling only; all statistical computations are deterministic.

---

## Requirements

- Python 3.8+
- Google Colab (recommended) or ≥12 GB RAM local environment
- ~3 GB disk space for GTEx v8 data files

See [`requirements.txt`](requirements.txt) for Python package dependencies.

---

## Citation

If you use this pipeline or build upon these findings, please cite:

```bibtex
@article{harbert2026sigma,
  title={Sigma-1 and Sigma-2 Receptors Exhibit Divergent Genome-Wide
         Co-Expression Networks in Human Brain Despite Shared
         Subcellular Localization},
  author={Harbert, Drake H.},
  journal={Journal of Neurochemistry},
  year={2026},
  note={Submitted},
  doi={PENDING}
}
```

---

## License

This project is licensed under the MIT License — see [`LICENSE`](LICENSE) for details.

---

## Contact

**Drake H. Harbert**  
Inner Architecture LLC  
Canton, OH 44721, USA  
ORCID: [0009-0007-7740-3616](https://orcid.org/0009-0007-7740-3616)  
GitHub: [@nwharbert8-ui](https://github.com/nwharbert8-ui)
