# Data Access

## GTEx v8

**Portal:** https://gtexportal.org/ | **Accession:** phs000424.v8.p2

### Required Files
1. `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz` (~1.5 GB)
2. `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt` (~15 MB)
3. `GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt` (~200 KB)

### Gene Filtering
Genes with **median TPM < 1.0** excluded → **16,212 genes** in BA9.
Expression transformed as log2(TPM + 1) before correlation analysis.

### Brain Regions
BA9 (n=209), Putamen (n=205), Hippocampus (n=197), Nuc. Acc. (n=246), BA24 (n=176)
