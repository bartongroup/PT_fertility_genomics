# PT_fertility_genomics

# Multi-omics identification of high-confidence testis-specific and sperm-retained fertility genes

## Overview

This repository contains a reproducible pipeline to prioritise **high-confidence male fertility candidate genes** by integrating:

- GTEx v11 tissue RNA-seq (testis specificity)
- Human sperm RNA-seq (GSE40181 TPM)
- ClinVar gene/variant evidence
- Human Phenotype Ontology (HPO) male infertility phenotypes
- Proteomics evidence (internal data) (optional)
- Functional annotation (MyGene.info; optional, requires internet)

The goal is to identify genes that are **testis-enriched**, **retained/expressed in sperm**, and **supported by clinical and/or phenotype evidence**.

All tabular outputs are **tab-separated (TSV)**.

---

## Inputs

You will need the following inputs available locally (downloaded by `download_data.sh` or manually where noted):

### GTEx v11
- Gene TPM GCT (RNASeQC) file (e.g. `GTEx_Analysis_*_gene_tpm.gct`)
- SampleAttributesDS metadata file

### HGNC
- `hgnc_complete_set.txt` (for mapping Ensembl and Entrez IDs to approved HGNC symbols)

### Sperm RNA-seq (GSE40181)
- TPM matrix file (gzipped TSV)

### ClinVar
- `tab_delimited/variant_summary.txt.gz`
- `tab_delimited/gene_condition_source_id` (uncompressed)

### HPO (for male infertility phenotype extraction)
- HPO ontology files and `genes_to_phenotype.txt` (exact filenames depend on your HPO download)

### Proteomics (optional)
- Proteomics Excel, containing a gene symbol column and sample presence columns

---

## Quick start

A typical run consists of:

1. Download inputs
2. Build ClinVar male infertility tiers (HPO + ClinVar)
3. Identify GTEx testis-specific genes and intersect with sperm expression
4. Add annotations (ClinVar / HPO / proteomics / MyGene)
5. Build final summary tables and plots



## Requirements

This pipeline combines Bash and Python scripts 

### System requirements


- Python ≥ 3.9 (tested with Python 3.9)

---

### Python dependencies

The following Python packages are required:

- pandas
- numpy
- requests
- openpyxl (for proteomics Excel files)
- tqdm (optional but useful for progress bars if added later)



```bash
pip install pandas numpy requests openpyxl

```

or conda:

```bash
conda create -n fertility_pipeline python=3.9 pandas numpy requests openpyxl
conda activate fertility_pipeline
```

---

## Pipeline 1: Testis specificity + sperm integration (GTEx + GSE40181)

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40181
https://www.gtexportal.org/home/downloads/adult-gtex/overview


This is driven by `find_testis_expression.sh` and calls the scripts below in order.

### Step 1. Compute testis specificity from GTEx

- **Script:** `scripts/gtex_testis_specificity_logged.py`
- **What it does:**
  - Loads GTEx gene TPM from a GCT file
  - Maps samples to tissues using SampleAttributesDS
  - Computes per-gene tissue median TPM values
  - Calculates tissue-specificity metrics per gene, including:
    - tau (0–1; higher = more tissue-specific)
    - log2 enrichment of testis median TPM vs maximum non-testis median TPM
    - fraction of testis samples where TPM passes a minimum threshold
    - top tissues by median TPM
  
Tau reference:  https://pubmed.ncbi.nlm.nih.gov/15388519/   

- **Key outputs:**
  - Ranked per-gene table (TSV)
  - Gene × tissue median TPM matrix (TSV)

### Step 2. Add HGNC-approved symbols (Ensembl → HGNC)

- **Script:** `scripts/add_hgnc_symbols.py`
- **What it does:**
  - Reads the ranked GTEx output (indexed by Ensembl gene IDs)
  - Strips Ensembl version suffixes (e.g. `.1`)
  - Joins against HGNC complete set to add a `hgnc_symbol` column
- **Key output:**
  - GTEx ranked table with HGNC symbols (TSV)

### external proteomics data:

wget https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/11/PXD014618/DBsearchpsm.csv

### Step 3. Annotate GTEx genes with ClinVar gene/condition evidence

- **Script:** `RNAseq/annotate_testis_genes_with_clinvar.py`
- **What it does:**
  - Adds ClinVar gene-condition associations from `gene_condition_source_id`
  - Optionally parses `variant_summary.txt.gz` and adds per-gene variant counts by clinical significance
  - Normalises gene symbols for robust joins
- **Key output:**
  - GTEx ranked table with ClinVar annotations (TSV)

### Step 4. Map sperm Entrez IDs to HGNC and intersect with GTEx

- **Script:** `RNAseq/map_sperm_entrez_to_hgnc_and_intersect.py`
- **What it does:**
  - Loads sperm TPM matrix (Entrez GeneID-based)
  - Uses HGNC complete set to map Entrez IDs → HGNC symbols
  - Computes sperm evidence columns such as:
    - `sperm_present_any` (TPM > threshold in at least one sample)
    - `sperm_present_frac`
    - `sperm_tpm_mean`, `sperm_tpm_median`
  - Merges sperm evidence onto the GTEx ranked table via HGNC symbol
- **Key outputs:**
  - Mapped sperm TPM table (TSV)
  - GTEx ranked table annotated with sperm evidence (TSV)

### Step 5. Extract a high-confidence testis+sperm gene set

- **Script:** `RNAseq/extract_high_confidence_testis_sperm_genes.py`
- **What it does:**
  - Filters the GTEx+sperm annotated table using configurable thresholds, typically:
    - `min_tau` (default 0.95)
    - `min_testis_tpm` (default 5)
    - `min_sperm_tpm` (default 0.1)
- **Key output:**
  - High-confidence gene set (TSV)

### Optional Step 6. Add functional annotation (MyGene.info)

- **Script:** `RNAseq/add_functional_annotation_mygene.py`
- **What it does:**
  - Queries MyGene.info in batches using HGNC symbols
  - Adds fields such as gene name, summary, gene type, GO terms
  - Supports a JSONL cache to avoid repeated queries
- **Key output:**
  - Functionally annotated GTEx+sperm table (TSV)

---

## Pipeline 2: Male infertility clinical tiers (HPO + ClinVar)

This is driven by `identify_genes_to_infertility.sh` and creates a ClinVar evidence set focused on male infertility phenotypes.

### Step 1. Extract male infertility genes from HPO

- **Script:** `extract_male_infertility_genes_hpo_ontology.py`
- **What it does:**
  - Traverses HPO ontology to collect genes associated with male infertility phenotypes
  - Can require “male” context and include disease terms depending on flags
- **Key outputs:**
  - HPO-derived gene summary tables (TSV)

### Step 2. Create a clean gene list from HPO outputs

- **Script:** `extract_gene_list_from_hpo_outputs.py`
- **What it does:**
  - Extracts unique gene symbols and Entrez IDs
  - Writes:
    - `male_infertility_gene_symbols.tsv`
    - `male_infertility_entrez_ids.tsv`
    - `male_infertility_gene_list.tsv` (combined)
- **Key outputs:**
  - Gene list TSVs for filtering ClinVar

### Step 3. Filter ClinVar variant_summary by the HPO-derived gene list

- **Script:** `filter_clinvar_variant_summary_by_gene_list.py`
- **What it does:**
  - Streams `variant_summary.txt.gz` and keeps only rows whose `GeneSymbol` matches the HPO-derived gene list
  - Writes counts by ClinicalSignificance and ReviewStatus
- **Key outputs:**
  - Filtered ClinVar subset (TSV)
  - Summary counts (TSV)

### Step 4. Restrict ClinVar variants to male infertility phenotype terms

- **Script:** `filter_clinvar_filtered_by_phenotype.py`
- **What it does:**
  - Filters ClinVar rows by keyword matching against `PhenotypeList`
  - Uses include and exclude keyword sets to avoid female-only contexts
  - Produces both full and minimal meeting-ready outputs
- **Key outputs:**
  - `male_infertility_clinvar_variants_by_phenotype.tsv`
  - `male_infertility_clinvar_variants_by_phenotype_minimal.tsv`
  - Summary counts TSV

### Step 5. Collapse to a single best record per variant

- **Script:** `collapse_clinvar_variants_to_best.py`
- **What it does:**
  - Groups by `AlleleID` (default) and selects one “best” record per variant using a rank:
    - Prefer GRCh38 over GRCh37
    - Prefer stronger ReviewStatus
    - Prefer more clinically relevant ClinicalSignificance
- **Key outputs:**
  - Best-per-variant TSV
  - Best-per-variant pathogenic-only TSV
  - Counts summary TSV

### Step 6. Define high-confidence ClinVar sets

- **Script:** `filter_clinvar_best_high_confidence.py`
- **What it does:**
  - Filters to high-confidence review categories, such as:
    - practice guideline
    - reviewed by expert panel
    - criteria provided, multiple submitters, no conflicts
  - Creates:
    - all high-confidence rows
    - high-confidence + pathogenic/likely pathogenic subset
  - Writes summary counts
- **Key outputs:**
  - `male_infertility_clinvar_variants_high_confidence.tsv`
  - `male_infertility_clinvar_variants_high_confidence_pathogenic.tsv`
  - Counts TSV

### Step 7. Produce concise ClinVar report tables

- **Script:** `make_clinvar_high_confidence_report.py`
- **What it does:**
  - Writes a boss-ready report table from the high-confidence pathogenic subset
  - Writes per-gene summaries with phenotype strings and review statuses
- **Key outputs:**
  - Concise report TSV
  - Per-gene summary TSV

---

## Pipeline 3: Overlap and tiering between testis+sperm genes and ClinVar infertility genes

- **Script:** `overlap_testis_specific_with_clinvar_tiers.py`
- **What it does:**
  - Intersects:
    - high-confidence testis+sperm gene set
    - ClinVar male infertility genes (best-per-variant tables)
  - Assigns tiers based on ReviewStatus and ClinicalSignificance (e.g. high-confidence pathogenic vs other)
  - Writes:
    - overlap table with per-gene evidence counts
    - testis-only gene list
    - ClinVar-only gene list
    - optional long-form evidence table

---

## Optional annotations

These can be applied to the GTEx+sperm table before extracting final gene sets.

### HPO gene-to-phenotype annotation
- **Script:** `add_hpo_annotation.py`
- **What it does:**
  - Joins HPO gene-to-phenotype associations onto the gene table
  - Summarises:
    - total HPO term count
    - reproductive-related term count (regex-based)
    - lists of HPO terms, IDs, and disease IDs

### Proteomics annotation
- **Script:** `add_proteomics_annotation.py`
- **What it does:**
  - Reads a proteomics Excel export and summarises evidence per gene
  - Supports filtering by:
    - FDR confidence (default “High”)
    - species (default “Homo sapiens”)
    - minimum unique peptides
  - Adds presence across “Found in Sample:” columns and accession summaries

### Functional annotation (MyGene.info)
- **Script:** `add_functional_annotation_mygene.py`
- **What it does:**
  - Adds names, summaries, gene types, GO terms via MyGene.info
  - Uses batching, retry logic, and an optional JSONL cache

---

## Outputs and results

At the end of the pipeline you should have:

### Core tables
- GTEx testis specificity ranked table (per-gene tau and enrichment metrics)
- GTEx ranked table with HGNC symbols
- GTEx ranked table annotated with:
  - ClinVar gene-condition summaries and optional variant counts
  - sperm TPM evidence columns (presence, mean/median)
  - optional HPO term summaries
  - optional proteomics evidence
  - optional functional annotation

### High-confidence gene sets
- `high_confidence_testis_sperm_genes.tsv`
- (optionally) `high_confidence_testis_sperm_genes_function.tsv` (if you re-filter after MyGene annotation)

### ClinVar male infertility tiers
- Best-per-variant and high-confidence subsets for male infertility phenotype-associated variants
- Gene-level summaries of high-confidence pathogenic variants

### Overlap summaries
- Intersection between the testis+sperm gene set and ClinVar male infertility tiers
- Testis-only and ClinVar-only gene lists
- Tiered overlap tables with per-gene evidence counts

---

##  notes

- All outputs are TSV (tab-separated).
- Scripts accept named CLI arguments; thresholds are configurable.
- Several scripts write logs to `logs/` when `--log_path` is provided.
- ClinVar filtering steps are designed to stream large inputs efficiently (no need to fully decompress variant_summary).

---

##  run order

1. `download_data.sh`
2. `identify_genes_to_infertility.sh`
3. `find_testis_expression.sh`
4. `mk_final_results.sh`

If you run things manually, follow the step order described above.

---

## Citation and licensing

See `LICENSE` for repository licensing terms.
