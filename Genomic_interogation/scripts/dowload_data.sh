#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# download.sh
#

###############################################################################


mkdir -p genome_resources
cd genome_resources

# 1) GENCODE GTF (choose a release; example placeholder name)
# Get the "comprehensive gene annotation (CHR)" GTF from the GENCODE release page:
#   https://www.gencodegenes.org/human/release_38.html
# Download the file ending in: comprehensive.chr.gtf.gz (or similar)
# (I’m not hardcoding the exact filename here because it changes by release.) :contentReference[oaicite:6]{index=6}

# 2) hg38 FASTA (UCSC)
wget -c -O hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/latest/hg38.fa.gz
wget -c -O hg38.fa.gz.fai https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/latest/hg38.fa.gz.fai || true

# 3) RepeatMasker table dump
wget -c -O rmsk.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz

# 4) Segmental duplications table dump
wget -c -O genomicSuperDups.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz

# 5) Recombination map (deCODE avg bigWig from UCSC recombRate track)
wget -c -O recombAvg.bw https://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/recombAvg.bw
