#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_gene_context_table.py

Build a per-gene genomic context feature table for a target gene list on GRCh38.

This script maps HGNC-style gene symbols (e.g., MNS1, ODAD2, PMFBP1) to GENCODE
GRCh38 gene models, then computes genomic-context features per gene.

Key outputs
-----------
- Gene body coordinates and gene-level TSS coordinates (sensitivity analysis)
- Local gene density around TSS (±50 kb, ±250 kb, ±1 Mb)
- Nearest neighbour distance between gene TSS positions (within chromosome)
- RepeatMasker repeat fraction in windows around TSS (±10 kb, ±50 kb, ±250 kb)
- Distance to nearest repeat (any and by class: LINE/SINE/LTR/DNA)
- Segmental duplication overlap fraction (gene body) and nearest distance
- Subtelomeric proximity (distance to nearest chromosome end + boolean flags)
- GC content for gene body and TSS flanks (±10 kb, ±50 kb, ±250 kb)
- Recombination mean/max in windows around TSS (±50 kb, ±250 kb, ±1 Mb)

All outputs are TSV.

Notes
-----
- Designed to work cleanly with UCSC hg38 resources (chr-prefixed chromosomes).
- Works with GENCODE primary assembly GTF by normalising chromosomes to chr*.
- Requires: pandas, numpy, pyranges, pyfaidx, pyBigWig

Example
-------
python build_gene_context_table.py \
  --gene_list_tsv sperm_genes.tsv \
  --gencode_gtf gencode.v49.primary_assembly.annotation.gtf.gz \
  --hg38_fasta hg38.fa.gz \
  --rmsk_tsv_gz rmsk.txt.gz \
  --segdup_tsv_gz genomicSuperDups.txt.gz \
  --recomb_bw recombAvg.bw \
  --out_tsv gene_context_features.tsv \
  --mapping_report_tsv gene_mapping_report.tsv \
  --verbose
"""

from __future__ import annotations

import argparse
import gzip
import logging
import math
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

import pyranges as pr
import pyBigWig
from pyfaidx import Fasta


WINDOWS_BP: Tuple[int, ...] = (50_000, 250_000, 1_000_000)
REPEAT_WINDOWS_BP: Tuple[int, ...] = (10_000, 50_000, 250_000)
SUBTELO_THRESHOLDS_BP: Tuple[int, ...] = (1_000_000, 5_000_000, 10_000_000)
REPEAT_CLASSES: Tuple[str, ...] = ("LINE", "SINE", "LTR", "DNA")

STD_CHROMS: Tuple[str, ...] = tuple([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])


@dataclass(frozen=True)
class GeneRecord:
    """
    Representation of a gene with coordinates on a chromosome.

    Attributes
    ----------
    gene_id
        Ensembl gene identifier (may include version suffix).
    gene_name
        Gene symbol from GENCODE (typically HGNC-like).
    chrom
        Chromosome name (expected to be chr-prefixed).
    start
        0-based inclusive start coordinate.
    end
        0-based exclusive end coordinate.
    strand
        '+' or '-'.
    """

    gene_id: str
    gene_name: str
    chrom: str
    start: int
    end: int
    strand: str

    @property
    def length(self) -> int:
        """Return gene length in base pairs."""
        return max(0, int(self.end) - int(self.start))

    @property
    def tss(self) -> int:
        """
        Return 0-based TSS coordinate.

        For '+' strand, TSS is start.
        For '-' strand, TSS is end - 1.
        """
        if self.strand == "+":
            return int(self.start)
        return max(int(self.start), int(self.end) - 1)


def setup_logger(verbose: bool) -> None:
    """
    Configure logging.

    Parameters
    ----------
    verbose
        If True, set log level to DEBUG; else INFO.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def open_maybe_gzip(path: str):
    """
    Open a text file that may be gzipped.

    Parameters
    ----------
    path
        Path to file, possibly ending in .gz

    Returns
    -------
    file handle
        Text mode handle.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt", encoding="utf-8")


def parse_gtf_attributes(attrs: str) -> Dict[str, str]:
    """
    Parse a GTF attributes column into a dictionary.

    Parameters
    ----------
    attrs
        The 9th column of a GTF row.

    Returns
    -------
    dict
        Mapping of attribute keys to values.
    """
    out: Dict[str, str] = {}
    for item in attrs.strip().split(";"):
        item = item.strip()
        if item == "":
            continue
        if " " not in item:
            continue
        key, value = item.split(" ", 1)
        out[key] = value.strip().strip('"')
    return out


def normalise_chrom(chrom: str) -> str:
    """
    Normalise chromosome names to UCSC-style chr-prefixed.

    Parameters
    ----------
    chrom
        Chromosome name from source.

    Returns
    -------
    str
        Chromosome name with chr prefix where appropriate.
    """
    c = str(chrom).strip()
    if c.startswith("chr"):
        return c
    return f"chr{c}"


def filter_standard_chromosomes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter to standard chromosomes chr1-22, chrX, chrY.

    Parameters
    ----------
    df
        DataFrame with a 'Chromosome' column.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame.
    """
    return df[df["Chromosome"].isin(STD_CHROMS)].copy()


def read_gene_list(gene_list_tsv: str) -> pd.DataFrame:
    """
    Read input gene list TSV.

    Parameters
    ----------
    gene_list_tsv
        Path to TSV containing column 'gene_key' (HGNC-like symbols).

    Returns
    -------
    pd.DataFrame
        Input genes, deduplicated and stripped.
    """
    df = pd.read_csv(gene_list_tsv, sep="\t", dtype=str).fillna("")
    if "gene_key" not in df.columns:
        raise ValueError("Input gene list must contain a 'gene_key' column.")
    df["gene_key"] = df["gene_key"].astype(str).str.strip()
    df = df[df["gene_key"] != ""].drop_duplicates(subset=["gene_key"]).copy()
    logging.info("Loaded %s unique gene_key entries", df.shape[0])
    return df


def parse_gtf_genes(gencode_gtf: str) -> pd.DataFrame:
    """
    Parse gene records from a GENCODE GTF.

    Parameters
    ----------
    gencode_gtf
        Path to a GTF (gz ok).

    Returns
    -------
    pd.DataFrame
        Gene table with 0-based half-open intervals and columns:
        Chromosome, Start, End, Strand, gene_id, gene_name.
    """
    records: List[Dict[str, object]] = []
    n_lines = 0
    n_genes = 0

    with open_maybe_gzip(gencode_gtf) as handle:
        for line in handle:
            n_lines += 1
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, _source, feature, start, end, _score, strand, _frame, attrs = parts
            if feature != "gene":
                continue
            if strand not in {"+", "-"}:
                continue

            try:
                start_1 = int(start)
                end_1 = int(end)
            except ValueError:
                continue

            attr_map = parse_gtf_attributes(attrs=attrs)
            gene_id = str(attr_map.get("gene_id", "")).strip()
            gene_name = str(attr_map.get("gene_name", "")).strip()
            if gene_id == "" or gene_name == "":
                continue

            start0 = start_1 - 1
            end0 = end_1

            chrom_norm = normalise_chrom(chrom=chrom)

            records.append(
                {
                    "Chromosome": chrom_norm,
                    "Start": int(start0),
                    "End": int(end0),
                    "Strand": strand,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                }
            )
            n_genes += 1

            if n_genes % 250000 == 0:
                logging.info("Parsed %s gene records so far...", n_genes)

    df = pd.DataFrame.from_records(records)
    if df.empty:
        raise ValueError("No gene records parsed from GTF. Check file.")
    df = filter_standard_chromosomes(df=df)
    logging.info(
        "Parsed %s gene records from GTF (%s lines scanned). After chr filter: %s",
        n_genes,
        n_lines,
        df.shape[0],
    )
    return df


def map_gene_symbols_to_gencode(
    gene_list: pd.DataFrame,
    gencode_genes: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Map input gene symbols to GENCODE genes via gene_name.

    Parameters
    ----------
    gene_list
        DataFrame containing 'gene_key' column.
    gencode_genes
        Parsed GENCODE gene records.

    Returns
    -------
    mapped_df
        Mapped gene table with unique mappings.
    report_df
        Mapping report with statuses: mapped, unmapped, ambiguous.
    """
    gencode = gencode_genes.copy()
    gencode["gene_name_upper"] = gencode["gene_name"].astype(str).str.upper()

    inp = gene_list.copy()
    inp["gene_key_upper"] = inp["gene_key"].astype(str).str.upper()

    merged = inp.merge(
        gencode,
        how="left",
        left_on="gene_key_upper",
        right_on="gene_name_upper",
    )

    report_rows: List[Dict[str, object]] = []
    mapped_rows: List[Dict[str, object]] = []

    for gene_key, sub in merged.groupby("gene_key", dropna=False):
        gene_key_str = str(gene_key)
        hits = sub.dropna(subset=["gene_id"])
        unique_gene_ids = sorted(set(hits["gene_id"].astype(str).tolist()))

        if len(unique_gene_ids) == 0:
            report_rows.append(
                {
                    "gene_key": gene_key_str,
                    "status": "unmapped",
                    "n_hits": 0,
                    "mapped_gene_id": "",
                    "mapped_gene_name": "",
                }
            )
            continue

        if len(unique_gene_ids) > 1:
            report_rows.append(
                {
                    "gene_key": gene_key_str,
                    "status": "ambiguous",
                    "n_hits": len(unique_gene_ids),
                    "mapped_gene_id": "|".join(unique_gene_ids),
                    "mapped_gene_name": "|".join(sorted(set(hits["gene_name"].tolist()))),
                }
            )
            continue

        one = hits.iloc[0]
        report_rows.append(
            {
                "gene_key": gene_key_str,
                "status": "mapped",
                "n_hits": 1,
                "mapped_gene_id": str(one["gene_id"]),
                "mapped_gene_name": str(one["gene_name"]),
            }
        )
        mapped_rows.append(
            {
                "gene_key": gene_key_str,
                "gene_id": str(one["gene_id"]),
                "gene_name": str(one["gene_name"]),
                "Chromosome": str(one["Chromosome"]),
                "Start": int(one["Start"]),
                "End": int(one["End"]),
                "Strand": str(one["Strand"]),
            }
        )

    mapped_df = pd.DataFrame(mapped_rows)
    report_df = pd.DataFrame(report_rows)

    logging.info(
        "Mapping summary: mapped=%s, unmapped=%s, ambiguous=%s",
        int((report_df["status"] == "mapped").sum()),
        int((report_df["status"] == "unmapped").sum()),
        int((report_df["status"] == "ambiguous").sum()),
    )
    return mapped_df, report_df


def build_gene_records(mapped_df: pd.DataFrame) -> List[GeneRecord]:
    """
    Build GeneRecord list from mapped gene table.

    Parameters
    ----------
    mapped_df
        DataFrame returned by map_gene_symbols_to_gencode().

    Returns
    -------
    list of GeneRecord
        Mapped gene records.
    """
    genes: List[GeneRecord] = []
    for row in mapped_df.itertuples(index=False):
        genes.append(
            GeneRecord(
                gene_id=str(row.gene_id),
                gene_name=str(row.gene_name),
                chrom=str(row.Chromosome),
                start=int(row.Start),
                end=int(row.End),
                strand=str(row.Strand),
            )
        )
    return genes


def make_tss_ranges(genes: Sequence[GeneRecord]) -> pr.PyRanges:
    """
    Create 1-bp PyRanges intervals at each gene TSS.

    Parameters
    ----------
    genes
        GeneRecord objects.

    Returns
    -------
    pr.PyRanges
        TSS intervals with 'gene_id' column.
    """
    rows = []
    for g in genes:
        tss = int(g.tss)
        rows.append({"Chromosome": g.chrom, "Start": tss, "End": tss + 1, "gene_id": g.gene_id})
    return pr.PyRanges(pd.DataFrame(rows))


def make_gene_body_ranges(genes: Sequence[GeneRecord]) -> pr.PyRanges:
    """
    Create PyRanges intervals spanning each gene body.

    Parameters
    ----------
    genes
        GeneRecord objects.

    Returns
    -------
    pr.PyRanges
        Gene body intervals with 'gene_id' column.
    """
    rows = []
    for g in genes:
        rows.append({"Chromosome": g.chrom, "Start": g.start, "End": g.end, "gene_id": g.gene_id})
    return pr.PyRanges(pd.DataFrame(rows))


def compute_nearest_gene_distance_tss(genes: Sequence[GeneRecord]) -> pd.Series:
    """
    Compute nearest neighbour distance between TSS positions (bp), per chromosome.

    Parameters
    ----------
    genes
        GeneRecord objects.

    Returns
    -------
    pd.Series
        Nearest neighbour distances, indexed by gene_id.
    """
    by_chr: Dict[str, List[Tuple[str, int]]] = {}
    for g in genes:
        by_chr.setdefault(g.chrom, []).append((g.gene_id, int(g.tss)))

    out: Dict[str, float] = {}
    for chrom, items in by_chr.items():
        items_sorted = sorted(items, key=lambda x: x[1])
        positions = np.array([p for _, p in items_sorted], dtype=int)
        gene_ids = [gid for gid, _ in items_sorted]

        if positions.size == 1:
            out[gene_ids[0]] = float("nan")
            continue

        diffs = np.diff(positions)
        left = np.concatenate(([np.iinfo(np.int64).max], diffs))
        right = np.concatenate((diffs, [np.iinfo(np.int64).max]))
        nearest = np.minimum(left, right).astype(float)

        for gid, d in zip(gene_ids, nearest):
            out[gid] = float(d)

    s = pd.Series(out, name="nn_dist_tss_bp")
    logging.info("Nearest neighbour distances computed for %s genes", s.shape[0])
    return s


def make_windows_around_points(points_pr: pr.PyRanges, window_bp: int) -> pr.PyRanges:
    """
    Create symmetric windows +/- window_bp around 1-bp point ranges.

    Parameters
    ----------
    points_pr
        1-bp point intervals (e.g., TSS positions).
    window_bp
        Window radius in base pairs.

    Returns
    -------
    pr.PyRanges
        Window intervals (clipped at 0).
    """
    df = points_pr.df.copy()
    df["Start"] = (df["Start"] - int(window_bp)).clip(lower=0).astype(int)
    df["End"] = (df["End"] + int(window_bp)).astype(int)
    return pr.PyRanges(df)


def compute_gene_density(
    tss_pr: pr.PyRanges,
    all_genes_pr: pr.PyRanges,
    windows_bp: Sequence[int],
) -> pd.DataFrame:
    """
    Compute counts of genes within windows around each TSS.

    Counts include the gene itself; we subtract 1 to exclude self.

    Parameters
    ----------
    tss_pr
        TSS PyRanges with 'gene_id'.
    all_genes_pr
        All genes PyRanges with 'gene_id'.
    windows_bp
        Sequence of window radii.

    Returns
    -------
    pd.DataFrame
        gene_id plus gene_count_tss_pm_* columns.
    """
    out = pd.DataFrame({"gene_id": tss_pr.df["gene_id"].astype(str)})
    for w in windows_bp:
        win = make_windows_around_points(points_pr=tss_pr, window_bp=int(w))
        counts = win.count_overlaps(all_genes_pr).df[["gene_id", "NumberOverlaps"]].copy()
        counts["NumberOverlaps"] = counts["NumberOverlaps"].astype(int) - 1
        counts = counts.rename(columns={"NumberOverlaps": f"gene_count_tss_pm_{w}"})
        out = out.merge(counts, how="left", on="gene_id")
        logging.debug("Gene density computed for window ±%s bp", w)
    return out


def parse_rmsk_table(rmsk_tsv_gz: str) -> pd.DataFrame:
    """
    Parse UCSC rmsk.txt.gz into a DataFrame with interval columns.

    Parameters
    ----------
    rmsk_tsv_gz
        Path to UCSC rmsk.txt.gz

    Returns
    -------
    pd.DataFrame
        Columns: Chromosome, Start, End, repClass
    """
    cols = [
        "bin",
        "swScore",
        "milliDiv",
        "milliDel",
        "milliIns",
        "genoName",
        "genoStart",
        "genoEnd",
        "genoLeft",
        "strand",
        "repName",
        "repClass",
        "repFamily",
        "repStart",
        "repEnd",
        "repLeft",
        "id",
    ]
    df = pd.read_csv(
        rmsk_tsv_gz,
        sep="\t",
        header=None,
        names=cols,
        compression="gzip",
        dtype={"genoName": str, "repClass": str},
        low_memory=False,
    )
    df = df.rename(columns={"genoName": "Chromosome", "genoStart": "Start", "genoEnd": "End"})
    df["Chromosome"] = df["Chromosome"].astype(str).map(normalise_chrom)
    df = filter_standard_chromosomes(df=df)
    df["repClass"] = df["repClass"].astype(str).str.split("/").str[0]
    out = df[["Chromosome", "Start", "End", "repClass"]].copy()
    logging.info("Loaded repeats: %s rows (std chroms)", out.shape[0])
    return out


def parse_segdup_table(segdup_tsv_gz: str) -> pd.DataFrame:
    """
    Parse UCSC genomicSuperDups.txt.gz into a DataFrame with interval columns.

    Parameters
    ----------
    segdup_tsv_gz
        Path to UCSC genomicSuperDups.txt.gz

    Returns
    -------
    pd.DataFrame
        Columns: Chromosome, Start, End
    """
    df = pd.read_csv(segdup_tsv_gz, sep="\t", header=None, compression="gzip", low_memory=False)
    if df.shape[1] < 4:
        raise ValueError("Segmental duplication table has unexpected format.")
    df = df.rename(columns={1: "Chromosome", 2: "Start", 3: "End"})
    df["Chromosome"] = df["Chromosome"].astype(str).map(normalise_chrom)
    df["Start"] = pd.to_numeric(df["Start"], errors="coerce").fillna(0).astype(int)
    df["End"] = pd.to_numeric(df["End"], errors="coerce").fillna(0).astype(int)
    df = filter_standard_chromosomes(df=df)
    out = df[["Chromosome", "Start", "End"]].copy()
    logging.info("Loaded segmental duplications: %s rows (std chroms)", out.shape[0])
    return out


def compute_repeat_fraction_in_windows(
    tss_pr: pr.PyRanges,
    repeats_pr: pr.PyRanges,
    windows_bp: Sequence[int],
) -> pd.DataFrame:
    """
    Compute fraction of bases covered by repeats in windows around each TSS.

    Parameters
    ----------
    tss_pr
        TSS PyRanges (1 bp intervals) with gene_id.
    repeats_pr
        Repeats PyRanges.
    windows_bp
        Window radii.

    Returns
    -------
    pd.DataFrame
        gene_id plus repeat_frac_tss_pm_* columns.
    """
    out = pd.DataFrame({"gene_id": tss_pr.df["gene_id"].astype(str)})

    for w in windows_bp:
        win = make_windows_around_points(points_pr=tss_pr, window_bp=int(w))
        joined = win.join(repeats_pr)
        if joined.df.empty:
            out[f"repeat_frac_tss_pm_{w}"] = 0.0
            logging.debug("No repeat overlaps found for window ±%s bp", w)
            continue

        jdf = joined.df.copy()
        overlap_start = np.maximum(jdf["Start"].to_numpy(), jdf["Start_b"].to_numpy())
        overlap_end = np.minimum(jdf["End"].to_numpy(), jdf["End_b"].to_numpy())
        overlap_len = np.maximum(0, overlap_end - overlap_start).astype(float)
        jdf["overlap_len"] = overlap_len

        sums = jdf.groupby("gene_id", as_index=False)["overlap_len"].sum()
        window_size = float(2 * int(w) + 1)
        sums[f"repeat_frac_tss_pm_{w}"] = (sums["overlap_len"] / window_size).clip(0.0, 1.0)
        sums = sums[["gene_id", f"repeat_frac_tss_pm_{w}"]]
        out = out.merge(sums, how="left", on="gene_id")
        out[f"repeat_frac_tss_pm_{w}"] = out[f"repeat_frac_tss_pm_{w}"].fillna(0.0)

        logging.debug("Repeat fraction computed for window ±%s bp", w)

    return out


def compute_distance_to_nearest(
    query_pr: pr.PyRanges,
    subject_pr: pr.PyRanges,
    out_name: str,
) -> pd.Series:
    """
    Compute distance from each query interval to the nearest subject interval.

    Parameters
    ----------
    query_pr
        Query PyRanges with 'gene_id'.
    subject_pr
        Subject PyRanges.
    out_name
        Name of the output Series.

    Returns
    -------
    pd.Series
        Distances in bp indexed by gene_id.
    """
    nearest = query_pr.nearest(subject_pr, how="any")
    ndf = nearest.df.copy()

    if ndf.empty:
        gene_ids = query_pr.df["gene_id"].astype(str).tolist()
        return pd.Series({gid: float("nan") for gid in gene_ids}, name=out_name, dtype=float)

    if "Distance" not in ndf.columns:
        ndf["Distance"] = (ndf["Start"] - ndf["Start_b"]).abs()

    dist = ndf.groupby("gene_id", as_index=False)["Distance"].min()
    s = dist.set_index("gene_id")["Distance"].astype(float)
    s.name = out_name
    return s


def compute_segdup_overlap_fraction(
    gene_body_pr: pr.PyRanges,
    segdup_pr: pr.PyRanges,
) -> pd.Series:
    """
    Compute fraction of each gene body overlapped by segmental duplications.

    Parameters
    ----------
    gene_body_pr
        Gene body intervals with gene_id.
    segdup_pr
        Segmental duplication intervals.

    Returns
    -------
    pd.Series
        Overlap fraction per gene_id.
    """
    joined = gene_body_pr.join(segdup_pr)
    gene_ids = gene_body_pr.df["gene_id"].astype(str).tolist()

    if joined.df.empty:
        return pd.Series({gid: 0.0 for gid in gene_ids}, name="segdup_overlap_frac_gene_body", dtype=float)

    jdf = joined.df.copy()
    overlap_start = np.maximum(jdf["Start"].to_numpy(), jdf["Start_b"].to_numpy())
    overlap_end = np.minimum(jdf["End"].to_numpy(), jdf["End_b"].to_numpy())
    overlap_len = np.maximum(0, overlap_end - overlap_start).astype(float)
    jdf["overlap_len"] = overlap_len

    sums = jdf.groupby("gene_id", as_index=False)["overlap_len"].sum()
    gdf = gene_body_pr.df[["gene_id", "Start", "End"]].copy()
    gdf["gene_len"] = (gdf["End"] - gdf["Start"]).astype(float).replace(0.0, np.nan)

    merged = gdf.merge(sums, how="left", on="gene_id").fillna({"overlap_len": 0.0})
    merged["segdup_overlap_frac_gene_body"] = (
        (merged["overlap_len"] / merged["gene_len"]).fillna(0.0).clip(0.0, 1.0)
    )
    s = merged.set_index("gene_id")["segdup_overlap_frac_gene_body"].astype(float)
    s.name = "segdup_overlap_frac_gene_body"
    return s


def load_fasta(hg38_fasta: str) -> Fasta:
    """
    Load FASTA with pyfaidx.

    Parameters
    ----------
    hg38_fasta
        Path to hg38 FASTA (gz ok).

    Returns
    -------
    Fasta
        pyfaidx Fasta object.
    """
    fasta = Fasta(hg38_fasta, as_raw=True, sequence_always_upper=True)
    logging.info("FASTA loaded with %s contigs", len(list(fasta.keys())))
    return fasta


def get_chrom_lengths(fasta: Fasta) -> Dict[str, int]:
    """
    Get chromosome lengths from FASTA.

    Parameters
    ----------
    fasta
        pyfaidx FASTA object.

    Returns
    -------
    dict
        Chromosome length mapping.
    """
    lengths: Dict[str, int] = {}
    for chrom in fasta.keys():
        lengths[str(chrom)] = len(fasta[str(chrom)])
    return lengths


def compute_gc_percent(seq: str) -> float:
    """
    Compute GC% over A/C/G/T bases, ignoring N/other.

    Parameters
    ----------
    seq
        DNA sequence string.

    Returns
    -------
    float
        GC percentage.
    """
    if not seq:
        return float("nan")
    s = seq.upper()
    a = s.count("A")
    c = s.count("C")
    g = s.count("G")
    t = s.count("T")
    denom = a + c + g + t
    if denom == 0:
        return float("nan")
    return 100.0 * float(c + g) / float(denom)


def safe_fetch_seq(fasta: Fasta, chrom: str, start: int, end: int) -> str:
    """
    Fetch sequence from FASTA, clipping to chromosome bounds.

    Parameters
    ----------
    fasta
        pyfaidx FASTA.
    chrom
        Chromosome name (must match FASTA keys).
    start
        0-based start.
    end
        0-based end (exclusive).

    Returns
    -------
    str
        Sequence (may be empty if unavailable).
    """
    if chrom not in fasta.keys():
        return ""
    chrom_len = len(fasta[chrom])
    s = max(0, int(start))
    e = min(int(end), int(chrom_len))
    if e <= s:
        return ""
    return str(fasta[chrom][s:e])


def compute_gc_features(
    genes: Sequence[GeneRecord],
    fasta: Fasta,
    flank_windows_bp: Sequence[int],
) -> pd.DataFrame:
    """
    Compute GC% for gene body and flanks around TSS.

    Parameters
    ----------
    genes
        GeneRecord objects.
    fasta
        FASTA object.
    flank_windows_bp
        Radii for TSS flanks.

    Returns
    -------
    pd.DataFrame
        GC features per gene_id.
    """
    rows: List[Dict[str, object]] = []
    for i, g in enumerate(genes, start=1):
        gene_seq = safe_fetch_seq(fasta=fasta, chrom=g.chrom, start=g.start, end=g.end)
        row: Dict[str, object] = {
            "gene_id": g.gene_id,
            "gc_gene_body_pct": compute_gc_percent(seq=gene_seq),
        }
        tss = int(g.tss)
        for w in flank_windows_bp:
            flank_seq = safe_fetch_seq(fasta=fasta, chrom=g.chrom, start=tss - int(w), end=tss + int(w) + 1)
            row[f"gc_tss_pm_{w}_pct"] = compute_gc_percent(seq=flank_seq)
        rows.append(row)

        if i % 500 == 0:
            logging.info("GC computed for %s/%s genes", i, len(genes))

    return pd.DataFrame(rows)


def compute_subtelomeric_features(
    genes: Sequence[GeneRecord],
    chrom_lengths: Dict[str, int],
) -> pd.DataFrame:
    """
    Compute distance to nearest chromosome end and subtelomeric flags.

    Parameters
    ----------
    genes
        GeneRecord objects.
    chrom_lengths
        Chromosome lengths.

    Returns
    -------
    pd.DataFrame
        Subtelomeric features per gene_id.
    """
    rows: List[Dict[str, object]] = []
    for g in genes:
        chrom_len = int(chrom_lengths.get(g.chrom, 0))
        tss = int(g.tss)
        if chrom_len <= 0:
            dist_end = float("nan")
        else:
            dist_left = float(tss)
            dist_right = float(max(0, chrom_len - (tss + 1)))
            dist_end = float(min(dist_left, dist_right))

        row: Dict[str, object] = {"gene_id": g.gene_id, "dist_to_chr_end_bp": dist_end}
        for thr in SUBTELO_THRESHOLDS_BP:
            row[f"is_within_{thr}_bp_of_chr_end"] = bool(not math.isnan(dist_end) and dist_end <= float(thr))
        rows.append(row)

    return pd.DataFrame(rows)


def open_bigwig(recomb_bw: str):
    """
    Open a bigWig file.

    Parameters
    ----------
    recomb_bw
        Path to bigWig.

    Returns
    -------
    pyBigWig.BigWigFile
        Open bigWig handle.
    """
    bw = pyBigWig.open(recomb_bw)
    if bw is None:
        raise ValueError(f"Failed to open bigWig: {recomb_bw}")
    return bw


def bw_stat(bw, chrom: str, start: int, end: int, stat: str) -> float:
    """
    Compute a bigWig summary statistic, returning NaN if missing.

    Parameters
    ----------
    bw
        Open bigWig handle.
    chrom
        Chromosome name.
    start
        Start coordinate.
    end
        End coordinate.
    stat
        Statistic type supported by pyBigWig (mean, max, min, etc.).

    Returns
    -------
    float
        Statistic value or NaN.
    """
    try:
        vals = bw.stats(chrom, int(start), int(end), type=stat)
    except RuntimeError:
        return float("nan")
    if not vals or vals[0] is None:
        return float("nan")
    return float(vals[0])


def compute_recombination_features(
    genes: Sequence[GeneRecord],
    recomb_bw: str,
    windows_bp: Sequence[int],
) -> pd.DataFrame:
    """
    Compute recombination statistics in windows around each TSS.

    Parameters
    ----------
    genes
        GeneRecord objects.
    recomb_bw
        Path to recombination bigWig (recombAvg.bw).
    windows_bp
        Radii for windows around TSS.

    Returns
    -------
    pd.DataFrame
        Recombination features per gene_id.
    """
    bw = open_bigwig(recomb_bw=recomb_bw)
    rows: List[Dict[str, object]] = []
    try:
        for i, g in enumerate(genes, start=1):
            tss = int(g.tss)
            row: Dict[str, object] = {"gene_id": g.gene_id}
            row["recomb_at_tss_mean"] = bw_stat(bw=bw, chrom=g.chrom, start=tss, end=tss + 1, stat="mean")

            for w in windows_bp:
                start = max(0, tss - int(w))
                end = tss + int(w) + 1
                row[f"recomb_tss_pm_{w}_mean"] = bw_stat(bw=bw, chrom=g.chrom, start=start, end=end, stat="mean")
                row[f"recomb_tss_pm_{w}_max"] = bw_stat(bw=bw, chrom=g.chrom, start=start, end=end, stat="max")

            rows.append(row)

            if i % 1000 == 0:
                logging.info("Recombination computed for %s/%s genes", i, len(genes))
    finally:
        bw.close()

    return pd.DataFrame(rows)


def write_tsv(df: pd.DataFrame, path: str) -> None:
    """
    Write a DataFrame to TSV.

    Parameters
    ----------
    df
        DataFrame to write.
    path
        Output path.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def build_feature_table(
    mapped_df: pd.DataFrame,
    all_genes_df: pd.DataFrame,
    rmsk_df: pd.DataFrame,
    segdup_df: pd.DataFrame,
    fasta: Fasta,
    recomb_bw: str,
) -> pd.DataFrame:
    """
    Compute and merge all genomic context features.

    Parameters
    ----------
    mapped_df
        Mapped genes (unique).
    all_genes_df
        All genes DataFrame for density computations.
    rmsk_df
        RepeatMasker intervals.
    segdup_df
        Segmental duplication intervals.
    fasta
        FASTA object.
    recomb_bw
        Recombination bigWig.

    Returns
    -------
    pd.DataFrame
        Feature table per gene.
    """
    genes = build_gene_records(mapped_df=mapped_df)
    logging.info("Computing features for %s mapped genes", len(genes))

    base = mapped_df.copy()
    base["gene_length_bp"] = (base["End"] - base["Start"]).astype(int)
    base["tss_0based"] = [int(g.tss) for g in genes]

    all_genes_pr = pr.PyRanges(all_genes_df[["Chromosome", "Start", "End", "gene_id"]].copy())
    tss_pr = make_tss_ranges(genes=genes)
    gene_body_pr = make_gene_body_ranges(genes=genes)

    logging.info("Computing nearest neighbour TSS distances")
    nn_dist = compute_nearest_gene_distance_tss(genes=genes).reset_index()
    nn_dist = nn_dist.rename(columns={"index": "gene_id"})
    base = base.merge(nn_dist, how="left", on="gene_id")

    logging.info("Computing gene density in windows around TSS")
    dens = compute_gene_density(tss_pr=tss_pr, all_genes_pr=all_genes_pr, windows_bp=WINDOWS_BP)
    base = base.merge(dens, how="left", on="gene_id")

    logging.info("Computing repeat features")
    repeats_pr = pr.PyRanges(rmsk_df)
    repeat_frac = compute_repeat_fraction_in_windows(tss_pr=tss_pr, repeats_pr=repeats_pr, windows_bp=REPEAT_WINDOWS_BP)
    base = base.merge(repeat_frac, how="left", on="gene_id")

    base["dist_nearest_repeat_any_bp"] = compute_distance_to_nearest(
        query_pr=tss_pr,
        subject_pr=repeats_pr,
        out_name="dist_nearest_repeat_any_bp",
    ).reindex(base["gene_id"].astype(str)).to_numpy()

    for rep_class in REPEAT_CLASSES:
        sub_df = rmsk_df[rmsk_df["repClass"] == rep_class].copy()
        col = f"dist_nearest_repeat_{rep_class.lower()}_bp"
        if sub_df.empty:
            logging.warning("No repeats found for class %s; setting %s to NaN", rep_class, col)
            base[col] = np.nan
            continue
        sub_pr = pr.PyRanges(sub_df[["Chromosome", "Start", "End"]].copy())
        dist = compute_distance_to_nearest(query_pr=tss_pr, subject_pr=sub_pr, out_name=col)
        base[col] = dist.reindex(base["gene_id"].astype(str)).to_numpy()

    logging.info("Computing segmental duplication features")
    segdup_pr = pr.PyRanges(segdup_df)
    seg_overlap = compute_segdup_overlap_fraction(gene_body_pr=gene_body_pr, segdup_pr=segdup_pr)
    base["segdup_overlap_frac_gene_body"] = seg_overlap.reindex(base["gene_id"].astype(str)).to_numpy()

    base["dist_nearest_segdup_bp"] = compute_distance_to_nearest(
        query_pr=gene_body_pr,
        subject_pr=segdup_pr,
        out_name="dist_nearest_segdup_bp",
    ).reindex(base["gene_id"].astype(str)).to_numpy()

    logging.info("Computing subtelomeric features")
    chrom_lengths = get_chrom_lengths(fasta=fasta)
    subtelo = compute_subtelomeric_features(genes=genes, chrom_lengths=chrom_lengths)
    base = base.merge(subtelo, how="left", on="gene_id")

    logging.info("Computing GC content features")
    gc_df = compute_gc_features(genes=genes, fasta=fasta, flank_windows_bp=REPEAT_WINDOWS_BP)
    base = base.merge(gc_df, how="left", on="gene_id")

    logging.info("Computing recombination features")
    recomb_df = compute_recombination_features(genes=genes, recomb_bw=recomb_bw, windows_bp=WINDOWS_BP)
    base = base.merge(recomb_df, how="left", on="gene_id")

    logging.info("Feature table complete: %s rows, %s columns", base.shape[0], base.shape[1])
    return base


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """
    Parse command-line arguments.

    Parameters
    ----------
    argv
        Argument list.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Build per-gene genomic context features (GRCh38).")
    parser.add_argument("--gene_list_tsv", required=True, help="Input TSV with column 'gene_key'.")
    parser.add_argument("--gencode_gtf", required=True, help="GENCODE GRCh38 GTF (primary assembly recommended).")
    parser.add_argument("--hg38_fasta", required=True, help="UCSC hg38 FASTA (gz ok).")
    parser.add_argument("--rmsk_tsv_gz", required=True, help="UCSC rmsk.txt.gz for hg38.")
    parser.add_argument("--segdup_tsv_gz", required=True, help="UCSC genomicSuperDups.txt.gz for hg38.")
    parser.add_argument("--recomb_bw", required=True, help="UCSC recombAvg.bw for hg38.")
    parser.add_argument("--out_tsv", required=True, help="Output TSV for gene context features.")
    parser.add_argument("--mapping_report_tsv", required=True, help="Output TSV mapping report.")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging.")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """
    Run the feature-building pipeline.

    Parameters
    ----------
    argv
        Optional list of arguments.

    Returns
    -------
    int
        Exit code.
    """
    args = parse_args(argv=argv)
    setup_logger(verbose=bool(args.verbose))

    logging.info("Loading input gene list: %s", args.gene_list_tsv)
    gene_list = read_gene_list(gene_list_tsv=args.gene_list_tsv)

    logging.info("Parsing GENCODE GTF (genes only): %s", args.gencode_gtf)
    gencode_genes = parse_gtf_genes(gencode_gtf=args.gencode_gtf)

    all_genes_df = gencode_genes[["Chromosome", "Start", "End", "gene_id"]].copy()
    logging.info("All genes for density: %s", all_genes_df.shape[0])

    logging.info("Mapping gene symbols to GENCODE gene_name")
    mapped_df, report_df = map_gene_symbols_to_gencode(gene_list=gene_list, gencode_genes=gencode_genes)

    write_tsv(df=report_df, path=args.mapping_report_tsv)
    if mapped_df.empty:
        logging.error("No genes mapped. See mapping report: %s", args.mapping_report_tsv)
        return 2

    logging.info("Loading repeats: %s", args.rmsk_tsv_gz)
    rmsk_df = parse_rmsk_table(rmsk_tsv_gz=args.rmsk_tsv_gz)

    logging.info("Loading segmental duplications: %s", args.segdup_tsv_gz)
    segdup_df = parse_segdup_table(segdup_tsv_gz=args.segdup_tsv_gz)

    logging.info("Loading FASTA (GC + chrom sizes): %s", args.hg38_fasta)
    fasta = load_fasta(hg38_fasta=args.hg38_fasta)

    logging.info("Building feature table")
    features = build_feature_table(
        mapped_df=mapped_df,
        all_genes_df=all_genes_df,
        rmsk_df=rmsk_df,
        segdup_df=segdup_df,
        fasta=fasta,
        recomb_bw=args.recomb_bw,
    )

    logging.info("Writing features TSV: %s", args.out_tsv)
    write_tsv(df=features, path=args.out_tsv)

    logging.info("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
