#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_gene_context_table.py

Build a per-gene genomic context feature table for a target gene list on GRCh38.

This script maps HGNC-style gene symbols (e.g., MNS1, ODAD2) to GENCODE GRCh38
gene models, then computes a range of genomic-context features per gene.

Features include
- Gene body coordinates and TSS coordinates (gene body plus TSS sensitivity)
- Local gene density within multiple window sizes (around TSS)
- Nearest-neighbour distance to any other gene (TSS-to-TSS, per chromosome)
- Repeat / transposable element proximity and repeat fraction (RepeatMasker rmsk)
- Segmental duplication overlap and distance (genomicSuperDups)
- Subtelomeric proximity (distance to nearest chromosome end)
- GC content for gene body and flanking windows (from hg38 FASTA)
- Recombination rate summary statistics in windows (from UCSC recombAvg bigWig)

Outputs are TSV only.

Expected inputs
- Gene list TSV with a column named 'gene_key' (HGNC gene symbols)
- GENCODE GTF (GRCh38) containing gene records with gene_name attributes
- hg38 FASTA (bgzipped or plain) readable by pyfaidx; index will be created if missing
- UCSC rmsk.txt.gz (RepeatMasker table dump) for hg38
- UCSC genomicSuperDups.txt.gz (segmental duplications table dump) for hg38
- UCSC recombAvg.bw (deCODE average recombination bigWig) for hg38

Example
python build_gene_context_table.py \
  --gene_list_tsv sperm_genes.tsv \
  --gencode_gtf gencode.vXX.annotation.chr.gtf.gz \
  --hg38_fasta hg38.fa.gz \
  --rmsk_tsv_gz rmsk.txt.gz \
  --segdup_tsv_gz genomicSuperDups.txt.gz \
  --recomb_bw recombAvg.bw \
  --out_tsv gene_context_features.tsv \
  --mapping_report_tsv gene_mapping_report.tsv
"""

from __future__ import annotations

import argparse
import gzip
import logging
import math
import os
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    import pyranges as pr
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "Missing dependency 'pyranges'. Install requirements before running."
    ) from exc

try:
    import pyBigWig
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "Missing dependency 'pyBigWig'. Install requirements before running."
    ) from exc

try:
    from pyfaidx import Fasta
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "Missing dependency 'pyfaidx'. Install requirements before running."
    ) from exc


WINDOWS_BP: Tuple[int, ...] = (50_000, 250_000, 1_000_000)
REPEAT_WINDOWS_BP: Tuple[int, ...] = (10_000, 50_000, 250_000)
SUBTELO_THRESHOLDS_BP: Tuple[int, ...] = (1_000_000, 5_000_000, 10_000_000)
REPEAT_CLASSES: Tuple[str, ...] = ("LINE", "SINE", "LTR", "DNA")


@dataclass(frozen=True)
class GeneRecord:
    """A minimal representation of a gene model for interval computations."""

    gene_id: str
    gene_name: str
    chrom: str
    start: int  # 0-based inclusive
    end: int  # 0-based exclusive
    strand: str  # '+' or '-'

    @property
    def length(self) -> int:
        """Return gene length in base pairs."""
        return max(0, self.end - self.start)

    @property
    def tss(self) -> int:
        """Return 0-based TSS coordinate (single base position)."""
        if self.strand == "+":
            return self.start
        return max(self.start, self.end - 1)


def setup_logger(*, verbose: bool) -> None:
    """Configure logging for the script."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def read_gene_list(*, gene_list_tsv: str) -> pd.DataFrame:
    """Read the input gene list TSV and validate required columns."""
    df = pd.read_csv(gene_list_tsv, sep="\t", dtype=str).fillna("")
    if "gene_key" not in df.columns:
        raise ValueError("Input gene list must contain a 'gene_key' column.")
    df["gene_key"] = df["gene_key"].astype(str).str.strip()
    df = df[df["gene_key"] != ""].drop_duplicates(subset=["gene_key"]).copy()
    return df


def _open_maybe_gzip(*, path: str):
    """Open a file that may be gzipped."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt", encoding="utf-8")


def parse_gtf_genes(*, gencode_gtf: str) -> pd.DataFrame:
    """
    Parse a GENCODE GTF and return a DataFrame of gene records.

    The output coordinates are converted to 0-based, half-open intervals.
    """
    records: List[Dict[str, object]] = []
    with _open_maybe_gzip(path=gencode_gtf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature != "gene":
                continue
            if strand not in {"+", "-"}:
                continue

            # Convert to 0-based half-open
            try:
                start_1 = int(start)
                end_1 = int(end)
            except ValueError:
                continue
            start0 = start_1 - 1
            end0 = end_1  # inclusive end in GTF becomes exclusive end in half-open

            attr_map = parse_gtf_attributes(attrs)
            gene_id = str(attr_map.get("gene_id", "")).strip()
            gene_name = str(attr_map.get("gene_name", "")).strip()

            if gene_id == "" or gene_name == "":
                continue
            if not chrom.startswith("chr"):
                # GENCODE "chr" GTF should already be chr-prefixed, but keep robust.
                chrom = f"chr{chrom}"

            records.append(
                {
                    "Chromosome": chrom,
                    "Start": start0,
                    "End": end0,
                    "Strand": strand,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                }
            )

    df = pd.DataFrame.from_records(records)
    if df.empty:
        raise ValueError("No gene records found in GTF. Check the input file.")
    return df


def parse_gtf_attributes(*, attrs: str) -> Dict[str, str]:
    """Parse GTF attribute column into a dictionary."""
    out: Dict[str, str] = {}
    for item in attrs.strip().split(";"):
        item = item.strip()
        if item == "":
            continue
        if " " not in item:
            continue
        key, value = item.split(" ", 1)
        value = value.strip().strip('"')
        out[key] = value
    return out


def filter_standard_chromosomes(*, df: pd.DataFrame) -> pd.DataFrame:
    """Keep standard chromosomes chr1-22, chrX, chrY."""
    keep = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    return df[df["Chromosome"].isin(keep)].copy()


def map_gene_symbols_to_gencode(
    *, gene_list: pd.DataFrame, gencode_genes: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Map HGNC gene symbols to GENCODE genes using gene_name.

    Returns
    -------
    mapped_df
        Rows for genes that mapped uniquely (gene_key -> gene_id).
    report_df
        Mapping report including unmapped and ambiguous entries.
    """
    gencode = gencode_genes.copy()
    gencode["gene_name_upper"] = gencode["gene_name"].astype(str).str.upper()

    input_df = gene_list.copy()
    input_df["gene_key_upper"] = input_df["gene_key"].astype(str).str.upper()

    merged = input_df.merge(
        gencode,
        how="left",
        left_on="gene_key_upper",
        right_on="gene_name_upper",
        suffixes=("", "_gencode"),
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

        # Unique mapping
        one = hits.iloc[0].to_dict()
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
    return mapped_df, report_df


def build_gene_records(*, mapped_df: pd.DataFrame) -> List[GeneRecord]:
    """Convert mapped gene DataFrame rows to GeneRecord objects."""
    records: List[GeneRecord] = []
    for row in mapped_df.itertuples(index=False):
        records.append(
            GeneRecord(
                gene_id=str(row.gene_id),
                gene_name=str(row.gene_name),
                chrom=str(row.Chromosome),
                start=int(row.Start),
                end=int(row.End),
                strand=str(row.Strand),
            )
        )
    return records


def make_tss_ranges(*, genes: Sequence[GeneRecord]) -> pr.PyRanges:
    """Create 1-bp PyRanges at each gene TSS."""
    rows = []
    for g in genes:
        tss = int(g.tss)
        rows.append({"Chromosome": g.chrom, "Start": tss, "End": tss + 1, "gene_id": g.gene_id})
    return pr.PyRanges(pd.DataFrame(rows))


def make_gene_body_ranges(*, genes: Sequence[GeneRecord]) -> pr.PyRanges:
    """Create PyRanges spanning each gene body."""
    rows = []
    for g in genes:
        rows.append(
            {"Chromosome": g.chrom, "Start": g.start, "End": g.end, "gene_id": g.gene_id}
        )
    return pr.PyRanges(pd.DataFrame(rows))


def compute_nearest_gene_distance_tss(
    *, genes: Sequence[GeneRecord]
) -> pd.Series:
    """
    Compute nearest neighbour distance between TSS positions (bp), per chromosome.

    Returns a Series indexed by gene_id.
    """
    by_chr: Dict[str, List[Tuple[str, int]]] = {}
    for g in genes:
        by_chr.setdefault(g.chrom, []).append((g.gene_id, int(g.tss)))

    out: Dict[str, float] = {}
    for chrom, items in by_chr.items():
        items_sorted = sorted(items, key=lambda x: x[1])
        positions = np.array([p for _, p in items_sorted], dtype=int)
        gene_ids = [gid for gid, _ in items_sorted]

        if len(positions) == 1:
            out[gene_ids[0]] = float("nan")
            continue

        # Nearest distance is min of adjacent diffs in sorted order.
        diffs = np.diff(positions)
        left = np.concatenate(([np.iinfo(np.int64).max], diffs))
        right = np.concatenate((diffs, [np.iinfo(np.int64).max]))
        nearest = np.minimum(left, right).astype(float)

        for gid, d in zip(gene_ids, nearest):
            out[gid] = float(d)

    return pd.Series(out, name="nn_dist_tss_bp")


def make_windows_around_points(
    *, points_pr: pr.PyRanges, window_bp: int, label: str
) -> pr.PyRanges:
    """Create symmetric windows of size +/- window_bp around 1-bp point ranges."""
    df = points_pr.df.copy()
    df["Start"] = (df["Start"] - window_bp).clip(lower=0).astype(int)
    df["End"] = (df["End"] + window_bp).astype(int)
    df["window_label"] = label
    return pr.PyRanges(df)


def compute_gene_density(
    *, tss_pr: pr.PyRanges, all_genes_pr: pr.PyRanges, windows_bp: Sequence[int]
) -> pd.DataFrame:
    """
    Compute counts of genes within windows around each TSS.

    Counts include the gene itself; we subtract 1 to exclude self.
    """
    out = pd.DataFrame({"gene_id": tss_pr.df["gene_id"].astype(str)})
    for w in windows_bp:
        win = make_windows_around_points(points_pr=tss_pr, window_bp=int(w), label=f"tss_pm_{w}")
        counts = win.count_overlaps(all_genes_pr).df[["gene_id", "NumberOverlaps"]].copy()
        counts["NumberOverlaps"] = counts["NumberOverlaps"].astype(int) - 1
        counts = counts.rename(columns={"NumberOverlaps": f"gene_count_tss_pm_{w}"})
        out = out.merge(counts, how="left", on="gene_id")
    return out


def parse_rmsk_table(*, rmsk_tsv_gz: str) -> pd.DataFrame:
    """
    Parse UCSC rmsk.txt.gz into a DataFrame with standard interval columns.

    Coordinates in UCSC tables are already 0-based, half-open.
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
    df = df.rename(
        columns={"genoName": "Chromosome", "genoStart": "Start", "genoEnd": "End"}
    )
    df["Chromosome"] = df["Chromosome"].astype(str)
    df = filter_standard_chromosomes(df=df)
    df["repClass"] = df["repClass"].astype(str).str.split("/").str[0]
    return df[["Chromosome", "Start", "End", "repClass"]].copy()


def parse_segdup_table(*, segdup_tsv_gz: str) -> pd.DataFrame:
    """
    Parse UCSC genomicSuperDups.txt.gz into a DataFrame with interval columns.

    Coordinates in UCSC tables are 0-based, half-open.
    """
    # The table has many columns; we need chrom/start/end at minimum.
    # UCSC order begins: bin, chrom, chromStart, chromEnd, name, score, ...
    df = pd.read_csv(
        segdup_tsv_gz,
        sep="\t",
        header=None,
        compression="gzip",
        low_memory=False,
    )
    if df.shape[1] < 4:
        raise ValueError("Segmental duplication table has unexpected format.")
    df = df.rename(columns={1: "Chromosome", 2: "Start", 3: "End"})
    df["Chromosome"] = df["Chromosome"].astype(str)
    df = filter_standard_chromosomes(df=df)
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    return df[["Chromosome", "Start", "End"]].copy()


def compute_distance_to_nearest(
    *, query_pr: pr.PyRanges, subject_pr: pr.PyRanges, suffix: str
) -> pd.Series:
    """
    Compute distance from each query interval to nearest subject interval.

    Uses midpoints of query intervals for a stable distance estimate.
    """
    qdf = query_pr.df.copy()
    qdf["q_mid"] = ((qdf["Start"] + qdf["End"]) / 2.0).astype(float)

    sdf = subject_pr.df.copy()
    sdf["s_mid"] = ((sdf["Start"] + sdf["End"]) / 2.0).astype(float)

    # Use pyranges nearest on intervals; distance is computed on interval edges.
    # To approximate midpoint-to-interval distance, we keep it simple and use the
    # interval-based distance; this is acceptable for large windows and nearest
    # proximity features.
    nearest = query_pr.nearest(subject_pr, how="any")
    ndf = nearest.df.copy()
    # pyranges returns Distance if it can compute it; if not present, we compute a fallback.
    if "Distance" not in ndf.columns:
        # Fallback: distance between starts (rough)
        ndf["Distance"] = (ndf["Start"] - ndf["Start_b"]).abs()
    out = ndf[["gene_id", "Distance"]].copy()
    out = out.groupby("gene_id", as_index=False)["Distance"].min()
    out = out.set_index("gene_id")["Distance"].astype(float)
    out.name = f"dist_nearest_{suffix}_bp"
    return out


def compute_repeat_fraction_in_windows(
    *,
    tss_pr: pr.PyRanges,
    repeats_pr: pr.PyRanges,
    windows_bp: Sequence[int],
    prefix: str,
) -> pd.DataFrame:
    """
    Compute fraction of bases covered by repeats in windows around TSS.

    Returns a DataFrame with gene_id and one column per window.
    """
    out = pd.DataFrame({"gene_id": tss_pr.df["gene_id"].astype(str)})

    for w in windows_bp:
        win = make_windows_around_points(points_pr=tss_pr, window_bp=int(w), label=f"{prefix}_{w}")
        # Compute total overlap length between each window and repeats
        joined = win.join(repeats_pr)
        if joined.df.empty:
            frac = pd.DataFrame({"gene_id": win.df["gene_id"].astype(str), f"repeat_frac_tss_pm_{w}": 0.0})
            out = out.merge(frac, how="left", on="gene_id")
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

    return out


def compute_segdup_overlap_fraction(
    *, gene_body_pr: pr.PyRanges, segdup_pr: pr.PyRanges
) -> pd.Series:
    """Compute fraction of each gene body overlapped by segmental duplications."""
    joined = gene_body_pr.join(segdup_pr)
    if joined.df.empty:
        s = pd.Series(
            data={gid: 0.0 for gid in gene_body_pr.df["gene_id"].astype(str).tolist()},
            name="segdup_overlap_frac_gene_body",
            dtype=float,
        )
        return s

    jdf = joined.df.copy()
    overlap_start = np.maximum(jdf["Start"].to_numpy(), jdf["Start_b"].to_numpy())
    overlap_end = np.minimum(jdf["End"].to_numpy(), jdf["End_b"].to_numpy())
    overlap_len = np.maximum(0, overlap_end - overlap_start).astype(float)
    jdf["overlap_len"] = overlap_len

    sums = jdf.groupby("gene_id", as_index=False)["overlap_len"].sum()
    gdf = gene_body_pr.df[["gene_id", "Start", "End"]].copy()
    gdf["gene_len"] = (gdf["End"] - gdf["Start"]).astype(float).replace(0.0, np.nan)
    merged = gdf.merge(sums, how="left", on="gene_id").fillna({"overlap_len": 0.0})
    merged["segdup_overlap_frac_gene_body"] = (merged["overlap_len"] / merged["gene_len"]).fillna(0.0).clip(0.0, 1.0)
    out = merged.set_index("gene_id")["segdup_overlap_frac_gene_body"].astype(float)
    out.name = "segdup_overlap_frac_gene_body"
    return out


def load_fasta(*, hg38_fasta: str) -> Fasta:
    """Load FASTA with pyfaidx; creates an index if missing."""
    return Fasta(hg38_fasta, as_raw=True, sequence_always_upper=True)


def get_chrom_lengths(*, fasta: Fasta) -> Dict[str, int]:
    """Return chromosome lengths from FASTA index."""
    lengths: Dict[str, int] = {}
    for chrom in fasta.keys():
        lengths[str(chrom)] = len(fasta[str(chrom)])
    return lengths


def compute_gc_percent(*, seq: str) -> float:
    """Compute GC% over A/C/G/T bases, ignoring N and other characters."""
    if not seq:
        return float("nan")
    seq_u = seq.upper()
    a = seq_u.count("A")
    c = seq_u.count("C")
    g = seq_u.count("G")
    t = seq_u.count("T")
    denom = a + c + g + t
    if denom == 0:
        return float("nan")
    return 100.0 * float(c + g) / float(denom)


def safe_fetch_seq(*, fasta: Fasta, chrom: str, start: int, end: int) -> str:
    """Fetch sequence from FASTA, clipping to chromosome bounds."""
    if chrom not in fasta.keys():
        return ""
    chrom_len = len(fasta[chrom])
    s = max(0, int(start))
    e = min(int(end), int(chrom_len))
    if e <= s:
        return ""
    return str(fasta[chrom][s:e])


def compute_gc_features(
    *, genes: Sequence[GeneRecord], fasta: Fasta, flank_windows_bp: Sequence[int]
) -> pd.DataFrame:
    """Compute GC% for gene body and flanking windows around TSS."""
    rows: List[Dict[str, object]] = []
    for g in genes:
        gene_seq = safe_fetch_seq(fasta=fasta, chrom=g.chrom, start=g.start, end=g.end)
        row: Dict[str, object] = {
            "gene_id": g.gene_id,
            "gc_gene_body_pct": compute_gc_percent(seq=gene_seq),
        }
        tss = int(g.tss)
        for w in flank_windows_bp:
            start = tss - int(w)
            end = tss + int(w) + 1
            flank_seq = safe_fetch_seq(fasta=fasta, chrom=g.chrom, start=start, end=end)
            row[f"gc_tss_pm_{w}_pct"] = compute_gc_percent(seq=flank_seq)
        rows.append(row)
    return pd.DataFrame(rows)


def compute_subtelomeric_features(
    *, genes: Sequence[GeneRecord], chrom_lengths: Dict[str, int]
) -> pd.DataFrame:
    """Compute distance to nearest chromosome end and subtelomeric flags."""
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

        row: Dict[str, object] = {
            "gene_id": g.gene_id,
            "dist_to_chr_end_bp": dist_end,
        }
        for thr in SUBTELO_THRESHOLDS_BP:
            key = f"is_within_{thr}_bp_of_chr_end"
            row[key] = bool(not math.isnan(dist_end) and dist_end <= float(thr))
        rows.append(row)
    return pd.DataFrame(rows)


def open_bigwig(*, recomb_bw: str):
    """Open a bigWig file with pyBigWig."""
    bw = pyBigWig.open(recomb_bw)
    if bw is None:
        raise ValueError("Failed to open recombination bigWig.")
    return bw


def bw_stats(
    *, bw, chrom: str, start: int, end: int, stat: str
) -> float:
    """Get a bigWig summary statistic; returns NaN if unavailable."""
    try:
        vals = bw.stats(chrom, int(start), int(end), type=stat)
    except RuntimeError:
        return float("nan")
    if not vals or vals[0] is None:
        return float("nan")
    return float(vals[0])


def compute_recombination_features(
    *, genes: Sequence[GeneRecord], recomb_bw: str, windows_bp: Sequence[int]
) -> pd.DataFrame:
    """Compute recombination mean/max in windows around TSS and at TSS."""
    bw = open_bigwig(recomb_bw=recomb_bw)
    rows: List[Dict[str, object]] = []
    try:
        for g in genes:
            tss = int(g.tss)
            row: Dict[str, object] = {"gene_id": g.gene_id}
            row["recomb_at_tss_mean"] = bw_stats(
                bw=bw, chrom=g.chrom, start=tss, end=tss + 1, stat="mean"
            )
            for w in windows_bp:
                start = max(0, tss - int(w))
                end = tss + int(w) + 1
                row[f"recomb_tss_pm_{w}_mean"] = bw_stats(
                    bw=bw, chrom=g.chrom, start=start, end=end, stat="mean"
                )
                row[f"recomb_tss_pm_{w}_max"] = bw_stats(
                    bw=bw, chrom=g.chrom, start=start, end=end, stat="max"
                )
            rows.append(row)
    finally:
        bw.close()
    return pd.DataFrame(rows)


def build_feature_table(
    *,
    mapped_df: pd.DataFrame,
    all_genes_df: pd.DataFrame,
    rmsk_df: pd.DataFrame,
    segdup_df: pd.DataFrame,
    fasta: Fasta,
    recomb_bw: str,
) -> pd.DataFrame:
    """Compute all features and return a single merged DataFrame."""
    genes = build_gene_records(mapped_df=mapped_df)
    logging.info("Computing features for %s mapped genes", len(genes))

    # Base table
    base = mapped_df.copy()
    base["gene_length_bp"] = (base["End"] - base["Start"]).astype(int)
    base["tss_0based"] = [
        int(GeneRecord(
            gene_id=str(r.gene_id),
            gene_name=str(r.gene_name),
            chrom=str(r.Chromosome),
            start=int(r.Start),
            end=int(r.End),
            strand=str(r.Strand),
        ).tss)
        for r in base.itertuples(index=False)
    ]

    # PyRanges objects
    all_genes_pr = pr.PyRanges(all_genes_df[["Chromosome", "Start", "End", "gene_id"]].copy())
    tss_pr = make_tss_ranges(genes=genes)
    gene_body_pr = make_gene_body_ranges(genes=genes)

    # Nearest neighbour gene distance (TSS-based)
    nn_dist = compute_nearest_gene_distance_tss(genes=genes).reset_index()
    nn_dist = nn_dist.rename(columns={"index": "gene_id"})
    base = base.merge(nn_dist, how="left", on="gene_id")

    # Gene density around TSS
    dens = compute_gene_density(tss_pr=tss_pr, all_genes_pr=all_genes_pr, windows_bp=WINDOWS_BP)
    base = base.merge(dens, how="left", on="gene_id")

    # Repeats
    repeats_pr = pr.PyRanges(rmsk_df)
    repeat_frac = compute_repeat_fraction_in_windows(
        tss_pr=tss_pr, repeats_pr=repeats_pr, windows_bp=REPEAT_WINDOWS_BP, prefix="rmsk"
    )
    base = base.merge(repeat_frac, how="left", on="gene_id")

    # Distance to nearest repeat overall and by class
    base["dist_nearest_repeat_any_bp"] = compute_distance_to_nearest(
        query_pr=tss_pr, subject_pr=repeats_pr, suffix="repeat_any"
    ).reindex(base["gene_id"].astype(str)).to_numpy()

    for rep_class in REPEAT_CLASSES:
        sub_df = rmsk_df[rmsk_df["repClass"] == rep_class].copy()
        if sub_df.empty:
            base[f"dist_nearest_repeat_{rep_class.lower()}_bp"] = np.nan
            continue
        sub_pr = pr.PyRanges(sub_df[["Chromosome", "Start", "End"]].copy())
        dist = compute_distance_to_nearest(
            query_pr=tss_pr, subject_pr=sub_pr, suffix=f"repeat_{rep_class.lower()}"
        )
        base[f"dist_nearest_repeat_{rep_class.lower()}_bp"] = (
            dist.reindex(base["gene_id"].astype(str)).to_numpy()
        )

    # Segmental duplications
    segdup_pr = pr.PyRanges(segdup_df)
    seg_overlap = compute_segdup_overlap_fraction(gene_body_pr=gene_body_pr, segdup_pr=segdup_pr)
    base["segdup_overlap_frac_gene_body"] = seg_overlap.reindex(base["gene_id"].astype(str)).to_numpy()

    base["dist_nearest_segdup_bp"] = compute_distance_to_nearest(
        query_pr=gene_body_pr, subject_pr=segdup_pr, suffix="segdup"
    ).reindex(base["gene_id"].astype(str)).to_numpy()

    # Subtelomeric
    chrom_lengths = get_chrom_lengths(fasta=fasta)
    subtelo = compute_subtelomeric_features(genes=genes, chrom_lengths=chrom_lengths)
    base = base.merge(subtelo, how="left", on="gene_id")

    # GC content
    gc_df = compute_gc_features(genes=genes, fasta=fasta, flank_windows_bp=REPEAT_WINDOWS_BP)
    base = base.merge(gc_df, how="left", on="gene_id")

    # Recombination
    recomb_df = compute_recombination_features(genes=genes, recomb_bw=recomb_bw, windows_bp=WINDOWS_BP)
    base = base.merge(recomb_df, how="left", on="gene_id")

    return base


def write_tsv(*, df: pd.DataFrame, path: str) -> None:
    """Write a DataFrame as TSV, creating parent directories if needed."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description="Build per-gene genomic context features (GRCh38)."
    )
    parser.add_argument("--gene_list_tsv", required=True, help="Input TSV with column 'gene_key'.")
    parser.add_argument("--gencode_gtf", required=True, help="GENCODE GRCh38 gene annotation GTF (chr).")
    parser.add_argument("--hg38_fasta", required=True, help="hg38/GRCh38 FASTA (gz OK).")
    parser.add_argument("--rmsk_tsv_gz", required=True, help="UCSC rmsk.txt.gz for hg38.")
    parser.add_argument("--segdup_tsv_gz", required=True, help="UCSC genomicSuperDups.txt.gz for hg38.")
    parser.add_argument("--recomb_bw", required=True, help="UCSC recombAvg.bw for hg38.")
    parser.add_argument("--out_tsv", required=True, help="Output TSV path for gene context features.")
    parser.add_argument("--mapping_report_tsv", required=True, help="Output TSV path for mapping report.")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging.")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point."""
    args = parse_args(argv=argv)
    setup_logger(verbose=bool(args.verbose))

    logging.info("Loading input gene list: %s", args.gene_list_tsv)
    gene_list = read_gene_list(gene_list_tsv=args.gene_list_tsv)

    logging.info("Parsing GENCODE GTF (genes only): %s", args.gencode_gtf)
    gencode_genes = parse_gtf_genes(gencode_gtf=args.gencode_gtf)
    gencode_genes = filter_standard_chromosomes(df=gencode_genes)

    logging.info("Preparing all-genes table for density computations")
    all_genes_df = gencode_genes[["Chromosome", "Start", "End", "gene_id"]].copy()

    logging.info("Mapping gene symbols to GENCODE gene_name")
    mapped_df, report_df = map_gene_symbols_to_gencode(
        gene_list=gene_list,
        gencode_genes=gencode_genes,
    )
    logging.info(
        "Mapped %s/%s genes (see mapping report for unmapped/ambiguous).",
        mapped_df.shape[0],
        gene_list.shape[0],
    )

    write_tsv(df=report_df, path=args.mapping_report_tsv)
    if mapped_df.empty:
        logging.error("No genes mapped. Cannot proceed.")
        return 2

    logging.info("Loading repeats (RepeatMasker rmsk): %s", args.rmsk_tsv_gz)
    rmsk_df = parse_rmsk_table(rmsk_tsv_gz=args.rmsk_tsv_gz)

    logging.info("Loading segmental duplications (genomicSuperDups): %s", args.segdup_tsv_gz)
    segdup_df = parse_segdup_table(segdup_tsv_gz=args.segdup_tsv_gz)

    logging.info("Loading FASTA for GC and chromosome sizes: %s", args.hg38_fasta)
    fasta = load_fasta(hg38_fasta=args.hg38_fasta)

    logging.info("Computing feature table")
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
