#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summarise_pxd014618_toppic_output_table.py

Summarise TopPIC '*_ms2.OUTPUT_TABLE' files from PRIDE PXD014618 into a
gene-level proteomics evidence table suitable for annotation of testis
candidate genes.

Why this is needed
------------------
The TopPIC OUTPUT_TABLE is proteoform/PRSM-level output. It does not provide
protein-level unique peptide counts or sequence coverage fields, so gene-level
evidence is derived instead from:
- Gene symbol parsed from 'Protein name' (GN=...)
- Spectral FDR q-values and proteoform FDR (where present)
- Feature intensity (as a proxy for abundance)

Outputs
-------
A tab-separated table containing, per gene:
- gene_key (normalised gene symbol)
- prot_present_any (True if present in >=1 replicate at the chosen FDR)
- prot_present_fraction (fraction of replicates where present at the chosen FDR)
- prot_n_replicates_present
- prot_n_replicates_total
- prot_n_prsms (rows supporting the gene across all replicates)
- prot_max_feature_intensity
- prot_min_q_value
- prot_min_proteoform_fdr
- plus per-replicate presence columns (optional)

All arguments are named. Output is TSV (not comma-separated).
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Summarise TopPIC '*_ms2.OUTPUT_TABLE' files from PXD014618 into a "
            "gene-level proteomics evidence TSV."
        )
    )
    parser.add_argument(
        "--input_tables",
        required=True,
        nargs="+",
        type=Path,
        help="One or more TopPIC '*_ms2.OUTPUT_TABLE' files to summarise.",
    )
    parser.add_argument(
        "--out_tsv",
        required=True,
        type=Path,
        help="Output TSV path.",
    )
    parser.add_argument(
        "--q_value_threshold",
        required=False,
        default=0.01,
        type=float,
        help="Spectral FDR q-value threshold for calling presence (default: 0.01).",
    )
    parser.add_argument(
        "--proteoform_fdr_threshold",
        required=False,
        default=0.01,
        type=float,
        help="Proteoform FDR threshold for calling presence (default: 0.01).",
    )
    parser.add_argument(
        "--require_both_thresholds",
        action="store_true",
        help=(
            "If set, require BOTH q-value and proteoform FDR thresholds to be met. "
            "If not set, presence requires q-value threshold only (proteoform FDR used if available)."
        ),
    )
    parser.add_argument(
        "--write_per_replicate_columns",
        action="store_true",
        help="If set, include per-replicate presence columns in the output.",
    )
    return parser.parse_args()


def _normalise_gene_symbol(value: str) -> str:
    """
    Normalise a gene symbol for joins.

    Parameters
    ----------
    value : str
        Gene symbol.

    Returns
    -------
    str
        Normalised gene symbol (uppercase, no whitespace).
    """
    v = str(value).strip()
    v = re.sub(r"\s+", "", v)
    return v.upper()




import re
from pathlib import Path
from typing import List

import pandas as pd


def _find_table_start_line(lines: List[str]) -> int:
    """
    Identify the line index where the tabular section begins.

    Parameters
    ----------
    lines : List[str]
        All lines in the file.

    Returns
    -------
    int
        Line index of the header row.

    Raises
    ------
    ValueError
        If the table header row cannot be found.
    """
    for i, line in enumerate(lines):
        if line.strip().startswith("Data file name"):
            return i
    raise ValueError("Could not find TopPIC table header line starting with 'Data file name'.")


def _read_toppic_output_table(path: Path) -> pd.DataFrame:
    """
    Read a TopPIC '*_ms2.OUTPUT_TABLE' file into a DataFrame.

    These files contain a parameter preamble followed by a tab-delimited table.
    Some rows may wrap across multiple lines (typically within the Proteoform
    field). This reader stitches wrapped lines back into a single record.

    Parameters
    ----------
    path : Path
        Path to a TopPIC OUTPUT_TABLE.

    Returns
    -------
    pd.DataFrame
        Parsed table.
    """
    lines = path.read_text(errors="replace").splitlines()
    start_idx = _find_table_start_line(lines)

    header = lines[start_idx].rstrip("\n")

    # Fix the known missing separator in the header
    header = header.replace("Feature intensityProtein name", "Feature intensity\tProtein name")

    # Data rows start with a file path ending in '.msalign' (often 'C:/...').
    row_start_re = re.compile(r"^(?:[A-Za-z]:[/\\]|/).+\.msalign\b")

    records: List[str] = []
    for raw in lines[start_idx + 1 :]:
        line = raw.rstrip("\n")
        if not line.strip():
            continue

        if row_start_re.match(line):
            records.append(line)
        else:
            # Continuation line (wrapped proteoform etc.) -> append to previous record
            if not records:
                continue
            records[-1] = records[-1] + " " + line.strip()

    table_text = header + "\n" + "\n".join(records)

    df = pd.read_csv(
        pd.io.common.StringIO(table_text),
        sep="\t",
        dtype=str,
        engine="python",
    )

    df.columns = [str(c).strip() for c in df.columns]
    return df




def _extract_gene_from_protein_name(protein_name: str) -> Optional[str]:
    """
    Extract gene symbol from a TopPIC 'Protein name' string containing 'GN=...'.

    Parameters
    ----------
    protein_name : str
        Protein name field, e.g. "... GN=PRM1 PE=1 SV=2".

    Returns
    -------
    Optional[str]
        Gene symbol, or None if not found.
    """
    if protein_name is None:
        return None
    m = re.search(r"GN=([A-Za-z0-9\-]+)", str(protein_name))
    if not m:
        return None
    return _normalise_gene_symbol(m.group(1))


def _to_float(series: pd.Series) -> pd.Series:
    """
    Convert a Series to numeric floats safely.

    Parameters
    ----------
    series : pd.Series
        Input series.

    Returns
    -------
    pd.Series
        Float series (NaN for non-numeric).
    """
    return pd.to_numeric(series, errors="coerce")


def summarise_tables(
    tables: List[Tuple[str, pd.DataFrame]],
    q_value_threshold: float,
    proteoform_fdr_threshold: float,
    require_both_thresholds: bool,
    write_per_replicate_columns: bool,
) -> pd.DataFrame:
    """
    Summarise multiple TopPIC tables into a combined gene-level evidence table.

    Parameters
    ----------
    tables : List[Tuple[str, pd.DataFrame]]
        List of (replicate_id, DataFrame) pairs.
    q_value_threshold : float
        Spectral q-value threshold.
    proteoform_fdr_threshold : float
        Proteoform FDR threshold.
    require_both_thresholds : bool
        Whether presence requires both thresholds.
    write_per_replicate_columns : bool
        Whether to include per-replicate presence columns.

    Returns
    -------
    pd.DataFrame
        Gene-level summary table.
    """
    per_rep_rows: List[pd.DataFrame] = []

    for rep_id, df in tables:
        if "Protein name" not in df.columns:
            raise ValueError(
                f"Expected column 'Protein name' not found in {rep_id}. "
                f"Columns seen: {list(df.columns)}"
            )

        work = df.copy()

        # Extract gene symbols
        work["gene_key"] = work["Protein name"].map(_extract_gene_from_protein_name)
        work = work[work["gene_key"].notna()].copy()
        work["gene_key"] = work["gene_key"].astype(str)

        # Extract key numeric columns (if present)
        q_col = "Q-value (spectral FDR)"
        pfdr_col = "Proteoform FDR"
        intensity_col = "Feature intensity"

        if q_col in work.columns:
            work["_q_value"] = _to_float(work[q_col])
        else:
            work["_q_value"] = np.nan

        if pfdr_col in work.columns:
            work["_proteoform_fdr"] = _to_float(work[pfdr_col])
        else:
            work["_proteoform_fdr"] = np.nan

        if intensity_col in work.columns:
            work["_feature_intensity"] = _to_float(work[intensity_col])
        else:
            work["_feature_intensity"] = np.nan

        # Presence call
        q_pass = work["_q_value"].le(float(q_value_threshold))
        if require_both_thresholds:
            pfdr_pass = work["_proteoform_fdr"].le(float(proteoform_fdr_threshold))
            present_mask = q_pass & pfdr_pass
        else:
            present_mask = q_pass

        work["_present"] = present_mask.fillna(False)

        rep_summary = (
            work.groupby("gene_key", as_index=False)
            .agg(
                prot_present=("._present".replace(".", "_"), "any")
                if False
                else ("_present", "any"),
                prot_n_prsms=("_present", "size"),
                prot_max_feature_intensity=("_feature_intensity", "max"),
                prot_min_q_value=("_q_value", "min"),
                prot_min_proteoform_fdr=("_proteoform_fdr", "min"),
            )
        )
        rep_summary["replicate_id"] = rep_id
        per_rep_rows.append(rep_summary)

    per_rep = pd.concat(per_rep_rows, ignore_index=True)

    # Combined summary across replicates
    n_reps_total = per_rep["replicate_id"].nunique()

    combined = (
        per_rep.groupby("gene_key", as_index=False)
        .agg(
            prot_n_prsms=("prot_n_prsms", "sum"),
            prot_max_feature_intensity=("prot_max_feature_intensity", "max"),
            prot_min_q_value=("prot_min_q_value", "min"),
            prot_min_proteoform_fdr=("prot_min_proteoform_fdr", "min"),
            prot_n_replicates_present=("prot_present", "sum"),
        )
    )
    combined["prot_n_replicates_total"] = int(n_reps_total)
    combined["prot_present_fraction"] = (
        combined["prot_n_replicates_present"] / combined["prot_n_replicates_total"]
    )
    combined["prot_present_any"] = combined["prot_n_replicates_present"].astype(int).ge(1)

    if write_per_replicate_columns:
        piv = per_rep.pivot_table(
            index="gene_key",
            columns="replicate_id",
            values="prot_present",
            aggfunc="max",
            fill_value=False,
        ).reset_index()
        piv.columns = [
            "gene_key"
            if c == "gene_key"
            else f"prot_present__{str(c)}"
            for c in piv.columns
        ]
        combined = combined.merge(piv, on="gene_key", how="left")

    # Sort for convenience
    combined = combined.sort_values(
        by=["prot_present_any", "prot_present_fraction", "prot_max_feature_intensity"],
        ascending=[False, False, False],
        na_position="last",
    )

    return combined


def main() -> None:
    """
    Entry point.
    """
    args = parse_args()

    tables: List[Tuple[str, pd.DataFrame]] = []
    for p in args.input_tables:
        rep_id = p.stem
        df = _read_toppic_output_table(path=p)
        tables.append((rep_id, df))

    out = summarise_tables(
        tables=tables,
        q_value_threshold=float(args.q_value_threshold),
        proteoform_fdr_threshold=float(args.proteoform_fdr_threshold),
        require_both_thresholds=bool(args.require_both_thresholds),
        write_per_replicate_columns=bool(args.write_per_replicate_columns),
    )

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)

    print(f"Replicates: {len(tables)}")
    print(f"Genes with proteomics evidence: {int(out['prot_present_any'].sum())}")
    print(f"Wrote: {args.out_tsv}")


if __name__ == "__main__":
    main()
