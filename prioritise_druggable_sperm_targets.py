#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prioritise_druggable_sperm_targets.py

Prioritise candidate sperm targets by combining:
- membership across multiple input gene lists (e.g. HPO, HC testis, literature)
- optional proteomics evidence tiers
- optional Open Targets tractability annotations
- optional genomic/regulatory context features (already computed per gene)

This script is intentionally conservative and transparent:
- it produces a ranked TSV with a simple, auditable scoring model
- it keeps a "score_components" column explaining where the score came from

Primary use-case
----------------
Identify novel, potentially druggable sperm protein targets by:
1) prioritising genes with proteomic support (protein in sperm)
2) prioritising genes with tractability evidence (small molecule / antibody / PROTAC)
3) prioritising genes with multiple independent lines of sperm relevance
4) optionally downweighting genes already present in a "literature" list (novelty)

Outputs
-------
Writes:
1) Main ranked table:
   - <out_prefix>.tsv
   - <out_prefix>.xlsx

2) Convenience subsets:
   - <out_prefix>__candidate_druggable_sperm_protein.tsv
   - <out_prefix>__candidate_druggable_sperm_protein.xlsx
   - (optional) <out_prefix>__topN.tsv and .xlsx if --top_n > 0

Excel workbook structure
------------------------
Each .xlsx contains multiple sheets (where applicable):
- ranked (always)
- candidate_druggable_sperm_protein (always; may be empty)
- topN (only if --top_n > 0; may be smaller than N after filtering)

Notes
-----
- All outputs are TSV and XLSX (tab-separated for TSV; Excel is workbook format).
- This script ranks candidates; it does not "prove" druggability.
"""

from __future__ import annotations

import argparse
import logging
import os
import re
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Weights:
    """Weights controlling the priority scoring model."""
    per_list_member: float = 1.0
    proteomics_detected: float = 1.0
    proteomics_strong: float = 3.0
    tract_small_molecule: float = 4.0
    tract_antibody: float = 2.0
    tract_protac: float = 2.0
    novelty_bonus: float = 1.5


def setup_logger(verbose: bool) -> None:
    """
    Configure logging.

    Parameters
    ----------
    verbose
        If True, enable DEBUG logging.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def read_tsv(path: str) -> pd.DataFrame:
    """
    Read a TSV with robust defaults.

    Parameters
    ----------
    path
        Input TSV path.

    Returns
    -------
    pd.DataFrame
        DataFrame (strings preserved where possible).
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = df.fillna("")
    return df


def normalise_gene_key(series: pd.Series) -> pd.Series:
    """
    Normalise gene symbols for reliable matching.

    Parameters
    ----------
    series
        Input gene key series.

    Returns
    -------
    pd.Series
        Upper-cased, stripped gene keys.
    """
    return series.astype(str).str.strip().str.upper()


def list_tsv_files(gene_lists_dir: str) -> List[str]:
    """
    List TSV files in a directory (non-recursive).

    Parameters
    ----------
    gene_lists_dir
        Directory containing gene list TSV files.

    Returns
    -------
    list of str
        TSV file paths.
    """
    files: List[str] = []
    for name in sorted(os.listdir(gene_lists_dir)):
        if name.startswith("."):
            continue
        if name.lower().endswith(".tsv"):
            files.append(os.path.join(gene_lists_dir, name))
    return files


def load_gene_lists(
    gene_lists_dir: str,
    gene_key_column: str,
) -> Tuple[pd.DataFrame, Dict[str, Set[str]]]:
    """
    Load all gene list TSVs from a directory and build a membership matrix.

    Parameters
    ----------
    gene_lists_dir
        Directory with TSV files.
    gene_key_column
        Column name containing HGNC-like gene symbols.

    Returns
    -------
    membership_df
        DataFrame with columns: gene_key, boolean columns per list, n_lists, list_names.
    list_to_genes
        Mapping from list stem name to set of gene_key values.
    """
    paths = list_tsv_files(gene_lists_dir=gene_lists_dir)
    if not paths:
        raise ValueError(f"No TSV files found in gene_lists_dir: {gene_lists_dir}")

    list_to_genes: Dict[str, Set[str]] = {}
    all_genes: Set[str] = set()

    logging.info("Discovered %s TSV files in gene_lists_dir", len(paths))

    for path in paths:
        stem = os.path.splitext(os.path.basename(path))[0]
        df = read_tsv(path)

        if gene_key_column not in df.columns:
            logging.warning("Skipping %s (missing column '%s')", path, gene_key_column)
            continue

        keys = normalise_gene_key(df[gene_key_column])
        keys = keys[keys != ""].drop_duplicates()
        gene_set = set(keys.tolist())

        list_to_genes[stem] = gene_set
        all_genes.update(gene_set)

        logging.info("Loaded list '%s': %s genes (%s)", stem, len(gene_set), path)

    if not list_to_genes:
        raise ValueError(
            f"No valid gene lists found (none had column '{gene_key_column}')."
        )

    all_genes_sorted = sorted(all_genes)
    membership = pd.DataFrame({"gene_key": all_genes_sorted})

    for stem, gene_set in list_to_genes.items():
        membership[stem] = membership["gene_key"].isin(gene_set)

    list_cols = [c for c in membership.columns if c != "gene_key"]
    membership["n_lists"] = membership[list_cols].sum(axis=1).astype(int)

    # Vectorised list_names (no row-wise apply)
    col_arr = np.array(list_cols, dtype=object)
    mask = membership[list_cols].to_numpy(dtype=bool)
    membership["list_names"] = [
        "|".join(col_arr[row_mask].tolist()) if row_mask.any() else ""
        for row_mask in mask
    ]

    logging.info(
        "Membership matrix: %s genes, %s list columns",
        membership.shape[0],
        len(list_cols),
    )

    return membership, list_to_genes


def derive_proteomics_class(df: pd.DataFrame) -> pd.DataFrame:
    """
    Derive a simple proteomics evidence class if possible.

    Logic (if numeric columns exist):
      Strong: unique_peptides >= 2 AND coverage_pct >= 10
      Detected: present_any == 1 OR (unique_peptides > 0 OR coverage_pct > 0)
      None: otherwise

    Parameters
    ----------
    df
        Proteomics DataFrame.

    Returns
    -------
    pd.DataFrame
        DataFrame with 'prot_class' column added if derivable.
    """
    out = df.copy()

    if "prot_class" in out.columns:
        non_empty = out["prot_class"].astype(str).str.strip().ne("")
        if non_empty.any():
            logging.info("Proteomics: using existing 'prot_class' column")
            return out

    def to_float(s: pd.Series) -> pd.Series:
        return pd.to_numeric(s, errors="coerce")

    present_any = to_float(out["prot_present_any"]) if "prot_present_any" in out.columns else None
    uniq = to_float(out["prot_unique_peptides_max"]) if "prot_unique_peptides_max" in out.columns else None
    cov = to_float(out["prot_coverage_pct_max"]) if "prot_coverage_pct_max" in out.columns else None

    if present_any is None and uniq is None and cov is None:
        logging.info("Proteomics: no prot_class or numeric evidence columns found; leaving as-is")
        return out

    prot_class: List[str] = []
    for i in range(out.shape[0]):
        pa = float(present_any.iloc[i]) if present_any is not None and not np.isnan(present_any.iloc[i]) else 0.0
        up = float(uniq.iloc[i]) if uniq is not None and not np.isnan(uniq.iloc[i]) else 0.0
        cv = float(cov.iloc[i]) if cov is not None and not np.isnan(cov.iloc[i]) else 0.0

        if up >= 2.0 and cv >= 10.0:
            prot_class.append("Strong")
        elif pa >= 1.0 or up > 0.0 or cv > 0.0:
            prot_class.append("Detected")
        else:
            prot_class.append("None")

    out["prot_class"] = prot_class
    logging.info("Proteomics: derived 'prot_class' for %s rows", out.shape[0])
    return out


def parse_bool_series(series: pd.Series) -> np.ndarray:
    """
    Parse a pandas Series into a boolean numpy array from common truthy tokens.

    Parameters
    ----------
    series
        Series to parse.

    Returns
    -------
    np.ndarray
        Boolean array.
    """
    vals = series.astype(str).str.strip().str.lower()
    return vals.isin({"true", "1", "t", "yes", "y"}).to_numpy()


def compute_priority_score(
    df: pd.DataFrame,
    weights: Weights,
    novelty_excluded_lists: Set[str],
    list_columns: List[str],
) -> pd.DataFrame:
    """
    Compute a transparent priority score and score components.

    Parameters
    ----------
    df
        DataFrame containing membership, optional proteomics, optional tractability.
    weights
        Scoring weights.
    novelty_excluded_lists
        Lists that, if present for a gene, remove novelty bonus.
    list_columns
        Boolean membership columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with added columns 'priority_score' and 'score_components'.
    """
    out = df.copy()

    components: List[List[str]] = [[] for _ in range(out.shape[0])]
    scores = np.zeros(out.shape[0], dtype=float)

    # List membership
    n_lists = pd.to_numeric(out.get("n_lists", 0), errors="coerce").fillna(0).astype(int).to_numpy()
    scores += n_lists.astype(float) * float(weights.per_list_member)
    for i, n in enumerate(n_lists.tolist()):
        if n > 0:
            components[i].append(f"lists:{n}*{weights.per_list_member:g}")

    # Proteomics
    if "prot_class" in out.columns:
        prot = out["prot_class"].fillna("None").astype(str).str.strip().str.title()
        for i, pc in enumerate(prot.tolist()):
            if pc == "Strong":
                scores[i] += float(weights.proteomics_strong)
                components[i].append(f"prot:Strong+{weights.proteomics_strong:g}")
            elif pc == "Detected":
                scores[i] += float(weights.proteomics_detected)
                components[i].append(f"prot:Detected+{weights.proteomics_detected:g}")

    # Tractability
    def bool_col(name: str) -> Optional[np.ndarray]:
        if name not in out.columns:
            return None
        return parse_bool_series(out[name])

    sm = bool_col("ot_any_small_molecule_tractable")
    ab = bool_col("ot_any_antibody_tractable")
    pr = bool_col("ot_any_protac_tractable")

    if sm is not None:
        idx = np.where(sm)[0]
        scores[idx] += float(weights.tract_small_molecule)
        for i in idx.tolist():
            components[i].append(f"tract:SM+{weights.tract_small_molecule:g}")

    if ab is not None:
        idx = np.where(ab)[0]
        scores[idx] += float(weights.tract_antibody)
        for i in idx.tolist():
            components[i].append(f"tract:Ab+{weights.tract_antibody:g}")

    if pr is not None:
        idx = np.where(pr)[0]
        scores[idx] += float(weights.tract_protac)
        for i in idx.tolist():
            components[i].append(f"tract:PROTAC+{weights.tract_protac:g}")

    # Novelty bonus: only if gene is not present in excluded lists
    novelty_mask = np.ones(out.shape[0], dtype=bool)
    if novelty_excluded_lists:
        for col in novelty_excluded_lists:
            if col in out.columns:
                # membership columns should be bool, but parse robustly anyway
                novelty_mask &= ~parse_bool_series(out[col])
        novelty_mask &= (n_lists > 0)

        scores[novelty_mask] += float(weights.novelty_bonus)
        for i in np.where(novelty_mask)[0].tolist():
            components[i].append(f"novelty+{weights.novelty_bonus:g}")

    out["priority_score"] = scores
    out["score_components"] = [";".join(c) for c in components]

    # Convenience booleans
    tract_any = np.zeros(out.shape[0], dtype=bool)
    for arr in [sm, ab, pr]:
        if arr is not None:
            tract_any |= arr
    out["ot_any_tractable"] = tract_any

    if "prot_class" in out.columns:
        out["prot_any_detected_or_strong"] = (
            out["prot_class"]
            .fillna("None")
            .astype(str)
            .str.strip()
            .str.title()
            .isin({"Detected", "Strong"})
        )
    else:
        out["prot_any_detected_or_strong"] = False

    out["candidate_druggable_sperm_protein"] = (
        out["ot_any_tractable"].astype(bool)
        & out["prot_any_detected_or_strong"].astype(bool)
    )

    return out


def ensure_out_prefix(path: str) -> str:
    """
    Convert an output path (possibly ending in .tsv/.xlsx) into a prefix.

    Parameters
    ----------
    path
        Output path or prefix.

    Returns
    -------
    str
        Prefix without extension.
    """
    if path.lower().endswith(".tsv"):
        return path[:-4]
    if path.lower().endswith(".xlsx"):
        return path[:-5]
    return path


def write_outputs(
    out_prefix: str,
    ranked: pd.DataFrame,
    candidate: pd.DataFrame,
    topn: Optional[pd.DataFrame],
) -> None:
    """
    Write TSV and Excel outputs.

    Parameters
    ----------
    out_prefix
        Output prefix (no extension).
    ranked
        Ranked full table.
    candidate
        Subset where candidate_druggable_sperm_protein is True.
    topn
        Optional top-N subset.
    """
    out_dir = os.path.dirname(os.path.abspath(out_prefix)) or "."
    os.makedirs(out_dir, exist_ok=True)

    ranked_tsv = f"{out_prefix}.tsv"
    ranked_xlsx = f"{out_prefix}.xlsx"

    ranked.to_csv(ranked_tsv, sep="\t", index=False)
    logging.info("Wrote ranked TSV: %s rows -> %s", ranked.shape[0], ranked_tsv)

    # Subset outputs
    cand_prefix = f"{out_prefix}__candidate_druggable_sperm_protein"
    cand_tsv = f"{cand_prefix}.tsv"
    cand_xlsx = f"{cand_prefix}.xlsx"

    candidate.to_csv(cand_tsv, sep="\t", index=False)
    logging.info("Wrote candidate TSV: %s rows -> %s", candidate.shape[0], cand_tsv)

    if topn is not None:
        top_prefix = f"{out_prefix}__top{topn.shape[0]}"
        top_tsv = f"{top_prefix}.tsv"
        top_xlsx = f"{top_prefix}.xlsx"
        topn.to_csv(top_tsv, sep="\t", index=False)
        logging.info("Wrote top-N TSV: %s rows -> %s", topn.shape[0], top_tsv)

    # Excel workbooks
    def write_excel(path: str, sheets: Dict[str, pd.DataFrame]) -> None:
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            for sheet_name, df in sheets.items():
                safe_name = sheet_name[:31]
                df.to_excel(writer, sheet_name=safe_name, index=False)

    write_excel(
        ranked_xlsx,
        sheets={
            "ranked": ranked,
            "candidate_druggable_sperm_protein": candidate,
            **({"topN": topn} if topn is not None else {}),
        },
    )
    logging.info("Wrote ranked Excel workbook -> %s", ranked_xlsx)

    write_excel(
        cand_xlsx,
        sheets={"candidate_druggable_sperm_protein": candidate},
    )
    logging.info("Wrote candidate Excel workbook -> %s", cand_xlsx)

    if topn is not None:
        top_prefix = f"{out_prefix}__top{topn.shape[0]}"
        top_xlsx = f"{top_prefix}.xlsx"
        write_excel(
            top_xlsx,
            sheets={"topN": topn},
        )
        logging.info("Wrote top-N Excel workbook -> %s", top_xlsx)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    p = argparse.ArgumentParser(description="Prioritise candidate druggable sperm targets.")
    p.add_argument("--gene_lists_dir", required=True, help="Directory of TSV gene lists.")
    p.add_argument("--gene_key_column", default="gene_key", help="Gene symbol column name in gene lists.")

    p.add_argument("--features_tsv", default="", help="Optional feature table TSV to carry through.")
    p.add_argument("--proteomics_tsv", default="", help="Optional proteomics TSV.")
    p.add_argument("--proteomics_gene_key_column", default="gene_key", help="Gene symbol column in proteomics TSV.")
    p.add_argument("--tractability_tsv", default="", help="Optional Open Targets tractability TSV.")
    p.add_argument("--tractability_gene_id_column", default="gene_id", help="Ensembl gene id column name in tractability TSV.")
    p.add_argument("--features_gene_id_column", default="gene_id", help="Ensembl gene id column name in features TSV.")
    p.add_argument(
        "--novelty_exclude_list_regex",
        default="",
        help="Regex; list(s) matching will be treated as 'known' and reduce novelty bonus (e.g. 'literature').",
    )

    # Weights
    p.add_argument("--w_per_list_member", type=float, default=1.0)
    p.add_argument("--w_prot_detected", type=float, default=1.0)
    p.add_argument("--w_prot_strong", type=float, default=3.0)
    p.add_argument("--w_tract_sm", type=float, default=4.0)
    p.add_argument("--w_tract_ab", type=float, default=2.0)
    p.add_argument("--w_tract_protac", type=float, default=2.0)
    p.add_argument("--w_novelty", type=float, default=1.5)

    p.add_argument("--min_lists", type=int, default=1, help="Filter: minimum number of lists a gene must appear in.")
    p.add_argument("--top_n", type=int, default=0, help="Optional: also write top N rows as separate outputs.")
    p.add_argument(
        "--out_prefix",
        required=True,
        help="Output prefix (or path ending .tsv/.xlsx). Writes both TSV and XLSX.",
    )
    p.add_argument("--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the prioritisation workflow."""
    args = parse_args(argv)
    setup_logger(verbose=args.verbose)

    start = datetime.now()

    if not os.path.isdir(args.gene_lists_dir):
        raise SystemExit(f"gene_lists_dir not found: {args.gene_lists_dir}")

    out_prefix = ensure_out_prefix(args.out_prefix)
    logging.info("Output prefix: %s", out_prefix)

    logging.info("Loading gene lists from: %s", args.gene_lists_dir)
    membership_df, _list_to_genes = load_gene_lists(
        gene_lists_dir=args.gene_lists_dir,
        gene_key_column=args.gene_key_column,
    )

    list_columns = [c for c in membership_df.columns if c not in {"gene_key", "n_lists", "list_names"}]
    logging.info("List columns: %s", ", ".join(list_columns))

    merged = membership_df.copy()
    merged["gene_key"] = normalise_gene_key(merged["gene_key"])

    # Optional: join to features for gene_id and context columns
    if args.features_tsv:
        logging.info("Loading features TSV: %s", args.features_tsv)
        feat = read_tsv(args.features_tsv)

        if "gene_key" not in feat.columns:
            if args.gene_key_column in feat.columns:
                feat = feat.rename(columns={args.gene_key_column: "gene_key"})
            else:
                raise SystemExit("features_tsv must contain a 'gene_key' column (or match --gene_key_column).")

        feat["gene_key"] = normalise_gene_key(feat["gene_key"])

        if feat["gene_key"].duplicated().any():
            dup_n = int(feat["gene_key"].duplicated().sum())
            logging.warning(
                "features_tsv has duplicated gene_key entries (%s); keeping first occurrence per gene_key",
                dup_n,
            )
            feat = feat.drop_duplicates(subset=["gene_key"], keep="first")

        merged = merged.merge(feat, on="gene_key", how="left")
        logging.info("After merging features: %s rows, %s columns", merged.shape[0], merged.shape[1])

        if args.features_gene_id_column in merged.columns and args.features_gene_id_column != "gene_id":
            merged = merged.rename(columns={args.features_gene_id_column: "gene_id"})
            logging.info("Standardised features gene id column to 'gene_id'")

        if "gene_id" in merged.columns:
            merged["gene_id"] = merged["gene_id"].astype(str).str.strip()
            missing_gene_id = int((merged["gene_id"] == "").sum())
            logging.info("gene_id present after features merge; missing gene_id for %s genes", missing_gene_id)
    else:
        logging.info("No features TSV provided; proceeding without genomic/regulatory context columns")

    # Optional: proteomics
    if args.proteomics_tsv:
        logging.info("Loading proteomics TSV: %s", args.proteomics_tsv)
        prot = read_tsv(args.proteomics_tsv)

        if args.proteomics_gene_key_column not in prot.columns:
            raise SystemExit(f"Proteomics TSV missing column: {args.proteomics_gene_key_column}")

        prot = prot.rename(columns={args.proteomics_gene_key_column: "gene_key"})
        prot["gene_key"] = normalise_gene_key(prot["gene_key"])
        prot = derive_proteomics_class(df=prot)

        if prot["gene_key"].duplicated().any():
            dup_n = int(prot["gene_key"].duplicated().sum())
            logging.warning(
                "proteomics_tsv has duplicated gene_key entries (%s); keeping best class per gene_key",
                dup_n,
            )
            rank = {"Strong": 2, "Detected": 1, "None": 0, "": -1}
            prot["_prot_rank"] = prot["prot_class"].map(rank).fillna(-1).astype(int)
            prot = (
                prot.sort_values("_prot_rank", ascending=False)
                .drop_duplicates("gene_key", keep="first")
                .drop(columns=["_prot_rank"])
            )

        keep_cols = ["gene_key", "prot_class"]
        for c in ["prot_present_any", "prot_unique_peptides_max", "prot_coverage_pct_max"]:
            if c in prot.columns:
                keep_cols.append(c)

        merged = merged.merge(prot[keep_cols], on="gene_key", how="left")

        # Ensure prot_class is present and does not contain NaN
        merged["prot_class"] = (
            merged.get("prot_class", "None")
            .fillna("None")
            .replace("", "None")
        )

        prot_counts = merged["prot_class"].value_counts(dropna=False).to_dict()
        logging.info("After merging proteomics: %s rows, %s columns", merged.shape[0], merged.shape[1])
        logging.info("Proteomics class counts (post-merge): %s", prot_counts)
    else:
        merged["prot_class"] = "None"
        logging.info("No proteomics TSV provided; prot_class set to 'None' for all genes")

    # Optional: tractability
    if args.tractability_tsv:
        logging.info("Loading tractability TSV: %s", args.tractability_tsv)
        tract = read_tsv(args.tractability_tsv)

        if args.tractability_gene_id_column not in tract.columns:
            raise SystemExit(f"Tractability TSV missing gene id column: {args.tractability_gene_id_column}")

        tract = tract.rename(columns={args.tractability_gene_id_column: "gene_id"})
        tract["gene_id"] = tract["gene_id"].astype(str).str.strip()

        tract_keep = ["gene_id"]
        for c in [
            "ot_any_small_molecule_tractable",
            "ot_any_antibody_tractable",
            "ot_any_protac_tractable",
            "ot_tractability_summary",
            "ot_approved_symbol",
        ]:
            if c in tract.columns:
                tract_keep.append(c)

        tract = tract[tract_keep].copy()

        if tract["gene_id"].duplicated().any():
            dup_n = int(tract["gene_id"].duplicated().sum())
            logging.warning("tractability_tsv has duplicated gene_id entries (%s); keeping first", dup_n)
            tract = tract.drop_duplicates(subset=["gene_id"], keep="first")

        if "gene_id" in merged.columns and merged["gene_id"].astype(str).str.strip().ne("").any():
            merged["gene_id"] = merged["gene_id"].astype(str).str.strip()
            merged = merged.merge(tract, on="gene_id", how="left")
            logging.info("After merging tractability by gene_id: %s rows, %s columns", merged.shape[0], merged.shape[1])

            tract_any = np.zeros(merged.shape[0], dtype=bool)
            for c in ["ot_any_small_molecule_tractable", "ot_any_antibody_tractable", "ot_any_protac_tractable"]:
                if c in merged.columns:
                    tract_any |= parse_bool_series(merged[c])
            logging.info("Genes with any tractability signal (post-merge): %s", int(tract_any.sum()))
        else:
            logging.warning(
                "No usable gene_id column available to merge tractability. "
                "Provide --features_tsv containing Ensembl gene_id to enable this step."
            )
            for c in [
                "ot_any_small_molecule_tractable",
                "ot_any_antibody_tractable",
                "ot_any_protac_tractable",
                "ot_tractability_summary",
                "ot_approved_symbol",
            ]:
                merged[c] = ""
    else:
        logging.info("No tractability TSV provided; tractability fields will be empty")
        for c in [
            "ot_any_small_molecule_tractable",
            "ot_any_antibody_tractable",
            "ot_any_protac_tractable",
            "ot_tractability_summary",
            "ot_approved_symbol",
        ]:
            if c not in merged.columns:
                merged[c] = ""

    # Novelty excluded lists
    novelty_excluded_lists: Set[str] = set()
    if args.novelty_exclude_list_regex:
        patt = re.compile(args.novelty_exclude_list_regex, flags=re.IGNORECASE)
        novelty_excluded_lists = {c for c in list_columns if patt.search(c)}
        logging.info(
            "Novelty exclude lists matched by regex '%s': %s",
            args.novelty_exclude_list_regex,
            ", ".join(sorted(novelty_excluded_lists)) if novelty_excluded_lists else "(none)",
        )
    else:
        logging.info("No novelty exclude regex provided; novelty bonus will not be applied")

    weights = Weights(
        per_list_member=args.w_per_list_member,
        proteomics_detected=args.w_prot_detected,
        proteomics_strong=args.w_prot_strong,
        tract_small_molecule=args.w_tract_sm,
        tract_antibody=args.w_tract_ab,
        tract_protac=args.w_tract_protac,
        novelty_bonus=args.w_novelty,
    )
    logging.info("Scoring weights: %s", weights)

    ranked = compute_priority_score(
        df=merged,
        weights=weights,
        novelty_excluded_lists=novelty_excluded_lists,
        list_columns=list_columns,
    )

    # Filter: minimum list membership
    ranked["n_lists"] = pd.to_numeric(ranked["n_lists"], errors="coerce").fillna(0).astype(int)
    before_filter = ranked.shape[0]
    ranked = ranked[ranked["n_lists"] >= int(args.min_lists)].copy()
    logging.info(
        "Filter min_lists >= %s: %s -> %s rows",
        args.min_lists,
        before_filter,
        ranked.shape[0],
    )

    # Sort: primary by priority_score, then tractability/proteomics/list count as tie-breakers
    ranked["prot_rank"] = ranked["prot_class"].fillna("None").map({"Strong": 2, "Detected": 1, "None": 0}).fillna(0).astype(int)
    ranked["tract_rank"] = ranked["ot_any_tractable"].astype(int)

    ranked = ranked.sort_values(
        by=["priority_score", "tract_rank", "prot_rank", "n_lists", "gene_key"],
        ascending=[False, False, False, False, True],
    ).drop(columns=["prot_rank", "tract_rank"], errors="ignore")

    # Column ordering
    front = [
        "gene_key",
        "gene_id",
        "priority_score",
        "score_components",
        "candidate_druggable_sperm_protein",
        "n_lists",
        "list_names",
        "prot_class",
        "prot_present_any",
        "prot_unique_peptides_max",
        "prot_coverage_pct_max",
        "ot_any_tractable",
        "ot_any_small_molecule_tractable",
        "ot_any_antibody_tractable",
        "ot_any_protac_tractable",
        "ot_approved_symbol",
        "ot_tractability_summary",
    ]
    cols = [c for c in front if c in ranked.columns] + [c for c in ranked.columns if c not in set(front)]
    ranked = ranked[cols]

    # Subsets
    candidate = ranked[ranked["candidate_druggable_sperm_protein"].astype(bool)].copy()
    logging.info(
        "Candidates with proteomics + any tractability: %s (of %s ranked)",
        candidate.shape[0],
        ranked.shape[0],
    )

    topn: Optional[pd.DataFrame] = None
    if int(args.top_n) > 0:
        topn = ranked.head(int(args.top_n)).copy()
        logging.info("Prepared top-N subset: %s rows", topn.shape[0])

    # Write outputs
    write_outputs(out_prefix=out_prefix, ranked=ranked, candidate=candidate, topn=topn)

    elapsed = datetime.now() - start
    logging.info("Done. Elapsed time: %s", str(elapsed).split(".")[0])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())