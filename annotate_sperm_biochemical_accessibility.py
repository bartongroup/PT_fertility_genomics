#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
annotate_sperm_biochemical_accessibility.py

Annotate per-gene tables with biochemical accessibility features relevant to
drug targeting (especially for drugs delivered via the female).

This script is intentionally pragmatic and works in two modes:

1) Offline (recommended on clusters):
   - You provide a UniProt-derived annotation TSV (pre-downloaded).
   - Optionally also provide a mapping TSV from gene_key -> UniProt accession.

2) Online (optional):
   - Query UniProt REST API for missing accessions (best-effort).
   - This is disabled by default because institutional proxies can be fiddly.

Input expectations
------------------
- Input is a tab-separated table containing at least a gene identifier column.
  By default: gene_key

- If your input already contains UniProt accessions (e.g. uniprot_acc),
  you can point --uniprot_acc_column at that column and skip mapping.

Outputs
-------
- A tab-separated TSV with new columns added:
  - uniprot_acc
  - uniprot_subcellular_location
  - uniprot_keywords
  - uniprot_go_cc
  - has_signal_peptide
  - has_transmembrane
  - is_secreted
  - is_membrane
  - is_extracellular
  - is_plasma_membrane
  - is_cell_surface_candidate
  - predicted_target_class (heuristic)
  - biochemical_accessibility_score

Notes on interpretation
-----------------------
- "Cell surface candidate" is a strong starting point for drug accessibility.
  We define it heuristically as plasma membrane + extracellular/exposed.

- This script does not attempt structural modelling. It is an annotation layer
  to help triage targets before deeper follow-up.

All inputs/outputs are TSV (no comma-separated outputs). All arguments are named.
"""

from __future__ import annotations

import argparse
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
from typing import Dict
import pandas as pd
import pandas as pd


LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class UniprotRecord:
    """
    Container for UniProt-derived fields used in downstream annotations.
    """

    uniprot_acc: str
    subcellular_location: str
    keywords: str
    go_cc: str
    signal_peptide: str
    transmembrane: str
    protein_name: str


def _configure_logging(verbose: bool) -> None:
    """
    Configure logging.

    Parameters
    ----------
    verbose : bool
        If True, log at DEBUG level; else INFO.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def _read_genes_master_from_excel(
    excel_path: Path,
    sheet_name: str = "Genes_Master",
) -> pd.DataFrame:
    """
    Read the Genes_Master sheet from Excel.

    Parameters
    ----------
    excel_path : Path
        Path to Excel workbook.
    sheet_name : str
        Sheet name to read.

    Returns
    -------
    pd.DataFrame
        Genes_Master table.
    """
    if not excel_path.exists():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")

    df = pd.read_excel(excel_path, sheet_name=sheet_name, dtype=str)
    return df


def _write_updated_excel(
    original_excel: Path,
    updated_genes_master: pd.DataFrame,
    out_excel: Path,
    sheet_name: str = "Genes_Master",
) -> None:
    """
    Copy Excel workbook and replace Genes_Master sheet.

    Parameters
    ----------
    original_excel : Path
        Source workbook.
    updated_genes_master : pd.DataFrame
        Annotated Genes_Master.
    out_excel : Path
        Output workbook.
    sheet_name : str
        Sheet to replace.
    """
    xls = pd.ExcelFile(original_excel)

    with pd.ExcelWriter(out_excel, engine="openpyxl") as writer:
        for sheet in xls.sheet_names:
            if sheet == sheet_name:
                updated_genes_master.to_excel(writer, sheet_name=sheet, index=False)
            else:
                df = pd.read_excel(original_excel, sheet_name=sheet, dtype=str)
                df.to_excel(writer, sheet_name=sheet, index=False)

def _read_tsv(path: Path) -> pd.DataFrame:
    """
    Read a TSV safely.

    Parameters
    ----------
    path : Path
        TSV path.

    Returns
    -------
    pd.DataFrame
        DataFrame with string dtype where possible.
    """
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)


def _write_tsv(df: pd.DataFrame, path: Path) -> None:
    """
    Write a TSV safely.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to write.
    path : Path
        Output path.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def _normalise_gene_key(value: object) -> str:
    """
    Normalise gene keys for joins.

    Parameters
    ----------
    value : object
        Any input value.

    Returns
    -------
    str
        Uppercase identifier with whitespace removed.
    """
    if value is None:
        return ""
    s = str(value).strip()
    s = re.sub(r"\s+", "", s)
    return s.upper()


def _normalise_uniprot_acc(value: object) -> str:
    """
    Normalise a UniProt accession.

    Parameters
    ----------
    value : object
        Raw accession.

    Returns
    -------
    str
        Trimmed accession string.
    """
    if value is None:
        return ""
    return str(value).strip()


def _as_bool_from_feature_text(value: object) -> bool:
    """
    Convert common UniProt feature texts into a boolean.

    Parameters
    ----------
    value : object
        Feature text (e.g. 'Signal peptide', 'Transmembrane', 'Yes', etc.)

    Returns
    -------
    bool
        True if value suggests feature presence.
    """
    if value is None:
        return False
    s = str(value).strip().lower()
    if s == "" or s in {"na", "none", "no", "false", "0"}:
        return False
    # UniProt TSV exports can contain strings like:
    # - "Signal peptide"
    # - "TRANSMEM"
    # - "Yes"
    return True


def _tokenise_text(value: object) -> str:
    """
    Normalise free-text annotation fields for simple substring matching.

    Parameters
    ----------
    value : object
        Any text.

    Returns
    -------
    str
        Lowercase, space-normalised string.
    """
    if value is None:
        return ""
    s = str(value).strip().lower()
    s = re.sub(r"\s+", " ", s)
    return s


def _contains_any(text: str, needles: Iterable[str]) -> bool:
    """
    Check whether any substring appears in the text.

    Parameters
    ----------
    text : str
        Haystack.
    needles : Iterable[str]
        Substrings to check.

    Returns
    -------
    bool
        True if any substring is present.
    """
    for n in needles:
        if n in text:
            return True
    return False


def _infer_target_class(
    protein_name: str,
    keywords: str,
    go_cc: str,
) -> str:
    """
    Heuristically infer a broad target class from UniProt text fields.

    Parameters
    ----------
    protein_name : str
        Protein name/description.
    keywords : str
        UniProt keywords.
    go_cc : str
        GO cellular component annotations.

    Returns
    -------
    str
        One of: GPCR, ion_channel, transporter, enzyme, receptor, secreted_factor,
        structural, unknown
    """
    pn = _tokenise_text(protein_name)
    kw = _tokenise_text(keywords)
    go = _tokenise_text(go_cc)
    text = " ".join([pn, kw, go]).strip()

    if _contains_any(text, ["g-protein coupled receptor", "gpcr"]):
        return "GPCR"
    if _contains_any(text, ["ion channel", "channel activity", "voltage-gated", "ligand-gated"]):
        return "ion_channel"
    if _contains_any(text, ["transporter", "solute carrier", "slc", "abc transporter", "pump", "atpase"]):
        return "transporter"
    if _contains_any(text, ["kinase", "phosphatase", "protease", "peptidase", "hydrolase", "transferase", "enzyme"]):
        return "enzyme"
    if _contains_any(text, ["receptor", "rtk", "tyrosine-protein kinase receptor"]):
        return "receptor"
    if _contains_any(text, ["secreted", "extracellular", "hormone", "cytokine"]):
        return "secreted_factor"
    if _contains_any(text, ["structural", "cytoskeleton", "axoneme", "flagellar", "tubulin", "dynein"]):
        return "structural"

    return "unknown"


def _build_accessibility_score(
    is_cell_surface_candidate: bool,
    is_secreted: bool,
    is_membrane: bool,
    target_class: str,
) -> int:
    """
    Build a simple integer score for biochemical accessibility.

    Parameters
    ----------
    is_cell_surface_candidate : bool
        Whether target is likely cell-surface exposed.
    is_secreted : bool
        Whether target is secreted/extracellular.
    is_membrane : bool
        Whether target is membrane-associated.
    target_class : str
        Broad inferred class.

    Returns
    -------
    int
        Integer score (higher is better for accessibility/druggability triage).
    """
    score = 0

    if is_cell_surface_candidate:
        score += 5
    elif is_secreted:
        score += 4
    elif is_membrane:
        score += 2

    if target_class in {"GPCR", "ion_channel", "transporter", "receptor"}:
        score += 3
    elif target_class == "enzyme":
        score += 2

    return score


def load_gene_to_uniprot_mapping(mapping_tsv: Path) -> Dict[str, str]:
    """
    Load a mapping TSV of gene_key -> uniprot_acc.

    Expected columns (case-insensitive):
    - gene_key (or gene / symbol)
    - uniprot_acc (or uniprot / entry)

    Parameters
    ----------
    mapping_tsv : Path
        Mapping TSV.

    Returns
    -------
    Dict[str, str]
        Mapping from normalised gene_key to UniProt accession.
    """
    df = _read_tsv(mapping_tsv)
    cols = {c.lower(): c for c in df.columns}

    gene_col = None
    for cand in ["gene_key", "gene", "symbol", "hgnc_symbol"]:
        if cand in cols:
            gene_col = cols[cand]
            break
    if gene_col is None:
        raise ValueError(
            f"Could not find gene column in mapping TSV: {mapping_tsv}. "
            f"Columns: {list(df.columns)}"
        )

    up_col = None
    for cand in ["uniprot_acc", "uniprot", "entry", "uniprotkb", "uniprot_id"]:
        if cand in cols:
            up_col = cols[cand]
            break
    if up_col is None:
        raise ValueError(
            f"Could not find UniProt accession column in mapping TSV: {mapping_tsv}. "
            f"Columns: {list(df.columns)}"
        )

    out: Dict[str, str] = {}
    for _, row in df.iterrows():
        g = _normalise_gene_key(row.get(gene_col))
        u = _normalise_uniprot_acc(row.get(up_col))
        if g and u:
            out[g] = u

    LOGGER.info("Loaded %s gene->UniProt mappings from %s", len(out), mapping_tsv)
    return out


def load_uniprot_annotation_table(uniprot_tsv: Path) -> Dict[str, UniprotRecord]:
    """
    Load a UniProt-derived TSV annotation table keyed by UniProt accession.

    This expects a TSV you have exported from UniProt (or elsewhere) that
    contains at least:
    - Entry (UniProt accession)
    And ideally:
    - Protein names
    - Subcellular location [CC]
    - Keywords
    - Gene Ontology (cellular component)
    - Signal peptide
    - Transmembrane

    The script is tolerant of missing optional fields; it will fill blanks.

    Parameters
    ----------
    uniprot_tsv : Path
        UniProt annotation TSV.

    Returns
    -------
    Dict[str, UniprotRecord]
        Mapping from UniProt accession to UniprotRecord.
    """
    df = _read_tsv(uniprot_tsv)
    cols = {c.lower(): c for c in df.columns}

    if "entry" not in cols:
        raise ValueError(
            f"UniProt TSV missing required 'Entry' column: {uniprot_tsv}. "
            f"Columns: {list(df.columns)}"
        )

    def _get_col(*candidates: str) -> Optional[str]:
        for cand in candidates:
            if cand.lower() in cols:
                return cols[cand.lower()]
        return None

    c_entry = cols["entry"]
    c_prot = _get_col("Protein names", "protein_name", "protein names")
    c_subc = _get_col("Subcellular location [CC]", "subcellular location [cc]", "subcellular_location")
    c_kw = _get_col("Keywords", "keyword")
    c_go = _get_col("Gene Ontology (cellular component)", "go [cc]", "go_cc", "gene ontology (cellular component)")
    c_sig = _get_col("Signal peptide", "signal_peptide")
    c_tm = _get_col("Transmembrane", "transmembrane")

    records: Dict[str, UniprotRecord] = {}
    for _, row in df.iterrows():
        acc = _normalise_uniprot_acc(row.get(c_entry))
        if not acc:
            continue

        rec = UniprotRecord(
            uniprot_acc=acc,
            subcellular_location=str(row.get(c_subc, "")).strip() if c_subc else "",
            keywords=str(row.get(c_kw, "")).strip() if c_kw else "",
            go_cc=str(row.get(c_go, "")).strip() if c_go else "",
            signal_peptide=str(row.get(c_sig, "")).strip() if c_sig else "",
            transmembrane=str(row.get(c_tm, "")).strip() if c_tm else "",
            protein_name=str(row.get(c_prot, "")).strip() if c_prot else "",
        )
        records[acc] = rec

    LOGGER.info("Loaded %s UniProt records from %s", len(records), uniprot_tsv)
    return records


def annotate_accessibility(
    df: pd.DataFrame,
    gene_col: str,
    uniprot_acc_col: str,
    gene_to_uniprot: Optional[Dict[str, str]],
    uniprot_records: Optional[Dict[str, UniprotRecord]],
) -> pd.DataFrame:
    """
    Annotate accessibility fields onto the input table.

    Parameters
    ----------
    df : pd.DataFrame
        Input table.
    gene_col : str
        Column containing gene identifiers.
    uniprot_acc_col : str
        Column name to use/store UniProt accessions.
    gene_to_uniprot : Optional[Dict[str, str]]
        Mapping from gene -> UniProt accession. If provided, used to fill missing accessions.
    uniprot_records : Optional[Dict[str, UniprotRecord]]
        UniProt annotation mapping by accession. If provided, used to add annotation fields.

    Returns
    -------
    pd.DataFrame
        Copy of df with added annotation columns.
    """
    out = df.copy()

    if gene_col not in out.columns:
        raise ValueError(
            f"Input table missing gene column '{gene_col}'. Columns: {list(out.columns)}"
        )

    out[gene_col] = out[gene_col].map(_normalise_gene_key)

    if uniprot_acc_col not in out.columns:
        out[uniprot_acc_col] = ""

    out[uniprot_acc_col] = out[uniprot_acc_col].fillna("").map(_normalise_uniprot_acc)

    if gene_to_uniprot is not None:
        missing_mask = out[uniprot_acc_col].astype(str).str.strip() == ""
        n_missing = int(missing_mask.sum())
        if n_missing > 0:
            LOGGER.info("Filling %s missing UniProt accessions using mapping", n_missing)
            out.loc[missing_mask, uniprot_acc_col] = (
                out.loc[missing_mask, gene_col].map(gene_to_uniprot).fillna("")
            )

    # Default empty annotation columns
    out["uniprot_subcellular_location"] = ""
    out["uniprot_keywords"] = ""
    out["uniprot_go_cc"] = ""
    out["uniprot_signal_peptide_raw"] = ""
    out["uniprot_transmembrane_raw"] = ""
    out["uniprot_protein_name"] = ""

    if uniprot_records is not None:
        LOGGER.info("Annotating UniProt fields from offline UniProt table")
        recs = out[uniprot_acc_col].map(uniprot_records)
        out["uniprot_subcellular_location"] = recs.map(
            lambda r: r.subcellular_location if r is not None else ""
        )
        out["uniprot_keywords"] = recs.map(lambda r: r.keywords if r is not None else "")
        out["uniprot_go_cc"] = recs.map(lambda r: r.go_cc if r is not None else "")
        out["uniprot_signal_peptide_raw"] = recs.map(
            lambda r: r.signal_peptide if r is not None else ""
        )
        out["uniprot_transmembrane_raw"] = recs.map(
            lambda r: r.transmembrane if r is not None else ""
        )
        out["uniprot_protein_name"] = recs.map(
            lambda r: r.protein_name if r is not None else ""
        )

    # Feature flags
    out["has_signal_peptide"] = out["uniprot_signal_peptide_raw"].map(_as_bool_from_feature_text)
    out["has_transmembrane"] = out["uniprot_transmembrane_raw"].map(_as_bool_from_feature_text)

    subc = out["uniprot_subcellular_location"].map(_tokenise_text)
    kw = out["uniprot_keywords"].map(_tokenise_text)
    go = out["uniprot_go_cc"].map(_tokenise_text)

    out["is_extracellular"] = subc.map(
        lambda s: _contains_any(s, ["extracellular", "secreted"])
    ) | go.map(
        lambda s: _contains_any(s, ["extracellular", "external side of plasma membrane"])
    )

    out["is_plasma_membrane"] = subc.map(
        lambda s: _contains_any(s, ["plasma membrane", "cell membrane"])
    ) | go.map(
        lambda s: _contains_any(s, ["plasma membrane", "cell surface", "external side of plasma membrane"])
    ) | kw.map(
        lambda s: _contains_any(s, ["membrane"])
    )

    out["is_membrane"] = out["has_transmembrane"] | out["is_plasma_membrane"] | kw.map(
        lambda s: _contains_any(s, ["membrane", "integral component of membrane"])
    )

    out["is_secreted"] = (
        out["has_signal_peptide"]
        & ~out["has_transmembrane"]
        & out["is_extracellular"]
    )

    # Stronger "surface" heuristic
    out["is_cell_surface_candidate"] = out["is_plasma_membrane"] & out["is_extracellular"]

    # Target class and score
    out["predicted_target_class"] = [
        _infer_target_class(
            protein_name=pn,
            keywords=kwd,
            go_cc=goc,
        )
        for pn, kwd, goc in zip(
            out["uniprot_protein_name"].fillna(""),
            out["uniprot_keywords"].fillna(""),
            out["uniprot_go_cc"].fillna(""),
            strict=False,
        )
    ]

    out["biochemical_accessibility_score"] = [
        _build_accessibility_score(
            is_cell_surface_candidate=bool(surf),
            is_secreted=bool(sec),
            is_membrane=bool(mem),
            target_class=str(tc),
        )
        for surf, sec, mem, tc in zip(
            out["is_cell_surface_candidate"].fillna(False),
            out["is_secreted"].fillna(False),
            out["is_membrane"].fillna(False),
            out["predicted_target_class"].fillna("unknown"),
            strict=False,
        )
    ]

    return out


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    p = argparse.ArgumentParser(
        description="Annotate biochemical accessibility features for sperm target prioritisation."
    )
    p.add_argument(
        "--in_tsv",
        required=True,
        type=Path,
        help="Input TSV with at least a gene identifier column (default: gene_key).",
    )
    p.add_argument(
        "--out_tsv",
        required=True,
        type=Path,
        help="Output TSV with accessibility annotations added.",
    )
    p.add_argument(
        "--gene_column",
        required=False,
        type=str,
        default="gene_key",
        help="Gene identifier column name (default: gene_key).",
    )
    p.add_argument(
        "--uniprot_acc_column",
        required=False,
        type=str,
        default="uniprot_acc",
        help="Column name to read/write UniProt accessions (default: uniprot_acc).",
    )
    p.add_argument(
        "--gene_to_uniprot_tsv",
        required=False,
        type=Path,
        default=None,
        help=(
            "Optional TSV mapping gene_key -> UniProt accession. "
            "Used to fill missing accessions."
        ),
    )
    p.add_argument(
        "--uniprot_annotation_tsv",
        required=False,
        type=Path,
        default=None,
        help=(
            "Optional UniProt annotation TSV exported offline, keyed by 'Entry'. "
            "If provided, used to add subcellular location, keywords, GO CC, "
            "signal peptide, transmembrane, and protein name."
        ),
    )
    p.add_argument(
        "--excel_in",
        required=True,
        type=Path,
        help="Input Excel workbook (e.g. SUMMARY_fertility_evidence.xlsx)",
    )

    p.add_argument(
        "--excel_out",
        required=True,
        type=Path,
        help="Output Excel workbook with annotated Genes_Master sheet.",
    )

    p.add_argument(
        "--sheet_name",
        required=False,
        default="Genes_Master",
        type=str,
        help="Sheet to annotate (default: Genes_Master).",
    )
    parser.add_argument(
        "--testis_annotated_override_tsv",
        required=False,
        type=Path,
        default=None,
        help=(
            "Optional override TSV for the testis annotated table used to seed Genes_Master. "
            "If set, this file is used instead of results/<testis_run_id>/08_proteomics_annotation/..."
        ),
    )
    p.add_argument(
        "--verbose",
        action="store_true",
        help="Verbose logging.",
    )
    return p.parse_args()



def main() -> None:
    """
    Entry point.
    """
    args = parse_args()
    _configure_logging(verbose=args.verbose)

    LOGGER.info("Loading Genes_Master from Excel")
    df = _read_genes_master_from_excel(
        excel_path=args.excel_in,
        sheet_name=args.sheet_name,
    )

    LOGGER.info("Loaded %s rows", df.shape[0])

    gene_to_uniprot = None
    if args.gene_to_uniprot_tsv is not None:
        if args.gene_to_uniprot_tsv.exists() and args.gene_to_uniprot_tsv.stat().st_size > 0:
            gene_to_uniprot = load_gene_to_uniprot_mapping(args.gene_to_uniprot_tsv)

    uniprot_records = None
    if args.uniprot_annotation_tsv is not None:
        if args.uniprot_annotation_tsv.exists() and args.uniprot_annotation_tsv.stat().st_size > 0:
            uniprot_records = load_uniprot_annotation_table(args.uniprot_annotation_tsv)

    LOGGER.info("Annotating biochemical accessibility")

    annotated = annotate_accessibility(
        df=df,
        gene_col=args.gene_column,
        uniprot_acc_col=args.uniprot_acc_column,
        gene_to_uniprot=gene_to_uniprot,
        uniprot_records=uniprot_records,
    )

    LOGGER.info("Writing updated Excel workbook")

    _write_updated_excel(
        original_excel=args.excel_in,
        updated_genes_master=annotated,
        out_excel=args.excel_out,
        sheet_name=args.sheet_name,
    )

    LOGGER.info("Done: %s", args.excel_out)



if __name__ == "__main__":
    main()
