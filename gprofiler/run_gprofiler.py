#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_gprofiler.py

Run g:Profiler enrichment (via gprofiler2 in R) for gene sets provided in either:

1) Excel workbook (.xlsx): each column is a gene set (original behaviour)
2) Single gene-list TSV/TXT file: one gene per row (new behaviour)

This script writes an R helper script into the output directory and runs it with Rscript.

Input modes
-----------
A) Excel mode (original)
   --excel <file.xlsx>

   Each column is treated as a separate gene set. Empty cells are allowed.
   By default, a custom background can be:
     - union of all genes across all columns (if --use_custom_bg is set and
       --background is not provided), or
     - a user-supplied background file (one gene per line), or
     - g:Profiler default background (if --use_custom_bg is not set)

B) Gene-list TSV/TXT mode (new)
   --gene_list_tsv <file.tsv>

   The file should contain one gene per row.
   Accepted formats:
     - A headered TSV with the gene IDs in a named column (default: gene_key),
       set via --gene_column
     - A single-column file without a header (set --no_header)

   The gene set name used for output filenames can be controlled with:
     --single_set_name <name>
   If not provided, it defaults to the input filename stem.

Outputs
-------
Per gene set:
- <set_name>_gprofiler_results.tsv
- <set_name>_gprofiler_plot.pdf   (if results exist)
- <set_name>_gprofiler_results.xlsx

Additionally (Excel mode, and TSV mode too for consistency):
- all_sets_gprofiler_results.tsv
- all_sets_gprofiler_results.xlsx

Examples
--------
Excel mode:
python run_gprofiler.py \
  --excel sperm_gene_sets.xlsx \
  --out_dir gprofiler_out \
  --organism hsapiens \
  --use_custom_bg

Single TSV mode:
python run_gprofiler.py \
  --gene_list_tsv sperm_only_genes.tsv \
  --out_dir gprofiler_out/sperm_only_genes \
  --organism hsapiens \
  --gene_column gene_key \
  --single_set_name sperm_only_genes

Notes
-----
- Requires R packages: gprofiler2, openxlsx, ggplot2
- All outputs are tab-separated (TSV) where applicable.
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import subprocess
from pathlib import Path
from typing import List, Optional


def setup_logging(*, out_dir: Path) -> None:
    """
    Configure logging to file and console.

    Parameters
    ----------
    out_dir
        Output directory where the log file will be written.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "gprofiler_run.log"

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(filename=str(log_path), mode="w", encoding="utf-8"),
            logging.StreamHandler(),
        ],
    )


def sanitise_set_name(*, name: str) -> str:
    """
    Convert an arbitrary name into a safe filename stem.

    Parameters
    ----------
    name
        Gene set name.

    Returns
    -------
    str
        Sanitised name suitable for filenames.
    """
    safe = str(name).strip()
    safe = re.sub(pattern=r"\s+", repl="_", string=safe)
    safe = re.sub(pattern=r"[^A-Za-z0-9_.-]+", repl="_", string=safe)
    safe = safe.strip("._-")
    return safe or "gene_set"


def build_r_script(*, out_dir: Path) -> Path:
    """
    Write an R script that reads either an Excel workbook or a gene list TSV and
    runs g:Profiler per gene set.

    Parameters
    ----------
    out_dir
        Output directory where the R script will be written.

    Returns
    -------
    Path
        Path to the written R script.
    """
    r_path = out_dir / "run_gprofiler.R"

    r_code = r"""#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(gprofiler2))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)

get_arg_value <- function(flag, default_value = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) {
    return(default_value)
  }
  if (idx == length(args)) {
    return(default_value)
  }
  return(args[idx + 1])
}

excel_path <- get_arg_value("--excel", NA)
gene_list_tsv <- get_arg_value("--gene_list_tsv", NA)

out_dir <- get_arg_value("--out_dir", NA)
organism <- get_arg_value("--organism", "hsapiens")

sheet_name <- get_arg_value("--sheet", NA)
use_custom_bg <- get_arg_value("--use_custom_bg", "false")
background_path <- get_arg_value("--background", NA)

gene_column <- get_arg_value("--gene_column", "gene_key")
no_header <- get_arg_value("--no_header", "false")
single_set_name <- get_arg_value("--single_set_name", NA)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

use_custom_bg_flag <- tolower(use_custom_bg) %in% c("true", "t", "1", "yes", "y")
no_header_flag <- tolower(no_header) %in% c("true", "t", "1", "yes", "y")

trim_gene_vector <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & x != ""]
  x <- unique(x)
  return(x)
}

make_safe_name <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\\\s+", "_", x)
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("^[._-]+|[._-]+$", "", x)
  if (nchar(x) == 0) {
    x <- "gene_set"
  }
  return(x)
}

flatten_list_columns <- function(df) {
  for (nm in names(df)) {
    if (is.list(df[[nm]])) {
      df[[nm]] <- vapply(
        df[[nm]],
        function(x) {
          if (is.null(x) || length(x) == 0) {
            return("")
          }
          paste(as.character(x), collapse = ";")
        },
        FUN.VALUE = character(1)
      )
    }
  }
  return(df)
}

read_excel_first_sheet <- function(path) {
  sheet_names <- openxlsx::getSheetNames(path)
  if (length(sheet_names) < 1) {
    stop("No sheets found in Excel file.")
  }
  df <- openxlsx::read.xlsx(path, sheet = sheet_names[1], colNames = TRUE)
  return(df)
}

read_excel_sheet <- function(path, sheet) {
  df <- openxlsx::read.xlsx(path, sheet = sheet, colNames = TRUE)
  return(df)
}

read_gene_list_file <- function(path, gene_column, no_header_flag) {
  if (no_header_flag) {
    df <- read.table(path, header = FALSE, sep = "\t", quote = "", comment.char = "",
                     stringsAsFactors = FALSE, fill = TRUE)
    colnames(df) <- c("gene")
    genes <- trim_gene_vector(df[["gene"]])
    return(list(df = data.frame(gene_set = genes), set_name = NA))
  }

  df <- read.table(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
                   stringsAsFactors = FALSE, fill = TRUE)

  if (!(gene_column %in% colnames(df))) {
    stop(paste0("Gene column not found in gene list TSV: ", gene_column))
  }

  genes <- trim_gene_vector(df[[gene_column]])
  return(list(df = data.frame(gene_set = genes), set_name = NA))
}

custom_bg <- NULL
if (use_custom_bg_flag) {
  if (!is.na(background_path)) {
    bg_lines <- readLines(background_path, warn = FALSE)
    bg_lines <- trimws(bg_lines)
    bg_lines <- bg_lines[bg_lines != ""]
    custom_bg <- unique(bg_lines)
  } else {
    custom_bg <- NULL
  }
}

combined_results <- data.frame()

run_one_set <- function(genes, set_name_raw, out_dir, organism, custom_bg) {
  set_name_safe <- make_safe_name(set_name_raw)

  message(sprintf("Processing set: %s (%d genes)", set_name_raw, length(genes)))

  if (length(genes) == 0) {
    message(sprintf("Skipping %s: no genes.", set_name_raw))
    return(invisible(NULL))
  }

  res <- NULL
  if (!is.null(custom_bg) && length(custom_bg) > 0) {
    res <- tryCatch(
      gost(query = genes, organism = organism, custom_bg = custom_bg),
      error = function(e) {
        message(sprintf("gost() failed for %s: %s", set_name_raw, e$message))
        return(NULL)
      }
    )
  } else {
    res <- tryCatch(
      gost(query = genes, organism = organism),
      error = function(e) {
        message(sprintf("gost() failed for %s: %s", set_name_raw, e$message))
        return(NULL)
      }
    )
  }

  out_tsv <- file.path(out_dir, paste0(set_name_safe, "_gprofiler_results.tsv"))
  out_pdf <- file.path(out_dir, paste0(set_name_safe, "_gprofiler_plot.pdf"))
  out_xlsx <- file.path(out_dir, paste0(set_name_safe, "_gprofiler_results.xlsx"))

  if (is.null(res) || is.null(res$result) || nrow(res$result) == 0) {
    message(sprintf("No enrichment results for %s. Writing empty TSV.", set_name_raw))
    write.table(
      data.frame(),
      file = out_tsv,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    write.xlsx(
      data.frame(),
      file = out_xlsx,
      sheetName = set_name_raw,
      overwrite = TRUE
    )
    return(invisible(NULL))
  }

  res_tbl <- res$result
  res_tbl$gene_set <- set_name_raw
  res_tbl <- flatten_list_columns(res_tbl)

  write.table(
    res_tbl,
    file = out_tsv,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, set_name_raw)
  openxlsx::writeData(wb, set_name_raw, res_tbl)
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

  p <- tryCatch(
    gostplot(res, interactive = FALSE),
    error = function(e) {
      message(sprintf("gostplot() failed for %s: %s", set_name_raw, e$message))
      return(NULL)
    }
  )

  if (!is.null(p)) {
    ggsave(filename = out_pdf, plot = p, width = 11, height = 8.5)
  }

  return(res_tbl)
}

# Decide mode
excel_mode <- !is.na(excel_path)
tsv_mode <- !is.na(gene_list_tsv)

if (!excel_mode && !tsv_mode) {
  stop("Missing input. Provide either --excel or --gene_list_tsv")
}
if (excel_mode && tsv_mode) {
  stop("Ambiguous input. Provide only one of --excel or --gene_list_tsv")
}

if (excel_mode) {
  if (is.na(out_dir)) {
    stop("Missing required arg: --out_dir")
  }

  if (is.na(sheet_name)) {
    input_df <- read_excel_first_sheet(excel_path)
  } else {
    input_df <- read_excel_sheet(excel_path, sheet_name)
  }

  if (ncol(input_df) < 1) {
    stop("Excel file has no columns to process.")
  }

  # If requested and this is a single-column Excel, rename that column for cleaner outputs
  if (!is.na(single_set_name) && ncol(input_df) == 1) {
    colnames(input_df) <- single_set_name
  }

  # If using custom background but no explicit background file, use union of all genes in Excel
  if (use_custom_bg_flag && is.null(custom_bg)) {
    all_genes <- c()
    for (col_i in seq_len(ncol(input_df))) {
      all_genes <- c(all_genes, trim_gene_vector(input_df[[col_i]]))
    }
    custom_bg <- unique(all_genes)
  }

  for (col_i in seq_len(ncol(input_df))) {
    set_name_raw <- colnames(input_df)[col_i]
    genes <- trim_gene_vector(input_df[[col_i]])

    res_tbl <- run_one_set(genes, set_name_raw, out_dir, organism, custom_bg)
    if (!is.null(res_tbl) && nrow(res_tbl) > 0) {
      combined_results <- rbind(combined_results, res_tbl)
    }
  }
}

if (tsv_mode) {
  if (is.na(out_dir)) {
    stop("Missing required arg: --out_dir")
  }

  gl <- read_gene_list_file(gene_list_tsv, gene_column, no_header_flag)
  genes <- trim_gene_vector(gl$df[["gene_set"]])

  # Choose set name
  set_name_raw <- NA
  if (!is.na(single_set_name)) {
    set_name_raw <- single_set_name
  } else {
    # default: filename stem
    base <- basename(gene_list_tsv)
    set_name_raw <- sub("\\\\.[^.]*$", "", base)
  }

  # If using custom background but no explicit background file:
  # In TSV mode we cannot infer a sensible "study universe" across many sets,
  # so we leave custom_bg as NULL unless user provided --background.
  res_tbl <- run_one_set(genes, set_name_raw, out_dir, organism, custom_bg)
  if (!is.null(res_tbl) && nrow(res_tbl) > 0) {
    combined_results <- rbind(combined_results, res_tbl)
  }
}

combined_out <- file.path(out_dir, "all_sets_gprofiler_results.tsv")
combined_xlsx <- file.path(out_dir, "all_sets_gprofiler_results.xlsx")

write.table(
  combined_results,
  file = combined_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

wb_all <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_all, "combined_results")
openxlsx::writeData(wb_all, "combined_results", combined_results)
openxlsx::saveWorkbook(wb_all, combined_xlsx, overwrite = TRUE)

message("Done.")
"""
    r_path.write_text(r_code, encoding="utf-8")
    os.chmod(r_path, 0o755)
    return r_path


def run_rscript(
    *,
    r_script_path: Path,
    out_dir: Path,
    organism: str,
    excel_path: Optional[Path],
    gene_list_tsv: Optional[Path],
    sheet: Optional[str],
    use_custom_bg: bool,
    background_path: Optional[Path],
    gene_column: str,
    no_header: bool,
    single_set_name: Optional[str],
) -> None:
    """
    Execute the R script with the provided arguments.

    Parameters
    ----------
    r_script_path
        Path to the R script.
    out_dir
        Output directory.
    organism
        g:Profiler organism name (e.g., hsapiens).
    excel_path
        Input Excel workbook path (if using Excel mode).
    gene_list_tsv
        Input gene-list TSV/TXT path (if using TSV mode).
    sheet
        Optional Excel sheet name (Excel mode only).
    use_custom_bg
        Whether to use a custom background.
    background_path
        Optional path to a background list file (one gene per line).
    gene_column
        Column name containing genes (TSV mode; ignored if --no_header is used).
    no_header
        Treat TSV as no-header single-column file (TSV mode).
    single_set_name
        Override the gene set name (TSV mode) or rename single-column Excel header (Excel mode).
    """
    cmd: List[str] = [
        "Rscript",
        str(r_script_path),
        "--out_dir",
        str(out_dir),
        "--organism",
        str(organism),
        "--use_custom_bg",
        "true" if use_custom_bg else "false",
        "--gene_column",
        str(gene_column),
        "--no_header",
        "true" if no_header else "false",
    ]

    if excel_path is not None:
        cmd.extend(["--excel", str(excel_path)])
        if sheet is not None:
            cmd.extend(["--sheet", str(sheet)])

    if gene_list_tsv is not None:
        cmd.extend(["--gene_list_tsv", str(gene_list_tsv)])

    if background_path is not None:
        cmd.extend(["--background", str(background_path)])

    if single_set_name is not None and str(single_set_name).strip() != "":
        cmd.extend(["--single_set_name", str(single_set_name)])

    logging.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Run g:Profiler enrichment from Excel (multi-set) or TSV (single set)."
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--excel",
        default=None,
        help="Path to input Excel (.xlsx) where each column is a gene set.",
    )
    group.add_argument(
        "--gene_list_tsv",
        default=None,
        help="Path to a single gene-list TSV/TXT (one gene per row).",
    )

    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for TSVs and PDFs.",
    )
    parser.add_argument(
        "--organism",
        default="hsapiens",
        help="g:Profiler organism code (default: hsapiens).",
    )

    parser.add_argument(
        "--sheet",
        default=None,
        help="Optional sheet name (Excel mode). If omitted, first sheet is used.",
    )

    parser.add_argument(
        "--use_custom_bg",
        action="store_true",
        help=(
            "Use custom background. In Excel mode, if --background is not provided, the "
            "union of all genes in the Excel workbook is used. In TSV mode, a custom "
            "background is only used if --background is provided."
        ),
    )
    parser.add_argument(
        "--background",
        default=None,
        help=(
            "Optional background list file (TSV/TXT), one gene per line. "
            "Only used if --use_custom_bg is set."
        ),
    )

    parser.add_argument(
        "--gene_column",
        default="gene_key",
        help=(
            "Gene column name in --gene_list_tsv (default: gene_key). "
            "Ignored if --no_header is set."
        ),
    )
    parser.add_argument(
        "--no_header",
        action="store_true",
        help="Treat --gene_list_tsv as a single-column file without a header.",
    )

    parser.add_argument(
        "--single_set_name",
        default=None,
        help=(
            "Override gene set name. In TSV mode, sets output prefixes to this value. "
            "In Excel mode, if the workbook has one column, renames that column for cleaner outputs."
        ),
    )

    return parser.parse_args()


def main() -> None:
    """
    Main entry point.
    """
    args = parse_args()

    out_dir = Path(args.out_dir).expanduser().resolve()
    setup_logging(out_dir=out_dir)

    excel_path: Optional[Path] = None
    gene_list_tsv: Optional[Path] = None

    if args.excel is not None:
        excel_path = Path(args.excel).expanduser().resolve()
        if not excel_path.exists():
            raise FileNotFoundError(f"Excel file not found: {excel_path}")

    if args.gene_list_tsv is not None:
        gene_list_tsv = Path(args.gene_list_tsv).expanduser().resolve()
        if not gene_list_tsv.exists():
            raise FileNotFoundError(f"Gene list TSV not found: {gene_list_tsv}")

    background_path: Optional[Path] = None
    if args.background is not None:
        background_path = Path(args.background).expanduser().resolve()
        if not background_path.exists():
            raise FileNotFoundError(f"Background file not found: {background_path}")

    r_script_path = build_r_script(out_dir=out_dir)

    try:
        run_rscript(
            r_script_path=r_script_path,
            out_dir=out_dir,
            organism=str(args.organism),
            excel_path=excel_path,
            gene_list_tsv=gene_list_tsv,
            sheet=args.sheet,
            use_custom_bg=bool(args.use_custom_bg),
            background_path=background_path,
            gene_column=str(args.gene_column),
            no_header=bool(args.no_header),
            single_set_name=args.single_set_name,
        )
    except subprocess.CalledProcessError as exc:
        logging.error("Rscript failed with return code %s", exc.returncode)
        raise


if __name__ == "__main__":
    main()