oxy).

You may also supply an explicit background TSV/TXT file (one gene per line),
or disable custom background to use g:Profiler defaults.
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
            logging.FileH#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_gprofiler_from_excel.py

Run g:Profiler enrichment (via gprofiler2 in R) for multiple gene sets stored as
columns in a single Excel workbook.

Input format
------------
An Excel file (.xlsx) where each column is a gene set and each row is a gene
symbol (e.g., ABCD1). Empty cells are allowed.

Outputs (TSV + PDF)
-------------------
For each gene set (column):
- <set_name>_gprofiler_results.tsv
- <set_name>_gprofiler_plot.pdf   (if results exist)

Additionally:
- all_sets_gprofiler_results.tsv  (combined results table)

Background choice
-----------------
By default, the script uses the union of all genes across all columns in the
Excel file as the custom background (a reasonable "study universe" prandler(filename=str(log_path), mode="w", encoding="utf-8"),
            logging.StreamHandler(),
        ],
    )


def sanitise_set_name(*, name: str) -> str:
    """
    Convert an arbitrary column header into a safe filename stem.

    Parameters
    ----------
    name
        Column header from the Excel file.

    Returns
    -------
    str
        Sanitised name suitable for filenames.
    """
    safe = name.strip()
    safe = re.sub(pattern=r"\s+", repl="_", string=safe)
    safe = re.sub(pattern=r"[^A-Za-z0-9_.-]+", repl="_", string=safe)
    safe = safe.strip("._-")
    return safe or "gene_set"


def write_text_lines(*, path: Path, lines: List[str]) -> None:
    """
    Write lines to a text file with newline separators.

    Parameters
    ----------
    path
        Output file path.
    lines
        Lines to write.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open(mode="w", encoding="utf-8") as handle:
        for line in lines:
            handle.write(f"{line}\n")


def build_r_script(*, out_dir: Path) -> Path:
    """
    Write an R script that reads the Excel and runs g:Profiler per column.

    Parameters
    ----------
    out_dir
        Output directory where the R script will be written.

    Returns
    -------
    Path
        Path to the written R script.
    """
    r_path = out_dir / "run_gprofiler_from_excel.R"

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
out_dir <- get_arg_value("--out_dir", NA)
organism <- get_arg_value("--organism", "hsapiens")
sheet_name <- get_arg_value("--sheet", NA)
use_custom_bg <- get_arg_value("--use_custom_bg", "true")
background_path <- get_arg_value("--background", NA)

if (is.na(excel_path) || is.na(out_dir)) {
  stop("Missing required args: --excel and --out_dir")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

if (is.na(sheet_name)) {
  input_df <- read_excel_first_sheet(excel_path)
} else {
  input_df <- read_excel_sheet(excel_path, sheet_name)
}

if (ncol(input_df) < 1) {
  stop("Excel file has no columns to process.")
}

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

all_genes <- c()
for (col_i in seq_len(ncol(input_df))) {
  all_genes <- c(all_genes, trim_gene_vector(input_df[[col_i]]))
}
all_genes <- unique(all_genes)

custom_bg <- NULL
use_custom_bg_flag <- tolower(use_custom_bg) %in% c("true", "t", "1", "yes", "y")

if (use_custom_bg_flag) {
  if (!is.na(background_path)) {
    bg_lines <- readLines(background_path, warn = FALSE)
    bg_lines <- trimws(bg_lines)
    bg_lines <- bg_lines[bg_lines != ""]
    custom_bg <- unique(bg_lines)
  } else {
    custom_bg <- all_genes
  }
}

combined_results <- data.frame()

for (col_i in seq_len(ncol(input_df))) {
  set_name_raw <- colnames(input_df)[col_i]
  set_name <- make_safe_name(set_name_raw)

  genes <- trim_gene_vector(input_df[[col_i]])

  message(sprintf("Processing set: %s (%d genes)", set_name_raw, length(genes)))

  if (length(genes) == 0) {
    message(sprintf("Skipping %s: no genes.", set_name_raw))
    next
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

  out_tsv <- file.path(out_dir, paste0(set_name, "_gprofiler_results.tsv"))
  out_pdf <- file.path(out_dir, paste0(set_name, "_gprofiler_plot.pdf"))

  if (is.null(res) || is.null(res$result) || nrow(res$result) == 0) {
    message(sprintf("No enrichment results for %s. Writing empty TSV.", set_name_raw))
    write.table(
      data.frame(),
      file = out_tsv,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    next
  }

  res_tbl <- res$result
  res_tbl$gene_set <- set_name_raw

  write.table(
    res_tbl,
    file = out_tsv,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  combined_results <- rbind(combined_results, res_tbl)

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
}

combined_out <- file.path(out_dir, "all_sets_gprofiler_results.tsv")
write.table(
  combined_results,
  file = combined_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Done.")
"""
    r_path.write_text(r_code, encoding="utf-8")
    os.chmod(r_path, 0o755)
    return r_path


def run_rscript(
    *,
    r_script_path: Path,
    excel_path: Path,
    out_dir: Path,
    organism: str,
    sheet: Optional[str],
    use_custom_bg: bool,
    background_path: Optional[Path],
) -> None:
    """
    Execute the R script with the provided arguments.

    Parameters
    ----------
    r_script_path
        Path to the R script.
    excel_path
        Input Excel workbook.
    out_dir
        Output directory.
    organism
        g:Profiler organism name (e.g., hsapiens).
    sheet
        Optional Excel sheet name. If None, first sheet is used.
    use_custom_bg
        Whether to use a custom background.
    background_path
        Optional path to a background list (one gene per line).
    """
    cmd = [
        "Rscript",
        str(r_script_path),
        "--excel",
        str(excel_path),
        "--out_dir",
        str(out_dir),
        "--organism",
        organism,
        "--use_custom_bg",
        "true" if use_custom_bg else "false",
    ]

    if sheet is not None:
        cmd.extend(["--sheet", sheet])

    if background_path is not None:
        cmd.extend(["--background", str(background_path)])

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
        description="Run g:Profiler enrichment per Excel column (gene set)."
    )
    parser.add_argument(
        "--excel",
        required=True,
        help="Path to input Excel (.xlsx) where each column is a gene set.",
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
        help="Optional sheet name. If omitted, first sheet is used.",
    )
    parser.add_argument(
        "--use_custom_bg",
        action="store_true",
        help=(
            "Use custom background. If --background is not provided, the "
            "union of all genes in the Excel workbook is used."
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
    return parser.parse_args()


def main() -> None:
    """
    Main entry point.
    """
    args = parse_args()

    excel_path = Path(args.excel).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()

    setup_logging(out_dir=out_dir)

    if not excel_path.exists():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")

    background_path = None
    if args.background is not None:
        background_path = Path(args.background).expanduser().resolve()
        if not background_path.exists():
            raise FileNotFoundError(f"Background file not found: {background_path}")

    r_script_path = build_r_script(out_dir=out_dir)

    try:
        run_rscript(
            r_script_path=r_script_path,
            excel_path=excel_path,
            out_dir=out_dir,
            organism=args.organism,
            sheet=args.sheet,
            use_custom_bg=bool(args.use_custom_bg),
            background_path=background_path,
        )
    except subprocess.CalledProcessError as exc:
        logging.error("Rscript failed with return code %s", exc.returncode)
        raise


if __name__ == "__main__":
    main()
