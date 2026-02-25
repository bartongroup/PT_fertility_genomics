#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
annotate_opentargets_tractability.py

Annotate genes with Open Targets tractability information via GraphQL.

Inputs
------
- genes TSV with a gene symbol column (default: gene_key)
- optional Ensembl gene id column (if present, preferred)

Outputs
-------
- TSV with tractability fields appended

Notes
-----
Open Targets provides tractability data for multiple modalities (small molecule,
antibody, PROTAC, other). :contentReference[oaicite:4]{index=4}

Example
-------
python annotate_opentargets_tractability.py \
  --in_tsv gene_context_features_universe.tsv \
  --gene_symbol_column gene_key \
  --ensembl_column gene_id \
  --out_tsv gene_context_features_universe_plus_tractability.tsv \
  --verbose
"""

from __future__ import annotations
import requests
import argparse
import json
import logging
import time
import urllib.request
from typing import Any, Dict, List, Optional, Sequence
from urllib.error import HTTPError, URLError

import pandas as pd


OT_API = "https://api.platform.opentargets.org/api/v4/graphql"


def setup_logger(verbose: bool) -> None:
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def gql_request(
    query: str,
    variables: Dict[str, Any],
    *,
    retries: int = 6,
    sleep_s: int = 5,
    timeout_s: int = 30,
) -> Dict[str, Any]:
    """
    Make a GraphQL request with retries (requests-based).

    This logs the HTTP response body on non-200 responses (e.g. 400), which is
    essential for debugging GraphQL schema/payload issues.
    """
    payload = {"query": query, "variables": variables}
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        "User-Agent": "PT_fertility_genomics/annotate_opentargets_tractability.py",
    }

    for attempt in range(1, retries + 1):
        try:
            r = requests.post(OT_API, json=payload, headers=headers, timeout=timeout_s)

            if r.status_code != 200:
                preview = (r.text or "")[:1000].replace("\n", " ").replace("\r", " ")
                raise RuntimeError(f"HTTP {r.status_code}: {preview}")

            out = r.json()
            if "errors" in out:
                raise RuntimeError(str(out["errors"])[:2000])

            return out

        except Exception as exc:
            logging.warning(
                "Open Targets request failed (%s/%s): %s",
                attempt,
                retries,
                str(exc)[:2000],
            )
            time.sleep(sleep_s)

    raise RuntimeError("Open Targets request failed after retries")


OT_QUERY = """
query TargetTractability($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    tractability {
      modality
      value
      label
    }
  }
}
"""


def summarise_tractability(tractability_list: Optional[List[Dict[str, Any]]]) -> Dict[str, Any]:
    """
    Convert tractability records into compact flags.

    Open Targets returns multiple records; we collapse into:
    - any_small_molecule_tractable
    - any_antibody_tractable
    - any_protac_tractable

    We also keep a compact text summary.
    """
    if not tractability_list:
        return {
            "ot_any_small_molecule_tractable": False,
            "ot_any_antibody_tractable": False,
            "ot_any_protac_tractable": False,
            "ot_tractability_summary": "",
        }

    sm = False
    ab = False
    pr = False
    bits: List[str] = []

    for rec in tractability_list:
        modality = str(rec.get("modality") or "").lower()
        label = str(rec.get("label") or "")
        value = rec.get("value")
        try:
            v = int(value)
        except Exception:
            v = 0

        if v == 1:
            if "small" in modality:
                sm = True
            if "antibody" in modality:
                ab = True
            if "protac" in modality:
                pr = True

        if label:
            bits.append(f"{rec.get('modality')}:{label}:{value}")

    return {
        "ot_any_small_molecule_tractable": sm,
        "ot_any_antibody_tractable": ab,
        "ot_any_protac_tractable": pr,
        "ot_tractability_summary": "|".join(bits)[:2000],
    }


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse CLI args."""
    p = argparse.ArgumentParser()
    p.add_argument("--in_tsv", required=True)
    p.add_argument("--gene_symbol_column", default="gene_key")
    p.add_argument("--ensembl_column", default="gene_id")
    p.add_argument("--out_tsv", required=True)
    p.add_argument("--max_genes", type=int, default=0)
    p.add_argument("--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run."""
    args = parse_args(argv)
    setup_logger(args.verbose)

    df = pd.read_csv(args.in_tsv, sep="\t")
    if args.ensembl_column not in df.columns:
        raise SystemExit(f"Missing Ensembl column: {args.ensembl_column}")

    ensembl_ids = df[args.ensembl_column].astype(str).dropna().drop_duplicates().tolist()
    logging.info("Unique Ensembl IDs to query: %s", len(ensembl_ids))
    
    if args.max_genes and args.max_genes > 0:
        ensembl_ids = ensembl_ids[: args.max_genes]
        logging.info("Limiting to max_genes: %s", len(ensembl_ids))

    out_rows: List[Dict[str, Any]] = []
    for i, ens in enumerate(ensembl_ids, start=1):
        ens_clean = ens.split(".")[0]
        if not ens_clean.startswith("ENSG"):
            out_rows.append(
                {
                    args.ensembl_column: ens,
                    "ot_any_small_molecule_tractable": False,
                    "ot_any_antibody_tractable": False,
                    "ot_any_protac_tractable": False,
                    "ot_tractability_summary": "",
                }
            )
            continue

        if i % 100 == 0:
            logging.info("Queried %s genes", i)

        resp = gql_request(OT_QUERY, {"ensemblId": ens_clean})
        target = (resp.get("data") or {}).get("target") or {}
        tract = target.get("tractability")
        summ = summarise_tractability(tract)

        summ[args.ensembl_column] = ens
        summ["ot_approved_symbol"] = target.get("approvedSymbol") or ""
        out_rows.append(summ)

    ann = pd.DataFrame(out_rows)
    merged = df.merge(ann, on=args.ensembl_column, how="left")

    merged.to_csv(args.out_tsv, sep="\t", index=False)
    logging.info("Wrote: %s", args.out_tsv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
