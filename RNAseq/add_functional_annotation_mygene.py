#!/usr/bin/env python3
"""
Add functional annotation to a TSV by querying MyGene.info for HGNC symbols.

This script reads an input TSV (for example your GTEx+ sperm annotated table),
uses a specified symbol column (default: hgnc_symbol), queries MyGene.info in
batches, and writes an output TSV with added functional annotation columns.

Outputs are tab-separated (TSV).

Notes
-----
- Requires internet access to MyGene.info.
- Uses conservative fields: symbol, name, summary, type_of_gene, go terms.
- Missing symbols are handled gracefully and left blank.

References
----------
MyGene.info documentation:
- https://docs.mygene.info/
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import requests


@dataclass(frozen=True)
class Config:
    """
    Configuration for MyGene.info functional annotation.

    Attributes
    ----------
    input_tsv
        Input TSV path.
    output_tsv
        Output TSV path.
    symbol_col
        Column containing HGNC symbols.
    taxid
        NCBI taxonomy ID (default 9606 for human).
    batch_size
        Number of symbols per API request.
    sleep_seconds
        Delay between API requests to be polite.
    timeout_seconds
        Request timeout.
    max_retries
        Retries per batch on transient failures.
    log_level
        Logging verbosity.
    log_path
        Optional log file.
    cache_jsonl
        Optional JSONL cache path of API results (one JSON per gene hit).
        If provided and file exists, cached entries will be reused.
    """

    input_tsv: str
    output_tsv: str
    symbol_col: str
    taxid: int
    batch_size: int
    sleep_seconds: float
    timeout_seconds: int
    max_retries: int
    log_level: str
    log_path: Optional[str]
    cache_jsonl: Optional[str]


def setup_logging(*, log_level: str, log_path: Optional[str]) -> logging.Logger:
    """
    Configure logging.

    Parameters
    ----------
    log_level
        Logging level (INFO, DEBUG).
    log_path
        Optional log file path.

    Returns
    -------
    logging.Logger
        Logger.
    """
    logger = logging.getLogger("add_functional_annotation_mygene")
    logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    logger.handlers = []
    logger.propagate = False

    fmt = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    sh = logging.StreamHandler(stream=sys.stdout)
    sh.setFormatter(fmt)
    sh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    logger.addHandler(sh)

    if log_path:
        fh = logging.FileHandler(log_path)
        fh.setFormatter(fmt)
        fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
        logger.addHandler(fh)

    return logger


def parse_args() -> Config:
    """
    Parse CLI arguments.

    Returns
    -------
    Config
        Parsed config.
    """
    p = argparse.ArgumentParser(description="Add functional annotation using MyGene.info.")
    p.add_argument("--input_tsv", required=True, help="Input TSV.")
    p.add_argument("--output_tsv", required=True, help="Output TSV.")
    p.add_argument("--symbol_col", default="hgnc_symbol", help="HGNC symbol column (default: hgnc_symbol).")
    p.add_argument("--taxid", type=int, default=9606, help="NCBI taxonomy ID (default: 9606).")
    p.add_argument("--batch_size", type=int, default=500, help="Symbols per API request (default: 500).")
    p.add_argument("--sleep_seconds", type=float, default=0.2, help="Sleep between requests (default: 0.2).")
    p.add_argument("--timeout_seconds", type=int, default=60, help="Request timeout seconds (default: 60).")
    p.add_argument("--max_retries", type=int, default=5, help="Max retries per batch (default: 5).")
    p.add_argument("--cache_jsonl", default=None, help="Optional JSONL cache path.")
    p.add_argument("--log_level", default="INFO", help="Log level.")
    p.add_argument("--log_path", default=None, help="Optional log file path.")
    a = p.parse_args()

    return Config(
        input_tsv=a.input_tsv,
        output_tsv=a.output_tsv,
        symbol_col=a.symbol_col,
        taxid=a.taxid,
        batch_size=a.batch_size,
        sleep_seconds=a.sleep_seconds,
        timeout_seconds=a.timeout_seconds,
        max_retries=a.max_retries,
        log_level=a.log_level,
        log_path=a.log_path,
        cache_jsonl=a.cache_jsonl,
    )


def load_cache(*, cache_jsonl: str, logger: logging.Logger) -> Dict[str, Dict[str, Any]]:
    """
    Load cached MyGene hits from a JSONL file.

    Parameters
    ----------
    cache_jsonl
        Path to JSONL cache.
    logger
        Logger.

    Returns
    -------
    dict
        Mapping from symbol (upper) to a MyGene hit dict.
    """
    cache: Dict[str, Dict[str, Any]] = {}
    try:
        with open(cache_jsonl, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                obj = json.loads(line)
                sym = str(obj.get("symbol", "")).strip().upper()
                if sym:
                    cache[sym] = obj
        logger.info("Loaded cache entries: %d", len(cache))
    except FileNotFoundError:
        logger.info("Cache file not found, starting without cache: %s", cache_jsonl)
    return cache


def append_cache(*, cache_jsonl: str, hits: List[Dict[str, Any]]) -> None:
    """
    Append MyGene hits to a JSONL cache file.

    Parameters
    ----------
    cache_jsonl
        Cache file path.
    hits
        List of hit dicts.
    """
    with open(cache_jsonl, "a", encoding="utf-8") as fh:
        for h in hits:
            fh.write(json.dumps(h, ensure_ascii=False) + "\n")


def chunk_list(*, items: List[str], chunk_size: int) -> List[List[str]]:
    """
    Split a list into chunks.

    Parameters
    ----------
    items
        Items to chunk.
    chunk_size
        Chunk size.

    Returns
    -------
    list
        List of chunks.
    """
    return [items[i : i + chunk_size] for i in range(0, len(items), chunk_size)]


def mygene_query_batch(
    *,
    symbols: List[str],
    taxid: int,
    timeout_seconds: int,
) -> List[Dict[str, Any]]:
    """
    Query MyGene.info for a batch of HGNC symbols via POST to /v3/query.

    This follows the documented "Batch queries via POST" behaviour where the
    response is a JSON list of per-query results. :contentReference[oaicite:2]{index=2}

    Parameters
    ----------
    symbols
        HGNC symbols.
    taxid
        Taxonomy ID (9606 for human).
    timeout_seconds
        HTTP timeout.

    Returns
    -------
    list
        List of hit dictionaries (excluding 'notfound' entries).
    """
    url = "https://mygene.info/v3/query"
    fields = "symbol,name,summary,type_of_gene,genomic_pos,go"

    # Batch POST expects q terms separated by comma/space/+ and uses scopes for the ID type. :contentReference[oaicite:3]{index=3}
    params = {
        "q": ",".join(symbols),
        "scopes": "symbol",
        "species": str(taxid),
        "fields": fields,
    }

    headers = {"content-type": "application/x-www-form-urlencoded"}
    resp = requests.post(url, data=params, headers=headers, timeout=timeout_seconds)
    resp.raise_for_status()

    payload = resp.json()
    if not isinstance(payload, list):
        raise ValueError(f"Unexpected response type from MyGene batch query: {type(payload)}")

    hits: List[Dict[str, Any]] = []
    for item in payload:
        if isinstance(item, dict) and not item.get("notfound", False):
            # Keep the best match per query term. If multiple matches occur, MyGene
            # may return multiple dicts with the same 'query'. We'll keep them all
            # here, and deduplicate later by symbol.
            hits.append(item)

    return hits



def safe_go_list(*, go_block: Any, aspect: str) -> str:
    """
    Convert a MyGene GO block into a semicolon-separated string for one aspect.

    Parameters
    ----------
    go_block
        GO block from MyGene, possibly missing or different shapes.
    aspect
        One of 'BP', 'MF', 'CC'.

    Returns
    -------
    str
        Semicolon-separated GO term names.
    """
    if not isinstance(go_block, dict):
        return ""
    key = {"BP": "BP", "MF": "MF", "CC": "CC"}.get(aspect)
    if key is None or key not in go_block:
        return ""

    terms = go_block.get(key)
    if terms is None:
        return ""
    # MyGene can return a dict for a single term or a list for multiple terms.
    if isinstance(terms, dict):
        name = str(terms.get("term", "")).strip()
        return name
    if isinstance(terms, list):
        names = []
        for t in terms:
            if isinstance(t, dict):
                name = str(t.get("term", "")).strip()
                if name:
                    names.append(name)
        # Deduplicate while preserving order
        seen = set()
        out = []
        for n in names:
            if n not in seen:
                seen.add(n)
                out.append(n)
        return ";".join(out)
    return ""


def build_annotation_table(
    *,
    symbols: List[str],
    cache: Dict[str, Dict[str, Any]],
    cfg: Config,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Build an annotation table for symbols.

    Parameters
    ----------
    symbols
        Unique HGNC symbols to annotate.
    cache
        Cache mapping symbol->hit.
    cfg
        Config.
    logger
        Logger.

    Returns
    -------
    pd.DataFrame
        Annotation table indexed by symbol upper with columns:
        gene_full_name, gene_summary, gene_type, go_bp, go_mf, go_cc.
    """
    wanted = [s.strip().upper() for s in symbols if str(s).strip()]
    wanted = sorted(set(wanted))

    hits_out: List[Dict[str, Any]] = []
    missing = [s for s in wanted if s not in cache]

    logger.info("Symbols requested: %d", len(wanted))
    logger.info("Symbols in cache: %d", len(cache))
    logger.info("Symbols to fetch from API: %d", len(missing))

    batches = chunk_list(items=missing, chunk_size=cfg.batch_size)

    for i, batch in enumerate(batches, start=1):
        logger.info("Fetching batch %d / %d (size=%d)", i, len(batches), len(batch))
        attempt = 0
        while True:
            attempt += 1
            try:
                hits = mygene_query_batch(
                    symbols=batch,
                    taxid=cfg.taxid,
                    timeout_seconds=cfg.timeout_seconds,
                )
                hits_out.extend(hits)
                if cfg.cache_jsonl:
                    append_cache(cache_jsonl=cfg.cache_jsonl, hits=hits)
                break
            except Exception as exc:
                if attempt >= cfg.max_retries:
                    raise
                logger.warning("Batch %d failed (attempt %d/%d): %s", i, attempt, cfg.max_retries, str(exc))
                time.sleep(min(5.0, cfg.sleep_seconds * attempt))

        time.sleep(cfg.sleep_seconds)

    # Update cache with fetched hits
    for h in hits_out:
        sym = str(h.get("symbol", "")).strip().upper()
        if sym:
            cache[sym] = h

    rows: List[Dict[str, Any]] = []
    for sym in wanted:
        h = cache.get(sym, {})
        rows.append(
            {
                "hgnc_symbol_norm": sym,
                "gene_full_name": str(h.get("name", "")).strip(),
                "gene_summary": str(h.get("summary", "")).strip(),
                "gene_type": str(h.get("type_of_gene", "")).strip(),
                "go_bp": safe_go_list(go_block=h.get("go"), aspect="BP"),
                "go_mf": safe_go_list(go_block=h.get("go"), aspect="MF"),
                "go_cc": safe_go_list(go_block=h.get("go"), aspect="CC"),
            }
        )

    ann = pd.DataFrame(rows).set_index("hgnc_symbol_norm", drop=True)
    return ann



def run(*, cfg: Config, logger: logging.Logger) -> None:
    """
    Run functional annotation and write output.

    Parameters
    ----------
    cfg
        Configuration.
    logger
        Logger.
    """
    df = pd.read_csv(cfg.input_tsv, sep="\t", index_col=0, low_memory=False)
    logger.info("Loaded input: %d rows x %d columns", df.shape[0], df.shape[1])

    if cfg.symbol_col not in df.columns:
        raise ValueError(
            f"Could not find symbol column '{cfg.symbol_col}'. "
            f"Available columns: {list(df.columns)[:30]}"
        )

    df = df.copy()
    df["hgnc_symbol_norm"] = df[cfg.symbol_col].astype(str).str.strip().str.upper()
    df.loc[df["hgnc_symbol_norm"].isin(["", "NAN"]), "hgnc_symbol_norm"] = ""


    # Build a fallback symbol from Description where HGNC is missing.
    # Many lncRNAs and some genes are present as symbols in 'Description'.
    fallback = df.get("Description", pd.Series([""] * df.shape[0], index=df.index)).astype(str).str.strip().str.upper()
    fallback = fallback.where(~fallback.str.startswith("ENSG"), "")

    df["query_symbol_norm"] = df["hgnc_symbol_norm"]
    df.loc[df["query_symbol_norm"] == "", "query_symbol_norm"] = fallback

    symbols = df.loc[df["query_symbol_norm"] != "", "query_symbol_norm"].tolist()


    cache: Dict[str, Dict[str, Any]] = {}
    if cfg.cache_jsonl:
        cache = load_cache(cache_jsonl=cfg.cache_jsonl, logger=logger)

    ann = build_annotation_table(symbols=symbols, cache=cache, cfg=cfg, logger=logger)
    logger.info("Annotation rows: %d", ann.shape[0])

    out = df.join(ann, on="query_symbol_norm")

    for c in ["gene_full_name", "gene_summary", "gene_type", "go_bp", "go_mf", "go_cc"]:
        if c in out.columns:
            out[c] = out[c].fillna("").astype(str)


    logger.info("After merge: %d rows x %d columns", out.shape[0], out.shape[1])

    logger.info("Writing output TSV: %s", cfg.output_tsv)
    out.to_csv(cfg.output_tsv, sep="\t", index=True)


def main() -> None:
    """
    Entry point.
    """
    cfg = parse_args()
    logger = setup_logging(log_level=cfg.log_level, log_path=cfg.log_path)
    run(cfg=cfg, logger=logger)


if __name__ == "__main__":
    main()
