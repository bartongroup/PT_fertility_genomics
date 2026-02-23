#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
download_encode_bigwig_tracks.py

Download ENCODE bigWig tracks for GRCh38 using the ENCODE REST API.

This script is intended to fetch genome-wide tracks that can be summarised
around gene bodies or TSS windows (e.g. replication timing, Lamin B1 signal).

Outputs
-------
- manifest TSV describing downloaded files
- downloaded bigWig files in the chosen output directory

Notes
-----
- ENCODE provides data mapped to GRCh38. :contentReference[oaicite:2]{index=2}
- You can select one or more biosamples and assays.

Example
-------
python download_encode_bigwig_tracks.py \
  --out_dir encode_tracks \
  --assembly GRCh38 \
  --assay_title "Repli-seq" \
  --assay_title "ChIP-seq" \
  --target "Lamin B1" \
  --biosample "H1-hESC" \
  --biosample "GM12878" \
  --max_files_per_query 3 \
  --verbose
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd


ENCODE_BASE = "https://www.encodeproject.org"
USER_AGENT = "PT_fertility_genomics/1.0 (contact: pthorpe001)"


@dataclass
class EncodeFile:
    """Container for ENCODE file metadata."""
    accession: str
    href: str
    file_format: str
    output_type: str
    biosample: str
    assay_title: str
    assembly: str
    target: str
    status: str
    file_size: int


def setup_logger(verbose: bool) -> None:
    """Set up logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def http_get_json(url: str) -> Dict:
    """GET JSON from a URL."""
    req = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": USER_AGENT,
        },
    )
    with urllib.request.urlopen(req) as resp:
        return json.loads(resp.read().decode("utf-8"))


def build_search_url(
    *,
    assembly: str,
    assay_titles: Sequence[str],
    biosamples: Sequence[str],
    targets: Sequence[str],
    file_format: str,
    output_type: str,
    limit: int
) -> str:
    """
    Build an ENCODE search URL for files.

    We query file objects directly for simplicity.

    Parameters
    ----------
    assembly
        Genome assembly, e.g. GRCh38.
    assay_titles
        Assay titles (e.g. "Repli-seq", "ChIP-seq").
    biosamples
        Biosample term names, e.g. "H1-hESC".
    targets
        Target labels (e.g. "Lamin B1") for ChIP-seq.
    file_format
        Usually "bigWig".
    output_type
        Common output type for signals is "signal p-value" or "fold change"
        depending on assay; here we default to "signal p-value" if present,
        but we also accept any signal-type bigWig below.
    limit
        ENCODE limit.

    Returns
    -------
    str
        Search URL.
    """
    params = [
        ("type", "File"),
        ("assembly", assembly),
        ("file_format", file_format),
        ("limit", str(limit)),
        ("format", "json"),
    ]
    for a in assay_titles:
        params.append(("assay_title", a))
    for b in biosamples:
        params.append(("biosample_ontology.term_name", b))
    for t in targets:
        params.append(("target.label", t))
    if output_type:
        params.append(("output_type", output_type))

    return f"{ENCODE_BASE}/search/?{urllib.parse.urlencode(params)}"


def parse_encode_files(payload: Dict) -> List[EncodeFile]:
    """Parse ENCODE file search results into a list."""
    files: List[EncodeFile] = []
    for item in payload.get("@graph", []):
        href = item.get("href")
        accession = item.get("accession")
        file_format = item.get("file_format", "")
        output_type = item.get("output_type", "")
        assembly = item.get("assembly", "")
        status = item.get("status", "")
        file_size = int(item.get("file_size") or 0)

        biosample = ""
        assay_title = ""
        target = ""
        dataset = item.get("dataset") or ""
        if dataset:
            try:
                ds = http_get_json(f"{ENCODE_BASE}{dataset}?format=json")
                assay_title = ds.get("assay_title") or ds.get("assay_term_name") or ""
                bios = ds.get("biosample_ontology", {})
                biosample = bios.get("term_name", "") if isinstance(bios, dict) else ""
                t = ds.get("target", {})
                target = t.get("label", "") if isinstance(t, dict) else ""
            except Exception:
                pass

        if not href or not accession:
            continue

        files.append(
            EncodeFile(
                accession=accession,
                href=href,
                file_format=file_format,
                output_type=output_type,
                biosample=biosample,
                assay_title=assay_title,
                assembly=assembly,
                target=target,
                status=status,
                file_size=file_size,
            )
        )
    return files


def download_file(url: str, out_path: str, *, retries: int = 8, sleep_s: int = 10) -> None:
    """
    Download a file with retries.

    Parameters
    ----------
    url
        Download URL.
    out_path
        Destination path.
    retries
        Number of retries.
    sleep_s
        Seconds between retries.
    """
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        logging.info("Already present, skipping: %s", out_path)
        return

    tmp_path = f"{out_path}.part"
    for attempt in range(1, retries + 1):
        try:
            logging.info("Downloading (%s/%s): %s", attempt, retries, url)
            req = urllib.request.Request(
                url,
                headers={"User-Agent": USER_AGENT},
            )
            with urllib.request.urlopen(req) as resp, open(tmp_path, "wb") as fh:
                while True:
                    chunk = resp.read(1024 * 1024 * 8)
                    if not chunk:
                        break
                    fh.write(chunk)

            os.replace(tmp_path, out_path)
            logging.info("Saved: %s", out_path)
            return
        except Exception as exc:
            logging.warning("Download failed: %s", exc)
            if os.path.exists(tmp_path):
                try:
                    os.remove(tmp_path)
                except OSError:
                    pass
            time.sleep(sleep_s)

    raise RuntimeError(f"Failed to download after {retries} attempts: {url}")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse CLI args."""
    p = argparse.ArgumentParser()
    p.add_argument("--out_dir", required=True)
    p.add_argument("--assembly", default="GRCh38")
    p.add_argument("--assay_title", action="append", default=[])
    p.add_argument("--biosample", action="append", default=[])
    p.add_argument("--target", action="append", default=[])
    p.add_argument("--file_format", default="bigWig")
    p.add_argument("--output_type", default="")
    p.add_argument("--max_files_per_query", type=int, default=3)
    p.add_argument("--manifest_tsv", default="encode_tracks_manifest.tsv")
    p.add_argument("--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run."""
    args = parse_args(argv)
    setup_logger(args.verbose)

    os.makedirs(args.out_dir, exist_ok=True)

    assay_titles = args.assay_title or ["Repli-seq"]
    biosamples = args.biosample or ["H1-hESC"]
    targets = args.target or []

    url = build_search_url(
        assembly=args.assembly,
        assay_titles=assay_titles,
        biosamples=biosamples,
        targets=targets,
        file_format=args.file_format,
        output_type=args.output_type,
        limit=200,
    )
    logging.info("ENCODE search URL: %s", url)

    payload = http_get_json(url)
    files = parse_encode_files(payload)

    if not files:
        raise SystemExit("No files returned. Try adjusting assay/biosample/target filters.")

    files = sorted(files, key=lambda x: (x.status != "released", x.file_size))
    files = files[: args.max_files_per_query]

    records: List[Dict] = []
    for f in files:
        dl_url = f"{ENCODE_BASE}{f.href}"
        safe_label = "_".join([x for x in [f.assay_title, f.biosample, f.target, f.accession] if x])
        safe_label = safe_label.replace(" ", "_").replace("/", "_")
        out_path = os.path.join(args.out_dir, f"{safe_label}.bigWig")

        download_file(dl_url, out_path)

        records.append(
            {
                "accession": f.accession,
                "assay_title": f.assay_title,
                "biosample": f.biosample,
                "target": f.target,
                "assembly": f.assembly,
                "output_type": f.output_type,
                "status": f.status,
                "file_size": f.file_size,
                "download_url": dl_url,
                "local_path": out_path,
            }
        )

    manifest_path = os.path.join(args.out_dir, args.manifest_tsv)
    pd.DataFrame(records).to_csv(manifest_path, sep="\t", index=False)
    logging.info("Wrote manifest: %s", manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
