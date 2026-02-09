#!/usr/bin/env python3
"""
Extract genes (and optionally diseases) associated with male infertility phenotypes using HPO ontology.

Inputs (from HPO downloads):
- hp.obo: HPO ontology (term names, parents/children, synonyms, alt_ids)
- genes_to_phenotype.txt: gene ↔ HPO mappings
- phenotype.hpoa: disease ↔ HPO annotations (optional enrichment output)

Method:
1) Start from a seed list of HPO IDs relevant to male infertility.
2) Expand to include all descendant terms using the ontology graph from hp.obo.
3) Map expanded HPO term set to genes via genes_to_phenotype.txt.
4) Optionally map expanded HPO term set to diseases via phenotype.hpoa.

Outputs are TSV (tab-separated).
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Deque, Dict, Iterable, List, Optional, Set, Tuple
import re


@dataclass(frozen=True)
class HpoTerm:
    """
    Representation of a single HPO term.

    Attributes:
        hpo_id: Primary HPO identifier (e.g. "HP:0000027").
        name: Term label.
        parents: Set of parent HPO IDs (is_a relationships).
        alt_ids: Set of alternative IDs that should be treated as synonyms of this term.
        synonyms: Set of synonym strings (informational; useful for reporting/search).
        is_obsolete: Whether term is marked obsolete in the ontology.
    """

    hpo_id: str
    name: str
    parents: Set[str]
    alt_ids: Set[str]
    synonyms: Set[str]
    is_obsolete: bool


DEFAULT_SEED_HPO_IDS: Tuple[str, ...] = (
    # Root concepts that are common entry points for male infertility work.
    # You can override/extend via --seed_hpo_ids or --seed_file.
    "HP:0000027",  # Azoospermia
    "HP:0000044",  # Hypogonadism
    "HP:0000789",  # Infertility
    "HP:0003251",  # Abnormal sperm morphology (term name may vary by version)
)


DEFAULT_SEED_PHRASES: Tuple[str, ...] = (
    "male infertility",
    "azoospermia",
    "non-obstructive azoospermia",
    "obstructive azoospermia",
    "oligozoospermia",
    "severe oligozoospermia",
    "asthenozoospermia",
    "teratozoospermia",
    "oligoasthenoteratozoospermia",
    "globozoospermia",
    "multiple morphological abnormalities of the sperm flagella",
    "sperm flagellum",
    "abnormal sperm motility",
    "abnormal sperm morphology",
    "spermatogenic failure",
    "maturation arrest",
    "sertoli cell-only",
    "hypospermatogenesis",
    "congenital bilateral absence of the vas deferens",
    "absent vas deferens",
    "ejaculatory dysfunction",
    "retrograde ejaculation",
    "hypospermia",
    "aspermia",
    "hypogonadism",
)


MALE_INCLUDE_REGEX = re.compile(
    r"(male infertility|male hypogonadism|sperm|semen|azoospermia|oligozoospermia|"
    r"asthenozoospermia|teratozoospermia|cryptozoospermia|globozoospermia|"
    r"macrozoospermia|vas deferens|ejaculat|spermat|testi|testicular|"
    r"flagell|axoneme|acrosom)",
    flags=re.IGNORECASE,
)

FEMALE_EXCLUDE_REGEX = re.compile(
    r"(female infertility|oocyte|ovary|uterus|uterine|endomet|menstrual|ovulat|"
    r"fallopian|cervix|vaginal|vagina)",
    flags=re.IGNORECASE,
)



def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Extract male infertility gene sets using HPO ontology expansion."
    )
    parser.add_argument(
        "--hpo_dir",
        required=True,
        help="Directory containing hp.obo, genes_to_phenotype.txt, phenotype.hpo",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Directory to write TSV outputs",
    )
    parser.add_argument(
        "--seed_hpo_ids",
        default="",
        help=(
            "Comma-separated list of seed HPO IDs to expand from. "
            "If omitted, a default seed set is used."
        ),
    )
    parser.add_argument(
        "--seed_file",
        default="",
        help=(
            "Optional path to a text file with one seed HPO ID per line. "
            "Lines starting with # are ignored."
        ),
    )
    parser.add_argument(
        "--include_diseases",
        action="store_true",
        help="If set, also extract diseases from phenotype.hpoa for the expanded HPO set.",
    )
    parser.add_argument(
        "--require_male",
        action="store_true",
        help=(
            "If set and phenotype.hpoa contains a Sex column, restrict disease annotations to MALE."
        ),
    )
    parser.add_argument(
        "--seed_phrases_file",
        default="",
        help=(
            "Optional path to a text file with one seed phrase per line. "
            "If provided, seeds will be discovered by matching these phrases against hp.obo "
            "term names/synonyms."
        ),
    )
    parser.add_argument(
        "--use_phrase_seeds",
        action="store_true",
        help="If set, discover seed HPO IDs by matching seed phrases against hp.obo.",
    )

    return parser.parse_args()


def load_seed_hpo_ids(seed_hpo_ids: str, seed_file: str) -> Set[str]:
    """
    Load seed HPO IDs from CLI and/or file.

    Args:
        seed_hpo_ids: Comma-separated list of HPO IDs.
        seed_file: Path to file containing one HPO ID per line.

    Returns:
        Set of seed HPO IDs.
    """
    seeds: Set[str] = set()

    if seed_hpo_ids.strip():
        for part in seed_hpo_ids.split(","):
            part = part.strip()
            if part:
                seeds.add(part)

    if seed_file.strip():
        path = Path(seed_file)
        with path.open(mode="r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                seeds.add(line)

    if not seeds:
        seeds.update(DEFAULT_SEED_HPO_IDS)

    return seeds


def load_seed_phrases(seed_phrases_file: str) -> List[str]:
    """
    Load seed phrases from file or return defaults.

    Args:
        seed_phrases_file: Optional path to phrase file.

    Returns:
        List of phrases.
    """
    if seed_phrases_file.strip() == "":
        return list(DEFAULT_SEED_PHRASES)

    phrases: List[str] = []
    with Path(seed_phrases_file).open(mode="r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith("#"):
                phrases.append(line)
    return phrases


def parse_hp_obo(hp_obo_path: Path) -> Dict[str, HpoTerm]:
    """
    Parse hp.obo and return a mapping of HPO ID to HpoTerm.

    This parser focuses on fields needed for expansion and reporting:
    - id
    - name
    - is_a
    - alt_id
    - synonym
    - is_obsolete

    Args:
        hp_obo_path: Path to hp.obo.

    Returns:
        Dictionary mapping HPO ID to HpoTerm.
    """
    terms_raw: Dict[str, Dict[str, object]] = {}
    current: Dict[str, object] = {}
    in_term = False

    def _commit_term(block: Dict[str, object]) -> None:
        term_id = str(block.get("id", "")).strip()
        if not term_id.startswith("HP:"):
            return
        terms_raw[term_id] = block

    with hp_obo_path.open(mode="r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("!"):
                continue

            if line == "[Term]":
                if in_term and current:
                    _commit_term(current)
                current = {
                    "id": "",
                    "name": "",
                    "parents": set(),
                    "alt_ids": set(),
                    "synonyms": set(),
                    "is_obsolete": False,
                }
                in_term = True
                continue

            if line.startswith("[") and line.endswith("]") and line != "[Term]":
                if in_term and current:
                    _commit_term(current)
                current = {}
                in_term = False
                continue

            if not in_term or ":" not in line:
                continue

            key, value = line.split(":", 1)
            key = key.strip()
            value = value.strip()

            if key == "id":
                current["id"] = value
            elif key == "name":
                current["name"] = value
            elif key == "is_a":
                parent_id = value.split("!")[0].strip()
                if parent_id.startswith("HP:"):
                    current["parents"].add(parent_id)
            elif key == "alt_id":
                if value.startswith("HP:"):
                    current["alt_ids"].add(value)
            elif key == "synonym":
                # synonym: "text" EXACT [source]
                match = re.search(r"\"(.*?)\"", value)
                if match:
                    current["synonyms"].add(match.group(1))
            elif key == "is_obsolete":
                current["is_obsolete"] = (value.lower() == "true")

    if in_term and current:
        _commit_term(current)

    terms: Dict[str, HpoTerm] = {}
    for hpo_id, block in terms_raw.items():
        terms[hpo_id] = HpoTerm(
            hpo_id=hpo_id,
            name=str(block.get("name", "")),
            parents=set(block.get("parents", set())),
            alt_ids=set(block.get("alt_ids", set())),
            synonyms=set(block.get("synonyms", set())),
            is_obsolete=bool(block.get("is_obsolete", False)),
        )

    return terms


def filter_terms_male_relevant(
    hpo_ids: Set[str],
    terms: Dict[str, HpoTerm],
) -> Set[str]:
    """
    Filter an expanded HPO term set to retain male infertility-relevant terms.

    The filter keeps terms whose name/synonyms match MALE_INCLUDE_REGEX and removes those
    whose name/synonyms match FEMALE_EXCLUDE_REGEX.

    Args:
        hpo_ids: Expanded HPO IDs.
        terms: HPO term dictionary.

    Returns:
        Filtered set of HPO IDs.
    """
    kept: Set[str] = set()

    for hpo_id in hpo_ids:
        term = terms.get(hpo_id)
        if term is None:
            continue

        text = " ".join([term.name] + sorted(term.synonyms))

        if FEMALE_EXCLUDE_REGEX.search(text):
            continue

        if MALE_INCLUDE_REGEX.search(text):
            kept.add(hpo_id)

    return kept


def build_children_index(terms: Dict[str, HpoTerm]) -> Dict[str, Set[str]]:
    """
    Build a parent -> children adjacency mapping.

    Args:
        terms: Mapping of HPO ID to HpoTerm.

    Returns:
        Dictionary where keys are parent IDs and values are sets of child IDs.
    """
    children: Dict[str, Set[str]] = defaultdict(set)
    for term in terms.values():
        for parent_id in term.parents:
            children[parent_id].add(term.hpo_id)
    return children


def find_hpo_ids_by_text(
    terms: Dict[str, HpoTerm],
    phrases: List[str],
    include_obsolete: bool = False,
) -> Set[str]:
    """
    Find HPO IDs whose name or synonyms match any of the provided phrases.

    Args:
        terms: Mapping of HPO ID to HpoTerm.
        phrases: List of phrases to search for (case-insensitive).
        include_obsolete: If True, include obsolete terms in matches.

    Returns:
        Set of matching HPO IDs (primary IDs).
    """
    if not phrases:
        return set()

    pattern = re.compile("|".join(re.escape(p) for p in phrases), flags=re.IGNORECASE)
    matched: Set[str] = set()

    for hpo_id, term in terms.items():
        if (not include_obsolete) and term.is_obsolete:
            continue

        haystack_parts = [term.name] + sorted(term.synonyms)
        haystack = " ".join(haystack_parts)

        if pattern.search(haystack):
            matched.add(hpo_id)

    return matched


def write_seed_terms(
    output_path: Path,
    seed_ids: Set[str],
    terms: Dict[str, HpoTerm],
) -> None:
    """
    Write the seed HPO terms used (IDs and names) for reproducibility.

    Args:
        output_path: TSV output path.
        seed_ids: Seed HPO IDs.
        terms: HPO terms dictionary.
    """
    with output_path.open(mode="w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["seed_hpo_id", "name", "n_synonyms", "synonyms"])
        for hpo_id in sorted(seed_ids):
            term = terms.get(hpo_id)
            if term is None:
                writer.writerow([hpo_id, "", "0", ""])
                continue
            synonyms = sorted(term.synonyms)
            writer.writerow([hpo_id, term.name, str(len(synonyms)), ";".join(synonyms)])


def expand_descendants(
    seed_ids: Set[str],
    terms: Dict[str, HpoTerm],
    children_index: Dict[str, Set[str]],
    include_obsolete: bool = False,
) -> Set[str]:
    """
    Expand a seed set of HPO IDs to include all descendant terms (ontology traversal).

    Args:
        seed_ids: Seed HPO IDs.
        terms: Mapping of HPO ID to HpoTerm.
        children_index: Parent -> children adjacency mapping.
        include_obsolete: If True, include obsolete terms in expansion.

    Returns:
        Set of expanded HPO IDs (primary IDs only).
    """
    expanded: Set[str] = set()
    queue: Deque[str] = deque()

    for hpo_id in seed_ids:
        queue.append(hpo_id)

    while queue:
        current = queue.popleft()
        if current in expanded:
            continue

        term = terms.get(current)
        if term is None:
            # Keep unknown seeds as-is (still useful for matching external mappings),
            # but do not expand further.
            expanded.add(current)
            continue

        if (not include_obsolete) and term.is_obsolete:
            continue

        expanded.add(current)

        for child_id in children_index.get(current, set()):
            if child_id not in expanded:
                queue.append(child_id)

    return expanded


def build_alt_id_index(terms: Dict[str, HpoTerm]) -> Dict[str, str]:
    """
    Build a mapping from alt_id -> primary_id.

    Args:
        terms: Mapping of HPO ID to HpoTerm.

    Returns:
        Dictionary mapping alternative IDs to their primary term ID.
    """
    alt_to_primary: Dict[str, str] = {}
    for term in terms.values():
        for alt_id in term.alt_ids:
            alt_to_primary[alt_id] = term.hpo_id
    return alt_to_primary


def write_hpo_term_set(
    output_path: Path,
    hpo_ids: Set[str],
    terms: Dict[str, HpoTerm],
) -> None:
    """
    Write the expanded HPO set with names and synonyms.

    Args:
        output_path: TSV output path.
        hpo_ids: Expanded HPO IDs (primary IDs).
        terms: HPO term dictionary.
    """
    with output_path.open(mode="w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["hpo_id", "name", "is_obsolete", "n_synonyms", "synonyms", "parents"])
        for hpo_id in sorted(hpo_ids):
            term = terms.get(hpo_id)
            if term is None:
                writer.writerow([hpo_id, "", "", "0", "", ""])
                continue
            synonyms = sorted(term.synonyms)
            parents = sorted(term.parents)
            writer.writerow(
                [
                    term.hpo_id,
                    term.name,
                    "true" if term.is_obsolete else "false",
                    str(len(synonyms)),
                    ";".join(synonyms),
                    ";".join(parents),
                ]
            )


def extract_genes_for_hpo_ids(
    genes_to_phenotype_path: Path,
    target_hpo_ids: Set[str],
    alt_to_primary: Dict[str, str],
) -> Tuple[Dict[Tuple[str, str], Set[str]], List[Tuple[str, str, str]]]:
    """
    Extract genes associated with target HPO IDs from genes_to_phenotype.txt.

    Args:
        genes_to_phenotype_path: Path to genes_to_phenotype.txt.
        target_hpo_ids: Set of primary HPO IDs to match.
        alt_to_primary: Mapping from alt HPO IDs to primary IDs.

    Returns:
        gene_to_hpos: {(gene_id, gene_symbol): {primary_hpo_ids}}
        gene_hpo_pairs: list of (gene_id, gene_symbol, primary_hpo_id)
    """
    gene_to_hpos: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    gene_hpo_pairs: List[Tuple[str, str, str]] = []

    with genes_to_phenotype_path.open(mode="r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            raw_hpo_id = (row.get("HPO_ID") or row.get("hpo_id") or row.get("phenotype_id") or "").strip()
            if not raw_hpo_id:
                continue

            primary_hpo_id = alt_to_primary.get(raw_hpo_id, raw_hpo_id)
            if primary_hpo_id not in target_hpo_ids:
                continue

            gene_id = (row.get("GeneID") or row.get("gene_id") or row.get("entrez_gene_id") or "").strip()
            gene_symbol = (row.get("GeneSymbol") or row.get("gene_symbol") or row.get("symbol") or "").strip()
            if not gene_id and not gene_symbol:
                continue

            key = (gene_id, gene_symbol)
            gene_to_hpos[key].add(primary_hpo_id)
            gene_hpo_pairs.append((gene_id, gene_symbol, primary_hpo_id))

    return gene_to_hpos, gene_hpo_pairs


def write_gene_summary(
    output_path: Path,
    gene_to_hpos: Dict[Tuple[str, str], Set[str]],
    terms: Dict[str, HpoTerm],
) -> None:
    """
    Write a per-gene summary TSV including HPO IDs and names.

    Args:
        output_path: Output TSV path.
        gene_to_hpos: Mapping of (gene_id, gene_symbol) to HPO IDs.
        terms: HPO terms dictionary (for names).
    """
    with output_path.open(mode="w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_id", "gene_symbol", "n_hpo_terms", "hpo_ids", "hpo_names"])
        for (gene_id, gene_symbol), hpo_ids in sorted(gene_to_hpos.items(), key=lambda x: (x[0][1], x[0][0])):
            hpo_ids_sorted = sorted(hpo_ids)
            hpo_names = [terms[h].name if h in terms else "" for h in hpo_ids_sorted]
            writer.writerow(
                [
                    gene_id,
                    gene_symbol,
                    str(len(hpo_ids_sorted)),
                    ";".join(hpo_ids_sorted),
                    ";".join(hpo_names),
                ]
            )


def write_gene_hpo_pairs(
    output_path: Path,
    gene_hpo_pairs: List[Tuple[str, str, str]],
    terms: Dict[str, HpoTerm],
) -> None:
    """
    Write long-format gene ↔ HPO mappings with term names.

    Args:
        output_path: Output TSV path.
        gene_hpo_pairs: List of (gene_id, gene_symbol, hpo_id).
        terms: HPO terms dictionary (for names).
    """
    with output_path.open(mode="w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_id", "gene_symbol", "hpo_id", "hpo_name"])
        for gene_id, gene_symbol, hpo_id in sorted(gene_hpo_pairs, key=lambda x: (x[1], x[0], x[2])):
            writer.writerow([gene_id, gene_symbol, hpo_id, terms[hpo_id].name if hpo_id in terms else ""])



def _normalise_header(header: str) -> str:
    """
    Normalise column headers for robust matching.

    Args:
        header: Original header.

    Returns:
        Normalised header string.
    """
    return re.sub(pattern=r"[^a-z0-9]+", repl="_", string=header.strip().lower()).strip("_")


def _find_column_index(headers: List[str], candidates: Iterable[str]) -> Optional[int]:
    """
    Find a column index by matching normalised header names.

    Args:
        headers: Column headers.
        candidates: Candidate header names.

    Returns:
        Index if found, otherwise None.
    """
    normalised = [_normalise_header(h) for h in headers]
    candidate_set = {_normalise_header(c) for c in candidates}
    for idx, h in enumerate(normalised):
        if h in candidate_set:
            return idx
    return None


def extract_diseases_for_hpo_ids(
    phenotype_hpoa_path: Path,
    target_hpo_ids: Set[str],
    alt_to_primary: Dict[str, str],
    require_male: bool,
) -> List[Tuple[str, str, str]]:
    """
    Extract disease annotations that use any of the target HPO IDs.

    This function is robust to HPOA files that start with comment lines beginning with '#'.
    It reads the first non-comment line as the header.

    Args:
        phenotype_hpoa_path: Path to phenotype.hpoa.
        target_hpo_ids: Target primary HPO IDs.
        alt_to_primary: alt_id -> primary_id mapping.
        require_male: If True, keep only rows annotated as male when Sex column exists.

    Returns:
        List of tuples (disease_id, disease_name, hpo_id_primary).
    """
    results: List[Tuple[str, str, str]] = []

    with phenotype_hpoa_path.open(mode="r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")

        headers: Optional[List[str]] = None
        for row in reader:
            if not row:
                continue
            if row[0].startswith("#"):
                continue
            headers = row
            break

        if headers is None:
            raise ValueError("phenotype.hpoa appears to contain no header row (only comments/empty lines).")

        disease_id_idx = _find_column_index(headers, candidates=["disease_id", "database_id"])
        disease_name_idx = _find_column_index(headers, candidates=["disease_name"])
        hpo_id_idx = _find_column_index(headers, candidates=["hpo_id", "phenotype_id"])
        sex_idx = _find_column_index(headers, candidates=["sex"])

        if disease_name_idx is None or hpo_id_idx is None:
            raise ValueError(
                "Could not find required columns in phenotype.hpoa. "
                "Expected at least 'disease_name' and 'hpo_id'. "
                f"Detected headers: {headers}"
            )

        for row in reader:
            if not row:
                continue
            if row[0].startswith("#"):
                continue

            raw_hpo_id = row[hpo_id_idx].strip() if len(row) > hpo_id_idx else ""
            if not raw_hpo_id:
                continue

            primary_hpo_id = alt_to_primary.get(raw_hpo_id, raw_hpo_id)
            if primary_hpo_id not in target_hpo_ids:
                continue

            if require_male and sex_idx is not None and len(row) > sex_idx:
                sex_val = row[sex_idx].strip().upper()
                if sex_val and sex_val not in {"MALE", "M"}:
                    continue

            disease_id = (
                row[disease_id_idx].strip()
                if (disease_id_idx is not None and len(row) > disease_id_idx)
                else ""
            )
            disease_name = row[disease_name_idx].strip()

            results.append((disease_id, disease_name, primary_hpo_id))

    return results





def write_disease_hpo_pairs(
    output_path: Path,
    disease_hpo_pairs: List[Tuple[str, str, str]],
    terms: Dict[str, HpoTerm],
) -> None:
    """
    Write long-format disease ↔ HPO mappings with term names.

    Args:
        output_path: Output TSV path.
        disease_hpo_pairs: List of (disease_id, disease_name, hpo_id).
        terms: HPO terms dictionary.
    """
    with output_path.open(mode="w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["disease_id", "disease_name", "hpo_id", "hpo_name"])
        for disease_id, disease_name, hpo_id in sorted(disease_hpo_pairs, key=lambda x: (x[1], x[0], x[2])):
            writer.writerow([disease_id, disease_name, hpo_id, terms[hpo_id].name if hpo_id in terms else ""])


def main() -> None:
    """
    Run ontology expansion and extract genes/diseases for expanded HPO term set.
    """
    args = parse_args()
    hpo_dir = Path(args.hpo_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    hp_obo_path = hpo_dir / "hp.obo"
    genes_to_phenotype_path = hpo_dir / "genes_to_phenotype.txt"
    phenotype_hpoa_path = hpo_dir / "phenotype.hpoa"

    if not hp_obo_path.exists():
        raise FileNotFoundError(f"Missing file: {hp_obo_path}")
    if not genes_to_phenotype_path.exists():
        raise FileNotFoundError(f"Missing file: {genes_to_phenotype_path}")
    if args.include_diseases and not phenotype_hpoa_path.exists():
        raise FileNotFoundError(f"Missing file: {phenotype_hpoa_path}")

    terms = parse_hp_obo(hp_obo_path=hp_obo_path)
    children_index = build_children_index(terms=terms)
    alt_to_primary = build_alt_id_index(terms=terms)

    if args.use_phrase_seeds:
        phrases = load_seed_phrases(seed_phrases_file=args.seed_phrases_file)
        seed_ids = find_hpo_ids_by_text(terms=terms, phrases=phrases, include_obsolete=False)
        if not seed_ids:
            raise ValueError(
                "No seed HPO IDs were found by phrase matching. "
                "Try broadening phrases or check hp.obo content."
            )
    else:
        seed_ids = load_seed_hpo_ids(seed_hpo_ids=args.seed_hpo_ids, seed_file=args.seed_file)

    write_seed_terms(
        output_path=out_dir / "male_infertility_seed_hpo_terms.tsv",
        seed_ids=seed_ids,
        terms=terms,
    )

    expanded_primary_ids = expand_descendants(
        seed_ids=seed_ids,
        terms=terms,
        children_index=children_index,
        include_obsolete=False,
    )


    # Also ensure seeds that are alt_ids map correctly if user supplies them.
    expanded_primary_ids = {alt_to_primary.get(h, h) for h in expanded_primary_ids}

    filtered_primary_ids = filter_terms_male_relevant(
    hpo_ids=expanded_primary_ids,
    terms=terms,
)

    # Replace downstream usage to use filtered_primary_ids
    expanded_primary_ids = filtered_primary_ids

    # female terms that creep in
    expanded_primary_ids.discard("HP:0000134")
    expanded_primary_ids.discard("HP:0008222")



    write_hpo_term_set(
        output_path=out_dir / "male_infertility_filtered_hpo_terms.tsv",
        hpo_ids=expanded_primary_ids,
        terms=terms,
    )


    write_hpo_term_set(
        output_path=out_dir / "male_infertility_expanded_hpo_terms.tsv",
        hpo_ids=expanded_primary_ids,
        terms=terms,
    )

    gene_to_hpos, gene_hpo_pairs = extract_genes_for_hpo_ids(
        genes_to_phenotype_path=genes_to_phenotype_path,
        target_hpo_ids=expanded_primary_ids,
        alt_to_primary=alt_to_primary,
    )

    write_gene_summary(
        output_path=out_dir / "male_infertility_genes_summary.tsv",
        gene_to_hpos=gene_to_hpos,
        terms=terms,
    )
    write_gene_hpo_pairs(
        output_path=out_dir / "male_infertility_gene_hpo_pairs.tsv",
        gene_hpo_pairs=gene_hpo_pairs,
        terms=terms,
    )

    if args.include_diseases:
        disease_hpo_pairs = extract_diseases_for_hpo_ids(
            phenotype_hpoa_path=phenotype_hpoa_path,
            target_hpo_ids=expanded_primary_ids,
            alt_to_primary=alt_to_primary,
            require_male=bool(args.require_male),
        )
        write_disease_hpo_pairs(
            output_path=out_dir / "male_infertility_disease_hpo_pairs.tsv",
            disease_hpo_pairs=disease_hpo_pairs,
            terms=terms,
        )


if __name__ == "__main__":
    main()
