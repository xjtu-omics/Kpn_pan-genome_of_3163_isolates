import argparse
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

DEFAULT_TARGET_GENES = [
    "group_7560",
    "group_12673",
    "group_2412",
    "bla_3~~~bla_1~~~bla_2~~~bla_6",
    "bla_7~~~bla_6~~~bla_3~~~bla_4",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Analyze gene order patterns per sample from CSV matrix and sample GFF files."
        )
    )
    parser.add_argument(
        "--csv",
        required=True,
        help="Path to the CSV file. The first column is read as row index.",
    )
    parser.add_argument(
        "--gff-dir",
        required=True,
        help="Directory containing sample GFF/GFF3 files.",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Output directory for result tables.",
    )
    parser.add_argument(
        "--target-genes",
        default=",".join(DEFAULT_TARGET_GENES),
        help=(
            "Comma-separated row names to extract from the CSV index. "
            "Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--sample-start-col",
        default=None,
        help=(
            "First sample column name, or 0-based column index after reading the CSV "
            "index. Columns before this are treated as annotations and skipped. "
            "Default: infer from the first column that matches a GFF file name."
        ),
    )
    parser.add_argument(
        "--min-genes-per-pattern",
        type=int,
        default=1,
        help="Minimum genes on a contig to form a pattern. Default: 1.",
    )
    return parser.parse_args()


def parse_attributes(attr_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    if not attr_text:
        return attrs
    for item in attr_text.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key.strip()] = value.strip()
        elif " " in item:
            key, value = item.split(" ", 1)
            attrs[key.strip()] = value.strip().strip('"')
    return attrs


def clean_id_variants(raw_id: str) -> List[str]:
    """
    Generate multiple key variants for robust lookup from CSV cell to GFF feature.
    """
    if not raw_id:
        return []
    rid = str(raw_id).strip()
    if not rid:
        return []
    keys = {rid}
    if ":" in rid:
        keys.add(rid.split(":", 1)[-1])
    for prefix in ("gene-", "cds-", "rna-", "mRNA-"):
        if rid.startswith(prefix):
            keys.add(rid[len(prefix) :])
    if "." in rid:
        keys.add(rid.split(".", 1)[0])
    return [k for k in keys if k]


def split_cell_gene_ids(cell_value: object) -> List[str]:
    """
    Split one CSV cell into one or more possible feature IDs.
    Supports common delimiters used in merged annotations.
    """
    if cell_value is None:
        return []
    if isinstance(cell_value, float) and math.isnan(cell_value):
        return []
    text = str(cell_value).strip()
    if not text or text.lower() == "nan":
        return []
    parts = re.split(r"[,\s;|~]+", text)
    return [p for p in parts if p]


def extract_numeric_part(text: str) -> Optional[int]:
    """
    Extract trailing numeric token from an ID-like string.
    Example: "KPN_00123" -> 123
    """
    if not text:
        return None
    m = re.search(r"(\d+)(?!.*\d)", str(text))
    if not m:
        return None
    try:
        return int(m.group(1))
    except ValueError:
        return None


def is_strictly_consecutive(nums: List[Optional[int]]) -> bool:
    """
    Check whether numeric IDs are consecutive with step +1 or -1.
    """
    if len(nums) <= 1:
        return True
    if any(n is None for n in nums):
        return False
    clean_nums = [int(n) for n in nums]
    diffs = [b - a for a, b in zip(clean_nums[:-1], clean_nums[1:])]
    return all(d == 1 for d in diffs) or all(d == -1 for d in diffs)


def load_gff_index(gff_path: Path) -> Tuple[Dict[str, dict], Dict[str, int], str]:
    """
    Return:
    - id_index: map from id-like key -> feature position info
    - contig_max_end: map contig -> max end coordinate
    - primary_contig: inferred chromosome/main contig
    """
    id_index: Dict[str, dict] = {}
    contig_max_end: Dict[str, int] = defaultdict(int)
    contig_declared_len: Dict[str, int] = {}
    chromosome_like_contigs: set = set()
    contig_cds_order_counter: Dict[str, int] = defaultdict(int)
    contig_feature_order_counter: Dict[str, int] = defaultdict(int)

    with gff_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith("##sequence-region"):
                # Example: ##sequence-region AXLF01000001.1 1 217478
                parts = line.strip().split()
                if len(parts) >= 4:
                    seqid = parts[1]
                    try:
                        end_pos = int(parts[3])
                        contig_declared_len[seqid] = end_pos
                    except ValueError:
                        pass
                continue
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs_raw = parts
            if ftype.lower() not in {"gene", "cds"}:
                continue
            try:
                s = int(start)
                e = int(end)
            except ValueError:
                continue

            attrs = parse_attributes(attrs_raw)
            contig_max_end[seqid] = max(contig_max_end[seqid], e)
            seqid_lower = seqid.lower()
            if "chromosome" in seqid_lower or "chrom" in seqid_lower:
                chromosome_like_contigs.add(seqid)

            contig_feature_order_counter[seqid] += 1
            feature_order = contig_feature_order_counter[seqid]
            cds_order: Optional[int] = None
            if ftype.lower() == "cds":
                contig_cds_order_counter[seqid] += 1
                cds_order = contig_cds_order_counter[seqid]

            candidates = set()
            for key_name in ("ID", "Name", "locus_tag", "gene", "old_locus_tag", "protein_id"):
                if key_name in attrs and attrs[key_name]:
                    for piece in str(attrs[key_name]).split(","):
                        candidates.update(clean_id_variants(piece))

            if not candidates:
                continue

            record = {
                "contig": seqid,
                "start": s,
                "end": e,
                "strand": strand if strand in {"+", "-"} else ".",
                "feature_type": ftype,
                "feature_order": feature_order,
                "cds_order": cds_order,
            }
            for cid in candidates:
                # Prefer CDS mapping when both gene and CDS share a similar identifier.
                if cid not in id_index:
                    id_index[cid] = record
                else:
                    old_ftype = str(id_index[cid].get("feature_type", "")).lower()
                    if old_ftype != "cds" and ftype.lower() == "cds":
                        id_index[cid] = record

    # Use declared sequence-region length when available; fallback to observed max CDS end.
    contig_lengths = dict(contig_max_end)
    for seqid, length in contig_declared_len.items():
        if length > contig_lengths.get(seqid, 0):
            contig_lengths[seqid] = length

    if not contig_lengths:
        primary = ""
    else:
        if chromosome_like_contigs:
            primary = max(chromosome_like_contigs, key=lambda c: contig_lengths.get(c, 0))
        else:
            primary = max(contig_lengths, key=contig_lengths.get)
    return id_index, contig_lengths, primary


def build_sample_gff_map(gff_dir: Path, samples: List[str]) -> Dict[str, Path]:
    gff_files = list(gff_dir.rglob("*.gff")) + list(gff_dir.rglob("*.gff3"))
    sample_to_path: Dict[str, Path] = {}
    lowered = [(p, p.name.lower(), p.stem.lower()) for p in gff_files]

    for sample in samples:
        token = sample.lower()
        hits = []
        for p, lname, lstem in lowered:
            if token in lname or token == lstem:
                hits.append(p)
        if len(hits) == 1:
            sample_to_path[sample] = hits[0]
        elif len(hits) > 1:
            # Prefer exact stem match if multiple candidates exist.
            exact = [p for p in hits if p.stem.lower() == token]
            if len(exact) == 1:
                sample_to_path[sample] = exact[0]
    return sample_to_path


def parse_target_genes(target_genes_text: str) -> List[str]:
    genes = [g.strip() for g in target_genes_text.split(",") if g.strip()]
    if not genes:
        raise ValueError("--target-genes did not contain any row names.")
    return genes


def resolve_target_index_rows(index_values: List[str], target_genes: List[str]) -> Dict[str, str]:
    target_to_index: Dict[str, str] = {}
    for target in target_genes:
        exact_matches = [idx for idx in index_values if idx == target]
        if exact_matches:
            target_to_index[target] = exact_matches[0]
            continue

        prefix_matches = [idx for idx in index_values if idx.startswith(target)]
        if len(prefix_matches) == 1:
            target_to_index[target] = prefix_matches[0]
        elif len(prefix_matches) > 1:
            raise ValueError(
                f"Target row prefix {target!r} matched multiple CSV index rows: "
                f"{prefix_matches}"
            )
        else:
            raise ValueError(
                f"CSV index is missing target row name/prefix: {target!r}"
            )
    return target_to_index


def infer_sample_columns(
    all_columns: List[str], gff_dir: Path, sample_start_col: Optional[str]
) -> List[str]:
    if not all_columns:
        raise ValueError("No columns found in the CSV after reading the first column as index.")

    if sample_start_col is not None:
        token = str(sample_start_col).strip()
        if token in all_columns:
            start_idx = all_columns.index(token)
        else:
            try:
                start_idx = int(token)
            except ValueError as exc:
                raise ValueError(
                    f"--sample-start-col must be an existing column name or a 0-based "
                    f"column index, got: {sample_start_col}"
                ) from exc
            if start_idx < 0 or start_idx >= len(all_columns):
                raise ValueError(
                    f"--sample-start-col index {start_idx} is out of range for "
                    f"{len(all_columns)} CSV columns."
                )
        return all_columns[start_idx:]

    for idx, col in enumerate(all_columns):
        if build_sample_gff_map(gff_dir, [col]):
            return all_columns[idx:]

    raise ValueError(
        "Could not infer the first sample column from GFF file names. "
        "Please set --sample-start-col to the first sample column name or 0-based index."
    )


def canonicalize_order(ordered_genes: List[str]) -> Tuple[List[str], str]:
    """
    Treat forward and reverse as equivalent:
    choose lexicographically smaller one as canonical order.
    """
    forward = list(ordered_genes)
    reverse = list(reversed(ordered_genes))
    if tuple(forward) <= tuple(reverse):
        return forward, "forward_as_canonical"
    return reverse, "reverse_as_canonical"


def main() -> None:
    args = parse_args()
    csv_path = Path(args.csv)
    gff_dir = Path(args.gff_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    raw_df = pd.read_csv(csv_path, index_col=0)
    raw_df.index = raw_df.index.astype(str)

    target_genes = parse_target_genes(args.target_genes)
    target_to_index = resolve_target_index_rows(raw_df.index.tolist(), target_genes)

    all_value_cols = [str(c) for c in raw_df.columns]
    raw_df.columns = all_value_cols
    sample_cols = infer_sample_columns(all_value_cols, gff_dir, args.sample_start_col)
    if not sample_cols:
        raise ValueError("No sample columns found.")

    matched_index_rows = [target_to_index[g] for g in target_genes]
    df = raw_df.loc[matched_index_rows, sample_cols].copy()
    df.index = target_genes
    df.insert(0, "csv_index_row", matched_index_rows)
    df.insert(0, "Gene", df.index.astype(str))

    genes_in_order = df["Gene"].astype(str).tolist()
    row_index_by_gene = {g: i for i, g in enumerate(genes_in_order)}

    sample_to_gff = build_sample_gff_map(gff_dir, sample_cols)
    missing_gff_samples = sorted(set(sample_cols) - set(sample_to_gff))

    detailed_rows: List[dict] = []
    sample_pattern_rows: List[dict] = []
    pattern_aggregate: Dict[str, dict] = {}

    for sample in sample_cols:
        gff_path = sample_to_gff.get(sample)
        if gff_path is None:
            continue

        id_index, contig_lengths, primary_contig = load_gff_index(gff_path)
        hits_by_contig: Dict[str, List[dict]] = defaultdict(list)

        for _, row in df.iterrows():
            gene_name = str(row["Gene"])
            cell = row[sample]
            raw_ids = split_cell_gene_ids(cell)
            if not raw_ids:
                continue

            matched = None
            matched_source_id = None
            matched_candidate_id = None
            for raw_id in raw_ids:
                for candidate in clean_id_variants(raw_id):
                    matched = id_index.get(candidate)
                    if matched:
                        matched_source_id = raw_id
                        matched_candidate_id = candidate
                        break
                if matched:
                    break

            detail = {
                "sample": sample,
                "gene": gene_name,
                "gene_row_index": row_index_by_gene.get(gene_name, -1),
                "csv_index_row": row["csv_index_row"],
                "gene_id_from_csv": str(cell).strip(),
                "matched_source_id": matched_source_id,
                "matched_candidate_id": matched_candidate_id,
                "gff_file": str(gff_path),
                "found_in_gff": bool(matched),
                "contig": matched["contig"] if matched else None,
                "start": matched["start"] if matched else None,
                "end": matched["end"] if matched else None,
                "strand": matched["strand"] if matched else None,
                "feature_order": matched.get("feature_order") if matched else None,
                "cds_order": matched.get("cds_order") if matched else None,
                "is_primary_contig": bool(matched and matched["contig"] == primary_contig),
            }
            detailed_rows.append(detail)

            if matched:
                hits_by_contig[matched["contig"]].append(
                    {
                        "gene": gene_name,
                        "gene_row_index": row_index_by_gene.get(gene_name, -1),
                        "start": matched["start"],
                        "end": matched["end"],
                        "strand": matched["strand"],
                        "matched_candidate_id": matched_candidate_id,
                        "feature_order": matched.get("feature_order"),
                        "cds_order": matched.get("cds_order"),
                    }
                )

        for contig, genes in hits_by_contig.items():
            if len(genes) < args.min_genes_per_pattern:
                continue
            genes_sorted = sorted(genes, key=lambda x: (x["start"], x["end"], x["gene_row_index"]))
            genomic_order = [g["gene"] for g in genes_sorted]
            row_order_index = [g["gene_row_index"] for g in genes_sorted]
            strands = Counter([g["strand"] for g in genes_sorted if g["strand"] in {"+", "-"}])

            # Keep genomic order by coordinate, then merge reverse-equivalent patterns.
            canonical_order, canonical_direction = canonicalize_order(genomic_order)
            normalized_order = canonical_order
            block_direction = "contig_coordinate_ascending"

            # Final filter rule: only count patterns whose matched CDS rows are tightly adjacent
            # in Panaroo GFF (same contig, no CDS inserted between them).
            cds_order_sequence = [g.get("cds_order") for g in genes_sorted]
            feature_order_sequence = [g.get("feature_order") for g in genes_sorted]
            adjacency_mode = "cds_order"
            adjacency_sequence = cds_order_sequence
            if any(v is None for v in adjacency_sequence):
                adjacency_mode = "feature_order_fallback"
                adjacency_sequence = feature_order_sequence

            if not is_strictly_consecutive(adjacency_sequence):
                continue

            # Keep numeric IDs as a diagnostic field only (not used as a hard filter).
            gff_numeric_order = [
                extract_numeric_part(g.get("matched_candidate_id")) for g in genes_sorted
            ]

            pattern_id = "|".join(normalized_order)
            pattern_key = f"{len(normalized_order)}::{pattern_id}"

            sample_pattern = {
                "sample": sample,
                "gff_file": str(gff_path),
                "contig": contig,
                "is_primary_contig": contig == primary_contig,
                "contig_length_estimate": contig_lengths.get(contig),
                "n_genes": len(genomic_order),
                "genomic_order": " -> ".join(genomic_order),
                "row_order_index": ",".join(map(str, row_order_index)),
                "strand_counts": json.dumps(strands, ensure_ascii=False),
                "block_direction": block_direction,
                "adjacency_mode": adjacency_mode,
                "adjacency_sequence": ",".join(map(str, adjacency_sequence)),
                "feature_order_sequence": ",".join(map(str, feature_order_sequence)),
                "cds_order_sequence": ",".join(map(str, cds_order_sequence)),
                "gff_numeric_order": ",".join(map(str, gff_numeric_order)),
                "canonical_direction": canonical_direction,
                "canonical_order": " -> ".join(canonical_order),
                "normalized_order": " -> ".join(normalized_order),
                "pattern_key": pattern_key,
            }
            sample_pattern_rows.append(sample_pattern)

            if pattern_key not in pattern_aggregate:
                pattern_aggregate[pattern_key] = {
                    "pattern_key": pattern_key,
                    "n_genes": len(normalized_order),
                    "canonical_order": " -> ".join(normalized_order),
                    "normalized_order": " -> ".join(normalized_order),
                    "total_occurrences": 0,
                    "primary_contig_occurrences": 0,
                    "samples": [],
                    "sample_contigs": [],
                }

            agg = pattern_aggregate[pattern_key]
            agg["total_occurrences"] += 1
            if contig == primary_contig:
                agg["primary_contig_occurrences"] += 1
            agg["samples"].append(sample)
            agg["sample_contigs"].append(f"{sample}:{contig}")

    detailed_df = pd.DataFrame(detailed_rows)
    sample_pattern_df = pd.DataFrame(sample_pattern_rows)

    agg_rows = []
    pattern_summary_cols = [
        "pattern_key",
        "n_genes",
        "canonical_order",
        "normalized_order",
        "total_occurrences",
        "primary_contig_occurrences",
        "samples",
        "sample_contigs",
    ]
    for _, v in pattern_aggregate.items():
        agg_rows.append(
            {
                "pattern_key": v["pattern_key"],
                "n_genes": v["n_genes"],
                "canonical_order": v["canonical_order"],
                "normalized_order": v["normalized_order"],
                "total_occurrences": v["total_occurrences"],
                "primary_contig_occurrences": v["primary_contig_occurrences"],
                "samples": ";".join(v["samples"]),
                "sample_contigs": ";".join(v["sample_contigs"]),
            }
        )
    pattern_summary_df = pd.DataFrame(agg_rows, columns=pattern_summary_cols)
    if not pattern_summary_df.empty:
        pattern_summary_df = pattern_summary_df.sort_values(
            by=["total_occurrences", "primary_contig_occurrences", "n_genes"],
            ascending=[False, False, False],
        )

    detailed_path = out_dir / "gene_position_details.tsv"
    sample_pattern_path = out_dir / "sample_contig_patterns.tsv"
    pattern_summary_path = out_dir / "pattern_summary.tsv"
    meta_path = out_dir / "run_metadata.json"

    detailed_df.to_csv(detailed_path, sep="\t", index=False)
    sample_pattern_df.to_csv(sample_pattern_path, sep="\t", index=False)
    pattern_summary_df.to_csv(pattern_summary_path, sep="\t", index=False)

    metadata = {
        "csv": str(csv_path),
        "gff_dir": str(gff_dir),
        "target_genes": target_genes,
        "target_gene_to_csv_index": target_to_index,
        "sample_start_col": args.sample_start_col,
        "sample_columns": sample_cols,
        "direction_equivalent_mode": True,
        "adjacent_gff_cds_order_rule": True,
        "adjacent_gff_numeric_rule": False,
        "n_total_rows_in_csv": int(raw_df.shape[0]),
        "n_selected_target_rows": int(df.shape[0]),
        "n_sample_columns": len(sample_cols),
        "n_samples_with_matched_gff": len(sample_to_gff),
        "missing_gff_samples": missing_gff_samples,
        "outputs": {
            "gene_position_details": str(detailed_path),
            "sample_contig_patterns": str(sample_pattern_path),
            "pattern_summary": str(pattern_summary_path),
        },
    }
    meta_path.write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print("Done.")
    print(f"Selected target rows: {len(target_genes)}")
    print(f"Sample columns: {len(sample_cols)}")
    print(f"Matched GFF files: {len(sample_to_gff)}")
    print(f"Patterns found: {len(pattern_summary_df)}")
    print(f"Outputs written to: {out_dir}")


if __name__ == "__main__":
    main()
