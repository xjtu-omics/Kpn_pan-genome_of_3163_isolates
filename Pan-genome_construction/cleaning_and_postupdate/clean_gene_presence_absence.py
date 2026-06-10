from pathlib import Path
import argparse
import csv
from collections import defaultdict

DEFAULT_INPUT_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT = "gene_contig_matrix.tsv"


def normalize_sample_id(sample_id):
    """
    Keep the same sample ID normalization rule as the original script.
    Example:
        _R_GNxxxx -> GNxxxx
    """
    if sample_id.startswith("_R_"):
        return sample_id[3:]
    return sample_id


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate a gene-by-sample matrix from .fas files. "
            "Rows are genes, columns are samples, and each cell contains "
            "the full contig/header names where the gene is present."
        )
    )
    parser.add_argument(
        "-i",
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"Directory containing .fas files (default: {DEFAULT_INPUT_DIR})",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path(DEFAULT_OUTPUT),
        help=f"Output matrix table path (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args()


def parse_fas_file(fas_path):
    """
    Parse one .fas file.

    Returns
    -------
    gene_name : str
        File name without suffix.
    sample_to_contigs : dict[str, list[str]]
        Mapping from normalized sample ID to all full contig/header names
        observed for this gene in that sample.

    Notes
    -----
    The full contig/header name is the complete FASTA header after removing
    the leading '>'.
    """
    gene_name = fas_path.stem
    sample_to_contigs = defaultdict(list)
    seen_pairs = set()

    with fas_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or not line.startswith(">"):
                continue

            full_contig_name = line[1:].strip()
            if not full_contig_name:
                continue

            sample_id = normalize_sample_id(
                full_contig_name.split(";", 1)[0].strip()
            )

            if not sample_id:
                continue

            # Avoid duplicated headers within the same sample/gene cell.
            pair = (sample_id, full_contig_name)
            if pair not in seen_pairs:
                sample_to_contigs[sample_id].append(full_contig_name)
                seen_pairs.add(pair)

    if not sample_to_contigs:
        raise ValueError(f"No FASTA headers found in {fas_path.name}")

    return gene_name, sample_to_contigs


def collect_gene_sample_matrix(input_dir):
    fas_files = sorted(input_dir.glob("*.fas"))
    if not fas_files:
        raise FileNotFoundError(f"No .fas files found in {input_dir}")

    gene_to_sample_contigs = {}
    all_samples = set()
    errors = []

    for fas_path in fas_files:
        try:
            gene_name, sample_to_contigs = parse_fas_file(fas_path)
            gene_to_sample_contigs[gene_name] = sample_to_contigs
            all_samples.update(sample_to_contigs.keys())
        except ValueError as exc:
            errors.append(f"{fas_path.name}: {exc}")

    if errors:
        error_text = "\n".join(errors)
        raise ValueError(
            "Some .fas files could not be parsed:\n"
            f"{error_text}"
        )

    return gene_to_sample_contigs, sorted(all_samples)


def write_gene_sample_matrix(gene_to_sample_contigs, samples, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)

    genes = sorted(gene_to_sample_contigs.keys())

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")

        # First column is the row index name.
        writer.writerow(["gene"] + samples)

        for gene in genes:
            row = [gene]
            sample_to_contigs = gene_to_sample_contigs[gene]

            for sample in samples:
                contigs = sample_to_contigs.get(sample, [])
                cell_value = "；".join(contigs) if contigs else ""
                row.append(cell_value)

            writer.writerow(row)


def main():
    args = parse_args()

    gene_to_sample_contigs, samples = collect_gene_sample_matrix(args.input_dir)
    write_gene_sample_matrix(gene_to_sample_contigs, samples, args.output)

    print(f"Output written to: {args.output}")
    print(f"Genes: {len(gene_to_sample_contigs)}")
    print(f"Samples: {len(samples)}")


if __name__ == "__main__":
    main()