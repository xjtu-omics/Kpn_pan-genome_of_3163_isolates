#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import os
import sqlite3
import tempfile
from pathlib import Path
from functools import lru_cache


def normalize_gff_sample_name(gff_file):
    """
    Generate possible sample names from the gff_file column in Panaroo gene_data.csv.

    Examples:
        /path/to/GCA_022067925.gff      -> GCA_022067925
        GCA_022067925.gff3             -> GCA_022067925
        GN191724.gff                   -> GN191724
        GCA_022067925                  -> GCA_022067925
    """
    raw = str(gff_file).strip()
    base = Path(raw).name

    candidates = set()
    candidates.add(raw)
    candidates.add(base)

    suffixes = [
        ".gff.gz",
        ".gff3.gz",
        ".gff",
        ".gff3",
        ".gbk",
        ".gbff",
    ]

    for suffix in suffixes:
        if base.endswith(suffix):
            candidates.add(base[: -len(suffix)])

    # General fallback
    candidates.add(Path(base).stem)

    return {x for x in candidates if x}


def normalize_input_sample_name(sample_name):
    """
    Normalize sample names from the input matrix cell before matching.

    Example:
        _R_GCA_000492295 -> GCA_000492295
    """
    sample_name = str(sample_name).strip()

    if sample_name.startswith("_R_"):
        sample_name = sample_name[3:]

    return sample_name


def connect_sqlite(db_path):
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA locking_mode=EXCLUSIVE")
    return conn


def init_mapping_db(conn):
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS mapping (
            full_key TEXT PRIMARY KEY,
            annotation_id TEXT NOT NULL
        )
        """
    )
    conn.commit()


def build_mapping_db(
    gene_data_csv,
    db_path,
    gene_data_delimiter=",",
    batch_size=100000,
    rebuild=False,
):
    """
    Stream gene_data.csv and build an on-disk SQLite mapping database.

    Mapping stored:
        sample_name;clustering_id -> annotation_id

    where sample_name is derived from the gff_file column.
    """
    if os.path.exists(db_path) and not rebuild:
        print(f"Using existing mapping database: {db_path}")
        return

    if os.path.exists(db_path) and rebuild:
        os.remove(db_path)

    print(f"Building mapping database: {db_path}")

    conn = connect_sqlite(db_path)
    init_mapping_db(conn)

    inserted_rows = 0
    processed_rows = 0
    batch = []

    with open(gene_data_csv, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter=gene_data_delimiter)

        if reader.fieldnames is None:
            raise ValueError("gene_data.csv appears to have no header.")

        fieldnames = [x.strip() for x in reader.fieldnames]
        fieldname_map = {x.strip(): x for x in reader.fieldnames}

        required_cols = ["gff_file", "clustering_id", "annotation_id"]
        missing = [col for col in required_cols if col not in fieldname_map]

        if missing:
            raise ValueError(
                f"gene_data.csv is missing required columns: {missing}\n"
                f"Observed columns: {fieldnames}"
            )

        gff_col = fieldname_map["gff_file"]
        clustering_col = fieldname_map["clustering_id"]
        annotation_col = fieldname_map["annotation_id"]

        for row in reader:
            processed_rows += 1

            gff_file = row.get(gff_col, "").strip()
            clustering_id = row.get(clustering_col, "").strip()
            annotation_id = row.get(annotation_col, "").strip()

            if not gff_file or not clustering_id or not annotation_id:
                continue

            sample_candidates = normalize_gff_sample_name(gff_file)

            for sample_name in sample_candidates:
                full_key = f"{sample_name};{clustering_id}"
                batch.append((full_key, annotation_id))

            if len(batch) >= batch_size:
                conn.executemany(
                    "INSERT OR IGNORE INTO mapping(full_key, annotation_id) VALUES (total, total)",
                    batch,
                )
                conn.commit()
                inserted_rows += len(batch)
                batch.clear()
                print(f"Processed gene_data rows: {processed_rows:,}")

        if batch:
            conn.executemany(
                "INSERT OR IGNORE INTO mapping(full_key, annotation_id) VALUES (total, total)",
                batch,
            )
            conn.commit()
            inserted_rows += len(batch)
            batch.clear()

    conn.execute("CREATE INDEX IF NOT EXISTS idx_mapping_key ON mapping(full_key)")
    conn.commit()

    unique_keys = conn.execute("SELECT COUNT(*) FROM mapping").fetchone()[0]
    conn.close()

    print(f"Finished building mapping database.")
    print(f"Processed gene_data rows: {processed_rows:,}")
    print(f"Candidate mapping entries inserted/ignored: {inserted_rows:,}")
    print(f"Unique mapping keys in database: {unique_keys:,}")


def make_lookup_function(conn, cache_size=500000):
    """
    Return a cached lookup function.

    It tries:
        1. exact full entry
        2. normalized sample name, especially removing leading _R_
        3. trimming trailing suffixes from clustering_id
    """

    @lru_cache(maxsize=cache_size)
    def query_db(full_key):
        cur = conn.execute(
            "SELECT annotation_id FROM mapping WHERE full_key = total LIMIT 1",
            (full_key,),
        )
        row = cur.fetchone()
        if row is None:
            return None
        return row[0]

    def lookup_annotation_id(full_entry):
        full_entry = str(full_entry).strip()

        if not full_entry:
            return ""

        # 1. Exact full entry
        hit = query_db(full_entry)
        if hit is not None:
            return hit

        if ";" not in full_entry:
            return full_entry

        # Split only by the first English semicolon.
        # Example:
        #   _R_GCA_000492295;2_refound_161
        sample_name, clustering_id = full_entry.split(";", 1)
        sample_name = sample_name.strip()
        clustering_id = clustering_id.strip()

        if not sample_name or not clustering_id:
            return full_entry

        normalized_sample_name = normalize_input_sample_name(sample_name)

        candidate_samples = []
        candidate_samples.append(sample_name)

        if normalized_sample_name != sample_name:
            candidate_samples.append(normalized_sample_name)

        # Avoid duplicates while preserving order
        seen_samples = set()
        candidate_samples = [
            x for x in candidate_samples
            if not (x in seen_samples or seen_samples.add(x))
        ]

        # 2. Try original clustering_id with candidate sample names
        for s in candidate_samples:
            candidate_key = f"{s};{clustering_id}"
            hit = query_db(candidate_key)
            if hit is not None:
                return hit

        # 3. Try trimming trailing suffixes from clustering_id
        #
        # Example:
        #   2825_0_4176_1 -> 2825_0_4176
        #
        # Exact match is always tried first, so this will not harm valid IDs like:
        #   2_refound_161
        tmp = clustering_id
        while "_" in tmp:
            tmp = tmp.rsplit("_", 1)[0]

            for s in candidate_samples:
                candidate_key = f"{s};{tmp}"
                hit = query_db(candidate_key)
                if hit is not None:
                    return hit

        return full_entry

    return lookup_annotation_id


def convert_cell(
    cell,
    lookup_annotation_id,
    input_entry_sep="；",
    output_entry_sep=";",
):
    """
    Convert one matrix cell.

    Input example:
        GCA_022067925;510_141_0；_R_GCA_000492295;2_refound_161

    Output example:
        IEIMDPGK_00002;XXXXXXX_00001
    """
    cell = str(cell).strip()

    if cell == "":
        return "", []

    entries = [x.strip() for x in cell.split(input_entry_sep) if x.strip()]

    converted_entries = []
    unmapped_entries = []

    for entry in entries:
        converted = lookup_annotation_id(entry)
        converted_entries.append(converted)

        if converted == entry:
            unmapped_entries.append(entry)

    return output_entry_sep.join(converted_entries), unmapped_entries


def convert_matrix_streaming(
    input_tsv,
    output_csv,
    db_path,
    input_delimiter="\t",
    output_delimiter=",",
    input_entry_sep="；",
    output_entry_sep=";",
    unmapped_output=None,
    cache_size=500000,
    progress_every=1000,
):
    """
    Stream input matrix row by row and write converted output row by row.
    """
    conn = sqlite3.connect(db_path)
    lookup_annotation_id = make_lookup_function(conn, cache_size=cache_size)

    total_rows = 0
    total_unmapped = 0

    unmapped_fh = None
    unmapped_writer = None

    if unmapped_output:
        unmapped_fh = open(unmapped_output, "w", encoding="utf-8", newline="")
        unmapped_writer = csv.writer(unmapped_fh, delimiter="\t")
        unmapped_writer.writerow(["row_id", "column_name", "unmapped_entry"])

    with open(input_tsv, "r", encoding="utf-8-sig", newline="") as fin, \
         open(output_csv, "w", encoding="utf-8", newline="") as fout:

        reader = csv.reader(fin, delimiter=input_delimiter)
        writer = csv.writer(fout, delimiter=output_delimiter)

        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("Input matrix is empty.")

        writer.writerow(header)

        if len(header) < 2:
            raise ValueError(
                "Input matrix should contain one index column plus sample columns."
            )

        sample_columns = header[1:]

        for row in reader:
            total_rows += 1

            if not row:
                continue

            row_id = row[0]

            # If some rows are shorter than the header, pad empty cells.
            if len(row) < len(header):
                row = row + [""] * (len(header) - len(row))

            # If some rows are longer than the header, keep extra columns but warn implicitly by processing them.
            out_row = [row_id]

            for col_idx, cell in enumerate(row[1:], start=1):
                column_name = (
                    header[col_idx]
                    if col_idx < len(header)
                    else f"extra_column_{col_idx}"
                )

                converted_cell, unmapped_entries = convert_cell(
                    cell,
                    lookup_annotation_id=lookup_annotation_id,
                    input_entry_sep=input_entry_sep,
                    output_entry_sep=output_entry_sep,
                )

                out_row.append(converted_cell)

                if unmapped_entries:
                    total_unmapped += len(unmapped_entries)

                    if unmapped_writer is not None:
                        for entry in unmapped_entries:
                            unmapped_writer.writerow([row_id, column_name, entry])

            writer.writerow(out_row)

            if progress_every and total_rows % progress_every == 0:
                print(f"Converted matrix rows: {total_rows:,}")

    if unmapped_fh is not None:
        unmapped_fh.close()

    conn.close()

    print("Finished matrix conversion.")
    print(f"Converted matrix rows: {total_rows:,}")
    print(f"Unmapped entries: {total_unmapped:,}")
    print(f"Output written to: {output_csv}")

    if unmapped_output:
        print(f"Unmapped entries written to: {unmapped_output}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Stream-convert a gene-by-sample TSV matrix from Panaroo internal IDs "
            "to Prokka annotation IDs using Panaroo gene_data.csv."
        )
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input TSV matrix. Rows are genes, columns are samples."
    )

    parser.add_argument(
        "-g", "--gene-data",
        required=True,
        help="Panaroo gene_data.csv."
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output CSV file."
    )

    parser.add_argument(
        "--db",
        default=None,
        help=(
            "SQLite mapping database path. "
            "If not provided, a database will be created next to the output file."
        )
    )

    parser.add_argument(
        "--rebuild-db",
        action="store_true",
        help="Force rebuilding the SQLite mapping database."
    )

    parser.add_argument(
        "--input-delimiter",
        default="\t",
        help="Input matrix delimiter. Default: tab."
    )

    parser.add_argument(
        "--output-delimiter",
        default=",",
        help="Output matrix delimiter. Default: comma."
    )

    parser.add_argument(
        "--gene-data-delimiter",
        default=",",
        help="Delimiter of gene_data.csv. Default: comma."
    )

    parser.add_argument(
        "--input-entry-sep",
        default="；",
        help="Delimiter between entries inside each input cell. Default: Chinese semicolon."
    )

    parser.add_argument(
        "--output-entry-sep",
        default=";",
        help="Delimiter between entries inside each output cell. Default: English semicolon."
    )

    parser.add_argument(
        "--unmapped-output",
        default=None,
        help="Optional TSV file recording unmapped entries."
    )

    parser.add_argument(
        "--batch-size",
        type=int,
        default=100000,
        help="Batch size for inserting mappings into SQLite. Default: 100000."
    )

    parser.add_argument(
        "--cache-size",
        type=int,
        default=500000,
        help="LRU cache size for mapping lookup. Default: 500000."
    )

    parser.add_argument(
        "--progress-every",
        type=int,
        default=1000,
        help="Print progress every N matrix rows. Default: 1000."
    )

    args = parser.parse_args()

    if args.db is None:
        output_path = Path(args.output)
        db_path = str(output_path.with_suffix(output_path.suffix + ".mapping.sqlite"))
    else:
        db_path = args.db

    build_mapping_db(
        gene_data_csv=args.gene_data,
        db_path=db_path,
        gene_data_delimiter=args.gene_data_delimiter,
        batch_size=args.batch_size,
        rebuild=args.rebuild_db,
    )

    convert_matrix_streaming(
        input_tsv=args.input,
        output_csv=args.output,
        db_path=db_path,
        input_delimiter=args.input_delimiter,
        output_delimiter=args.output_delimiter,
        input_entry_sep=args.input_entry_sep,
        output_entry_sep=args.output_entry_sep,
        unmapped_output=args.unmapped_output,
        cache_size=args.cache_size,
        progress_every=args.progress_every,
    )


if __name__ == "__main__":
    main()