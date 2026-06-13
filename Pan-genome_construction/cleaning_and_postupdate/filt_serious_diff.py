#!/usr/bin/env python3
import sys
import pandas as pd

def filter_serious_discrepancies(input_file, output_file):
    # Read the table
    df = pd.read_csv(input_file, sep="\t")

    # Convert to numeric values; keep NA as NaN
    df["Reference_length"] = pd.to_numeric(df["Reference_length"], errors="coerce")
    df["Mode_length"] = pd.to_numeric(df["Mode_length"], errors="coerce")

    # Compute differences
    df["Abs_diff"] = (df["Mode_length"] - df["Reference_length"]).abs()
    df["Rel_diff"] = df["Abs_diff"] / df["Reference_length"]

    # Filter entries with serious differences: >100 bp and >20%
    serious = df[(df["Abs_diff"] > 100) & (df["Rel_diff"] > 0.2)]

    # Save results
    serious.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <msa_length_summary.tsv> <serious_diff.tsv>")
        sys.exit(1)
    filter_serious_discrepancies(sys.argv[1], sys.argv[2])
