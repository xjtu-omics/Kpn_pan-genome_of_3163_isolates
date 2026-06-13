#!/usr/bin/env python3
import os
import sys
import shutil
import pandas as pd

def move_serious_msa(serious_file, msa_dir, out_dir):
    # Read the table of serious differences
    df = pd.read_csv(serious_file, sep="\t")

    # Create the output directory
    os.makedirs(out_dir, exist_ok=True)

    # Iterate over MSA file names in the table
    for msa_file in df["MSA_file"]:
        src = os.path.join(msa_dir, msa_file)
        dst = os.path.join(out_dir, msa_file)

        if os.path.exists(src):
            shutil.move(src, dst)
            print(f"Moved: {msa_file}")
        else:
            print(f"Warning: {msa_file} not found in {msa_dir}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <serious_diff.tsv> <msa_dir> <output_dir>")
        sys.exit(1)
    move_serious_msa(sys.argv[1], sys.argv[2], sys.argv[3])