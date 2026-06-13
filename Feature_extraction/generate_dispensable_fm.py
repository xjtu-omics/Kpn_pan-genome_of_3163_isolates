import os
import re
import json
import pandas as pd
from Bio import SeqIO

# Input path
input_dir = "./both-align-results-strict-adv/dispensable_gene_sequences/"   # Change this to your folder path
mapping_json = "./both-align-results-strict-adv/amr_panaroo_dict.json"  # JSON file storing the column-name mapping
output_csv = "./Panaroo-DownStream-both/final_dispensable_fm.csv"

# Use regex to extract the group name: from G to the first semicolon
sample_pattern = re.compile(r"^(G[^;]+)")

gene_to_samples = {}

# Iterate over all .aln.fas files
for filename in os.listdir(input_dir):
    if filename.endswith(".aln.fas"):
        gene_name = filename.replace(".aln.fas", "")
        file_path = os.path.join(input_dir, filename)

        samples = set()
        for record in SeqIO.parse(file_path, "fasta"):
            match = sample_pattern.search(record.id)
            if match:
                samples.add(match.group(1))

        gene_to_samples[gene_name] = samples

# Collect all samples
all_samples = sorted({s for samples in gene_to_samples.values() for s in samples})

# Build the 0/1 matrix
df = pd.DataFrame(0, index=all_samples, columns=gene_to_samples.keys())

for gene, samples in gene_to_samples.items():
    df.loc[list(samples), gene] = 1

# Read the mapping and replace column names
if os.path.exists(mapping_json):
    with open(mapping_json, "r") as f:
        mapping = json.load(f)
    df.rename(columns=mapping, inplace=True)

# Check that each column has at least four 1s or at least four 0s, and filter rare features
rare_threshold = 3
bad_cols = [col for col in df.columns if df[col].sum() <= rare_threshold or df[col].sum() >= 3163-rare_threshold]
if bad_cols:
    print("The following columns have fewer than four 1s or fewer than four 0s; please check:")
    for col in bad_cols:
        print(f"{col}: {df[col].sum()}  ones")

# Save as CSV
df.to_csv(output_csv)

print(f"Results saved to {output_csv}")
