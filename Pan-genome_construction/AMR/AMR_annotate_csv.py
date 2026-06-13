import pandas as pd

# Input file
panaroo_csv = "./both-align-results-strict-adv/gene_presence_absence.csv"
abricate_files = {
    "CARD": "./AMR_annotation/both_card_results_longCentroidID.tab",
    "NCBI": "./AMR_annotation/both_ncbi_results_longCentroidID.tab",
    "ResFinder": "./AMR_annotation/both_resfinder_results_longCentroidID.tab"
}
#output_all = "./both-align-results-strict/panaroo_with_resistance_all.csv"
output_filtered = "./both-align-results-strict-adv/amr_gene_presence_absence_longCentroidID.csv"

# Read the panaroo output and use low_memory=False to avoid dtype warnings
df = pd.read_csv(panaroo_csv, low_memory=False)

# Process each abricate file in turn
for db_name, ab_file in abricate_files.items():
    print(f"Processing {db_name} database: {ab_file}")

    # Read the abricate file
    ab_df = pd.read_csv(ab_file, sep="\t", header=None)

    # Take column 2 (gene name) and column 6 (resistance gene name)
    ab_map = dict(zip(ab_df.iloc[:, 1], ab_df.iloc[:, 5]))
    print(ab_map)

    # Add a new column to the panaroo table
    df[db_name] = df["Gene"].map(ab_map)

# -------- Adjust column order: move the three AMR columns to the front --------
amr_cols = list(abricate_files.keys())  # ["CARD", "NCBI", "ResFinder"]
other_cols = [c for c in df.columns if c not in amr_cols]
df = df[amr_cols + other_cols]

# Save the full table with resistance annotations
#df.to_csv(output_all, index=False)

# Keep rows with a hit in at least one database
filtered_df = df[df[list(abricate_files.keys())].notna().any(axis=1)]

# Save the filtered table
filtered_df.to_csv(output_filtered, index=False)

print(f"Processing complete.\nFiltered table: {output_filtered}")
