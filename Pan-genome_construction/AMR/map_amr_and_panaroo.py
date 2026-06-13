import pandas as pd
import json

# Read CSV
df = pd.read_csv("./both-align-results-strict-adv/amr_gene_presence_absence.csv")

# Add a prefix before AMR name
df["AMR name"] = "amr_" + df["AMR name"].astype(str)

# Check whether AMR names are duplicated
duplicates = df["AMR name"][df["AMR name"].duplicated()].unique()
if len(duplicates) > 0:
    print("⚠️ Found duplicated AMR names:")
    for name in duplicates:
        print("  ", name)
else:
    print("✅ No duplicated AMR names found")

# Build the dictionary
gene_dict = dict(zip(df["Gene"], df["AMR name"]))

# Save as a JSON file
with open("./both-align-results-strict-adv/amr_panaroo_dict.json", "w") as f:
    json.dump(gene_dict, f, indent=4)

print("Dictionary saved to ./both-align-results-strict-adv/amr_panaroo_dict.json")
