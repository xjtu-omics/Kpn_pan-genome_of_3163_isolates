from Bio import SeqIO
import os
import re

# Input and output folder paths
input_folder = "./both-align-results-strict-adv/aligned_gene_sequences/"   # Change this to the folder path for your original MSA FASTA files
output_folder = "./both-align-results-strict-adv/gene_sequences/" # Change this to the folder path where you want to save results

# Create the output folder if it does not exist
os.makedirs(output_folder, exist_ok=True)

# Iterate over all FASTA files under the folder
for filename in os.listdir(input_folder):
    if filename.endswith(".fasta") or filename.endswith(".fas"):
        input_path = os.path.join(input_folder, filename)
        output_path = os.path.join(output_folder, filename)

        with open(output_path, "w") as out_handle:
            for record in SeqIO.parse(input_path, "fasta"):
                # Keep only A, T, C, and G
                clean_seq = re.sub(r"[^ATCGatcg]", "", str(record.seq))
                record.seq = clean_seq
                SeqIO.write(record, out_handle, "fasta")

print("All FASTA files processed; results saved in:", output_folder)
