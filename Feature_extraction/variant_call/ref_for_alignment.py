import os
import shutil

# Input source folder and target folder paths
src_dir = "./both-align-results-strict-adv/core_gene_references/"           # Source folder
dst_dir = "./both-align-results-strict-adv/alignment/"           # Target folder

# Ensure the target folder exists
os.makedirs(dst_dir, exist_ok=True)

for filename in os.listdir(src_dir):
    if filename.endswith(".fa"):
        # Remove the suffix
        name_without_suffix = filename[:-len(".fa")]

        # Create the corresponding subfolder
        subfolder = os.path.join(dst_dir, name_without_suffix)
        os.makedirs(subfolder, exist_ok=True)

        # Source file path
        src_path = os.path.join(src_dir, filename)

        # Target file path, renamed to ref.fa
        dst_path = os.path.join(subfolder, "ref.fa")

        # Copy file
        shutil.copy(src_path, dst_path)

print("Reference sequences for alignment have been distributed.")
