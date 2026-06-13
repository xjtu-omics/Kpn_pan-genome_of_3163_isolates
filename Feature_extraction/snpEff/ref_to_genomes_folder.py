import os
import shutil

def extract_and_rename_files(source_dir, target_dir):
    try:
        # Ensure the target folder exists
        os.makedirs(target_dir, exist_ok=True)

        # Iterate over all subfolders under source_dir
        for root, dirs, files in os.walk(source_dir):
            for dir_name in dirs:
                folder_path = os.path.join(root, dir_name)
                ref_file = os.path.join(folder_path, "ref.fa")

                # Check whether ref.fa exists
                if os.path.exists(ref_file):
                    # Build the target file path using the folder name as the new file name
                    target_file = os.path.join(target_dir, f"{dir_name}.fa")

                    # Copy the file to the target folder
                    shutil.copy(ref_file, target_file)
                    print(f"Copy successful: {ref_file} -> {target_file}")
                else:
                    print(f"ref.fa file not found: {folder_path}")
    except Exception as e:
        print(f"Operation failed:{e}")

# Example usage
source_directory = "./both-align-results-strict-adv/alignment/"  # Replace with the source directory containing subfolders
target_directory = "./snpEff/data/genomes/"        # Replace with the target directory path
extract_and_rename_files(source_directory, target_directory)
