import os
import json

def replace_contig_name(folder_path, output_file, name_dict):
    with open(output_file, 'w') as outfile:
        # Iterate over all files in the folder
        for filename in os.listdir(folder_path):
            if filename.endswith('.fa'):
                file_path = os.path.join(folder_path, filename)
                # Remove the file-name suffix, such as .fa
                contig_name = filename.replace('.fa', '')

                # Use the dictionary for the second replacement
                for key, value in name_dict.items():
                    if key==contig_name:
                        contig_name = value
                        break

                with open(file_path, 'r') as infile:
                    # Read the file content
                    lines = infile.readlines()

                    # Replace the contig name in the first row of the FASTA file
                    if lines[0].startswith('>'):
                        # Replace the contig name
                        lines[0] = f'>{contig_name}\n'

                    # Write the updated content
                    outfile.write(''.join(lines))  # Write the contig name and DNA sequence

def read_json_file(json_file_path):
    with open(json_file_path, 'r') as json_file:
        return json.load(json_file)

folder_path = "./both-align-results-strict-adv/core_gene_references/"  # Folder path
output_file = "./both-align-results-strict-adv/all_core_align_references.fasta"  # Merged FASTA file
json_dict = read_json_file("./both-align-results-strict-adv/amr_panaroo_dict.json")

replace_contig_name(folder_path, output_file, json_dict)

###Merge all FASTA files into one FASTA file