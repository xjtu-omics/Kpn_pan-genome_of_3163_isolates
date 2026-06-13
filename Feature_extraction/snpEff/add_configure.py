import os

def get_file_names_without_extension(folder_path):
    file_names = []
    # Iterate over all files in the folder
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            # Separate file name and extension
            name, _ = os.path.splitext(file)
            file_names.append(name)
    return file_names

if __name__=='__main__':
    folder_path = "./both-align-results-strict-adv/core_gene_references/"
    file_names = get_file_names_without_extension(folder_path) # Get the list of file names without extensions (all gene names)
    print(len(file_names))

    snp_configure_path="./snpEff/snpEff.config"
    with open(snp_configure_path,'a') as wfile:
        for gene in file_names:
            wfile.write('\n'+gene+'.genome : '+gene)
            wfile.write('\n'+gene+'.reference : ./snpEff/data/genomes/')
            wfile.write('\n'+gene+'.retrieval_date : 2025.10.05')
            wfile.write('\n')
