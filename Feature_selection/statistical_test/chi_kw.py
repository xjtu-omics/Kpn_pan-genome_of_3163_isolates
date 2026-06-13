# Find feature-antibiotic pairs significant in both chi-square and K-W tests
# Data storage: tentatively use a dictionary with features as keys and antibiotic lists as values
# General approach: iterate through chi and KW result files together, read two dataframes for the same gene each time, find feature-antibiotic pairs significant in both chi and KW, and store them in a dictionary
import pandas as pd
import numpy as np
import json
import os

result_dict={}

pheno_path="./Panaroo-DownStream-both/phenotypes.csv"
phenos=pd.read_csv(pheno_path,index_col=0)
med_list=phenos.columns.tolist()
for med in med_list:
    result_dict[med]=[]

def select_significance(df_chi,df_kw,gene):

    for index,row in df_chi.iterrows():
        for column,value in row.items():
            if index in df_kw.index and column in df_kw.columns:
                if not np.isnan(df_kw.at[index,column]):
                    # Take the intersection/union of chi and K-W results
                    if value < 0.05 and df_kw.at[index,column] < 0.05:
                        print(gene+index)
                        result_dict[column].append(gene+index)

if __name__ == "__main__":

    chi_file_path="./Panaroo-DownStream-both/statistics-test/dispensable_chisquare.csv"
    kw_file_path="./Panaroo-DownStream-both/statistics-test/dispensable_kw.csv"

    if os.path.exists(chi_file_path) and os.path.exists(kw_file_path):
        df_chi=pd.read_csv(chi_file_path,index_col=0)
        df_kw=pd.read_csv(kw_file_path,index_col=0)
        select_significance(df_chi,df_kw,'')
    with open("./Panaroo-DownStream-both/statistics-test/chi_and_kw_dispensable.json", 'w') as f:
        json.dump(result_dict, f, indent=4)