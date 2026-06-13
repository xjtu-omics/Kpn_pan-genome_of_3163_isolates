###Assemble feature matrices for the newly selected feature sets
import os
import list2fm as lf
import pandas as pd

# Read the phenotype matrix

pheno_path="./PATRIC_data/Drug_phenotypes.csv"
phenos=pd.read_csv(pheno_path,index_col=0)
'''
phenos=pd.read_csv("./Origin_data/Kpn_antibiotics_AST_phenotypes.csv",index_col=0)
'''

# Loop over each antibiotic feature list, extract the corresponding feature columns from dis_fm, and combine them with labels to assemble the training matrix
feature_folder_path="./Panaroo-DownStream-patric/feature_set_determine/5f_avgauc_train/result/"
folder_list = [f for f in os.listdir(feature_folder_path)] # Antibiotic list for the Anhui hospital samples
#folder_list = ['ETP'] # Run ETP only

for f in folder_list:

    # Read the original full feature matrix (100% samples)
    #whole_core_fm=pd.read_csv('./Core/filted_whole_fm.csv',index_col=0)
    #whole_dis_fm=pd.read_csv('./Dispensable/filted_whole_fm.csv',index_col=0)

    # Read the core and dispensable feature matrices for antibiotic f (30%/70% samples)
    dataset="train"
    whole_fm_path=f"./Panaroo-DownStream-patric/ml_sort/result/{f}/"
    whole_core_fm=pd.read_csv(whole_fm_path+'core_'+dataset+'_set.csv',index_col=0)
    whole_dis_fm=pd.read_csv(whole_fm_path+'dis_'+dataset+'_set.csv',index_col=0)

    core_list=[]
    dis_list=[]
    # Read the feature file
    with open(f'./Panaroo-DownStream-patric/feature_set_determine/5f_avgauc_train/result/{f}/feature_select_result.txt','r') as file:
        for row in file:
            if 'pvalue' in row:
                row_info=row.split(',')
                gene=row_info[0]
                if '->' not in row:
                    dis_list.append(gene) # Dispensable-gene list
                else:
                    core_list.append(gene) # Core-variant list
    core_fm = whole_core_fm[core_list].copy()
    dispensable_fm = whole_dis_fm[dis_list].copy()
    total_fm = pd.concat([core_fm, dispensable_fm], axis=1) # Merge
    total_fm['label'] = phenos[f].reindex(total_fm.index)
    total_fm = total_fm[total_fm['label'].isin(['R', 'S'])]
    total_fm['label'] = total_fm['label'].replace({'R': 1, 'S': 0})

    total_fm.to_csv(f"./Panaroo-DownStream-patric/feature_set_determine/5f_avgauc_train/result_{dataset}_matrix/{f}.csv")