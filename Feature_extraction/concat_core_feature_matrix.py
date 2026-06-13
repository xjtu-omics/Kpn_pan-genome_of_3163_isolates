import pandas as pd
import os

# 读耐药表型
pheno_csvpath="./Panaroo-DownStream-both/phenotypes.csv"
phenos=pd.read_csv(pheno_csvpath,index_col=0)
GN_list=phenos.index.tolist()
bigFM=pd.DataFrame(index=GN_list)

def add_column_genename(fm,gene):
    c=fm.columns.tolist()
    for i in range(0,len(c)):
        c[i]=gene+'_'+c[i]
    fm.columns=c
    return(fm)

if __name__=='__main__':

    folder_path="./both-align-results-strict-adv/feature_matrix/"

    file_names=os.listdir(folder_path)
    for f in file_names:
        gene=(f.split('.'))[0]
        FM=pd.read_csv(folder_path+f, index_col=0)
        print(gene,len(FM.columns))
        if not FM.empty:
            nFM=add_column_genename(FM,gene)
            bigFM=pd.concat([bigFM,nFM],axis=1)
            print(len(bigFM.columns))

    bigFM.to_csv('./Panaroo-DownStream-both/final_core_feature_matrix.csv')