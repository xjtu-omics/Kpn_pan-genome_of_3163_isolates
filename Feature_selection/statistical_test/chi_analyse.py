import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency

def cal_chi(pheno_list,feature_list):

    R0=0
    R1=0
    S0=0
    S1=0

    for index,value in pheno_list.items():
        if value=="R" and feature_list[index]==0:
            R0+=1
        elif value=="R" and feature_list[index]==1:
            R1+=1
        elif value=="S" and feature_list[index]==0:
            S0+=1
        elif value=="S" and feature_list[index]==1:
            S1+=1

    return R0,R1,S0,S1

pheno_path="./Panaroo-DownStream-both/phenotypes.csv"
pheno=pd.read_csv(pheno_path,index_col=0)
med_class=pheno.columns.tolist()

fm_path="./Panaroo-DownStream-both/final_core_feature_matrix.csv"
fm=pd.read_csv(fm_path,index_col=0)
variants=fm.columns.tolist()
df=pd.DataFrame(columns=med_class,index=variants)
for index,value in df.iterrows(): #Iterate over each variant
    for column,item in value.items(): #Iterate over each antibiotic
        pheno_list=pheno[column]
        feature_list=fm[index]
        R0,R1,S0,S1=cal_chi(pheno_list,feature_list)
        if 0 not in [R0,R1,S0,S1]:
            chi_arr=np.array([[R0,R1],[S0,S1]])
            df.at[index,column]=chi2_contingency(chi_arr).pvalue
        elif R0==0 or R1==0 or S0==0 or S1==0 or (R0+S1)==0 or (R1+S0)==0:
            R0+=0.5
            R1+=0.5
            S0+=0.5
            S1+=0.5
            chi_arr=np.array([[R0,R1],[S0,S1]])
            df.at[index,column]=chi2_contingency(chi_arr).pvalue
        else:
            df.at[index,column]=float('nan')
df.to_csv("./Panaroo-DownStream-both/statistics-test/core_chisquare.csv")