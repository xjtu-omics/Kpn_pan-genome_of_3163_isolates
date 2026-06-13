# 找卡方检验和k-w检验同时具有显著性的（特征，药物）对
# 数据存储：初步拟定为字典，key为特征，value为药物list
# 大概思路：同时遍历chi和kw的结果文件，每次读两个df（同一基因），找到既在chi中显著又在kw中显著的特征-药物对，存入字典
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
                    # 取chi和k-w交集/并集
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