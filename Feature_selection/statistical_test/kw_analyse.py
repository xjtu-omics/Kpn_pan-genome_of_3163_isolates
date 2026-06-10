import pandas as pd
from scipy.stats import kruskal
from toolkit import merge_label_to_fm

# 读RS耐药表型数据
pheno_path="/data/home/sfwang/kpn/Panaroo-DownStream-both/phenotypes.csv"
phenos=pd.read_csv(pheno_path,index_col=0)

fm_path="/data/home/sfwang/kpn/Panaroo-DownStream-both/final_core_feature_matrix.csv"
fm=pd.read_csv(fm_path,index_col=0)

# 移除方差为0的特征（所有值相同）
fm = fm.loc[:, fm.nunique() > 1]

med_list=phenos.columns.tolist()
variants=fm.columns.tolist()
result=pd.DataFrame(columns=med_list,index=variants)

# 获取标签
#med='AMP'
for med in med_list:

    RS_data=phenos[med]
    # 合并标签与大变异矩阵
    df=merge_label_to_fm(fm,RS_data)

    # 假设 df 是包含特征矩阵的 DataFrame，y 是目标变量
    X = df.drop(columns='label')  # 删除目标变量列，得到特征矩阵
    y = df['label']  # 目标变量

    # 移除所有列值相同的特征
    X = X.loc[:, X.nunique() > 1]  # 只保留有多个不同值的特征

    # 存储每个特征的 p 值
    p_values = []

    for feature in X.columns:
        
        # 按表型分组
        grouped = [X[feature][y == c] for c in y.unique()]

        # 只保留样本数至少 2 且组内有变化（非只有一种特征值）的组
        valid_groups = [g for g in grouped if len(g) > 1 and g.nunique() > 1]

        if len(valid_groups) >= 2:
            # 进行 Kruskal-Wallis 检验
            _, p_value = kruskal(*valid_groups)
        else:
            p_value = float('nan')

        p_values.append(p_value)
        result.loc[feature, med] = p_value

    # 将 p 值与特征名组合成 DataFrame
    p_values_df = pd.DataFrame({
        'Feature': X.columns,
        'P-Value': p_values
    })

result.to_csv("/data/home/sfwang/kpn/Panaroo-DownStream-both/statistics-test/core_kw.csv")