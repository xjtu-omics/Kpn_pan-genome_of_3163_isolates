import matplotlib.pyplot as plt
import scipy.stats as stats
import matplotlib as mpl

mpl.rcParams['font.family'] = 'Arial'

INPUT_FILE = "C:/Users/yxLig/Nutstore/1/我的坚果云 (2)/Klebsiella_pneumoniae/figures-new/fig2-中间文件/anhui_sample_gene_counts.txt"
OUTPUT_FILE = "gene_count_trapezoid_plot.pdf"

def read_data(filename):
    # 返回字典，key为样本名，value为N50值
    with open(filename, 'r') as f:
        lines = f.readlines()
    data = {}
    for line in lines:
        parts = line.strip().split()
        if len(parts) == 2:
            strain, n50 = parts[0], float(parts[1])
            data[strain]=n50
    return data

data1 = read_data("C:/Users/yxLig/Nutstore/1/我的坚果云 (2)/Klebsiella_pneumoniae/figures-new/fig2-中间文件/anhui_contig_n50_1174.txt")
data2 = read_data("C:/Users/yxLig/Nutstore/1/我的坚果云 (2)/Klebsiella_pneumoniae/figures-new/fig2-中间文件/anhui_contig_n50_491.txt")

# 读取并排序
data = []
with open(INPUT_FILE) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) != 2:
            continue
        sample, count = parts
        if sample in data1:
            try:
                data.append((sample, int(count),'Pure short-read',data1[sample]))
            except ValueError:
                continue
        elif sample in data2:
            try:
                data.append((sample, int(count),'Hybrid long+short-read',data2[sample]))
            except ValueError:
                continue

data_sorted = sorted(data, key=lambda x: x[3], reverse=True)

# 提取排序后数据
gene_count_1174 = [x[1] for x in data_sorted if x[2]=='Pure short-read']
gene_count_491 = [x[1] for x in data_sorted if x[2]=='Hybrid long+short-read']
print(len(gene_count_1174), len(gene_count_491))

import seaborn as sns
import pandas as pd
df = pd.DataFrame({
    'Value': gene_count_1174 + gene_count_491,  # 合并两个列表数据
    'List': ['Pure short-read'] * len(gene_count_1174) + ['Hybrid long+short-read'] * len(gene_count_491)  # 标记数据来源
})

# 计算 Wilcoxon rank-sum (Mann-Whitney U) 检验的 p-value
u_stat, p_value = stats.mannwhitneyu(gene_count_491, gene_count_1174, alternative='greater')
print(f'Mann-Whitney U Statistic: {u_stat}')
print(f'Wilcoxon rank-sum test p-value: {p_value}')

# 创建图形
plt.figure(figsize=(4, 6))

# 使用seaborn绘制箱线图
sns.boxplot(x='List', 
            y='Value', 
            data=df,
            palette=["#91CAE8", "#F48892"],
            width=0.4,  # 设置箱体宽度
            flierprops=dict(marker='o', color='r', markersize=5, markeredgewidth=1))  # 设置散点样式

# 设置纵轴的上限
plt.ylim(4500, 7500) 

# 设置标题和标签
plt.title('Gene Numbers from Two Assembly Methods', fontsize=12)
plt.xlabel('Assembly Method', fontsize=10)
plt.ylabel('Gene Number', fontsize=10)

# 显示图形
plt.tight_layout()
plt.show()