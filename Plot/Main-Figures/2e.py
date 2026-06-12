import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl

mpl.rcParams['font.family'] = 'Arial'

# 1. 整理数据
data = {
    'Category': ['SNP', 'Insertion', 'Deletion'],
    'Shared': [279061, 316, 266],
    'PATRIC Only': [144078, 180, 286],
    'Hui-net Only': [232212, 634, 1142],
}
df = pd.DataFrame(data).set_index('Category')

# 2. 创建画布：1行2列，宽度比例 1:2 (左边SNP，右边两个Indel)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.5, 6), gridspec_kw={'width_ratios': [1, 2]})
colors = ['#464555', '#aba9bc', '#d89b99'] # 共有, GN, GC

# 绘制左侧：SNP (大刻度)
df.loc[['SNP']].plot(kind='bar', stacked=True, ax=ax1, color=colors, legend=False, width=0.4)
ax1.set_title('SNP Count', fontsize=12)
ax1.set_ylabel('Number of Variants')
ax1.set_xticklabels(['SNP'], rotation=0)

# 绘制右侧：Insertion & Deletion (小刻度)
df.loc[['Insertion', 'Deletion']].plot(kind='bar', stacked=True, ax=ax2, color=colors, width=0.4)
ax2.set_title('Indel Counts', fontsize=12)
ax2.set_xticklabels(['Insertion', 'Deletion'], rotation=0)
ax2.legend(title="Type", loc='upper left')

# 在柱子顶部添加总数标签
for ax, cats in zip([ax1, ax2], [['SNP'], ['Insertion', 'Deletion']]):
    for i, cat in enumerate(cats):
        total = df.loc[cat].sum()
        ax.text(i, total + (total*0.005), f'{total}', ha='center', va='bottom')

plt.tight_layout()
plt.show()