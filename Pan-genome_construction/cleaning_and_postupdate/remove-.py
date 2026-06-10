from Bio import SeqIO
import os
import re

# 输入和输出文件夹路径
input_folder = "/data/home/sfwang/kpn/Panaroo/both-align-results-strict-adv/aligned_gene_sequences/"   # 改成你的原始msa fasta文件夹路径
output_folder = "/data/home/sfwang/kpn/Panaroo/both-align-results-strict-adv/gene_sequences/" # 改成你想保存结果的文件夹路径

# 如果输出文件夹不存在就创建
os.makedirs(output_folder, exist_ok=True)

# 遍历文件夹下所有 fasta 文件
for filename in os.listdir(input_folder):
    if filename.endswith(".fasta") or filename.endswith(".fas"):
        input_path = os.path.join(input_folder, filename)
        output_path = os.path.join(output_folder, filename)

        with open(output_path, "w") as out_handle:
            for record in SeqIO.parse(input_path, "fasta"):
                # 只保留 A, T, C, G
                clean_seq = re.sub(r"[^ATCGatcg]", "", str(record.seq))
                record.seq = clean_seq
                SeqIO.write(record, out_handle, "fasta")

print("所有 fasta 文件处理完成，结果保存在:", output_folder)
