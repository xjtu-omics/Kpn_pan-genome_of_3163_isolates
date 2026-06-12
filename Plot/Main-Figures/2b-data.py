import os

def count_genes_in_gff(gff_file):
    count = 0
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            feature_type = parts[2]
            # 这里统计CDS数量，按需可修改
            if feature_type == "CDS":
                count += 1
    return count

def main(input_dir, output_file):
    gene_counts = []
    for filename in os.listdir(input_dir):
        if filename.endswith(".gff"):
            filepath = os.path.join(input_dir, filename)
            sample_id = os.path.splitext(filename)[0]
            count = count_genes_in_gff(filepath)
            gene_counts.append((sample_id, count))

    # 按基因数倒序排序
    gene_counts.sort(key=lambda x: x[1], reverse=True)

    # 写入输出文件
    with open(output_file, 'w') as out_f:
        for sample_id, count in gene_counts:
            out_f.write(f"{sample_id} {count}\n")

    print(f"统计完成，结果写入 {output_file}")

if __name__ == "__main__":
    # 你可以修改下面的路径
    input_directory = "/data/home/sfwang/kpn/GFFs/"  # GFF文件所在目录
    output_txt = "./sample_gene_counts.txt"
    main(input_directory, output_txt)