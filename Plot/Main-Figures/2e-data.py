import os
import gzip
import sys

def is_vcf(filename):
    return filename.endswith(".vcf") or filename.endswith(".vcf.gz")

def open_vcf(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')

def classify_variant(ref, alt):
    # 结构变异，如 <DEL>、<INS>
    if alt.startswith('<') and alt.endswith('>'):
        return "SV"
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    elif len(ref) < len(alt):
        return "Insertion"
    elif len(ref) > len(alt):
        return "Deletion"
    elif len(ref) == len(alt) and len(ref) > 1:
        return "Complex"
    else:
        return None

def count_variants_in_directory(directory):

    # check
    outliers = []
    outliers_ids = []

    # 初始化变异类型计数器：每种变异类型都有 both, GN_only, GC_only 三种分类
    variant_counts = {
        "SNP": {"both": 0, "GN_only": 0, "GC_only": 0},
        "Insertion": {"both": 0, "GN_only": 0, "GC_only": 0},
        "Deletion": {"both": 0, "GN_only": 0, "GC_only": 0},
        "Complex": {"both": 0, "GN_only": 0, "GC_only": 0},
        "SV": {"both": 0, "GN_only": 0, "GC_only": 0}
    }
    
    for filename in os.listdir(directory):
        print(f"Processing file: {filename}")
        if not is_vcf(filename):
            continue
        filepath = os.path.join(directory, filename)
        with open_vcf(filepath) as f:
            gn_indices = []  # GN开头样本的列索引（基于完整fields）
            gc_indices = []  # GC开头样本的列索引（基于完整fields）
            header_parsed = False
            
            for line in f:
                # 解析header行，获取样本名称和对应的列索引
                if line.startswith('#CHROM'):
                    fields = line.strip().split('\t')
                    # VCF标准格式：#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [样本1 样本2 ...]
                    # 样本从第10列开始（0-indexed为第9列）
                    for col_idx in range(9, len(fields)):
                        sample_name = fields[col_idx]
                        # 从完整路径中提取样本ID前缀（GCF_xxxxx 或 GN开头）
                        # 例如: ./aaaT~~~aaaT_2~~~aaaT_1/GCF_001701365-1464_8_79.fa.bam
                        # 提取: GCF_001701365
                        sample_basename = sample_name.split('/')[-1]  # 获取文件名
                        sample_id = sample_basename.split('-')[0]  # 获取 GCF_xxxxx 或 GN
                        #print(sample_id)
                        
                        if sample_id.startswith('GN') or sample_id.startswith('_R_GN'):
                            gn_indices.append(col_idx)
                        elif sample_id.startswith('GC') or sample_id.startswith('_R_GC'):
                            gc_indices.append(col_idx)
                        else:
                            outliers_ids.append(sample_id)
                    
                    header_parsed = True
                    continue
                
                if header_parsed and len(gn_indices) == 0 and len(gc_indices) == 0:
                    print(f"Warning: No GN or GC samples found in file {filename}")                
                    sys.exit(0)
                    sys.exit("error message")
                
                if line.startswith('#'):
                    continue
                
                # 必须先有header信息才能处理数据行
                if not header_parsed:
                    print("Error: VCF header not found before data lines.")
                    continue
                
                fields = line.strip().split('\t')
                ref = fields[3]
                alt = fields[4]
                
                # 判断某个样本是否有该变异（不是 .:.）
                def has_variant(col_idx):
                    if col_idx >= len(fields):
                        return False
                    return fields[col_idx] != ".:."
                
                vt = classify_variant(ref, alt)
                if vt is None:
                    continue
                
                # 检查GN和GC样本是否有该变异
                gn_has = any(has_variant(col_idx) for col_idx in gn_indices)
                gc_has = any(has_variant(col_idx) for col_idx in gc_indices)
                
                # 统计到相应的分类
                if gn_has and gc_has:
                    variant_counts[vt]["both"] += 1
                elif gn_has:
                    variant_counts[vt]["GN_only"] += 1
                elif gc_has:
                    variant_counts[vt]["GC_only"] += 1
                else:
                    outliers.append([filename,fields[1],fields[3],fields[4]])  # 既不在GN也不在GC中出现的变异，视为异常情况

            print(len(gn_indices), len(gc_indices))

    for i in outliers:
        print(i)

    print(outliers_ids)
    
    return variant_counts

# 修改路径为你的 VCF 文件夹路径
vcf_dir = "./both-align-results-strict-adv/ann_vcf/"
variant_counts = count_variants_in_directory(vcf_dir)

print("🧬 Total variant counts across all VCF files:")
print("\n" + "="*60)
for vtype in ["SNP", "Insertion", "Deletion", "Complex", "SV"]:
    counts = variant_counts[vtype]
    total = counts["both"] + counts["GN_only"] + counts["GC_only"]
    print(f"\n{vtype}:")
    print(f"  总数: {total}")
    print(f"  两种样本都有: {counts['both']}")
    print(f"  仅GN样本有: {counts['GN_only']}")
    print(f"  仅GC样本有: {counts['GC_only']}")