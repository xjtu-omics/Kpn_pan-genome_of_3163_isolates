import os
import shutil

def extract_and_rename_files(source_dir, target_dir):
    try:
        # 确保目标文件夹存在
        os.makedirs(target_dir, exist_ok=True)

        # 遍历 source_dir 下的所有子文件夹
        for root, dirs, files in os.walk(source_dir):
            for dir_name in dirs:
                folder_path = os.path.join(root, dir_name)
                ref_file = os.path.join(folder_path, "ref.fa")
                
                # 检查 ref.fa 是否存在
                if os.path.exists(ref_file):
                    # 构造目标文件路径，以文件夹名为新文件名
                    target_file = os.path.join(target_dir, f"{dir_name}.fa")
                    
                    # 复制文件到目标文件夹
                    shutil.copy(ref_file, target_file)
                    print(f"复制成功: {ref_file} -> {target_file}")
                else:
                    print(f"未找到 ref.fa 文件: {folder_path}")
    except Exception as e:
        print(f"操作失败：{e}")

# 示例使用
source_directory = "/data/home/sfwang/kpn/Panaroo/both-align-results-strict-adv/alignment/"  # 替换为包含子文件夹的源目录
target_directory = "/data/home/sfwang/software/snpEff/data/genomes/"        # 替换为目标目录路径
extract_and_rename_files(source_directory, target_directory)
