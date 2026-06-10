from pathlib import Path

# ====== 你需要修改的参数 ======
file_list = "empty_files.txt"   # 文件列表
NEW_LINE = "	.	+	0	ID=UNKNOWNG_00001;Name=genex;db_xref=COG:COG1289;gene=genex;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P46481;locus_tag=UNKNOWNG_00001;product=unknown"
# =================================

second_lines = {}  # 用来保存每个文件的第二行

with open(file_list, "r") as f:
    for line in f:
        filepath = Path(line.strip())
        filename = (str(filepath).split('/')[-1]).split('.')[0]
        new_filepath = './data/' + filename + '/genes.gff'
        if not filepath.exists():
            continue

        # 读取原文件内容
        with open(new_filepath, "r") as rf:
            lines = rf.readlines()

        # 至少需要有两行
        if len(lines) < 2:
            continue

        # 1️⃣ 读取第二行（索引 1）
        second_lines[str(new_filepath)] = lines[1].rstrip("\n")
        second_lines_info=second_lines[str(new_filepath)].split(' ')

        complete_new_line = second_lines_info[1] + '	Prodigal:002006	CDS	1	' + second_lines_info[3] + NEW_LINE

        # 2️⃣ 在第三行前插入新行（索引 2）
        if len(lines) >= 2:
            lines.insert(2, complete_new_line+"\n")

        # 3️⃣ 写回文件
        with open(new_filepath, "w") as wf:
            wf.writelines(lines)