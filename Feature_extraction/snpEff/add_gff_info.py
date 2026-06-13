from pathlib import Path

# ====== Parameters you need to modify ======
file_list = "empty_files.txt"   # File list
NEW_LINE = "	.	+	0	ID=UNKNOWNG_00001;Name=genex;db_xref=COG:COG1289;gene=genex;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P46481;locus_tag=UNKNOWNG_00001;product=unknown"
# =================================

second_lines = {}  # Used to save the second row of each file

with open(file_list, "r") as f:
    for line in f:
        filepath = Path(line.strip())
        filename = (str(filepath).split('/')[-1]).split('.')[0]
        new_filepath = './data/' + filename + '/genes.gff'
        if not filepath.exists():
            continue

        # Read the original file content
        with open(new_filepath, "r") as rf:
            lines = rf.readlines()

        # At least two rows are required
        if len(lines) < 2:
            continue

        # 1️⃣ Read the second row (index 1)
        second_lines[str(new_filepath)] = lines[1].rstrip("\n")
        second_lines_info=second_lines[str(new_filepath)].split(' ')

        complete_new_line = second_lines_info[1] + '	Prodigal:002006	CDS	1	' + second_lines_info[3] + NEW_LINE

        # 2️⃣ Insert a new row before the third row (index 2)
        if len(lines) >= 2:
            lines.insert(2, complete_new_line+"\n")

        # 3️⃣ Write back to file
        with open(new_filepath, "w") as wf:
            wf.writelines(lines)