#!/usr/bin/env python3
import os
import sys
from collections import Counter
from Bio import SeqIO

def read_reference(ref_fa):
    """Read reference sequence lengths"""
    ref_lengths = {}
    for record in SeqIO.parse(ref_fa, "fasta"):
        ref_lengths[record.id] = len(record.seq)
    return ref_lengths

def analyze_msa(msa_file):
    """Summarize the contig length distribution in MSA files"""
    lengths = [len(record.seq) for record in SeqIO.parse(msa_file, "fasta")]
    if not lengths:
        return None, None, None, None

    # Mode
    counter = Counter(lengths)
    mode_length = counter.most_common(1)[0][0]

    # Number of sequences outside the mode +/- 20%
    lower = mode_length * 0.8
    upper = mode_length * 1.2
    outliers = sum(1 for l in lengths if l < lower or l > upper)

    return mode_length, outliers, len(lengths), lengths

def main(ref_fa, msa_dir, output_file):
    ref_lengths = read_reference(ref_fa)

    with open(output_file, "w") as out:
        out.write("MSA_file\tMode_length\tReference_length\tMSA_total\tOutlier_count\n")

        for msa_file in sorted(os.listdir(msa_dir)):
            if not msa_file.endswith(".aln.fas"):
                continue

            msa_path = os.path.join(msa_dir, msa_file)
            mode_length, outliers, total, _ = analyze_msa(msa_path)
            if mode_length is None:
                continue

            # Remove the .aln.fas suffix to obtain the contig ID
            contig_id = msa_file.replace(".aln.fas", "")
            ref_length = ref_lengths.get(contig_id, "NA")

            out.write(f"{msa_file}\t{mode_length}\t{ref_length}\t{total}\t{outliers}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <pan_genome_reference.fa> <msa_dir> <output.tsv>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])