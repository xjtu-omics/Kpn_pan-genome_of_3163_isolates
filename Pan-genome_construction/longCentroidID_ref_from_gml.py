import networkx as nx
import csv
import ast

gml_file = "./both-align-results-strict-adv/final_graph.gml"
out_file = "./both-align-results-strict-adv/cluster_centroid_summary.tsv"

'''
# Read GML; networkx loads node attributes lazily and will not print everything at once
G = nx.read_gml(gml_file)

with open(out_file, "w") as f:
    f.write("cluster\tlongCentroidID\tcentroid\n")

    for node, data in G.nodes(data=True):
        cluster_name = data.get("name", node)

        long_centroid = data.get("longCentroidID", "NA")
        centroid = data.get("centroid", "NA")

        # Handle cases where centroid is a list or set
        if isinstance(centroid, (list, set, tuple)):
            centroid = ";".join(map(str, centroid))

        f.write(f"{cluster_name}\t{long_centroid}\t{centroid}\n")
'''

csv_file = "./both-align-results-strict-adv/gene_data.csv"
out_fasta = "./both-align-results-strict-adv/complete_pan_genome_reference.fasta"
tsv_file = out_file

# --------------------------------------------------
# 1. Read clusters.tsv
#    Extract: cluster_name -> clustering_id
# --------------------------------------------------
cluster_to_cid = {}

with open(tsv_file) as f:
    next(f)  # skip header
    for line in f:
        fields = line.rstrip("\n").split("\t")
        cluster = fields[0].strip()
        long_centroid_raw = fields[1].strip()

        try:
            # Parse the string "[300, '0_1_9']" into a Python list
            parsed = ast.literal_eval(long_centroid_raw)

            # Convention: the second element is the real clustering_id
            cid = str(parsed[1]).strip()

            cluster_to_cid[cluster] = cid

        except Exception as e:
            # If parsing fails, skip directly; raising is also possible
            continue

print(f"[INFO] clusters parsed: {len(cluster_to_cid)}")

# --------------------------------------------------
# 2. Read sequences.csv
#    Build: clustering_id -> dna_sequence
# --------------------------------------------------
seq_dict = {}

with open(csv_file) as f:
    reader = csv.DictReader(f)
    for row in reader:
        cid = row["clustering_id"].strip()
        dna = row["dna_sequence"].strip()

        # If a clustering_id appears multiple times, keep the first occurrence here
        if cid not in seq_dict:
            seq_dict[cid] = dna

print(f"[INFO] unique clustering_id in CSV: {len(seq_dict)}")

# --------------------------------------------------
# 3. Write FASTA: one sequence per cluster
# --------------------------------------------------
written = 0
with open(out_fasta, "w") as out:
    for cluster, cid in cluster_to_cid.items():
        seq = seq_dict.get(cid)
        if seq:
            out.write(f">{cluster}\n{seq}\n")
            written += 1

print(f"[INFO] FASTA written: {written}")