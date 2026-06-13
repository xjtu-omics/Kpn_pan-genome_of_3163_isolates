import os
import gzip
import sys

def is_vcf(filename):
    return filename.endswith(".vcf") or filename.endswith(".vcf.gz")

def open_vcf(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')

def classify_variant(ref, alt):
    # Structural variants, such as <DEL> and <INS>
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

    # Initialize variant-type counters: each variant type has three categories: both, GN_only, and GC_only
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
            gn_indices = []  # Column indices for samples starting with GN, based on the full fields
            gc_indices = []  # Column indices for samples starting with GC, based on the full fields
            header_parsed = False

            for line in f:
                # Parse the header line to get sample names and corresponding column indices
                if line.startswith('#CHROM'):
                    fields = line.strip().split('\t')
                    # VCF standard format: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [sample1 sample2 ...]
                    # Samples start from column 10; zero-indexed column 9
                    for col_idx in range(9, len(fields)):
                        sample_name = fields[col_idx]
                        # Extract the sample ID prefix from the full path, either GCF_xxxxx or starting with GN
                        # Example: ./aaaT~~~aaaT_2~~~aaaT_1/GCF_001701365-1464_8_79.fa.bam
                        # Extract: GCF_001701365
                        sample_basename = sample_name.split('/')[-1]  # Get the file name
                        sample_id = sample_basename.split('-')[0]  # Get GCF_xxxxx or GN
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

                # Header information is required before data rows can be processed
                if not header_parsed:
                    print("Error: VCF header not found before data lines.")
                    continue

                fields = line.strip().split('\t')
                ref = fields[3]
                alt = fields[4]

                # Determine whether a sample has this variant, meaning it is not .:.
                def has_variant(col_idx):
                    if col_idx >= len(fields):
                        return False
                    return fields[col_idx] != ".:."

                vt = classify_variant(ref, alt)
                if vt is None:
                    continue

                # Check whether GN and GC samples have this variant
                gn_has = any(has_variant(col_idx) for col_idx in gn_indices)
                gc_has = any(has_variant(col_idx) for col_idx in gc_indices)

                # Count into the corresponding category
                if gn_has and gc_has:
                    variant_counts[vt]["both"] += 1
                elif gn_has:
                    variant_counts[vt]["GN_only"] += 1
                elif gc_has:
                    variant_counts[vt]["GC_only"] += 1
                else:
                    outliers.append([filename,fields[1],fields[3],fields[4]])  # Variants appearing in neither GN nor GC are treated as outliers

            print(len(gn_indices), len(gc_indices))

    for i in outliers:
        print(i)

    print(outliers_ids)

    return variant_counts

# Change the path to your VCF folder path
vcf_dir = "./both-align-results-strict-adv/ann_vcf/"
variant_counts = count_variants_in_directory(vcf_dir)

print("🧬 Total variant counts across all VCF files:")
print("\n" + "="*60)
for vtype in ["SNP", "Insertion", "Deletion", "Complex", "SV"]:
    counts = variant_counts[vtype]
    total = counts["both"] + counts["GN_only"] + counts["GC_only"]
    print(f"\n{vtype}:")
    print(f"  Total: {total}")
    print(f"  Both sample groups: {counts['both']}")
    print(f"  GN samples only: {counts['GN_only']}")
    print(f"  GC samples only: {counts['GC_only']}")