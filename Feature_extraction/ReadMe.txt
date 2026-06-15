Feature extraction scripts
==========================

This folder contains scripts for converting core-gene variants and
dispensable-gene presence/absence into machine-learning feature matrices.

General usage
-------------

These scripts mostly use fixed relative paths. Check the path variables near
the top or in the __main__ block before running.

Scripts
-------

vcf2matrix.py
    Purpose:
        Convert snpEff-annotated VCF files into binary feature matrices after
        removing synonymous variants. It also merges contigs from the same
        isolate, fills missing isolates with zero rows, filters rare features,
        and replaces AMR-related cluster names using amr_panaroo_dict.json.
    Main inputs:
        ./both-align-results-strict-adv/ann_vcf/
        ./Panaroo-DownStream-both/phenotypes.csv
        ./both-align-results-strict-adv/amr_panaroo_dict.json
    Main output:
        Per-gene CSV files under ./both-align-results-strict-adv/feature_matrix/
    Usage:

            python vcf2matrix.py

concat_core_feature_matrix.py
    Purpose:
        Concatenate per-gene core-variant feature matrices into one final core
        feature matrix. Feature names are prefixed with the gene name.
    Main inputs:
        ./both-align-results-strict-adv/feature_matrix/
        ./Panaroo-DownStream-both/phenotypes.csv
    Main output:
        ./Panaroo-DownStream-both/final_core_feature_matrix.csv
    Usage:

            python concat_core_feature_matrix.py

generate_dispensable_fm.py
    Purpose:
        Generate a binary dispensable-gene presence/absence matrix from
        dispensable gene MSA files. Columns are renamed with AMR names when a
        mapping exists, and features with too few 1s or too few 0s are reported
        as rare/uninformative.
    Main inputs:
        ./both-align-results-strict-adv/dispensable_gene_sequences/
        ./both-align-results-strict-adv/amr_panaroo_dict.json
    Main output:
        ./Panaroo-DownStream-both/final_dispensable_fm.csv
    Usage:

            python generate_dispensable_fm.py
