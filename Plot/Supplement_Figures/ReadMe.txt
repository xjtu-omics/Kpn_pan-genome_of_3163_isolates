Supplementary figure scripts
============================

This folder contains scripts for supplementary figures. Most scripts use fixed
local paths to prepared data files. Edit those paths before running.

General usage
-------------

Python scripts:

    python <script>.py

R scripts:

    Rscript <script>.R

Scripts
-------

sfig5.py
    Purpose:
        Generate supplementary figure panels related to selected variants,
        feature matrices, log-odds ratios, and contingency-table style plots.
        The script includes data extraction from annotated VCF-derived feature
        matrices and plotting code for specific candidate examples.
    Usage:
        Edit the input paths, target_gene, selected_antibiotics, and output
        path variables in the script, then run:

            python sfig5.py

sfig11.R
    Purpose:
        Generate Supplementary Fig. 11 from prepared R-readable data tables.
    Usage:
        Edit input/output paths in the R script if needed, then run:

            Rscript sfig11.R
