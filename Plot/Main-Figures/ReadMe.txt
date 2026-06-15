Main figure scripts
===================

This folder contains scripts used to generate main-text figure panels and
intermediate plotting data. Most scripts use fixed local paths to prepared
tables, VCF files, Excel sheets, or figure-intermediate folders. Edit those
paths before running.

General usage
-------------

Python figure scripts:

    python <script>.py

R figure scripts:

    Rscript <script>.R

Scripts
-------

2a.py
    Purpose:
        Plot contig N50 values for pure short-read and hybrid long+short-read
        assemblies.

2b-data.py
    Purpose:
        Count CDS entries in GFF files and write per-sample gene counts.

2b.py
    Purpose:
        Plot gene-number distributions for different assembly methods and run
        a Mann-Whitney U test.

2c.py
    Purpose:
        Draw pan-genome accumulation curves from gene_contig_matrix.tsv and
        classify genes as core, dispensable, or rare by prevalence.

2c.R
    Purpose:
        R plotting script for the corresponding Fig. 2c panel.

2d-data.py
    Purpose:
        Count AMR genes across core, dispensable, rare, and invalid gene
        groups using amr_panaroo_dict.json and classified gene folders.

2e-data.py
    Purpose:
        Count variant types from VCF files and split counts by sample source
        groups such as Hui-net/GN and PATRIC/GC.

2e.py
    Purpose:
        Draw stacked bar plots for SNP, insertion, and deletion counts.

2f-data.py
    Purpose:
        Count snpEff annotation effects from annotated VCF files.

2f-bar.py
    Purpose:
        Draw a bar plot of snpEff variant-effect categories.

2f-pie.py
    Purpose:
        Draw a pie chart grouping variant effects into synonymous, missense,
        and other categories.

2ghi.py
    Purpose:
        Generate additional Fig. 2 panels from prepared pangenome or feature
        summary tables.

3b-data.py
    Purpose:
        Summarize which selected features are supported by RF, LR, SVM, and
        XGBoost rankings, and write model_evidence_support_summary.csv.

3b.py
    Purpose:
        Plot model-support overlap for selected features, including UpSet-style
        summaries.

3c.py
    Purpose:
        Plot selected candidate feature summaries for Fig. 3c.

3d.py
    Purpose:
        Plot candidate-feature or antibiotic-level summary data for Fig. 3d.

3de-u_test.py
    Purpose:
        Run statistical comparisons used for Fig. 3d/3e panels.

3e.py
    Purpose:
        Plot feature or model-performance summaries used in Fig. 3e.

4a-data.py
    Purpose:
        Prepare candidate-feature distribution data for Fig. 4a.

4a.py
    Purpose:
        Plot known/unknown selected candidate distributions by antibiotic.

4bc-graph-data.py
    Purpose:
        Prepare graph or network input data for Fig. 4b/4c.

4bc-lor-data.py
    Purpose:
        Compute log-odds-ratio related data for Fig. 4b/4c.

4bc-lor.py
    Purpose:
        Plot log-odds-ratio results for selected candidate features.

4d-data.R
    Purpose:
        Prepare R-side data for Fig. 4d.

4d.R
    Purpose:
        Draw Fig. 4d from prepared data.

4e.R
    Purpose:
        Draw Fig. 4e from prepared data.

4f.R
    Purpose:
        Draw Fig. 4f from prepared data.

4g-data.R
    Purpose:
        Prepare R-side data for Fig. 4g.

4g.R
    Purpose:
        Draw Fig. 4g from prepared data.
