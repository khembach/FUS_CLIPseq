# Analysis of mouse FUS CLIP-seq data

This repository contains the code for the analysis of FUS CLIP-seq data from WT mouse brain homogenate and synaptoneurosome preparations. Each sample consists of two file for the two different smears that we see on the blot: One at 70k where FUS runs and another smear at 115k. In our analyses, we focused on the 70k samples.

The FASTQ files were processed with a modified version of the [ARMOR workflow](https://github.com/csoneson/ARMOR). The repository includes all code that was used to explore the data. For the final analysis, we used peaks called by [CLIPper](https://github.com/YeoLab/clipper) and a custom filtering based on MA plots (see `Rmd/target_gene_selection_analysis.Rmd`).
