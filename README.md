# thp-1-hpylori-rna-seq  
## Infection of THP-1 macrophages with Helicobacter pylori B8 strain

This repository contains an RNA-seq analysis workflow for investigating transcriptional responses of THP-1 macrophages infected with the *Helicobacter pylori* B8 strain at 6 h and 24 h post infection.

After read trimming and alignment to the human reference genome (hg38) using the RNA-STAR aligner, read counts were generated using featureCounts. These preprocessing steps were performed on the US Galaxy server (usegalaxy.org) [1]. The downstream differential expression analysis was conducted in R and RStudio.

---

## Pipeline overview

The workflow includes:

* preprocessing and merging of featureCounts output files  
* differential expression analysis using DESeq2  
* identification of shared and unique differentially expressed genes (DEGs)  
* heatmap visualization of selected gene sets  
* Gene Ontology (GO) enrichment analysis  
* Venn diagram analysis of DEG overlap  
* chord diagram visualization of GO–gene relationships  

---

## Input data

* Raw count tables generated using featureCounts in Galaxy
* Samples: WT and Mock conditions at 6 h and 24 h post infection
* Reference genome: hg38 (Homo sapiens)

---

## Software and tools

* R (version 4.4.1)
* RStudio (version 2023.09.1)
* DESeq2
* clusterProfiler
* org.Hs.eg.db
* ggplot2
* pheatmap
* circlize
* eulerr
* enrichplot
* DOSE

---

## Reproducibility

All analyses were performed in R. A full record of the R session, including package versions and system information, is provided in `sessionInfo.txt`.

---

## Outputs

The workflow generates the following results:

* Differential expression tables (full and filtered DEGs, significant up-/downregulated saved as tsv files)
* Annotated gene lists (Ensembl IDs, gene symbols, Entrez IDs)
* Venn diagrams of shared and unique DEGs for each timepoint
* Heatmaps of top 50 DEGs
* GO enrichment plots (upregulated and downregulated genes)
* Chord diagrams linking GO terms and genes

---

## Author

Dr. Angelika Lahnsteiner

---

## References:
[1] The Galaxy Community. The Galaxy platform for accessible, reproducible, and collaborative data analyses: 2024 update, Nucleic Acids Research, 2024, 52(W1):W83-W94. doi:10.1093/nar/gkae410
