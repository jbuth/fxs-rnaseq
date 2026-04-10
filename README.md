# fxs-rnaseq
Co-expression network analysis (WGCNA) of RNA-seq data from a mouse model of Fragile X syndrome (Fmr1 KO versus WT).

## This repository contains R code for WGCNA implemented in:  
Suresh A, Kourdougli N, Buth JE, Sánchez-León CA, Wall LT, Tran AT, Miranda-Rottmann S, Araya R, Gandal MJ, Portera-Cailliau C. Translatome profiling reveals opposing alterations in inhibitory and excitatory neurons of Fragile X mice and identifies EPAC2 as a therapeutic target. bioRxiv 2025.04.21.649817; doi: https://doi.org/10.1101/2025.04.21.649817

## Dataset Description (69 samples x 17,634 genes):  
The dataset includes RNA-seq from 69 samples across genotype, brain region, and cell type.
- Genotypes: Fmr1 KO and WT mice  
- Brain Regions: S1 and V1  
- Cell Types: Pvalb (PV) and Camk2a (CAMK2A)

## Data Availability:
- The raw and processed dataset will be available at GEO after publication
  
## Analysis Overview:
Normalized expression matrix → WGCNA consensus network → module detection → module-trait relationships → functional enrichments

## Code Overview:
- 01_consensus_network.R
    - consensus network construction across cell types (PV and CAMK2A), module assignments, module-trait correlation, module eigengene expression, and module enrichments (protein-protein interactions, gene ontology, and Fisher's Exact test to test overlap with external gene lists)
- 02_biorxiv_figure_3A_and_3B.R
    - code to calculate and plot the dendrogram and module effect size in Figure 3A-B.
