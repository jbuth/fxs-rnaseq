# fxs-rnaseq
Co-expression network analysis (WGCNA) of RNA-seq data from a mouse model of Fragile X syndrome.

## This repository contains R code for weighted gene co-expression network analysis (WGCNA) applied to:  
Suresh A, Kourdougli N, Buth JE, Sánchez-León CA, Wall LT, Tran AT, Miranda-Rottmann S, Araya R, Gandal MJ, Portera-Cailliau C. Translatome profiling reveals opposing alterations in inhibitory and excitatory neurons of Fragile X mice and identifies EPAC2 as a therapeutic target. bioRxiv 2025.04.21.649817; doi: https://doi.org/10.1101/2025.04.21.649817

## Dataset Description (69 samples):  
The dataset includes RNA-seq from 69 samples across genotype, brain region, and cell type.
- Genotypes: Fmr1 KO and WT mice  
- Brain Regions: S1 and V1  
- Cell Types: PV and CAMK2A

## The included RData object contains the following processed data:  
- datExpr: counts (17,634 genes x 69 samples)  
- datExpr.vst: VST-normalized counts (17,634 genes x 69 samples) 
- datMeta: sample metadata including CRE_line, Genotype, Cell_type, Region, Sex, Group, and Mouse  
- geneAnno: gene annotation for the 17,634 genes in datExpr
