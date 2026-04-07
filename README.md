# fxs-rnaseq

## This repository contains R code used for the WGCNA analysis in:  
Suresh A, Kourdougli N, Buth JE, Sánchez-León CA, Wall LT, Tran AT, Miranda-Rottmann S, Araya R, Gandal MJ, Portera-Cailliau C. Translatome profiling reveals opposing alterations in inhibitory and excitatory neurons of Fragile X mice and identifies EPAC2 as a therapeutic target. bioRxiv 2025.04.21.649817; doi: https://doi.org/10.1101/2025.04.21.649817

## Dataset Description (69 samples):  
### Genotypes = Fmr1 KO and WT mice  
### Brain Regions = S1 and V1  
### Cell Types = PV and CAMK2A

## The included RData object contains:  
datExpr = counts (17634 genes as rows, 69 samples as columns)  
datExpr.vst = vst normalized counts (17634 genes as rows, 69 samples as columns)  
datMeta (sample metadata including CRE_line, Genotype, Cell_type, Brain_Region, Sex, Group, and Mouse)  
geneAnno (gene annotation for 17634 genes in datExpr)
