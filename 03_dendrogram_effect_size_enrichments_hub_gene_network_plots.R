
#################################
## --- 01: Setup workspace --- ##
#################################

# clear workspace
rm(list=ls()) 

# set the default for strings as characters not factors
options(stringsAsFactors = FALSE)

# set working directory as R (above data, output

library(DESeq2);  library(nlme); library(gtools); library(grid); 
library(gridExtra); library(ggplot2); library(gtable); library(gprofiler2); 
library(WGCNA); library(igraph); library(ggrepel); library(ggpubr); 
library(gProfileR); library(ggsignif); library(rstatix); library(STRINGdb); library(gridGraphics); 

# output folder and network name
root.dir <- getwd()
networkdir = paste0("output/Consensus_network_signed_minModsize50")
network.name = "Consensus_network_signed_minModsize50"

# Create a folder for this network (change name for each network)
if (networkdir %in% list.files() == FALSE) {dir.create(file.path(networkdir), showWarnings = FALSE)}

###########################
## --- 02: Load data --- ##
###########################

load(file=paste0("data/Consensus_signed_pwr_12.RData"))

ls()
    #  [1] "consensuskME" "consExpr"     "consMEs"     
    #  [4] "consTree"     "geneAnno"     "labels"      
    #  [7] "MEs"          "moduleColors" "moduleLabels"
    # [10] "myLabels"     "net"          "network.name"
    # [13] "networkdir"   "root.dir"      "softPower" 
rm(rootdir, networkdir, labels)

kME = consensuskME
colors = moduleColors

#############################################################################################
## --- 03: Panel A. WGCNA dendrogram with module colors & other co-variates underneath --- ##
#############################################################################################

# Get additional gene info to put under dendrogram from biomaRt ####
library(biomaRt);
mart.mouse <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                dataset="mmusculus_gene_ensembl")
attributes <- listAttributes(mart.mouse)
  # transcript_start  Transcript start (bp) feature_page
  # transcript_end  Transcript end (bp) feature_page
  # transcription_start_site  Transcription start site (TSS)  feature_page
  # transcript_length Transcript length (including UTRs and CDS)  feature_page
  # transcript_count  Transcript count  feature_page
  # percentage_gene_gc_content  Gene % GC content feature_page
  # gene_biotype  Gene type feature_page

getinfo <- c("ensembl_gene_id", "external_gene_name", 
             #"hsapiens_homolog_ensembl_gene", 
             #"hsapiens_homolog_associated_gene_name",
             "percentage_gene_gc_content", "gene_biotype",
             "transcript_count")

bm.mouse <- getBM(attributes = getinfo, mart=mart.mouse)
dim(bm.mouse)
    # [1] 57010     5
head(bm.mouse)
# there are many gene_biotype categories, will reduce to TRUE/FALSE columns for protein_coding, and maybe miRNA, lncRNA, rRNA, snoRNA?
bm.mouse$protein_coding = rep("FALSE", nrow(bm.mouse))
bm.mouse$protein_coding[bm.mouse$gene_biotype=="protein_coding"] <- TRUE
    # table(bm.mouse$protein_coding)
    # FALSE  TRUE 
    # 35162 21848 
bm.mouse$miRNA = rep("FALSE", nrow(bm.mouse))
bm.mouse$miRNA[bm.mouse$gene_biotype=="miRNA"] <- TRUE
    # table(bm.mouse$miRNA)
    # FALSE  TRUE 
    # 54804  2206
bm.mouse$lncRNA = rep("FALSE", nrow(bm.mouse))
bm.mouse$lncRNA[bm.mouse$gene_biotype=="lncRNA"] <- TRUE
    # table(bm.mouse$lncRNA)
    # FALSE  TRUE 
    # 45411 11599

# Subset biomaRt results (bm.mouse) to genes in dataset ####
idx = match(geneAnno$ensembl_gene_id, bm.mouse$ensembl_gene_id)
bm.mouse.sub = bm.mouse[idx, ]
    # dim(geneAnno) # [1] 17633    10
    # dim(bm.mouse) # [1] 57010     8
    # dim(bm.mouse.sub) # [1] 17633     8

# Add the new info to geneAnno ####
geneAnno$percentage_gene_gc_content = bm.mouse.sub$percentage_gene_gc_content
geneAnno$gene_biotype = bm.mouse.sub$gene_biotype
geneAnno$transcript_count = bm.mouse.sub$transcript_count
geneAnno$protein_coding = bm.mouse.sub$protein_coding
table(geneAnno$protein_coding)
    # FALSE  TRUE 
    #  2977 14652 
geneAnno$miRNA = bm.mouse.sub$miRNA
    # table(geneAnno$miRNA)
    # FALSE 
    # 17629
geneAnno$lncRNA = bm.mouse.sub$lncRNA
    # table(geneAnno$lncRNA)
    # FALSE  TRUE 
    # 15869  1760
# Load and format FMRP targets & SFARI gene lists
load(file="data/External_gene_lists/formatted_gene_lists.RData")  
    # gene.list.1 = FMRP targets (842 genes)
      # FMRP targets in dataset = 716

idx = match(geneAnno$external_gene_name, gene.list.1)
geneAnno$FMRP_targets = !is.na(idx)
    # table(geneAnno$FMRP_targets)
    # FALSE  TRUE 
    # 16916   717
den.fmrp.target = rep(NA, nrow(geneAnno))
den.fmrp.target[geneAnno$FMRP_targets=="FALSE"] <- "grey90"
den.fmrp.target[geneAnno$FMRP_targets=="TRUE"] <- "black" # "#D82243"

# den.gene.length = colorRampPalette(c("white", "#D80C32"))(nrow(geneAnno))[rank(geneAnno$length)]
den.gene.length = colorRampPalette(c("white", "black"))(nrow(geneAnno))[rank(geneAnno$length)]

den.num.transcript = colorRampPalette(c("white", "blue"))(nrow(geneAnno))[rank(geneAnno$transcript_count)]

den.gc = colorRampPalette(c("white", "#117733"))(nrow(geneAnno))[rank(geneAnno$percentage_gene_gc_content)]

den.protein.coding = rep(NA, nrow(geneAnno))
den.protein.coding[geneAnno$protein_coding=="FALSE"] <- "grey90"
den.protein.coding[geneAnno$protein_coding=="TRUE"] <- "black" # "#D82243"

den.lncRNA = rep(NA, nrow(geneAnno))
den.lncRNA[geneAnno$lncRNA=="FALSE"] <- "grey90"
den.lncRNA[geneAnno$lncRNA=="TRUE"] <- "black" # "#D82243"

# plot main fig dendrogram ####
pdf(file=paste0("output/Consensus_network_signed_minModsize50/01a_Consensus_network_signed_minModsize50_sft12_Dendrogram_Main_fig_version.pdf"), w=4, h=3)
plotDendroAndColors(consTree, dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05, colorHeightMax = 0.3,
                    colorHeightBase = 0.4, main=paste0( 
                      "CAMK2/PV Signed Consensus\nWGCNA (",
                      length(unique(net$colors))," Modules)"),
                    colors=cbind(moduleColors, den.gene.length, 
                                 den.protein.coding, den.fmrp.target),
                    groupLabels=c("Module Colors", "Gene Length", "Protein Coding",
                                  "FMRP targets"), 
                    marAll = c(1, 5, 3, 1), saveMar = TRUE,
                    cex.lab = 1, cex.axis = 1, cex.main = 1)
dev.off()



# plot supplemental fig dendrogram ####

# add in avg gene expression for diff subsets

library(dplyr);

# consExpr has samps as rows & genes as columns
df = rbind(consExpr[[1]]$data, consExpr[[2]]$data)
match(rownames(df), datMeta$Samples)

df.S1 = df[datMeta$Region=="S1",]
S1.avg = data.frame("S1"=scale(rowMeans(t(df.S1))))
den.S1 = colorRampPalette(c("white", "#D80C32"))(nrow(S1.avg))[rank(S1.avg$S1)]

df.V1 = df[datMeta$Region=="V1",]
V1.avg = data.frame("V1"=scale(rowMeans(t(df.V1))))
den.V1 = colorRampPalette(c("white", "#D80C32"))(nrow(V1.avg))[rank(V1.avg$V1)]

df.WT = df[datMeta$Genotype=="WT",]
WT.avg = data.frame("WT"=rowMeans(t(df.WT)))
den.WT = colorRampPalette(c("lightblue", "#D80C32"))(nrow(WT.avg))[rank(WT.avg$WT)]

df.KO = df[datMeta$Genotype=="KO",]
KO.avg = data.frame("KO"=rowMeans(t(df.KO)))
den.KO = colorRampPalette(c("lightblue", "#D80C32"))(nrow(KO.avg))[rank(KO.avg$KO)]

df.region = data.frame(S1.avg, V1.avg)
df.region$subtract = df.region$S1 - df.region$V1
den.region = colorRampPalette(c("blue", "white", "#D80C32"))(nrow(df.region))[rank(df.region$subtract)]

df.geno = data.frame(WT.avg, KO.avg)
df.geno$subtract = df.geno$KO - df.geno$WT
den.genotype = colorRampPalette(c("blue", "white", "#D80C32"))(nrow(df.geno))[rank(df.geno$subtract)]

df.CAMK2 = df[datMeta$Cell_type=="CAMK2",]
CAMK2.avg = data.frame("CAMK2"=rowMeans(t(df.CAMK2)))

df.PV = df[datMeta$Cell_type=="PV",]
PV.avg = data.frame("PV"=rowMeans(t(df.PV)))

df.cell = data.frame(CAMK2.avg, PV.avg)
df.cell$subtract = df.cell$CAMK2 - df.cell$PV
den.cell = colorRampPalette(c("blue", "white", "#D80C32"))(nrow(df.cell))[rank(df.cell$subtract)]

pdf(file=paste0("output/Consensus_network_signed_minModsize50/01b_Consensus_network_signed_minModsize50_sft12_Dendrogram_Supplemental_fig_version.pdf"), w=4, h=4)

  plotDendroAndColors(consTree, 
                      dendroLabels = FALSE,
                      hang = 0.03, 
                      addGuide = TRUE, 
                      guideHang = 0.05,
                      colorHeightMax = 0.3, 
                      colorHeightBase = 0.4,
                      main=paste0("CAMK2/PV Signed Consensus\nWGCNA (", length(unique(net$colors))," Modules)"),
                      colors=cbind(moduleColors, den.gc, den.lncRNA, den.genotype, den.region, den.cell),
                      groupLabels=c("Module Colors", "% GC", "lncRNA", "KO - WT", "S1 - V1", "CAMK2 - PV"),
                      marAll = c(1, 5, 3, 1), saveMar = TRUE,
                      cex.lab = 1, cex.axis = 1, cex.main = 1)
dev.off()

##################################################
## --- 04: Panel B. Module fold change plot --- ##
##################################################

# Graph is a heatmap, x is the module names, y is log2(FC), color is bars is -log10(pval), and then asterisks are added for significance
library(nlme); library(emmeans);

MEs.all.export$Cell_Genotype = paste0(MEs.all.export$Cell_type, 
                                      "_", MEs.all.export$Genotype)

# swap order of factor levels so contrast is KO - WT below
datMeta$Genotype = factor(datMeta$Genotype, levels=c("KO", "WT"))
res = list() # empty list for results

# Get effect sizes for each module ####
for(i in 1:ncol(MEs)) {
  
  # Select Module #
  me = MEs[,i]; m = colnames(MEs)[i]; # extract ME values & name
  moduleColor = c = substr(m,3,nchar(m)) # remove "ME" from name

  # extract out estimate/pvals #
  # linear model with adjustment for subject (same mouse has S1 and V1 samples)
  lm.test = lme(me ~ Genotype*Cell_type*Region + Sex, random=~1|Mouse,
                data=datMeta)
  # pairwise comparisons - fdr adjusted for the 4 tests
  em.test = emmeans(lm.test, pairwise ~ Genotype*Cell_type*Region, 
                  by=c("Cell_type", "Region"), cross.adjust="fdr")
  # format results 
  res[[i]] = data.frame(moduleColor = moduleColor, 
                        comparison=c("CAMK2_S1", "PV_S1", "CAMK2_V1", "PV_V1"),
                        em.test$contrasts %>% rbind(adjust="fdr"))
    # res looks like this for each module
      #   moduleColor comparison contrast Cell_type Region
      # 1       black   CAMK2_S1  KO - WT     CAMK2     S1
      # 2       black      PV_S1  KO - WT        PV     S1
      # 3       black   CAMK2_V1  KO - WT     CAMK2     V1
      # 4       black      PV_V1  KO - WT        PV     V1
      #      estimate         SE df   t.ratio      p.value
      # 1 0.266569873 0.07091670 32 3.7589155 0.0006859072
      # 2 0.063664041 0.06103245 32 1.0431180 0.3047103776
      # 3 0.186762929 0.07331749 28 2.5473177 0.0166338126
      # 4 0.008092855 0.06103245 28 0.1325992 0.8954586882
}

# convert "res" list into data.frame
all.mod.res = do.call(rbind.data.frame, res)
# add column for -log10(pval)
all.mod.res$neg.log10.pval = -log10(all.mod.res$p.value)
table(all.mod.res$p.value < 0.05) # these are fdr corrected pvals for 4 tests
    # FALSE  TRUE 
    #    93    23 
# add column with *, **, *** for significance
all.mod.res$asterisks = stars.pval(all.mod.res$p.value)

# export 
write.csv(all.mod.res, file=paste0("output/Consensus_network_signed_minModsize50/05a_Consensus_network_signed_minModsize50_sft12_module_fold_change_statistics.csv"))

# Plot effect size results ####
# select just the needed columns
plot.dat = all.mod.res[, c("moduleColor", "comparison", 
                           "estimate", "neg.log10.pval")]

# change moduleColor to a factor and order mods by effect size
plot.dat$moduleColor = factor(plot.dat$moduleColor,
                          levels=unique(plot.dat$moduleColor[order(plot.dat$estimate,
                                                                decreasing=FALSE)]))
# change comparison to a factor 
plot.dat$comparison = factor(plot.dat$comparison, 
                             levels = c("CAMK2_S1", "CAMK2_V1", "PV_S1", "PV_V1"))

# create a column for the bar color assignments (for the scale)
plot.dat$pval.color = rep("NA", nrow(plot.dat))
plot.dat$pval.color[plot.dat$neg.log10.pval < 1.3] <- 1
plot.dat$pval.color[plot.dat$neg.log10.pval > 1.3 & plot.dat$neg.log10.pval < 2.5] <- 2
plot.dat$pval.color[plot.dat$neg.log10.pval > 2.5 & plot.dat$neg.log10.pval < 3.5] <- 3
plot.dat$pval.color[plot.dat$neg.log10.pval > 3.5 & plot.dat$neg.log10.pval < 4.5] <- 4
plot.dat$pval.color[plot.dat$neg.log10.pval > 4.5 & plot.dat$neg.log10.pval < 5.5] <- 5

# ggplot of effect sizes (4 columns) ####
FC.plot <- ggplot(data=plot.dat, aes(x=estimate, y=moduleColor, fill=pval.color)) +
  geom_bar(aes(x=estimate, y=moduleColor, fill=pval.color),
           stat="summary", col="black") +
  facet_wrap(~comparison, ncol=4) + # plot 1x4 
  theme_bw() +
  xlab("Effect Size (estimate)") +
  ylab("Module") +
  scale_fill_manual(values = c("1"="grey90", "2"="#EA9C9C",
                               "3"="#D47767", "4"="#E24935",
                               "5"="#C50A0A"), # assign colors
                    labels=c("n.s.", "1.3-2.5", "2.5-3.5",
                             "3.5-4.5", "4.5-5.5"), # range labels
                    name="-log10(p.value)") # legend label

pdf(file=paste0("output/Consensus_network_signed_minModsize50/02_Consensus_network_signed_minModsize50_sft12_module_fold_change_plot_4columns.pdf"),h=5, w=8)
print(FC.plot)
dev.off()

############################################################################################
## --- 05: Panel C. Update SFARI gene list enrichments (split categories / syndromic) --- ##
############################################################################################

SFARI.list = read.csv(file="data/External_gene_lists/SFARI_2023-07-12/SFARI-Gene_genes_04-03-2023release_07-12-2023export.csv")
table(SFARI.list$gene.score, SFARI.list$syndromic)
  #     0   1
  # 1 117 115
  # 2 646  53
  # 3  84  19
# Keep syndromic and highly ranked (1 and 2) genes
# separate as syndromic & rest of 1 & 2?
idx.syn = na.omit(rownames(SFARI.list)[SFARI.list$syndromic=="1"])
syndromic = SFARI.list[idx.syn,]
table(syndromic$gene.score)
    #   1   2   3 
    # 115  53  19 

idx = na.omit(rownames(SFARI.list)[SFARI.list$gene.score=="1"])
cat1 = SFARI.list[idx, ]
dim(cat1)
    # [1] 232  10 
table(cat1$syndromic)
    #   0   1 
    # 117 115

idx2 = na.omit(rownames(SFARI.list)[SFARI.list$gene.score=="2"])
cat2 = SFARI.list[idx2, ]
dim(cat2)
    # [1] 699  10
table(cat2$syndromic)
    #   0   1 
    # 646  53 

cat1.cat2 = data.frame(rbind(cat1, cat2))
dim(cat1.cat2)
    # [1] 931  10

# fisher's exact test w/ multiple testing correction based on num modules (29)
cat1.cat2.in.data = na.omit(intersect(geneAnno$hsapiens_homolog_ensembl_gene, cat1.cat2$ensembl.id))
cat1.cat2.in.data = cat1.cat2.in.data[!cat1.cat2.in.data==""]
length(cat1.cat2.in.data)
    # [1] 781

enrich.res = list()
for(i in 1:ncol(MEs)) {
  
  # Select Module ####
  me = MEs[,i]; m = colnames(MEs)[i]; # extract ME values & name
  moduleColor = c = substr(m,3,nchar(m)) # remove "ME" from name
  # get mod genes (mouse ensembl_ids)
  mod.genes = names(moduleColors[moduleColors==moduleColor]) 
  idx = na.omit(match(mod.genes, geneAnno$ensembl_gene_id))
  mod.geneAnno = geneAnno[idx, ]
  
  write.csv(mod.geneAnno, file=paste0("output/Consensus_network_signed_minModsize50/", moduleColor, "_module_gene_annotation.csv"))

  # fisher's exact test ####
  # contingency table
  mod.size = unique(mod.geneAnno$hsapiens_homolog_ensembl_gene) 
  
  mod.size = mod.size[!mod.size==""]
      # 298 genes in module
  
  list.size = unique(intersect(cat1.cat2$ensembl.id,
                           geneAnno$hsapiens_homolog_ensembl_gene)) 
  
  list.size = list.size[!list.size==""]
      # 782/931 in our dataset
 
   mod.list.overlap = unique(intersect(mod.size, list.size))
  
   all.genes = geneAnno$hsapiens_homolog_ensembl_gene
      # 55 overlap
  
   list.not.mod = unique(setdiff(list.size, mod.size)) 
  
  contingencyTable <- data.frame( 
    modGenes.No=c(
      length(list.size)-length(mod.list.overlap), # mod.No/ext.Yes
      length(all.genes)-length(unique(c(mod.size, list.size)))), # mod.No/ext.No
    mod.Genes.Yes=c(
      length(mod.list.overlap), # mod.Yes/ext.Yes
      length(mod.size)-length(mod.list.overlap))) # mod.Yes/ext.No
  
  row.names(contingencyTable) <- c("externalGenes.Yes", "externalGenes.No")
          # black mod as example
          #                   modGenes.No mod.Genes.Yes
          # externalGenes.Yes         727            54    781 external.list.size
          # externalGenes.No        16608           244
          #                                         298    17633 total genes
          #                                     mod.size
  
  # Fisher's exact test
  fisher.res = fisher.test(contingencyTable)

  # Format fisher results 
  enrich.res[[i]] = data.frame(
    moduleColor = moduleColor,
    gene.list = "SFARI.cat1.and.cat2",
    method = fisher.res[["method"]],
    alternative = fisher.res[["alternative"]],
    conf.int.lower = fisher.res[["conf.int"]][[1]],
    conf.int.higher = fisher.res[["conf.int"]][[2]],
    odds.ratio = fisher.res[["estimate"]][["odds ratio"]],
    pval = fisher.res[["p.value"]],
    adj.pval = p.adjust(fisher.res[["p.value"]], method="fdr", n=29)
    )
}

# convert "res" list into data.frame with results for all modules
all.mod.enrich.res = do.call(rbind.data.frame, enrich.res)

# add column for -log10(pval)
all.mod.enrich.res$neg.log10.adj.pval = -log10(all.mod.enrich.res$adj.pval)
table(all.mod.enrich.res$adj.pval < 0.05) 
    # FALSE  TRUE 
    #    21     8 

all.mod.enrich.res$moduleColor[all.mod.enrich.res$adj.pval < 0.05]
    # [1] "turquoise"  "black"      "brown"      "darkorange" "grey60"    
    # [6] "blue"       "pink"       "cyan" 

write.csv(all.mod.enrich.res, file="output/Consensus_network_signed_minModsize50/SFARI_categories_1_and_2_module_list_enrichments.csv")

########################################################
## --- Hub genes for black and lightgreen modules --- ##
########################################################

# Black module hubs ####
kME.black = kME[,c("MEblack", "Gene", "hsapiens_homolog_ensembl_gene",
                   "hsapiens_homolog_associated_gene_name",
                   "FMRP_targets", "SFARI", "Fmr1_KO_CA1_DE_genes",
                   "Fmr1_KO_Granule_DE_genes")]
kME.black.ordered = kME.black[order(kME.black$MEblack, decreasing=TRUE),]
kME.black.ordered$Gene[1:20]
    #  [1] "Smg1"    "Cplane1" "Scn3a"   "Xrn1"    "Tasor"   "Vps13a"  "Kmt5b"  
    #  [8] "Srrm2"   "Nktr"    "Csmd1"   "Kmt2b"   "Birc6"   "Nexmif"  "Med23"  
    # [15] "Ncan"    "Pogz"    "Otud4"   "Tjp1"    "Trrap"   "Atxn2l" 

black.mod.genes = names(moduleColors[moduleColors=="black"]) 
idx = na.omit(match(black.mod.genes, geneAnno$ensembl_gene_id))
black.mod.geneAnno = geneAnno[idx, ]
dim(black.mod.geneAnno)
    # [1] 328  20

idx = na.omit(match(black.mod.geneAnno$external_gene_name, kME.black.ordered$Gene))
kME.black.ordered = kME.black.ordered[idx, ]
kME.black.ordered = kME.black.ordered[order(kME.black.ordered$MEblack, decreasing=TRUE),]
dim(kME.black.ordered)
    # [1] 328   8
table(kME.black.ordered$FMRP_targets) # 41 genes
table(kME.black.ordered$SFARI) # 54 genes
table(kME.black.ordered$Fmr1_KO_CA1_DE_genes) # none
table(kME.black.ordered$Fmr1_KO_Granule_DE_genes) # 4 genes

kME.black.hubs = kME.black.ordered[1:20, ]
kME.black.hubs[kME.black.hubs$FMRP_targets=="TRUE"]

table(kME.black.hubs$FMRP_targets) # 4 genes
table(kME.black.hubs$SFARI) # 8 genes
table(kME.black.hubs$Fmr1_KO_CA1_DE_genes) # none
table(kME.black.hubs$Fmr1_KO_Granule_DE_genes) # none
kME.black.hubs$Gene[kME.black.hubs$FMRP_targets]
    # [1] "Srrm2" "Birc6" "Ncan"  "Trrap"
kME.black.hubs$Gene[kME.black.hubs$SFARI]
    # [1] "Kmt5b"  "Srrm2"  "Csmd1"  "Birc6"  
    # "Nexmif" "Med23"  "Pogz"   "Trrap"

signif(kME.black.hubs$MEblack, 3)
#  [1] 0.919 0.915 0.910 0.894 0.891 0.889 0.888 0.887 0.886 0.885 0.885 0.885
# [13] 0.883 0.881 0.878 0.877 0.875 0.875 0.874 0.870


# Lightgreen module hubs ####
kME.lightgreen = kME[,c("MElightgreen", "Gene", "hsapiens_homolog_ensembl_gene",
                   "hsapiens_homolog_associated_gene_name",
                   "FMRP_targets", "SFARI", "Fmr1_KO_CA1_DE_genes",
                   "Fmr1_KO_Granule_DE_genes")]
kME.lightgreen.ordered = kME.lightgreen[order(kME.lightgreen$MElightgreen, decreasing=TRUE),]
kME.lightgreen.ordered$Gene[1:20]
    #  [1] "Slc1a1"   "Sema3e"   "Grm1"     "Wfs1"     "Smap2"    "Ehd3"    
    #  [7] "Cds1"     "Ube2ql1"  "Galnt17"  "Agfg2"    "Bcl6"     "Tmem150c"
    # [13] "Fgfr1"    "Fam110b"  "Rps6ka2"  "Slc25a25" "Pim1"     "Coro1c"  
    # [19] "Bace1"    "Mrps30" 

lightgreen.mod.genes = names(moduleColors[moduleColors=="lightgreen"]) 
idx = na.omit(match(lightgreen.mod.genes, geneAnno$ensembl_gene_id))
lightgreen.mod.geneAnno = geneAnno[idx, ]
dim(lightgreen.mod.geneAnno)
    # [1] 151  20

idx = na.omit(match(lightgreen.mod.geneAnno$external_gene_name, kME.lightgreen.ordered$Gene))
kME.lightgreen.ordered = kME.lightgreen.ordered[idx, ]
kME.lightgreen.ordered = kME.lightgreen.ordered[order(kME.lightgreen.ordered$MElightgreen, decreasing=TRUE),]
dim(kME.lightgreen.ordered)
    # [1] 151   8
table(kME.lightgreen.ordered$FMRP_targets) # 13
table(kME.lightgreen.ordered$SFARI) # 10
table(kME.lightgreen.ordered$Fmr1_KO_CA1_DE_genes) # 4
table(kME.lightgreen.ordered$Fmr1_KO_Granule_DE_genes) # 3

kME.lightgreen.hubs = kME.lightgreen.ordered[1:20, ]
kME.lightgreen.hubs[kME.lightgreen.hubs$FMRP_targets=="TRUE"]

table(kME.lightgreen.hubs$FMRP_targets) # none
table(kME.lightgreen.hubs$SFARI) # 3 genes
table(kME.lightgreen.hubs$Fmr1_KO_CA1_DE_genes) # 1 gene
table(kME.lightgreen.hubs$Fmr1_KO_Granule_DE_genes) # none
kME.lightgreen.hubs$Gene[kME.lightgreen.hubs$FMRP_targets]

kME.lightgreen.hubs$Gene[kME.lightgreen.hubs$SFARI]
    # [1] "Slc1a1"  "Fgfr1"   "Rps6ka2"
kME.lightgreen.hubs$Gene[kME.lightgreen.hubs$Fmr1_KO_CA1_DE_genes]
    # [1] "Sema3e"
kME.lightgreen.hubs$Gene[kME.lightgreen.hubs$Fmr1_KO_Granule_DE_genes]

signif(kME.lightgreen.hubs$MElightgreen, 3)
#  [1] 0.878 0.876 0.865 0.853 0.849 0.834 0.821 0.814 0.810 0.805 0.803 0.796
# [13] 0.795 0.794 0.792 0.781 0.781 0.780 0.777 0.772









