##################################
## --- 01: Set-up workspace --- ##
##################################

# clear workspace
rm(list=ls()) 
# set the default for strings as characters not factors
options(stringsAsFactors = FALSE)

library(DESeq2);  library(nlme); library(gtools); library(grid); library(gridExtra); library(ggplot2); library(gtable); library(gprofiler2); library(WGCNA); library(igraph); library(ggrepel); library(ggpubr); library(gProfileR); library(ggsignif); library(rstatix); library(STRINGdb); library(gridGraphics); 

# create shortcut for output folder and network name
rootdir <- getwd()
networkdir = paste0(rootdir, "/Network_Analysis/Networks/Signed_consensus_minModsize50_network")
network.name = "Signed_consensus_minModsize50_network_"

# Create a folder for this network (change name for each network)
if (networkdir %in% list.files() == FALSE) {dir.create(file.path(networkdir))}

########################################################################################
## --- 02: Calculate consensus network with minModsize = 50 (Set1=PV, Set2=CAMK2) --- ##
########################################################################################

load(file=paste0(rootdir, "/Network_Analysis/multiExpr_vst_norm_softthresh.RData"))

setwd(networkdir); # change to networdir for rest

# create data list
consExpr <- list(PV=multiExpr[[2]], CAMK2=multiExpr[[3]])
summary(consExpr)
names(consExpr)

# Check samples & genes
exp1 = goodSamplesGenes(consExpr[[1]]$data)
consExpr[[1]]$data = consExpr[[1]]$data[,exp1$goodGenes] # Removed ENSMUSG00000069014 from both
exp2 = goodSamplesGenes(consExpr[[2]]$data)
consExpr[[2]]$data = consExpr[[2]]$data[,exp1$goodGenes]

# create consensus network
softPower = 12; TOMfilebase = paste0(networkdir, "/", network.name, "consensus"); 
net = blockwiseConsensusModules(consExpr, maxBlockSize = 18000, corType = "bicor", power = softPower,
                                networkType = "signed", calibrationQuantile = 0.8, 
                                networkCalibration  = 'single quantile', consensusQuantile = 0.5, 
                                useMean = T, saveConsensusTOMs = TRUE, deepSplit = 2, 
                                minModuleSize = 50, # <-- increased to 50 minimum
                                pamStage = FALSE, reassignThresholdPS = 1e-10, mergeCutHeight = 0.25, verbose = 5)

consMEs = net$multiMEs; moduleLabels = net$colors; moduleColors = (moduleLabels); labels = ""; 
consTree = net$dendrograms[[1]]; MEs = rbind(net$multiMEs$PV$data, net$multiMEs$CAMK2$data); 
myLabels = data.frame(as.list(as.data.frame(moduleLabels)));
rownames(myLabels) = colnames(consExpr$PV$data);

kME = consensusKME(consExpr, moduleLabels=t(myLabels), multiEigengenes = net$multiMEs,  
                   consensusQuantile = 0.5, signed = TRUE, corAndPvalueFnc = bicorAndPvalue);
consensuskME = kME[,grep("consensus.kME", colnames(kME))];
idx = match(rownames(consensuskME), geneAnno$ensembl_gene_id);
geneAnno = geneAnno[idx,]; consensuskME$Gene = geneAnno$external_gene_name; 
colnames(consensuskME) = gsub("consensus.k", "", colnames(consensuskME))

# plot dendrogram
pdf(file=paste0(networkdir, "/01_", network.name, "sft", softPower, "_dendrogram.pdf"), w=8, h=5)
  plotDendroAndColors(consTree, moduleColors, "Module colors", dendroLabels = FALSE, hang = 0.03,
                  addGuide = TRUE, guideHang = 0.05, main = paste0("Consensus dendrogram sft ",
                  softPower, "\nnumModules=", length(unique(net$colors))))
dev.off()

## --- save network --- ##
save(net, geneAnno, softPower, consExpr, 
     consMEs, moduleLabels, moduleColors, 
     labels, consTree, MEs, myLabels,
     consensuskME, rootdir, networkdir, network.name, 
     file = paste0(networkdir, "/Consensus_signed_pwr_", softPower, ".RData"))

kME = consensuskME
colors = moduleColors

################################################
## --- 03: Export gene lists from all mods -- ##
################################################

list_of_genes = list()
for (i in 1:ncol(MEs)) {

  m = colnames(MEs)[i]; # extract out ME values and name
  moduleColor = c = substr(m,3,nchar(m)) 
  idx = which(colors==c) # in order by kME
  query = rownames(kME)[idx[order(kME[idx,paste0('ME',c)],decreasing = T)]]
  mod.kME = kME[idx[order(kME[idx,paste0('ME',c)],decreasing = T)],c(paste0('ME',c),"Gene")]
    
  idx = match(query, geneAnno$ensembl_gene_id)
  moduleGenes = geneAnno$ensembl_gene_id[idx] # ensembl_ids
  module.symbols = geneAnno$external_gene_name[idx] # gene names
  module.hsapiens_homolog_ensembl_gene = geneAnno$hsapiens_homolog_ensembl_gene[idx] #homolog ids
  module.hsapiens_homolog_associated_gene_name = geneAnno$hsapiens_homolog_associated_gene_name[idx] #homolog symbols

  list_of_genes[[i]] = list(moduleName = m, moduleColor = c, 
              moduleGenes = moduleGenes, 
              module.symbols = module.symbols, 
              module.hsapiens_homolog_ensembl_gene = module.hsapiens_homolog_ensembl_gene,
              module.hsapiens_homolog_associated_gene_name = module.hsapiens_homolog_associated_gene_name,
              kME = mod.kME[,1])
}
names(list_of_genes) = colnames(MEs)
save(list_of_genes, file=paste0(networkdir, "/", "All_module_gene_lists.RData"))

########################################
## --- 04: Module PPI enrichments --- ##
########################################

# BiocManager::install("STRINGdb", version="3.12")
library(STRINGdb);
    # Score threshold is not specified. We will be using medium stringency cut-off of 400.
PPI.enrich.res = list();
for (i in 1:length(list_of_genes)) {
  m = list_of_genes[[i]]$moduleName; moduleColor = c = list_of_genes[[i]]$moduleColor;
  moduleSymbols = list_of_genes[[i]]$module.symbols
  string_db <- STRINGdb$new( version="11", species=10090, input_directory="")
  mygenes = data.frame(gene=moduleSymbols); mapped <- string_db$map(mygenes, "gene", removeUnmappedRows = TRUE)
  PPI.enrich.res[[i]] = string_db$get_ppi_enrichment(mapped$STRING_id)
}
names(PPI.enrich.res) = colnames(MEs)

df=data.frame(matrix(NA, nrow = length(PPI.enrich.res), ncol = 3))
rownames(df) = names(PPI.enrich.res)
colnames(df) = c("enrichment", "edges", "lambda")

for (j in 1:length(PPI.enrich.res)) {
  if ((j==29)==TRUE) {
    df[j,] = c(0, 0, 0)
  } else {
  df[j,] = c(PPI.enrich.res[[j]]$enrichment, PPI.enrich.res[[j]]$edges, PPI.enrich.res[[j]]$lambda)
  }
}
df$stars = stars.pval(df$enrichment)

if (paste0(networkdir,"/PPI_enrichment") %in% list.files() == FALSE) {dir.create(file.path(paste0(networkdir,"/PPI_enrichment")))}

mod.sizes = data.frame(table(moduleColors)); # get module sizes
mod.sizes$moduleColors = paste0("ME", mod.sizes$moduleColors); # add "ME" to front to match df
idx = match(rownames(df), mod.sizes$moduleColors); mod.sizes = mod.sizes[idx, ]; # match row order
module.size = mod.sizes$Freq;

df.final = data.frame(module = rownames(df), module.size = module.size, expected.interactions = df$lambda, interactions = df$edges, p.value = df$enrichment, significance=df$stars);
df.final$module = paste0(rownames(df.final), "_", df.final$module)
g.df.final = tableGrob(df.final, rows = NULL); grid.arrange(g.df.final)
write.csv(df.final, file=paste0(networkdir, "/PPI_enrichment/module_PPI_enrichments.csv"))

##########################################
## --- 05: Module Trait Correlation --- ##
##########################################

nSets = checkSets(consExpr)$nSets
exprSize = checkSets(consExpr);
nSets = exprSize$nSets;

# check if order same
datMeta = rbind(consExpr[["PV"]]$meta, consExpr[["CAMK2"]]$meta)

idx1 = match(rownames(MEs), datMeta$Samples[datMeta$Cell_type=="PV"])
idx2 = match(rownames(MEs), datMeta$Samples[datMeta$Cell_type=="CAMK2"]) 
# looks good, already in order go to next step

datMeta.Traits = datMeta[,c("Genotype","Cell_type","Region","Sex","Mouse")]
colnames(datMeta.Traits) = c("Genotype", "Cell_type", "Region", "Sex", "Mouse")
rownames(datMeta.Traits) = datMeta$Samples
datTraits = data.frame(
  Genotype = as.numeric(factor(datMeta.Traits$Genotype, levels=c("WT", "KO"))),
  Cell_type = as.numeric(factor(datMeta.Traits$Cell_type, levels=c("PV", "CAMK2"))),
  Region = as.numeric(factor(datMeta.Traits$Region, levels=c("S1", "V1"))),
  Sex = as.numeric(factor(datMeta.Traits$Sex, levels=c("F", "M"))),
  Mouse = as.numeric(factor(datMeta.Traits$Mouse)))

MEs.all = rbind(net$multiMEs$PV$data, net$multiMEs$CAMK2$data);

modTraitCor = cor(MEs.all, datTraits, use = "p");
modTraitCor[1:3, 1:3]
modTraitP = corPvalueStudent(modTraitCor, nrow(datTraits))
modTraitP[1:3, 1:3]

adjusted.p = data.frame(matrix(NA, nrow = nrow(modTraitP), ncol = ncol(modTraitP)))
for (col in 1:ncol(modTraitP)) {
adjusted.p[ ,col] = p.adjust(modTraitP[ ,col], method="fdr")
}
colnames(adjusted.p) = colnames(modTraitP)
rownames(adjusted.p) = rownames(modTraitP)
adjusted.p.stars = data.frame(matrix(NA, nrow = nrow(adjusted.p), ncol = ncol(adjusted.p)))
for (col in 1:ncol(modTraitP)) {
adjusted.p.stars[,col] = stars.pval(adjusted.p[,col])
}
colnames(adjusted.p.stars) = colnames(adjusted.p)
rownames(adjusted.p.stars) = rownames(adjusted.p)

idx = match(rownames(adjusted.p.stars), rownames(modTraitCor))
modTraitCor = modTraitCor[idx,]

save(modTraitCor, adjusted.p, adjusted.p.stars, 
     file=paste0(networkdir, "/module_trait_relationship_results.RData"))

# Display the correlation values within a heatmap plot
pdf(file=paste0(networkdir,"/03_Module_Trait_Relationship.pdf"),h=10, w=9)
par(mar = c(8, 12.5, 3, 3))
labeledHeatmap(Matrix = modTraitCor, xLabels = names(datTraits), yLabels = rownames(modTraitCor), 
    ySymbols = rownames(modTraitCor), colorLabels = FALSE, colors = blueWhiteRed(50), 
    textMatrix = adjusted.p.stars, setStdMargins = FALSE, cex.text = 1, zlim = c(-1, 
        1), main = paste("Consensus Module-trait relationships with fdr correction\npval<0.001=*** pval<0.01=** pval<0.05=* pval<0.01=."),
    verticalSeparator.x=seq(1:ncol(modTraitCor)), horizontalSeparator.y=seq(1:nrow(modTraitCor)))
dev.off()

##############################################################################
## --- 06: Calculate module gene-list enrichments (Fisher's Exact Test) --- ##
##############################################################################

if (paste0(networkdir, "/module_gene_list_enrichments") %in% list.files() == FALSE) {dir.create(file.path(paste0(networkdir, "/module_gene_list_enrichments")))}

load(file=paste0(rootdir, "/gene_list_enrichments/formatted_gene_lists.RData"))
gene.list.for.enrichment = list(gene.list.1, gene.list.4, gene.list.5, gene.list.7) 
names(gene.list.for.enrichment) = c("FMRP mRNA targets", "Fmr1 KO CA1 DEGs", "Fmr1 KO Granule Cells DEGs", 
                                    "SFARI genes")

all.mod.padj = data.frame(matrix(NA, nrow = length(gene.list.for.enrichment), ncol = ncol(MEs)))
rownames(all.mod.padj) = names(gene.list.for.enrichment)
colnames(all.mod.padj) = colnames(MEs)

all.mod.odds = data.frame(matrix(NA, nrow = length(gene.list.for.enrichment), ncol = ncol(MEs)))
rownames(all.mod.odds) = names(gene.list.for.enrichment)
colnames(all.mod.odds) = colnames(MEs)

for (i in (1:ncol(MEs))) {
    
    me = MEs[,i]; m = colnames(MEs)[i]; # extract out ME values and name
    moduleColor = c = substr(m,3,nchar(m)) # remove "ME" from beginning, just color now
    name =  paste0(i, "_", c) 
    moduleGenes = list_of_genes[[i]]$moduleGenes
    moduleSymbols = list_of_genes[[i]]$module.symbols
    mod.genes = data.frame(Gene=moduleSymbols)

    enrichment.res = data.frame(external.list.size=c(length(gene.list.1), length(gene.list.4),
                                                length(gene.list.5),length(gene.list.7)),
                                mod.list.size=length(mod.genes$Gene), 
                                mod.overlap.size=rep(NA,4),
                                mod.list.odds=rep(NA,4), 
                                mod.list.pval=rep(NA,4), mod.list.padj=rep(NA,4))
    colnames(enrichment.res) = c("external.list.size", paste0(name,".mod.list.size"),
                                 paste0(name,".mod.overlap.size"), paste0(name,".mod.list.odds"),
                                 paste0(name, ".mod.list.pval"), paste0(name,".mod.list.padj"))
    rownames(enrichment.res) = names(gene.list.for.enrichment)
    
  for (k in 1:length(gene.list.for.enrichment)) {
        
          gene.list.test = gene.list.for.enrichment[[k]]
          
          # overlap with TEST list
          all_genes = length(rownames(geneAnno))
          TEST.in.exp.genes = intersect(gene.list.test,geneAnno$external_gene_name)
          TEST.pos = length(intersect(mod.genes$Gene,TEST.in.exp.genes))
          top.pos = length(mod.genes$Gene)
          TEST.in.exp.length = length(TEST.in.exp.genes)
        
        # Fisher Exact for up genes
        fisher.result <- fisher.test(data.frame(
          TEST.gene=c(TEST.pos, 
                      TEST.in.exp.length-TEST.pos),
          Not.TEST.gene=c(top.pos-TEST.pos,
                          all_genes-TEST.pos-(TEST.in.exp.length-TEST.pos)-(top.pos-TEST.pos))))
        
        enrichment.res[k,paste0(name,".mod.list.odds")] <- fisher.result$estimate
        enrichment.res[k,paste0(name,".mod.list.pval")] <- fisher.result$p.value
        enrichment.res[k,paste0(name,".mod.overlap.size")] <- TEST.pos
  }     
    enrichment.res[,paste0(name,".mod.list.padj")] = p.adjust(enrichment.res[,paste0(name,".mod.list.pval")], method = "bonferroni")
    
    all.mod.padj[,i] = enrichment.res[,paste0(name,".mod.list.padj")]
    all.mod.odds[,i] = enrichment.res[,paste0(name,".mod.list.odds")]
    
    write.csv(enrichment.res, file=paste0(networkdir, "/module_gene_list_enrichments/", name,
                                          "_fishers_exact_gene_list_enrichments.csv"))

    save(mod.genes, gene.list.for.enrichment, enrichment.res,
       file=paste0(networkdir, "/module_gene_list_enrichments/", name, "_module_gene_list_enrichments.RData"))
 }

# plot enrichments
all.mod.padj = all.mod.padj[,!colnames(all.mod.padj)=="MEgrey"]
all.mod.odds = all.mod.odds[,!colnames(all.mod.odds)=="MEgrey"]
test = reshape2::melt(all.mod.padj)
test$gene.list = rep(rownames(all.mod.padj), ncol(all.mod.padj))
test$gene.list = factor(test$gene.list, levels=c("FMRP mRNA targets", "SFARI genes", "Fmr1 KO Granule Cells DEGs", "Fmr1 KO CA1 DEGs"))
test$variable = factor(test$variable, levels = unique(test$variable))
colnames(test) = c("Module", "p.adj", "gene.list"); #levels(test$Module) = rev(levels(test$Module));

mod.padj <- ggplot(test, aes(x=Module, y=-log10(p.adj), fill=gene.list)) +
  geom_bar(position="dodge", stat="identity") + coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype="dotted", color = "black", size=0.5) +
  theme_bw() + facet_wrap(~gene.list, ncol=4) + theme(legend.position="none") + xlab("") 

pdf(paste0(networkdir, "/module_gene_list_enrichments/Module_list_enrichment_results.pdf"),w=12,h=8)   
  print(mod.padj)
dev.off()

######################################################
## --- 07: External Gene list overlap with DEGs --- ##
######################################################

      load(file=paste0(rootdir, "/All_DEG_lists.RData"))
      overlaps = list()
        for (i in 1:length(list_of_genes)) {
          overlaps[[i]] = list(intersect(list_of_genes[[i]]$module.symbols, gene.list.for.enrichment[[1]]),
                                intersect(list_of_genes[[i]]$module.symbols, gene.list.for.enrichment[[2]]),
                                intersect(list_of_genes[[i]]$module.symbols, gene.list.for.enrichment[[3]]),
                                intersect(list_of_genes[[i]]$module.symbols, gene.list.for.enrichment[[4]]),
                                intersect(list_of_genes[[i]]$module.symbols, DEGs[,1]),
                                intersect(list_of_genes[[i]]$module.symbols, DEGs[,2]),
                                intersect(list_of_genes[[i]]$module.symbols, DEGs[,3]),
                                intersect(list_of_genes[[i]]$module.symbols, DEGs[,4]),
                                intersect(list_of_genes[[i]]$module.symbols, DEGs[,5]),
                                intersect(list_of_genes[[i]]$module.symbols, DEGs[,6]),
                                intersect(list_of_genes[[i]]$module.symbols, DEGs[,7]),
                                intersect(list_of_genes[[i]]$module.symbols, DEGs[,8]))
        }
        names(overlaps) = names(list_of_genes)
        for (i in 1:length(overlaps)) { 
            names(overlaps[[i]]) = c(names(gene.list.for.enrichment), colnames(DEGs))
        }

##################################################
## --- 08: setup datMeta & data for figures --- ##
##################################################

datMeta = rbind(consExpr[[1]]$meta, consExpr[[2]]$meta) # check if this has final swap fixed
datMeta$Condition = paste0(datMeta$Cell_type,"_",datMeta$Genotype, "_",datMeta$Region)
datMeta$Condition = factor(datMeta$Condition)
datMeta$Condition = factor(datMeta$Condition)

datMeta$Cell_Geno = paste0(datMeta$Cell_type,"_",datMeta$Genotype)
datMeta$Cell_Geno = factor(datMeta$Cell_Geno)

datMeta$Geno_Region = paste0(datMeta$Genotype, "_", datMeta$Region)
datMeta$Geno_Region = factor(datMeta$Geno_Region)

idx = match(rownames(MEs),datMeta$Samples) 
datMeta = datMeta[idx, ] # make sure has same order as MEs

idx = match(colnames(consExpr[[1]]$data), geneAnno$ensembl_gene_id)
geneAnno = geneAnno[idx,] # make geneAnno match

##########################################################################
## --- 09: Loop through each mod and plot figure to inspect results --- ##
##########################################################################

for(i in 1:ncol(MEs)) {
  
  # Select Module ####
  me = MEs[,i]; m = colnames(MEs)[i]; # extract out ME values and name
  moduleColor = c = substr(m,3,nchar(m)) # remove "ME" from beginning, just color now
  
  # ANOVA ####
  a = anova(lme(me ~Genotype*Cell_type*Region + Sex, data=datMeta, random=~1|Mouse))

  p.Genotype = a$`p-value`[2]; p.GenotypeXCell_type=a$`p-value`[6]; p.GenotypeXRegion=a$`p-value`[7];
  p.GenotypeXCell_typeXRegion=a$`p-value`[9];
  
  if (( p.Genotype < 0.05 | p.GenotypeXCell_type < 0.05 | p.GenotypeXRegion < 0.05 | p.GenotypeXCell_typeXRegion < 0.05)==TRUE) { # if anova sig then plot the rest
      print(i)

      # BOXPLOT + pairwise.test - mod eigengenes ####
      dat.sub <- data.frame(Module=MEs[,i], Cell_Geno=datMeta$Cell_Geno, 
                            Region=datMeta$Region, Condition=datMeta$Condition)

      stat.test <- dat.sub %>% group_by(Region) %>% t_test(Module ~ Cell_Geno) %>%
        adjust_pvalue(method = "fdr") %>% add_significance()
      stat.test <- stat.test %>% add_xy_position(fun="max", x = "Cell_Geno")
      
      g.boxplot <- ggplot(dat.sub, aes(x=Cell_Geno, y=as.numeric(Module))) +
        geom_boxplot() +
        geom_point(aes(color=Cell_Geno), position=position_jitter(.1), size=3, alpha=0.5) + 
        facet_wrap(~Region) + xlab("") + ylab("eigengenes") + 
        theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none",
              text = element_text(size=20), plot.title  = element_text(hjust=0), 
              plot.subtitle = element_text(hjust=0, size=18)) +
        labs(title=paste0(m, " (nGenes=",length(list_of_genes[[i]]$moduleGenes),")"), 
             subtitle = paste0('Genotype p=', signif(p.Genotype,2), " ", stars.pval(p.Genotype), 
                        '\nGenotypeXCell_type p=', signif(p.GenotypeXCell_type,2), " ", 
                        stars.pval(p.GenotypeXCell_type),
                        '\nGenotypeXRegion p=', signif(p.GenotypeXRegion,2),
                        " ", stars.pval(p.GenotypeXRegion),
                        '\nGenotypeXCell_typeXRegion p=', signif(p.GenotypeXCell_typeXRegion,2),
                        " ", stars.pval(p.GenotypeXCell_typeXRegion)),
             caption="fdr adjusted p<0.001=***  p<0.01=** p<0.05=* p<0.1=.") + 
        stat_pvalue_manual(stat.test, hide.ns = TRUE, size = 6, bracket.size = 0.5) 

      # GO Terms gprofiler1 fdr ####
      query = as.vector(list_of_genes[[i]]$module.symbols)
      go.g2 = gprofiler(query, organism="mmusculus", custom_bg = geneAnno$external_gene_name, 
                 correction_method = "fdr", hier_filtering = "strong", ordered_query = T, 
                 significant = T, exclude_iea = F, region_query = F, max_p_value = 0.05,
                 max_set_size=1000, numeric_ns = "", include_graph = F, src_filter = c("GO", "KEGG","REACTOME"))
      go.g2 = go.g2[order(go.g2$p.value),]
      write.csv(go.g2,file=paste0(networkdir, "/Mod", i, "_", c, "_GO_ms_gprofiler1_setsize1000_fdr.csv"))
      
      plotme = go.g2[1:12, ]; plotme$term.name <- factor(plotme$term.name, levels = plotme$term.name)
      
      g.go.ms.res <- ggplot(plotme[1:12,], aes(x=term.name,y=-log10(p.value))) + 
            theme_bw(base_size = 16) +
            geom_bar(stat="identity", position = position_dodge(0.1), col="black", fill="grey80") + 
            coord_flip() + xlab("") + scale_x_discrete(limits = rev(levels(plotme$term.name))) +
            geom_hline(yintercept=-log10(0.05), lty=2, color="black") + 
            theme(axis.text.x = element_text(colour = "black"), 
                  axis.text.y = element_text(colour = "black"),
                  plot.title = element_text(size = 16)) +
        ggtitle(paste0("Gene Ontology Terms - ",c)) + theme(legend.position = "none")
      
      # HUB genes ####
      hubGenes = data.frame(gene = list_of_genes[[i]]$module.symbols[1:20], 
                            kME = list_of_genes[[i]]$kME[1:20])
      hubGenes$gene = factor(hubGenes$gene, levels=unique(hubGenes$gene))
      g.kME <- ggplot(hubGenes, aes(x=gene, y=kME)) +
        geom_bar(stat="identity", position = position_dodge(0.1), col="black", fill="grey80") +
        coord_flip() + ggtitle("Top 20 Hub Genes") +
        scale_x_discrete(limits = rev(levels(hubGenes$gene))) + theme_bw(base_size = 16) +
        theme(axis.text.y = element_text(colour="black"))
      
      # PPI enrichments ####
      PPI = df.final[gsub(".*_ME","",df.final$module)==c,]
      PPI2 = t(PPI[,-1])
      colnames(PPI2) = PPI$module
      g.PPI = tableGrob(PPI2)
      title3 <- textGrob("PPI enrichment", gp=gpar(fontsize=14))
      padding3 <- unit(5,"mm")
      g.PPI <- gtable_add_rows(g.PPI, heights = grobHeight(title3) + padding3, pos = 0)
      g.PPI <- gtable_add_grob(g.PPI, title3, 1, 1, 1, ncol(g.PPI))
      
      # Gene list overlaps ####
      overlap.sub = overlaps[[i]]
      overlap.table = data.frame(List=names(overlap.sub))
      for (j in 1:length(overlap.sub)) {
        overlap.table$ModGenes_Overlap[j] = length(overlap.sub[[j]])
      }
      overlap.table$List = factor(overlap.table$List, levels=unique(overlap.table$List))
      cbPalette = c("#4531AB", "#117733", "#44AA99", "#88CCEE", "#CFDCC5", "#DDCC77", "#FC98E4", "#CC6677",
                    "#AA4499", "#882255", "#847C7D", "#333232")
      overlap.plot = overlap.table[5:12,]
      overlap.plot$List = factor(overlap.plot$List, levels=unique(overlap.plot$List));
      
      overlap.size <- ggplot(overlap.plot, aes(x=List, y=ModGenes_Overlap, fill=List)) +
        geom_bar(stat="identity", col="black") +
        theme_bw(base_size = 16) +
        geom_text(aes(label = ModGenes_Overlap), position = position_dodge(width= 1),
                  vjust= 0.5, hjust = -0.3, size=3) +
        coord_flip() + ylab("Number of Genes") + xlab("") +
        theme(legend.position = "none", axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black")) + ggtitle("Module DEG List Overlap") +
        scale_fill_manual(values=cbPalette) +
        scale_x_discrete(limits = rev(levels(overlap.plot$List)))  +
        #ylim(0, max(overlap.table$ModGenes_Overlap)+3) +
        scale_y_continuous(expand = c(0.1,0.01))
      
      # Gene list enrichments ####
      mod.gene.list.enrich = data.frame(gene.list=rownames(all.mod.odds), odds.ratio=signif(all.mod.odds[,m],2),
                                        padj=signif(all.mod.padj[,m],2), 
                                        significance = stars.pval(all.mod.padj[,m]))
      g.mod.gene.list.enrich = tableGrob(mod.gene.list.enrich, rows=NULL)
      title4 <- textGrob("External gene list enrichments\np<0.001=***  p<0.01=** p<0.05=* p<0.1=.",
                         gp=gpar(fontsize=14)); padding4 <- unit(5,"mm");
      g.mod.gene.list.enrich <- gtable_add_rows(g.mod.gene.list.enrich, heights = grobHeight(title4) + padding4,
                                                pos = 0);
      g.mod.gene.list.enrich <- gtable_add_grob(g.mod.gene.list.enrich, title4, 1, 1, 1,
                                                ncol(g.mod.gene.list.enrich));
      # as barplot
      gene.list.plot = mod.gene.list.enrich; 
      gene.list.plot$gene.list =
        c(paste0(gene.list.plot$gene.list[1], " (", overlap.table$ModGenes_Overlap[1], "/842)"),
          paste0(gene.list.plot$gene.list[2], " (", overlap.table$ModGenes_Overlap[2], "/110)"),
          paste0(gene.list.plot$gene.list[3], " (", overlap.table$ModGenes_Overlap[3], "/110)"), 
          paste0(gene.list.plot$gene.list[4], " (", overlap.table$ModGenes_Overlap[4], "/402)"))
      gene.list.plot$gene.list <- factor(gene.list.plot$gene.list, levels = unique(gene.list.plot$gene.list))

      gg.gene.list.plot <- ggplot(gene.list.plot, aes(x=gene.list, y=-log10(padj), fill=gene.list)) +
            theme_bw(base_size = 16) +
            geom_bar(stat="identity", position = position_dodge(0.1), col="black") +
            geom_text(aes(label = significance), colour = "black", hjust = 0.5, vjust = -0.02,
                      size = 8, position = position_dodge(width = 0.5), angle=-90) +
            coord_flip() + xlab("") +
            geom_hline(yintercept=-log10(0.05), lty=2, color="black") + 
            theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
                  plot.title = element_text(size = 16)) + ggtitle(paste0("Gene List Enrichment - ",c)) +
            theme(legend.position = "none") + ylim(0, max(-log10(gene.list.plot$padj))+0.5) +
            scale_fill_manual(values = c("#4531AB", "#117733", "#44AA99", "#88CCEE")) +
            scale_x_discrete(limits = rev(levels(gene.list.plot$gene.list)))

      # Mod~Trait Cor ####
          load(file=paste0(networkdir,"/module_trait_relationship_results.RData"))
            adjusted.p.sub = data.frame(adjusted.p[m, ])
            adjusted.p.stars.sub = data.frame(adjusted.p.stars[m, ])
            cor.sub = data.frame(modTraitCor[m, ])
            modTraitCor.sub = data.frame(cor.sub, 
                                         adjusted.p=t(adjusted.p.sub),
                                         adjusted.p.stars=t(adjusted.p.stars.sub))
            modTraitCor.sub = cbind(rownames(modTraitCor.sub), modTraitCor.sub)
            colnames(modTraitCor.sub) = c("Condition", "Cor", "fdr.adjusted.p", "significance")


            par(mar=c(6, 11, 3, 6));
            labeledHeatmap(Matrix = as.matrix(modTraitCor.sub$Cor), 
                           xLabels = "",  yLabels = rownames(modTraitCor.sub), 
                           ySymbols = rownames(modTraitCor.sub), 
                           colorLabels = FALSE, colors = blueWhiteRed(50), 
                           textMatrix = as.matrix(modTraitCor.sub$significance),
                           setStdMargins = FALSE, zlim = c(-1, 1), 
                           main = paste0("ME ~ Covariate Association"), cex.main=1.5,
                           #verticalSeparator.x = NULL,
                           cex.lab.x=1, cex.lab.y=1.5, cex.lab=0.9,
                           horizontalSeparator.col = "black", horizontalSeparator.y = c(0,1,2,3,4,5), 
                           verticalSeparator.col = "black")
            grid.echo()
            g.cor <- grid.grab()
            grid.arrange(g.cor) 
      # Save figure ####
      pdf(file = paste0(networkdir, "/Mod", i, "_",c,".pdf"), w=18,h=20)
        grid.arrange(grobs=list(g.cor, g.go.ms.res, g.kME, g.PPI, gg.gene.list.plot, overlap.size, g.boxplot),
                     layout_matrix=rbind(c(7,7,7,4,1,1), 
                                         c(7,7,7,4,1,1), 
                                         c(7,7,7,4,1,1), 
                                         c(5,5,5,6,6,6),
                                         c(5,5,5,6,6,6),
                                         c(2,2,2,2,3,3),
                                         c(2,2,2,2,3,3)), 
                     padding=unit(1,"line"), vp=viewport(width=0.95, height=0.95))
      dev.off()
      
      # save gene list & kME
      mod_genes = data.frame(moduleName = m, moduleColor = c, 
              moduleGenes = list_of_genes[[i]]$module.symbols, module.symbols = list_of_genes[[i]]$module.symbols, 
              module.hsapiens_homolog_ensembl_gene = list_of_genes[[i]]$module.hsapiens_homolog_ensembl_gene,
              module.hsapiens_homolog_associated_gene_name = list_of_genes[[i]]$module.hsapiens_homolog_associated_gene_name, kME = list_of_genes[[i]]$kME)
      write.csv(mod_genes, file = paste0(networkdir, "/Mod", i, "_",c,"_genelist.csv"))

    #rm(me, m, c, a, res.table, g.boxplot, dat.sub, g.pvalTable, g.anovaTable)
      
  } else {
  }
}
