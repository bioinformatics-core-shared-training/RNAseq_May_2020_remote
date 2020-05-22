#------------------------------------
# Aim: include biological knowledge in differential expression analysis
#------------------------------------

# standard DE considers genes independently, ignores correlation between genes.
# little overlap between the top ranked genes for the two data sets
# due to difference in sample collection, preparation, assay, etc

# small differences (20%) in separate genes can combine into strong correlation
# between gene set and sample group.
# less senstitive to technical and sampling differences

# over-representation
#-----------------------
# get list of DE genes (DEGs)
# for a given gene set,
# - make 2x2 table
# - to compare the number of DEGs in that gene set to those that are not
# - eg test significance with hypergeometric

# cons: depends on cut-off used to define DEG list.

# aggregate score
#-----------------------

# assign scores to each geneset based on all the gene-specific scores for that gene set
# (KS test)

# leading edge: genes in the gene set with the highest correlation to the trait
# "leading-edge subset to be those genes in the gene set S that appear in the ranked list L
# at, or before, the point where the running sum reaches its maximum deviation from zero."

# normalizes the ES based on the number of gene members in the data set,
# for comparison across gene sets of different sizes

# greater weight is placed on genes with stronger correlation with the phenotype

# (FDR) that is based upon the distribution of results during repetitive,
# random assignments of class designations

# gene-wise permutations
# “This approach is not strictly accurate because it ignores gene-gene correlations
# and will overestimate the significance levels and may lead to false positives.”

#------------------------------------
# GSEA - gene set enrichment analysis
#------------------------------------

# decide on the trait to analyse
# contrast LvV

# 1 - rank genes using a metric measuring
# the correlation of expression with the trait
# shrunken LFC
# 2 - identify rank of the genes in the gene set analysed
# 3 - calculate enrichment score

library(tidyverse)
library(fgsea)

# load DE results
#-----------------------
load("Robjects/Annotated_Results_LvV.RData")
# keeps annotLvV and shrinkLvV

# rank the genes on LFC
#-----------------------
gseaDat <- filter(shrinkLvV, !is.na(Entrez))

rankData <- gseaDat$logFC
rankData
names(rankData) <- gseaDat$Entrez
rankData

# load pathways
#-----------------------
load("Robjects/mouse_H_v5.RData")
pathwaysH <- Mm.H

# conduct analysis
#-----------------------
fgseaRes <- fgsea(pathwaysH, 
                  rankData, 
                  minSize = 15, 
                  maxSize = 500, 
                  nperm = 1000)
# check output
fgseaRes
fgseaRes %>%
  arrange(desc(abs(NES))) %>%
  as_tibble()

# plot
plotEnrichment(pathwaysH[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]], rankData)

# table plot 
topPathways <- fgseaRes %>% 
  top_n(20, wt=-padj) %>% 
  arrange(-NES) %>% 
  pull(pathway)
topPathways

plotGseaTable(pathwaysH[topPathways], 
              rankData, 
              fgseaRes, 
              gseaParam = 0.5)

# challenge 1
#-----------------------
# rank the genes by significance,
# and sign of LFC,
# to have upregulated genes first
# and downregulated genes last 

# hint -log10( PVAL ) * sign( LFC )

library(tidyverse)
library(fgsea)

# DE results
load("Robjects/Annotated_Results_LvV.RData")
# keeps annotLvV and shrinkLvV

# rank the genes on LFC
gseaDat <- filter(shrinkLvV, !is.na(Entrez))

ranks <- -log10(gseaDat$pvalue) * sign(gseaDat$logFC) # or gseaDat$stat
names(ranks) <- gseaDat$Entrez

# load pathways
load("Robjects/mouse_H_v5.RData")
pathwaysH <- Mm.H

# run
fgseaResH <- fgsea(pathwaysH,
                   ranks,
                   minSize = 15,
                   maxSize = 500,
                   nperm = 1000)
# view, or:
fgseaResH %>% 
  top_n(20, wt=-padj) %>% 
  arrange(-NES) %>%
  as_tibble() %>%
  pull(pathway)

# plot NES
plot(fgseaRes$NES, fgseaResH$NES)
abline(a=0, b=1, lty=2)

# plot padj
plot(-log10(fgseaRes$padj), -log10(fgseaResH$padj))
abline(a=0, b=1, lty=2)
abline(v=-log10(0.05), lty=3)
abline(h=-log10(0.05), lty=3)

#------------------------------------
# skipping goseq
#------------------------------------
# see material

#------------------------------------
# clusterProfiler
#------------------------------------
library(clusterProfiler)

# over-representation test
# mention other tests, eg GSEA

# select differentially expressed genes DEGs
#----------------------
# shrinkLvV + filter() or results(, lfcThreshold = xx) then shrink?

sigGenes <- shrinkLvV %>% 
  filter(FDR < 0.05 & !is.na(FDR) & 
           abs(logFC) > 1 & 
           !is.na(Entrez)) %>% 
  pull(Entrez)
sigGenes

# test KEGG pathways
#----------------------
# explain KEGG
kk <- enrichKEGG(gene = sigGenes, organism = 'mmu')
head(kk, n=10) %>% as_tibble()

browseKEGG(kk, 'mmu03320')

# pathview
library(pathview)

# write plot to file in working directory,
# check with getwd()
getwd()
pathview(gene.data = logFC, 
         pathway.id = "mmu03320", 
         species = "mmu", 
         limit = list(gene=5, cpd=1))

# challenge
#----------------------
# use genes that are statistically significant at FDR < 0.01
# use 'mmu04060' pathway

library(clusterProfiler)
library(pathview)
sigGenes <- shrinkLvV$FDR < 0.01 & !is.na(shrinkLvV$FDR)

logFC <- shrinkLvV$logFC[sigGenes]
names(logFC) <- shrinkLvV$Entrez[sigGenes]

# write plot to file in working directory,
# check with getwd()
getwd()
pathview(gene.data = logFC, 
         pathway.id = "mmu04060", 
         species = "mmu", 
         limit = list(gene=5, cpd=1))

browseKEGG(kk, 'mm04060')
# check the clusterProfiler documentation at
# https://yulab-smu.github.io/clusterProfiler-book
