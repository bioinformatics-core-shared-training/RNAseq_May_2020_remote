library(EnsDb.Mmusculus.v79)
library(DESeq2)
library(tidyverse)

# Annotation

load("Robjects/DE.RData")

columns(EnsDb.Mmusculus.v79)
keytypes(EnsDb.Mmusculus.v79)

ourCols <- c("SYMBOL", "GENEID", "ENTREZID")
ourKeys <- rownames(resLvV)[1:1000]

annot <- AnnotationDbi::select(EnsDb.Mmusculus.v79,
                               keys = ourKeys,
                               columns = ourCols,
                               keytype = "GENEID")
head(annot)

length(unique(annot$ENTREZID))
dim(annot)

annot %>%
  add_count(GENEID) %>%
  dplyr::filter(n>1)

# Challenge 1

ourCols <- c("SYMBOL", "GENEID", "ENTREZID", "GENEBIOTYPE")
ourKeys <- rownames(resLvV)

annot <- AnnotationDbi::select(EnsDb.Mmusculus.v79,
                               keys = ourKeys,
                               columns = ourCols,
                               keytype = "GENEID")
head(annot)
multiple <- annot %>%
  add_count(GENEID) %>%
  dplyr::filter(n>1)
length(unique(multiple$SYMBOL))

load("Robjects/Ensembl_annotations.RData")
colnames(ensemblAnnot)

annotLvV <- as.data.frame(resLvV) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(ensemblAnnot, "GeneID") %>% 
  rename(logFC=log2FoldChange, FDR=padj)

write_tsv(annotLvV, "results/VirginVsLactating_Results_Annotated.txt")

annotLvV %>%
  arrange(FDR) %>%
  head(10)

# Visualisation

ddsShrink <- lfcShrink(ddsObj, coef = "Status_lactate_vs_virgin")

shrinkLvV <- as.data.frame(ddsShrink) %>%
  rownames_to_column("GeneID") %>%
  left_join(ensemblAnnot, "GeneID") %>%
  rename(logFC=log2FoldChange, FDR=padj)

#pvalue histogram
hist(shrinkLvV$pvalue)

#ma plots
plotMA(ddsShrink, alpha = 0.05)

cutoff <- sort(shrinkLvV$pvalue)[10]
shrinkLvV <- shrinkLvV %>%
  mutate(TopGeneLabel=ifelse(pvalue<=cutoff, Symbol, ""))

ggplot(shrinkLvV, aes(x = log2(baseMean), y = logFC)) +
  geom_point(aes(colour = FDR < 0.05), shape=20, size = 0.5) +
  geom_text(aes(label=TopGeneLabel)) +
  labs(x = "mean of normalised counts", y= "log fold change")

# Volcano plot (challenge 2)

filtTab <- shrinkLvV %>%
  filter(!is.na(FDR)) %>%
  mutate(`-log10(FDR)`= -log10(FDR))

ggplot(filtTab, aes(x = logFC, y = `-log10(FDR)`)) +
  geom_point(aes(colour=FDR < 0.05), size=1)

# Heatmaps

library(ComplexHeatmap)
library(circlize)

sigGenes <- as.data.frame(shrinkLvV) %>%
  top_n(150, wt=-FDR) %>%
  pull("GeneID")

plotDat <- vst(ddsObj)[sigGenes,] %>%
  assay()
z.mat <- t(scale(t(plotDat), center = TRUE, scale = TRUE))

myPalette <- c("red3", "ivory", "blue3")
myRamp <- colorRamp2(c(-2,0,2), myPalette)

Heatmap(z.mat, name = "z-score",
        col = myRamp,
        show_row_names = FALSE,
        cluster_columns = FALSE)


# cluster the data and split the tree
hcDat <- hclust(dist(z.mat))
cutGroups <- cutree(hcDat, h=4)

ha1 = HeatmapAnnotation(df = colData(ddsObj)[,c("CellType", "Status")])

Heatmap(z.mat, name = "z-score",
        col = myRamp,            
        show_row_name = FALSE,
        cluster_columns = FALSE,
        split=cutGroups,
        rect_gp = gpar(col = "darkgrey", lwd=0.5),
        top_annotation = ha1)

# Don't forget to check out more plots in the extended materials section on the website!


