# pre-processing

# use chat for quick questions?
# copy/paste own code?

#-----------------------------------------
# a RNA-seq data set, its parts
# QC, read count distribution, profiles
# DESeq2 and normalisation
#-----------------------------------------

# Mouse mammary gland
#-----------------------------------------
# Cell proliferation, differentiation and death in mammary epithelium
# MCL-1 (mitochondria function, targetted by some drug against cancer)
# https://pubmed.ncbi.nlm.nih.gov/25730472/
# https://pubmed.ncbi.nlm.nih.gov/30425521/

# two types of cell:
# - basal stem-cell enriched cells (B)
# - committed luminal cells (L)

# three types of mice:
# - virgin
# - pregnant
# - lactating

# so 2 x 3 == 6 classes
# 2 replicates per class.
# 12 samples in total

library(DESeq2)
library(tidyverse)

# sample info
#-----------------------------------------
getwd()
dir()
dir("data")

sampleinfo <- read_tsv("data/SampleInfo.txt") 
sampleinfo

# read counts
#-----------------------------------------
# HISAT2
# CRGm38
# featureCounts

seqdata <- read_tsv("data/GSE60450_Lactation.featureCounts", comment="#")
seqdata

#-----------------------------------------
# dplyr for data manipulation
#-----------------------------------------
# simpler, more readable code

# only keep basal cells
# only keep CellType and Group
# rename CellType to Cell

newTable <- sampleinfo
# only keep basal cells
newTable <- filter(newTable, CellType=="basal")
newTable
# only keep CellType and Group
newTable <- select(newTable, CellType, Group)
# rename CellType to Cell
newTable <- rename(newTable, Cell=CellType)
newTable

# with pipe %>%
newTable <- sampleinfo %>% filter(CellType=="basal") %>%
  select(CellType, Group) %>%
  rename(Cell=CellType)
newTable

#-----------------------------------------
# format data to use DESeq2
#-----------------------------------------
# count matrix
# sample info

# count matrix
#-----------------------------------------
head(seqdata)

# have row named after genes
# keep counts only and rename headers
# see stringr

countdata <- seqdata %>%
  column_to_rownames("Geneid") %>% 
  rename_all(str_remove, ".bam") %>%
  select(sampleinfo$Sample) %>%
  as.matrix()
countdata
dim(countdata)

# filter genes
#-----------------------------------------
# aim: exclude uninformative genes, multiple correction, size
# which genes? with few reads
# DESeq2 does its own: independent filtering

# keep genes with at least 6 reads
dim(countdata)
# count total number of reads per gene
# have vector to store test
keep <- rowSums(countdata) > 5
table(keep, useNA="always")
summary(rowSums(countdata))
# subset
countdata <- countdata[keep,]
dim(countdata)
head(countdata)

#-----------------------------------------
# Quality check
#-----------------------------------------
# library size
# count distribution
# expression profiles

# library size
#-----------------------------------------
librarySizes <- colSums(countdata)
librarySizes
summary(librarySizes)

barplot(librarySizes,
        names=names(librarySizes),
        las=2,
        main="library sizes"
        )
abline(h=20e6,lty=2)

# count distribution - boxplots
#-----------------------------------------
# log2 scale
# color

logcounts <- log2(countdata + 1)
class(logcounts)
summary(logcounts[,1])
summary(logcounts)

# colors
statusCol <- match(sampleinfo$Status, c("virgin", "pregnant", "lactate")) + 1

boxplot(logcounts,
        xlab="",
        ylab="log2(counts)",
        las=2,
        col=statusCol
        )
abline(h=median(logcounts), col="blue")

# Challenge 1
#-----------------------------------------
# use DESeq2 rlog() to transform the raw counts
# plot the distribution, with boxplots
# compare with plot above
# 10 min

?rlog
rlogcounts <- rlog(countdata)
boxplot(rlogcounts,
        xlab="",
        ylab="log2(counts), rlog",
        las=2,
        col=statusCol
)
abline(h=median(rlogcounts), col="blue")

#-----------------------------------------
# expression profiles - PCA
#-----------------------------------------
# unsupervised analysis to check similarity between samples here
# visually compare observation to expectation
# spot outliers, batch effect

# PCA: identify 'metagenes' that capture the largest amounts of variation in the data set
# replicates within a sample group should cluster closely
# sample groups may/should form distinct clusters

library(ggfortify)

rlogcounts <- rlog(countdata)

pcDat <- prcomp(t(rlogcounts))
class(pcDat)
autoplot(pcDat)

# color and shape
autoplot(pcDat,
  data = sampleinfo,
  colour = "CellType",
  shape = "Status",
  size = 5
)

# discussion:
# www.menti.com code 62 87 48

# label samples
library(ggrepel)
autoplot(pcDat,
         data = sampleinfo,
         colour = "CellType",
         shape = "Status",
         size = 5
) +
geom_text_repel(aes(x=PC1, y=PC2, label=Sample), box.padding=0.8)

# fix sample type swap
# MCL1.DG should be basal
# MCL1.LA should be luminal

sampleinfo <- sampleinfo %>%
  mutate(CellType=ifelse(Sample=="MCL1.DG", "basal", CellType)) %>%
  mutate(CellType=ifelse(Sample=="MCL1.LA", "luminal", CellType)) %>%
  mutate(Group=str_c(CellType, ".", Status))
sampleinfo

autoplot(pcDat,
         data = sampleinfo,
         colour = "CellType",
         shape = "Status",
         size = 5
)
