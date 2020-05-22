library(tidyverse)
library(DESeq2)

load("Robjects/preprocessing.RData")

sampleinfo
dim(countdata)

# Sample meta data & Count data
# What we need now - linear model
#
# y = mx + c
# y ~ mx + c

# Let's apply a simple model iwth just Status

design <- as.formula(~ Status)
design
class(design)

# model matrix

model.matrix(design, data = sampleinfo)

# change status to factor

sampleinfo$Status <- factor(sampleinfo$Status, level = c("virgin", "pregnant", "lactate"))
model.matrix(design, data = sampleinfo)

# Build a DESeq2 object

ddsObj.raw <- DESeqDataSetFromMatrix(countData = countdata,
                                     colData = sampleinfo,
                                     design = design)


# The DESeq2 workflow

ddsObj <- estimateSizeFactors(ddsObj.raw)

colData(ddsObj.raw)
colData(ddsObj)

## have a look at the scaling factors

logcounts <- log2(countdata + 1)
limma::plotMA(logcounts)
abline(h=0, col="red")

normalizedCounts <- counts(ddsObj, normalized = TRUE)
logNorm <- log2(normalizedCounts + 1)

limma::plotMA(logNorm)
abline(h=0, col="blue")

# Estimate dispersion

ddsObj <- estimateDispersions(ddsObj)
plotDispEsts(ddsObj)

# Apply negative binomial model and test

ddsObj <- nbinomWaldTest(ddsObj)


## The DESeq command

ddsObj <- DESeq(ddsObj.raw)

## Generating a results table

res <- results(ddsObj, alpha = 0.05)
res

model.matrix(design, data=sampleinfo)


# look at the coefficients

resultsNames(ddsObj)

resPvV <- results(ddsObj,
                  name = "Status_pregnant_vs_virgin",
                  alpha = 0.05)
resPvV

# select top 100 genes

topGenes <- as.data.frame(resPvV) %>%
  rownames_to_column("GeneID") %>%
  arrange(padj) %>% 
  head(100)
topGenes

# Exercise 1

design <- as.formula(~ CellType + Status)
ddsObj.raw <- DESeqDataSetFromMatrix(countData = countdata, 
                                     colData = sampleinfo,
                                     design = design)
ddsObj <- DESeq(ddsObj.raw)
resLvV <- results(ddsObj, alpha=0.05)
resLvV

## Challenge 1

resLvB <- results(ddsObj, name = "CellType_luminal_vs_basal", alpha = 0.05)
resLvB

sum(resLvB$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)


# Challenge 2 - contrasts

# as ratios P/V  and L/V
# if we want is P/L = (P/V) / (L/V)

resPvL <- results(ddsObj,
                  contrast = c("Status", "pregnant", "lactate"),
                  alpha = 0.05)
resPvL

## PCA

vstcounts <- vst(ddsObj.raw, blind=TRUE)
plotPCA(vstcounts, intgroup=c("Status", "CellType"))


## I think that only CellType is important

#### create the simpler (reduced) model

design.reduced <- as.formula(~CellType)

# ~CellType + Status
ddsObjC <- DESeq(ddsObj, test = "LRT", reduced = design.reduced)

resCvCS <- results(ddsObjC)
resCvCS

sum(resCvCS$padj < 0.05, na.rm = TRUE)

# Exercise 2

# 1. Create a new DESeq2 object using a model with an interaction between CellType and Status (~ CellType * Status)

design <- as.formula(~ CellType + Status + CellType:Status)

ddsObj2.raw <- DESeqDataSetFromMatrix(countData = countdata, 
                                      colData = sampleinfo,
                                      design = design)

# 2. Use the LRT to compare this to the simpler additive model (~CellType + Status)

# create the simpler model
design.reduced <- as.formula(~ CellType + Status)
# Compare the two designs
ddsObjC2 <- DESeq(ddsObj2.raw, test="LRT", reduced=design.reduced)
# get the results table
resCSvCxS <- results(ddsObjC2)
# For how many genes is interaction model a better fit?
sum(resCSvCxS$padj < 0.05, na.rm=TRUE)





