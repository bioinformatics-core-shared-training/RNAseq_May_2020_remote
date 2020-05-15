# Introduction to RNA-seq data analysis - Extended Materials
### 19th - 21st May 2020
#### Taught remotely
#### Bioinformatics Training, Craik-Marshall Building, Downing Site, University of Cambridge

![](images/CRUK_Cambridge_Major Centre logo.jpg)

## Outline

In this workshop, you will be learning how to analyse RNA-seq data.  This will
include read alignment, quality control, quantification against a reference,
reading the count data into R, performing differential expression analysis, and
gene set testing, with a focus on the DESeq2 analysis workflow. You will learn
how to generate common plots for analysis and visualisation of gene expression
data, such as boxplots and heatmaps. 

This workshop is aimed at biologists interested in learning how to perform
differential expression analysis of RNA-seq data. 

### Day 1

1. Raw read file format and QC

    - [Introductory slides](html/B_FastQC.html)  
    - [Practical](extended_html/B_FastQC_practical.html)  
    - [Practical solutions](Extended_html/B_FastQC_practical.Solutions.html)  
    - [Example of using Trimmomatic to remove adapters](extended_html/Trimming.html)

2. Short read alignment with HISAT2  

    - [Introductory slides](html/C_Alignment_with_HISAT2.html)  
    - [Practical](extended_html/C_Alignment_with_HISAT2_practical.html)  
    - [Practical solutions](extended_html/C_Alignment_with_HISAT2.Solutions.html)  

3. QC of alignment  

    - [Introductory slides](html/D_QC_of_aligned_reads.html)  
    - [Practical](extended_html/D_QC_of_aligned_reads_practical.html)  
    - [Practical solutions](extended_html/D_QC_of_aligned_reads.Solutions.html)  

4. Quantification with SubRead - Abbi Edwards  
    - [Introductory slides](html/E_Read_Counts_with_Subread.html)  
    - [Practical](extended_html/E_Read_Counts_with_Subread.practical.html)  
    - [Practical solutions](extended_html/E_Read_Counts_with_Subread.Solutions.html)  

### Day 2

1. [RNA-seq Pre-processing](extended_html/02_Preprocessing_Data.html)

    - [Practical solutions](extended_html/02_Preprocessing_Data.Solutions.html)  

### Day 3

1. [Differential Expression for RNA-seq](extended_html/04_DE_analysis_with_DESeq2.html)

    - [Practical solutions](extended_html/04_DE_analysis.Solutions.html)

2. [Annotation and Visualisation of RNA-seq results](extended_html/05_Annotation_and_Visualisation.html)

    - [practical 
solutions](extended_html/05_Annotation_and_Visualisation.Solutions.html)

3, [Gene-set testing](extended_html/06_Gene_set_testing.html)

    - [practical solutions](extended_html/06_Gene_set_testing.Solutions.html)

## Source Materials for Practicals

The all of the lecture slides and other source materials, including R code and 
practical solutions, can be found in the course's [Github 
repository](https://github.com/bioinformatics-core-shared-training/RNAseq_September_2019)

### Supplementary lessons

Introductory R materials:

- [Introduction to R](https://bioinformatics-core-shared-training.github.io/r-intro/)

Additional RNAseq materials:

- [Downloading files from SRA and aligning](Supplementary_Materials/S1_Getting_raw_reads_from_SRA.html)
- [Additional annotation and plotting](Supplementary_Materials/S3_Annotation_and_Visualisation.nb.html)

Data: Example Mouse mammary data (fastq files): 
	[https://figshare.com/s/f5d63d8c265a05618137](https://figshare.com/s/f5d63d8c265a05618137)

### Additional resources

[Bioconductor help](https://www.bioconductor.org/help/)  
[Biostars](https://www.biostars.org/)  
[SEQanswers](http://seqanswers.com/)  

## Acknowledgements

This course is based on the course [RNAseq analysis in R](http://combine-australia.github.io/2016-05-11-RNAseq/) prepared by [Combine Australia](https://combine.org.au/) and delivered on May 11/12th 2016 in Carlton. We are extremely grateful to the authors for making their materials available; Maria Doyle, Belinda Phipson, Matt Ritchie, Anna Trigos, Harriet Dashnow, Charity Law.

![](images/combine_banner_small.png)
