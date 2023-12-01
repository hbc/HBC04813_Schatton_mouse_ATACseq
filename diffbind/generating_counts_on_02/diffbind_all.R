# Load libraries
library(SummarizedExperiment)
library(DESeq2)
library(tidyverse)
library(Category)
library(ChIPpeakAnno)
library(DiffBind)


samples <- read.csv("/n/data1/cores/bcbio/PIs/tobias_schatton/schatton_atacseq_hbc04813/diffbind/bowtie2/diffbind_all.csv")

brca <- dba(sampleSheet=samples)


brca_count <- dba.count(brca, bUseSummarizeOverlaps=TRUE,bParallel=T,summits=75)
write_rds(brca_count, "/n/data1/cores/bcbio/PIs/tobias_schatton/schatton_atacseq_hbc04813/diffbind/bowtie2/all_counts.RDS")


