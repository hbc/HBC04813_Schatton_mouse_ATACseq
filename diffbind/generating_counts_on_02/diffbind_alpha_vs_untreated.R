# Load libraries
library(SummarizedExperiment)
library(DESeq2)
library(tidyverse)
library(Category)
library(ChIPpeakAnno)
library(DiffBind)




samples <- read.csv("/n/data1/cores/bcbio/PIs/tobias_schatton/schatton_atacseq_hbc04813/diffbind/bowtie2/alpha_vs_untreated_mergesplit.csv")

brca <- dba(sampleSheet=samples)


brca_count <- dba.count(brca, bUseSummarizeOverlaps=TRUE,bParallel=T,summits=75)
write_rds(brca_count, "/n/data1/cores/bcbio/PIs/tobias_schatton/schatton_atacseq_hbc04813/diffbind/bowtie2/alpha_vs_untreated_mergesplit_counts.RDS")


samples2 <- read.csv("/n/data1/cores/bcbio/PIs/tobias_schatton/schatton_atacseq_hbc04813/diffbind/bowtie2/alpha_vs_untreated_singlesplit.csv")

brca2 <- dba(sampleSheet=samples2)


brca_count2 <- dba.count(brca2, bUseSummarizeOverlaps=TRUE,bParallel=T,summits=75)
write_rds(brca_count2, "/n/data1/cores/bcbio/PIs/tobias_schatton/schatton_atacseq_hbc04813/diffbind/bowtie2/alpha_vs_untreated_singlesplit_counts.RDS")

