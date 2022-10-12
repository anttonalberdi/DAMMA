gene_annotations <- read.table("data/source/gene_annotations.tsv",header=TRUE,sep="\t")
library(gtools)
gene_annotations <- gene_annotations[mixedorder(gene_annotations$genome),]
gene_expression <- read.table("data/source/gene_expression.tsv",header=TRUE,sep="\t",row.names=1)
genome_counts <- read.table("data/source/genome_counts.tsv",header=TRUE,sep="\t",row.names=1)
genome_quality <- read.table("data/source/genome_quality.tsv",header=TRUE,sep="\t")

pathway_table_FD1 <- read.table("data/source/DAMMA_db1.tsv",header=TRUE,sep="\t")
pathway_table_FD2 <- read.table("data/source/DAMMA_db2.tsv",header=TRUE,sep="\t")
pathway_table_FD3 <- read.table("data/source/DAMMA_db3.tsv",header=TRUE,sep="\t")
pathway_table_FD4 <- read.table("data/source/DAMMA_db4.tsv",header=TRUE,sep="\t")
pathway_table_FD5 <- read.table("data/source/DAMMA_db5.tsv",header=TRUE,sep="\t")
pathway_table_FD6 <- read.table("data/source/DAMMA_db6.tsv",header=TRUE,sep="\t")
pathway_table_FD7 <- read.table("data/source/DAMMA_db7.tsv",header=TRUE,sep="\t")

#Use latest version as default
pathway_table <- pathway_table_FD7
save(pathway_table,gene_annotations,gene_expression,genome_counts,genome_quality,
  pathway_table_FD1,
  pathway_table_FD2,
  pathway_table_FD3,
  pathway_table_FD4,
  pathway_table_FD5,
  pathway_table_FD6,
  pathway_table_FD7,
  file="data/damma_data.RData")

#Remove and update DAMMA
remove.packages("DAMMA")

#
library(roxygen2)
roxygenize()

library(devtools)
load_all('/Users/anttonalberdi/github/DAMMA')
