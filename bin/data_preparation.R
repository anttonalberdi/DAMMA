annotations_example <- fread("data/annotations.tsv")
functions_table_FD1 <- read.table("data/DAMMA_functions_FD1.tsv",header=TRUE,sep="\t")
functions_table_FD2 <- read.table("data/DAMMA_functions_FD2.tsv",header=TRUE,sep="\t")
functions_table_FD3 <- read.table("data/DAMMA_functions_FD3.tsv",header=TRUE,sep="\t")

#Use latest version as default
functions_table <- functions_table_FD3
save(functions_table,annotations_example,
  functions_table_FD1,
  functions_table_FD2,
  functions_table_FD3,
  file="data/damma_data.RData")
