#' Generate community-level MCIs from microbial gene expression data
#'
#' @param annotations Table containing Genome identifiers and annotation codes
#' @param functions_table Table containing definitions and metadata of metabolic functions (included in DAMMA package)
#' @param expression Gene expression data
#' @param fullness_table_list A list of pathway-expression matrices (one table per genome)
#' @param genecol Column index (number) of the annotations table containing the gene identifiers
#' @param genomecol Column index (number) of the annotations table containing the genome identifiers
#' @param keggcol Column index(es) of the annotations table in which to search for KEGG KO annotations
#' @param eccol Column index(es) of the annotations table in which to search for Enzyme Commision (EC) annotations
#' @param pepcol Column index(es) of the annotations table in which to search for Peptidase annotations
#' @importFrom stringr str_extract str_match_all
#' @return A pathway fullness vector (if no abundance data are provided) or matrix (if abundance data are provided)
#' @examples
#' damma_community_expression(annotations,functions_table,expression,fullness_table_list,genecol,genomecol,keggcol,eccol,pepcol)
#' damma_community_expression(annotations=gene_annotations,functions_table,expression=gene_expression,genecol=1,genomecol=2,keggcol=9,eccol=c(10,19),pepcol=12)
#' @export

damma_community_expression <- function(annotations,functions_table,expression,fullness_table_list,genecol,genomecol,keggcol,eccol,pepcol){

annotations=gene_annotations
functions_table
expression=gene_expression[,-1]
rownames(expression)=gene_expression[,1]
genecol=1
genomecol=2
keggcol=9
eccol=c(10,19)
pepcol=12

  #Prepare tables
  annotations <- as.data.frame(annotations)
  annotations2 <- annotations[,c(genecol,genomecol,keggcol,eccol,pepcol)]
  colnames(annotations2) <- c("Genes","Genomes",paste0("K",c(1:length(keggcol))),paste0("E",c(1:length(eccol))),paste0("P",c(1:length(pepcol))))
  functions_table2 <- as.data.frame(functions_table)

  #####
  # Generate Genome- and Community-level fullness scores
  #####

  if(missing(fullness_table_list)){
    # If fullness table is not provided, run damma_expression() function to create it
    cat("Running DAMMA expression script for creating the eMCI table...\n")
    cat("\tNote: if you have already run DAMMA expression on this data set,\n")
    cat("\tyou can include the list of fullness tables to avoid this step.\n")
    fullness_table_list2 <- damma_expression(expression,annotations2,functions_table,genecol=1,genomecol=2,keggcol=3,eccol=c(4,5),pepcol=6)
  }else{
    cat("Fullness table has been provided.\n")
    fullness_table_list2 <- fullness_table_list
  }

  #Compute overall community-level fullness table
  cat("Computing overall community-level eMCIs...\n")
  cat("\tExtracting annotations...\n")

  #Validate expression table
  sharedgenes <- intersect(rownames(expression),annotations2$Genes)

  #Filter expression table
  expression2 <- expression[sharedgenes,]

  #Filter annotations table
  annotations3 <- annotations2[annotations2$Genes %in% sharedgenes,]

  expression_table <- data.frame()
  #K00000
  annotations3 <- annotations3[order(annotations3$K1),]
  kegg <- str_extract(annotations3$K1, "K[0-9]+")
  kegg <- sort(unique(kegg[!is.na(kegg)]))
  for(k in kegg){
    genes <- annotations3[grep(k, annotations3$K1),"Genes"]
    expression3 <- expression2[genes,]
    if(dim(expression3)[1]>1){
      expression3 <- colSums(expression3,na.rm=TRUE)
      expression3 <- t(expression3)
      rownames(expression3) <- k
      expression_table <- rbind(expression_table,expression3)
    }
  }

  cat("\t\tProcessing EC annotations...\n", sep = "")
  annotations3 <- annotations3[order(annotations3$E1,annotations3$E2),]
  #[EC:0.0.0.0]
  EC1 <- unlist(str_match_all(annotations3$E1, "(?<=\\[EC:).+?(?=\\])"))
  EC1 <- unique(unlist(strsplit(EC1, " ")))
  EC1 <- EC1[!grepl("-", EC1, fixed = TRUE)]
  #(EC 0.0.0.0)
  EC2 <- unlist(str_match_all(annotations3$E2, "(?<=\\(EC ).+?(?=\\))"))
  EC2 <- unique(unlist(strsplit(EC2, " ")))
  EC2 <- EC2[!grepl("-", EC2, fixed = TRUE)]
  EC <- unique(EC1,EC2)
  EC <- sort(EC[!is.na(EC)])
  for(e in EC){
    genes1 <- annotations3[(grep(e, annotations3$E1)),"Genes"]
    genes2 <- annotations3[(grep(e, annotations3$E2)),"Genes"]
    genes <- unique(c(genes1,genes2))
    expression3 <- expression2[genes,]
    if(dim(expression3)[1]>1){
      expression3 <- colSums(expression3,na.rm=TRUE)
      expression3 <- t(expression3)
      rownames(expression3) <- e
      expression_table <- rbind(expression_table,expression3)
    }
  }

  cat("\t\tProcessing peptidase annotations...\n", sep = "")
  #Peptidases
  pep <- unique(annotations3$P1)
  pep <- pep[(pep != "") & (!is.na(pep))]
  annotations3 <- annotations3[order(annotations3$P1),]
  for(p in pep){
    genes <- annotations3[grep(k, annotations3$P1),"Genes"]
    expression3 <- expression2[genes,]
    if(dim(expression3)[1]>1){
      expression3 <- colSums(expression3,na.rm=TRUE)
      expression3 <- t(expression3)
      rownames(expression3) <- k
      expression_table <- rbind(expression_table,expression3)
    }
  }

  #Compute expression scores
  cat("\t\tCalculating community expression metrics...\n")
  suppressWarnings(
    for(f in c(1:nrow(functions_table))){
      definition=functions_table[f,"Definition"]
      #cat("\tFunction ",paste0(f,"/",nrow(functions)),"\n")
      expression_fullness <- compute_fullness_expression(definition,expression_table)
      if(f == 1){
        #Create list if it is the first function
        expression_fullness_list <- expression_fullness
      }else{
        #Append to list if it is not the first function
        expression_fullness_list <- Map(c, expression_fullness_list, expression_fullness)
      }
    }
  )

  #Convert sample list to matrix
  expression_fullness_list <- lapply(expression_fullness_list,function(x) as.numeric(x))
  expression_fullness_table <- do.call(rbind, expression_fullness_list)
  colnames(expression_fullness_table) <- functions_table$Code


  ########
  # FROM HERE ONWARDS TO BE UPDATED
  ########


  #####
  # Generate Genome- and Community-level fullness scores
  #####

    cat("Weighing fullness scores by relative abundances...\n")
    samples <- colnames(expression)

    community_table <- c()
    for(s in samples){
      relabun <- abundance_table2[,s]
      weigedfullness <- sweep(fullness_table2, 1, relabun, FUN = "*")
      community_row <- colSums(weigedfullness) * community_fullness_vector
      community_table <- rbind(community_table,community_row)
    }
    rownames(community_table) <- samples

  return(community_table)
}
