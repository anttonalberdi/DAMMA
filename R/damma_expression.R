#' Generate the fullness table of metabolic pathways/modules from a Genome annotation table
#'
#' @param expression Table containing normalised gene expression data with genes in rows and samples in columns
#' @param annotations Table containing gene and genome identifiers, and annotation codes
#' @param functions Table containing definitions and metadata of metabolic functions (provided by GAMMA)
#' @param genecol Column index (number) of the annotations table containing the gene identifiers
#' @param genomecol Column index (number) of the annotations table containing the genome identifiers
#' @param keggcol Column index(es) of the annotations table in which to search for KEGG KO annotations
#' @param eccol Column index(es) of the annotations table in which to search for Enzyme Commision (EC) annotations
#' @param pepcol Column index(es) of the annotations table in which to search for Peptidase annotations
#' @importFrom stringr str_extract str_match_all
#' @return A list of pathway-expression matrices (one table per genome)
#' @examples
#' damma_expression(expression,annotations,functions,genecol,genomecol,keggcol,eccol,pepcol)
#' damma_expression(expression,annotations,functions,genecol=1,genomecol=2,keggcol=9,eccol=c(10,19),pepcol=12)
#' @export

#UNDER DEVELOPMENT
damma_expression <- function(expression,annotations,functions,genecol,genomecol,keggcol,eccol,pepcol){

  cat("STARTING dammaE ANALYSIS\n(Note this may take a while)...\n")

  #Simplify annotations table
  annotations <- as.data.frame(annotations)
  annotations2 <- annotations[,c(genecol,genomecol,keggcol,eccol,pepcol)]
  colnames(annotations2) <- c("Genes","Genomes",paste0("K",c(1:length(keggcol))),paste0("E",c(1:length(eccol))),paste0("P",c(1:length(pepcol))))

  #Validate expression table
  sharedgenes <- intersect(rownames(expression),annotations2$Genes)

  #Filter expression table
  expression2 <- as.data.frame(expression[sharedgenes,])

  #Filter annotations table
  annotations3 <- annotations2[annotations2$Genes %in% sharedgenes,]

  #List Genomes
  Genomes <- unique(annotations3$Genomes)

  #Calculate expression values for each Genome
  cat("Calculating function expression values for Genome:\n")
  m=0
  expression_fullness_table_list <- list()
  for(Genome in Genomes){
    m=m+1
    cat("\t",Genome," (",m,"/",length(Genomes),")\n", sep = "")
    cat("\t\tProcessing KEGG annotations...\n", sep = "")
    #Fetch Genome annotations

    expression_table <- data.frame()
    annotations_Genome <- annotations3[annotations3$Genomes == Genome,]
    #K00000
    annotations_Genome2 <- annotations_Genome[order(annotations_Genome$K1),]
    annotations_Genome2 <- annotations_Genome[annotations_Genome$K1 != "",]
    kegg <- str_extract(annotations_Genome2$K1, "K[0-9]+")
    kegg <- sort(unique(kegg[!is.na(kegg)]))
    for(k in kegg){
      genes <- annotations_Genome2[grep(k, annotations_Genome2$K1),"Genes"]
      expression3 <- expression2[genes,]
      if(dim(expression3)[1]>1){
        expression3 <- colSums(expression3,na.rm=TRUE)
        expression3 <- t(expression3)
        rownames(expression3) <- k
        expression_table <- rbind(expression_table,expression3)
      }
    }

    cat("\t\tProcessing EC annotations...\n", sep = "")
    annotations_Genome2 <- annotations_Genome[order(annotations_Genome$E1,annotations_Genome$E2),]
    annotations_Genome2 <- annotations_Genome2[(annotations_Genome2$E1 != "") | (annotations_Genome2$E2 != ""),]
    #[EC:0.0.0.0]
    EC1 <- unlist(str_match_all(annotations_Genome2$E1, "(?<=\\[EC:).+?(?=\\])"))
    EC1 <- unique(unlist(strsplit(EC1, " ")))
    EC1 <- EC1[!grepl("-", EC1, fixed = TRUE)]
    #(EC 0.0.0.0)
    EC2 <- unlist(str_match_all(annotations_Genome2$E2, "(?<=\\(EC ).+?(?=\\))"))
    EC2 <- unique(unlist(strsplit(EC2, " ")))
    EC2 <- EC2[!grepl("-", EC2, fixed = TRUE)]
    EC <- unique(EC1,EC2)
    EC <- sort(EC[!is.na(EC)])
    for(e in EC){
      genes1 <- annotations_Genome2[(grep(e, annotations_Genome2$E1)),"Genes"]
      genes2 <- annotations_Genome2[(grep(e, annotations_Genome2$E2)),"Genes"]
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
    pep <- unique(annotations_Genome$P1)
    pep <- pep[(pep != "") & (!is.na(pep))]
    annotations_Genome2 <- annotations_Genome[order(annotations_Genome$P1),]
    annotations_Genome2 <- annotations_Genome2[annotations_Genome2$P1 != "",]
    for(p in pep){
      genes <- annotations_Genome2[grep(p, annotations_Genome2$P1),"Genes"]
      expression3 <- expression2[genes,]
      if(dim(expression3)[1]>1){
        expression3 <- colSums(expression3,na.rm=TRUE)
        expression3 <- t(expression3)
        rownames(expression3) <- p
        expression_table <- rbind(expression_table,expression3)
      }
    }

    #Compute expression scores
    cat("\t\tCalculating expression metrics...\n")
    suppressWarnings(
      for(f in c(1:nrow(functions))){
        definition=functions[f,"Definition"]
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
    colnames(expression_fullness_table) <- functions$Code

    #Append to Genome list
    expression_fullness_table_list[[Genome]] <- expression_fullness_table
  }

  return(expression_fullness_table_list)

}
