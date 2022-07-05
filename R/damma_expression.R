#' Generate the fullness table of metabolic pathways/modules from a MAG annotation table
#'
#' @param expression Table containing normalised gene expression data with genes in rows and samples in columns
#' @param annotations Table containing gene and MAG identifiers, and annotation codes
#' @param functions Table containing definitions and metadata of metabolic functions (provided by GAMMA)
#' @param genecol Column index (number) of the annotations table containing the gene identifiers
#' @param magcol Column index (number) of the annotations table containing the MAG identifiers
#' @param keggcol Column index(es) of the annotations table in which to search for KEGG KO annotations
#' @param eccol Column index(es) of the annotations table in which to search for Enzyme Commision (EC) annotations
#' @param pepcol Column index(es) of the annotations table in which to search for Peptidase annotations
#' @import data.table
#' @importFrom stringr str_extract str_match_all
#' @importFrom dplyr pull
#' @return A fullness matrix
#' @examples
#' damma(annotations,functions,magcol,keggcol,eccol,pepcol)
#' damma(annotations,functions,magcol=2,keggcol=9,eccol=c(10,19),pepcol=12)
#' @export

#UNDER DEVELOPMENT
damma_expression <- function(expression,annotations,functions,genecol,magcol,keggcol,eccol,pepcol){

#### UNDER DEVELOPMENT ####

  #Simplify annotations table
  setDT(annotations)
  annotations2 <- annotations[,c(genecol,magcol,keggcol,eccol,pepcol), with=FALSE]
  colnames(annotations2) <- c("Genes","MAGs",paste0("K",c(1:length(keggcol))),paste0("E",c(1:length(eccol))),paste0("P",c(1:length(pepcol))))

  #Validate expression table
  sharedgenes <- intersect(rownames(expression),annotations2$Genes)

  #Filter expression table
  expression2 <- expression[sharedgenes,]

  #Filter annotations table
  annotations3 <- annotations2[Genes %in% sharedgenes,]

  #List MAGs
  MAGs <- unique(dplyr::pull(annotations3, MAGs))

  #Calculate fullness values
  expression_table <- c()
  cat("Calculating fulness values for MAG:\n")
  m=0
  for(MAG in MAGs){
    m=m+1
    cat("\t",MAG," (",m,"/",length(MAGs),")\n", sep = "")
    #Fetch MAG annotations
    annotations_MAG <- annotations3[annotations3$MAGs == MAG]
    #K00000
    kegg <- str_extract(annotations_MAG$K1, "K[0-9]+")
    kegg <- unique(kegg[!is.na(kegg)])
    for(k in kegg){
      genes <- annotations_MAG[grep(k, annotations_MAG$K1),"Genes"]$Genes
      expression3 <- expression2[genes,]
      if(dim(expression3)[1]>1){
        expression3 <- colSums(expression3,na.rm=TRUE)
        expression3 <- t(expression3)
      }
      rownames(expression3) <- k
      expression_table <- rbind(expression_table,expression3)
    }
    #[EC:0.0.0.0]
    EC1 <- unlist(str_match_all(annotations_MAG$E1, "(?<=\\[EC:).+?(?=\\])"))
    EC1 <- unique(unlist(strsplit(EC1, " ")))
    EC1 <- EC1[!grepl("-", EC1, fixed = TRUE)]
    #(EC 0.0.0.0)
    EC2 <- unlist(str_match_all(annotations_MAG$E2, "(?<=\\(EC ).+?(?=\\))"))
    EC2 <- unique(unlist(strsplit(EC2, " ")))
    EC2 <- EC2[!grepl("-", EC2, fixed = TRUE)]
    EC <- unique(EC1,EC2)
    EC <- EC[!is.na(EC)]
    for(e in EC){
      genes1 <- annotations_MAG[(grep(e, annotations_MAG$E1)),"Genes"]$Genes
      genes2 <- annotations_MAG[(grep(e, annotations_MAG$E2)),"Genes"]$Genes
      genes <- unique(c(genes1,genes2))
      expression3 <- expression2[genes,]
      if(dim(expression3)[1]>1){
        expression3 <- colSums(expression3,na.rm=TRUE)
        expression3 <- t(expression3)
      }
      rownames(expression3) <- e
      expression_table <- rbind(expression_table,expression3)
    }
    #Peptidases
    pep <- unique(annotations_MAG$P1)
    pep <- pep[pep != ""]
    for(p in pep){
      genes <- annotations_MAG[grep(k, annotations_MAG$P1),"Genes"]$Genes
      expression3 <- expression2[genes,]
      if(dim(expression3)[1]>1){
        expression3 <- colSums(expression3,na.rm=TRUE)
        expression3 <- t(expression3)
      }
      rownames(expression3) <- k
      expression_table <- rbind(expression_table,expression3)
    }

    #Compute completeness scores TO BE UPDATED
    #OUTPUT A LIST OF TABLES RATHER THAN A TABLE
    for(f in c(1:nrow(functions))){
      definition=functions[f,"Definition"]
      cat("\tFunction ",paste0(f,"/",nrow(functions)),"\n")
      expression_fullness <- compute_fullness_expression(definition,expression_table)
      #Create list object if it is the first iteration
      if(f == 1){expression_fullness_list <- expression_fullness}
      expression_fullness_list <- Map(c, expression_fullness_list, expression_fullness)
    }
  }
  rownames(fullness_table) <- MAGs
  colnames(fullness_table) <- functions$Code
  fullness_table[is.na(fullness_table)] <- 0
  return(fullness_table)


}
