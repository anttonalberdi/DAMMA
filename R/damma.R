#' Generate the fullness table of metabolic pathways/modules from a MAG annotation table
#'
#' @param annotations Table containing MAG identifiers and annotation codes
#' @param functions Table containing definitions and metadata of metabolic functions (provided by GAMMA)
#' @param magcol Column index (number) of the annotations table containing the MAG identifiers
#' @param keggcol Column index(es) of the annotations table in which to search for KEGG KO annotations
#' @param eccol Column index(es) of the annotations table in which to search for Enzyme Commision (EC) annotations
#' @param pepcol Column index(es) of the annotations table in which to search for Peptidase annotations
#' @import stringr
#' @import dplyr
#' @return A fullness matrix
#' @examples
#' damma(annotations,functions,magcol,keggcol,eccol,pepcol)
#' damma(annotations,functions,magcol=2,keggcol=9,eccol=c(10,19),pepcol=12)
#' @export

damma <- function(annotations,functions,magcol,keggcol,eccol,pepcol){

  #Simplify annotations table
  annotations2 <- annotations[,c(magcol,keggcol,eccol,pepcol), with=FALSE]
  colnames(annotations2) <- c("MAGs",paste0("K",c(1:length(keggcol))),paste0("E",c(1:length(eccol))),paste0("P",c(1:length(pepcol))))

  #annotations2[ , Annotations := do.call(paste, c(.SD, sep = " ")), .SDcols=-1]
  #annotations2 <- annotations2[, .(MAGs, Annotations)]

  #List MAGs
  MAGs <- unique(dplyr::pull(annotations2, MAGs))

  #Calculate fullness values
  fullness_table <- c()
  cat("Calculating fulness values for MAG:\n")
  m=0
  for(MAG in MAGs){
    m=m+1
    cat("\t",MAG," (",m,"/",length(MAGs),")\n", sep = "")
    #Fetch MAG annotations
    annotations_MAG <- annotations2[annotations2$MAGs == MAG]
    #K00000
    kegg <- str_extract(annotations_MAG$K1, "K[0-9]+")
    kegg <- unique(kegg[!is.na(kegg)])
    #[EC:0.0.0.0]
    EC <- unlist(str_match_all(annotations_MAG$E1, "(?<=\\[EC:).+?(?=\\])"))
    EC <- unique(unlist(strsplit(EC, " ")))
    EC <- EC[!grepl("-", EC, fixed = TRUE)]
    #(EC 0.0.0.0)
    EC2 <- unlist(str_match_all(annotations_MAG$E2, "(?<=\\(EC ).+?(?=\\))"))
    EC2 <- unique(unlist(strsplit(EC2, " ")))
    EC2 <- EC2[!grepl("-", EC2, fixed = TRUE)]
    #Peptidases
    pep <- unique(annotations_MAG$P1)
    pep <- pep[pep != ""]
    #Concatenate all annotations
    present <- unique(c(kegg,EC,EC2,pep))
    #Compute completeness scores
    fullness_vector <- c()
    for(f in c(1:nrow(functions))){
      definition=functions[f,5]
      fullness <- compute_fullness(definition,present)
      fullness_vector <- c(fullness_vector,fullness)
    }
    fullness_table <- rbind(fullness_table,fullness_vector)
  }
  rownames(fullness_table) <- MAGs
  colnames(fullness_table) <- rownames(functions)
  return(fullness_table)
}
