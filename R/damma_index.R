#' Aggregates raw fullness values into meaningful indices at the compound level
#'
#' @param MCI_table Function distillation table outputted by damma() or damma_correction
#' @param functions_table Table containing definitions and metadata of metabolic functions (provided by GAMMA)
#' @param weighed Whether to weigh average MCI values according to their hierarchical classification
#' @return A MCI value for each genome or community
#' @import tidyverse
#' @import dplyr
#' @import tibble
#' @examples
#' damma_compounds(compound_table)
#' @export

damma_index <- function(MCI_table,pathway_table){

  #Detect input type
  if(any(colnames(MCI_table) %in% pathway_table$Code)){level="P"}
  if(any(colnames(MCI_table) %in% pathway_table$Compound)){level="C"}
  if(any(colnames(MCI_table) %in% pathway_table$Function)){level="F"}

  if(level == "P"){
    MCI_table <- damma_functions(damma_compounds(MCI_table,pathway_table),pathway_table)
  }else if(level == "C"){
    MCI_table <- damma_functions(MCI_table,pathway_table)
  }else if(level == "F"){
    MCI_table <- MCI_table
  }else{
    stop("MCI table and pathway table do not coincide")
  }

  MCI <- rowMeans(MCI_table)

  return(MCI)

}
