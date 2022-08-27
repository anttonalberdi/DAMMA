#' Aggregates compound-level MCIs into function-level MCIs
#'
#' @param MCI_table Function distillation table outputted by damma() or damma_correction
#' @param pathway_table Table containing definitions and metadata of metabolic functions (provided by GAMMA)
#' @param transform Whether to transform to 0-1 scale. Default=TRUE
#' @import tidyverse
#' @import dplyr
#' @return A MCI table aggregated at the compound level
#' @examples
#' dama_functions(compound_table)
#' @export


damma_functions <- function(MCI_table,pathway_table){

  #Detect input type
  if(any(colnames(MCI_table) %in% pathway_table$Code)){level="P"}
  if(any(colnames(MCI_table) %in% pathway_table$Compound)){level="C"}

  if(level == "P"){
    MCI_table <- damma_compounds(MCI_table,pathway_table)
  }else if(level == "C"){
    MCI_table <- MCI_table
  }else{
    stop("MCI table and pathway table do not coincide")
  }

  functions_table <- t(as.data.frame(t(MCI_table)) %>%
    rownames_to_column('Compound') %>%
    left_join(pathway_table[,c('Compound', 'Function')], by = 'Compound') %>%
    group_by(Function) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    arrange(factor(Function, levels = unique(pathway_table$Function))) %>%
    column_to_rownames('Function'))

  return(functions_table)
}
