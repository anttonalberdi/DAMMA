#' Aggregates compound-level MCIs into function-level MCIs
#'
#' @param compounds_table Compound table outputted by aggregate_compounds()
#' @param pathway_table Table containing definitions and metadata of metabolic functions (provided by GAMMA)
#' @param transform Whether to transform to 0-1 scale. Default=TRUE
#' @import tidyverse
#' @import dplyr
#' @return A MCI table aggregated at the compound level
#' @examples
#' dama_functions(compound_table)
#' @export


damma_functions <- function(compounds_table,pathway_table,transform=TRUE){

  if(transform == TRUE){
  functions_table <- t(as.data.frame(t(compounds_table)) %>%
    rownames_to_column('Compound') %>%
    left_join(pathway_table[,c('Compound', 'Function')], by = 'Compound') %>%
    group_by(Function) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    arrange(factor(Function, levels = unique(pathway_table$Function))) %>%
    column_to_rownames('Function'))
  }

  if(transform == FALSE){
  functions_table <- t(as.data.frame(t(compounds_table)) %>%
    rownames_to_column('Compound') %>%
    left_join(pathway_table[,c('Compound', 'Function')], by = 'Compound') %>%
    group_by(Function) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>%
    arrange(factor(Function, levels = unique(pathway_table$Function))) %>%
    column_to_rownames('Function'))
  }

  return(functions_table)
}
