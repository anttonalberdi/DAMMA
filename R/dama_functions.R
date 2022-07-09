#' Aggregates compound values into indices at the function level
#'
#' @param compound_table Compound table outputted by aggregate_compounds()
#' @param functions Table containing definitions and metadata of metabolic functions (provided by GAMMA)
#' @param transform Whether to transform to 0-1 scale. Default=TRUE
#' @return A function table aggregated at the compound level
#' @examples
#' dama_functions(compound_table)
#' @export


damma_functions <- function(compounds_table,functions_table,transform=TRUE){

  if(transform == TRUE){
  functions_matrix <- t(as.data.frame(t(compounds_table)) %>%
    rownames_to_column('Compound') %>%
    left_join(functions_table[,c('Compound', 'Function')], by = 'Compound') %>%
    group_by(Function) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    arrange(factor(Function, levels = unique(functions_table$Function))) %>%
    column_to_rownames('Function'))
  }

  if(transform == FALSE){
  functions_matrix <- t(as.data.frame(t(compounds_table)) %>%
    rownames_to_column('Compound') %>%
    left_join(functions_table[,c('Compound', 'Function')], by = 'Compound') %>%
    group_by(Function) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>%
    arrange(factor(Function, levels = unique(functions_table$Function))) %>%
    column_to_rownames('Function'))
  }

  return(functions_matrix)
}
