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


damma_functions <- function(compounds_table,pathway_table){

  functions_table <- t(as.data.frame(t(compounds_table)) %>%
    rownames_to_column('Code_compound') %>%
    left_join(pathway_table[,c('Code_compound', 'Code_function')], by = 'Code_compound') %>%
    group_by(Code_function) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    arrange(factor(Code_function, levels = unique(pathway_table$Code_function))) %>%
    column_to_rownames('Code_function'))

  return(functions_table)
}
