#' Aggregates pathway-level MCIs into compound-level MCIs
#'
#' @param MCI_table Function distillation table outputted by damma() or damma_correction
#' @param pathway_table Table containing definitions and metadata of metabolic functions (provided by GAMMA)
#' @return A MCI table aggregated at the compound level
#' @import tidyverse
#' @import dplyr
#' @import tibble
#' @examples
#' damma_compounds(compound_table,pathway_table)
#' @export

damma_compounds <- function(MCI_table,pathway_table){

   compounds_table <- t(as.data.frame(t(MCI_table)) %>%
     rownames_to_column('Code_pathway') %>%
     left_join(pathway_table[,c('Code_pathway', 'Code_compound')], by = 'Code_pathway') %>%
     group_by(Code_compound) %>%
     summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE))) %>%
     arrange(factor(Code_compound, levels = unique(pathway_table$Code_compound))) %>%
     column_to_rownames('Code_compound'))

   return(compounds_table)

}
