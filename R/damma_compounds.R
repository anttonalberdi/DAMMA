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
     rownames_to_column('Code') %>%
     left_join(pathway_table[,c('Code', 'Compound')], by = 'Code') %>%
     group_by(Compound) %>%
     summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE))) %>%
     arrange(factor(Compound, levels = unique(pathway_table$Compound))) %>%
     column_to_rownames('Compound'))

   return(compounds_table)

}
