#' Aggregates raw fullness values into meaningful indices at the compound level
#'
#' @param distilled_table Function distillation table outputted by damma() or damma_correction
#' @param functions_table Table containing definitions and metadata of metabolic functions (provided by GAMMA)
#' @return A function table aggregated at the compound level
#' @import tidyverse
#' @import dplyr
#' @import tibble
#' @examples
#' damma_compounds(compound_table)
#' @export

damma_compounds <- function(distilled_table,functions_table){

   compounds_table <- t(as.data.frame(t(distilled_table)) %>%
     rownames_to_column('Code') %>%
     left_join(functions_table[,c('Code', 'Compound')], by = 'Code') %>%
     group_by(Compound) %>%
     summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE))) %>%
     arrange(factor(Compound, levels = unique(functions_table$Compound))) %>%
     column_to_rownames('Compound'))

   return(compounds_table)

}
