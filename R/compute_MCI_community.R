#' Calculates the MCI of a metabolic pathway at the community level
#'
#' @param definition Definition string of a given metabolic pathway/module
#' @param abundance_vector Vector containing relative abundance data per identifier
#' @return A MCI value
#' @examples
#' compute_MCI_community(definition,abundance_vector)
#' @export

compute_MCI_community <- function(definition,abundance_vector){
  #If using EC codes
  if (grepl(".", definition, fixed = TRUE)){
    names(abundance_vector) <- gsub(".","_",names(abundance_vector),fixed=TRUE)
    definition <- gsub(".","_",definition,fixed=TRUE)
  }
  #Decompose definition
  def_decomp <- unlist(strsplit(definition, "(?=[ ( ),+]+)", perl=TRUE))
  #Set levels
  def_level <- set_levels(def_decomp)
  #Definition-level table
  def_table <- create_step_matrix(def_decomp,def_level)
  #List levels
  levels <- colnames(def_table[,c(3:ncol(def_table))])
  #Iterate calculation across levels
  for(level in rev(levels)){
    definition <- distillate_definition_community(definition, def_table, level, abundance_vector)
    if(level != "L0_group"){
      def_decomp <- unlist(strsplit(definition, "(?=[ ( ),+]+)", perl=TRUE))
      def_level <- set_levels(def_decomp)
      def_table <- create_step_matrix(def_decomp,def_level)
    }
  }
  #Return value
  return(as.numeric(definition))
}
