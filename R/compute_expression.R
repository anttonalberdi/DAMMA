#' Calculates the overall fullness-expression of a metabolic pathway/module
#'
#' @param definition Definition string of a given metabolic pathway/module
#' @param expression_table Table containing expression data per functon
#' @return A fullness-expression value
#' @examples
#' compute_expression(definition,expression_table)
#' @export

#UNDER DEVELOPMENT
compute_expression <- function(definition,expression_table){

  #TO BE UPDATED
  #Expression * fullness of each hierarchical value

  #If using EC codes
  if (grepl(".", definition, fixed = TRUE)){
    present <- gsub(".","_",present,fixed=TRUE)
    definition <- gsub(".","_",definition,fixed=TRUE)
  }
  #Lis present annotations units
  present <- nrow(expression_table)
  #Decompose definition
  def_decomp <- decompose_definition(definition)
  #Set levels
  def_level <- set_levels(def_decomp)
  #Definition-level table
  def_table <- create_step_matrix(def_decomp,def_level)
  #Calculate number of levels
  levels <- names(colSums(def_table[,c(3:8)],na.rm=TRUE)[colSums(def_table[,c(3:8)],na.rm=TRUE)>0])
  #Iterate calculation across levels
  for(level in rev(levels)){
    definition <- distillate_definition_expression(definition, def_table, level, expression_vector)
    if(level != "L0_group"){
      def_decomp <- decompose_definition(definition)
      def_level <- set_levels(def_decomp)
      def_table <- create_step_matrix(def_decomp,def_level)
    }
    }
  #Return value
  return(as.numeric(definition))
}
