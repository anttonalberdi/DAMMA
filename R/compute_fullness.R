#' Calculates the overall fullness of a metabolic pathway/module
#'
#' @param definition Definition string of a given metabolic pathway/module
#' @param present Vector of functional units present in the genome
#' @return A fullness value
#' @examples
#' compute_fullness(definition,present)
#' compute_fullness("K01580 (K13524,K07250,K00823,K16871) (K00135,K00139,K17761)",c("K01580","K00823","K16871"))
#' @export

compute_fullness <- function(definition,present){
  #If using EC codes
  if (grepl(".", definition, fixed = TRUE)){
    present <- gsub(".","_",present,fixed=TRUE)
    definition <- gsub(".","_",definition,fixed=TRUE)
  }
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
    definition <- distillate_definition(definition, def_table, level, present)
    if(level != "L0_group"){
      def_decomp <- decompose_definition(definition)
      def_level <- set_levels(def_decomp)
      def_table <- create_step_matrix(def_decomp,def_level)
    }
    }
  #Return value
  return(as.numeric(definition))
}
