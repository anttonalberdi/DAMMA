#' Generates the scores of each hierarchical level of a metabolic pathway required to calculate community-level MCIs
#'
#' @param definition Definition string of a given metabolic pathway/module
#' @param def_table Decomposed hierarchy matrix produced by create_step_matrix.R
#' @param level Hierarchical level
#' @param abundance_vector Vector containing relative abundance data per identifier
#' @importFrom stringr str_sub
#' @import dplyr
#' @return A (partially) distilled definition string
#' @examples
#' distillate_definition(definition, def_table, level, present)
#' @export

#UNDER DEVELOPMENT
distillate_definition_community <- function(definition, def_table, level, abundance_vector){

  if(level == "L0_group"){
    def_table$clusters <- 0
    def_table_sub <- def_table[complete.cases(def_table[,level]),]
    clusters <- unique(def_table_sub$clusters)
  }else if(level == "L1_group"){
    def_table$clusters <- def_table[,"L0_group"]
    def_table_sub <- def_table[complete.cases(def_table[,level]),]
    clusters <- unique(def_table_sub$clusters)
  }else{
    def_table$clusters <- unite(def_table[,c(3:(ncol(def_table)-1))],col='clusters',colnames(def_table[,c(3:(ncol(def_table)-1))]), sep='-')
    def_table_sub <- def_table[complete.cases(def_table[,level]),]
    clusters <- unique(def_table_sub$clusters)[,1]
  }

  for (c in clusters){
    subdef <- def_table_sub[def_table_sub$clusters == c,"def_decomp"]
    if(" " %in% subdef | "+" %in% subdef){
      subdef2 <- subdef[(subdef != " ") & (subdef != "+")]
      subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)] <- abundance_vector[match(subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)], names(abundance_vector))]
      subdef2[is.na(subdef2)] <- 0
      value=round(mean(as.numeric(subdef2)),2)
    } else if("," %in% subdef){
      subdef2 <- subdef[subdef != ","]
      subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)] <- abundance_vector[match(subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)], names(abundance_vector))]
      subdef2[is.na(subdef2)] <- 0
      value=round(max(as.numeric(subdef2)),2)
    } else {
      subdef2 <- subdef
      subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)] <- abundance_vector[match(subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)], names(abundance_vector))]
      subdef2[is.na(subdef2)] <- 0
      value=round(max(as.numeric(subdef2)),2)
    }
    if(level == "L0_group"){
      definition <- gsub(paste(subdef,collapse=""),value,definition, fixed = TRUE)
    }else{
      definition <- gsub(paste(c("(",subdef,")"),collapse=""),value,definition, fixed = TRUE)
    }
  }
    return(definition)
}
