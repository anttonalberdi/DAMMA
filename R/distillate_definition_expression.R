#' Generates the scores of each hierarchical level of a metabolic pathway required to calculate gene expression/abundance-based MCIs
#'
#' @param sample Code of the sample to be analysed
#' @param definition_expression Definition-expression string
#' @param def_table Decomposed hierarchy matrix produced by create_step_matrix.R
#' @param level Hierarchical level of the pathway definition
#' @param expression_table Table of expression values of functional units present in the genome
#' @importFrom stringr str_sub
#' @return A (partially) distilled definition string
#' @examples
#' distillate_definition_expression(definition_expression, def_table, level, expression_table)
#' @export

#UNDER DEVELOPMENT
distillate_definition_expression <- function(sample, definition_expression, def_table, level, expression_table){

  if(level == "L0_group"){
    def_table$clusters <- 0
    def_table_sub <- def_table[complete.cases(def_table[,level]),]
    clusters <- unique(def_table_sub$clusters)
  }else if(level == "L1_group"){
    def_table$clusters <- def_table[,"L0_group"]
    def_table_sub <- def_table[complete.cases(def_table[,level]),]
    clusters <- unique(def_table_sub$clusters)
  }else{
    def_table$clusters <- do.call(paste, c(def_table[,c(3:(ncol(def_table)-1))], sep="-"))
    def_table_sub <- def_table[complete.cases(def_table[,level]),]
    clusters <- unique(def_table_sub$clusters)
  }


    for (c in clusters){
      subdef <- def_table_sub[def_table_sub$clusters == c,"def_decomp"]
      if(" " %in% subdef | "+" %in% subdef){
        subdef2 <- subdef[(subdef != " ") & (subdef != "+")]
        subdef2_code <- subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)]
        subdef2_number <- as.numeric(subdef2[!subdef2 %in% subdef2_code])
        if(length(subdef2_code) > 0){
          subdef2_expression <- expression_table[rownames(expression_table) %in% subdef2_code,sample]
          if(length(subdef2_expression) == 0){
            subdef2_expression <- 0
          }
        }else{
          subdef2_expression <- NA
        }
        if(length(subdef2_expression) > 0 | is.na(subdef2_expression)){
          value <- mean(c(subdef2_number,subdef2_expression),na.rm=TRUE)
        }else{
          value <- 0
        }
      }else if("," %in% subdef){
        subdef2 <- subdef[subdef != ","]
        subdef2_code <- subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)]
        subdef2_number <- as.numeric(subdef2[!subdef2 %in% subdef2_code])
        if(length(subdef2_code) > 0){
          subdef2_expression <- expression_table[rownames(expression_table) %in% subdef2_code,sample]
          if(length(subdef2_expression) == 0){
            subdef2_expression <- 0
          }
        }else{
          subdef2_expression <- NA
        }
        if(length(subdef2_expression) > 0 | is.na(subdef2_expression)){
          value <- max(c(subdef2_number,subdef2_expression),na.rm=TRUE)
        }else{
          value <- 0
        }
      } else {
        subdef2 <- subdef
        subdef2_code <- subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)]
        subdef2_number <- as.numeric(subdef2[!subdef2 %in% subdef2_code])
        if(length(subdef2_code) > 0){
          subdef2_expression <- expression_table[rownames(expression_table) %in% subdef2_code,sample]
          if(length(subdef2_expression) == 0){
            subdef2_expression <- 0
          }
        }else{
          subdef2_expression <- NA
        }
        if(length(subdef2_expression) > 0 | is.na(subdef2_expression)){
          value <- max(c(subdef2_number,subdef2_expression),na.rm=TRUE)
        }else{
          value <- 0
        }
      }
      if(level == "L0_group"){
        definition_expression <- gsub(paste(subdef,collapse=""),value,definition_expression, fixed = TRUE)
      }else{
        definition_expression <- gsub(paste(c("(",subdef,")"),collapse=""),value,definition_expression, fixed = TRUE)
      }
    }

    return(definition_expression)
}
