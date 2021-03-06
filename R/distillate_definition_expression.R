#' Calculates the expression-fullness of each hierarchical level of a metabolic pathway/module
#'
#' @param definition_expression Definition-expression string
#' @param def_table Decomposed hierarchy matrix produced by create_step_matrix.R
#' @param level Hierarchical level
#' @param expression_vector Vector of expression values of functional units present in the genome
#' @importFrom stringr str_sub
#' @return A (partially) distilled definition string
#' @examples
#' distillate_definition_expression(definition_expression, def_table, level, present)
#' @export

#UNDER DEVELOPMENT
distillate_definition_expression <- function(sample,definition_expression, def_table, level, expression_table){

  if (level == "L5_group"){
    #L5
    def_table_sub <- def_table[,c("L0_group","L1_group","L2_group","L3_group","L4_group","L5_group")][complete.cases(def_table[,c("L0_group","L1_group","L2_group","L3_group","L4_group","L5_group")]),]
    L5_clusters <- unique(paste(def_table_sub$L0_group,def_table_sub$L1_group,def_table_sub$L2_group,def_table_sub$L3_group,def_table_sub$L4_group,sep="-"))
    for (c in L5_clusters){
      subdef <- def_table[(def_table$L0_group == str_sub(c,1,1)) & (def_table$L1_group == str_sub(c,3,3)) & (def_table$L2_group == str_sub(c,5,5)) & (def_table$L3_group == str_sub(c,7,7)) & (def_table$L4_group == str_sub(c,9,9)) & (!is.na(def_table$L5_group)),"def_decomp"]
      if(" " %in% subdef | "+" %in% subdef){
        subdef2 <- subdef[(subdef != " ") & (subdef != "+")]
        subdef2_code <- subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)]
        subdef2_number <- as.numeric(subdef2[!subdef2 %in% subdef2_code])
        subdef2_expression <- expression_table[rownames(expression_table) %in% subdef2_code,sample]
        if(length(subdef2_expression) > 0){
          value <- mean(subdef2_expression)
        }else{
          value <- 0
        }
      } else if("," %in% subdef){
        subdef2 <- subdef[subdef != ","]
        subdef2_code <- subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)]
        subdef2_number <- as.numeric(subdef2[!subdef2 %in% subdef2_code])
        subdef2_expression <- expression_table[rownames(expression_table) %in% subdef2_code,sample]
        if(length(subdef2_expression) > 0){
          value <- max(subdef2_expression)
        }else{
          value <- 0
        }
      } else {
        value <- 0
      }
      #Replace definition by value
      definition_expression <- gsub(paste(c("(",subdef,")"),collapse=""),value,definition_expression, fixed = TRUE)
    }
  }

  if (level == "L4_group"){
    #L4
    def_table_sub <- def_table[,c("L0_group","L1_group","L2_group","L3_group","L4_group")][complete.cases(def_table[,c("L0_group","L1_group","L2_group","L3_group","L4_group")]),]
    L4_clusters <- unique(paste(def_table_sub$L0_group,def_table_sub$L1_group,def_table_sub$L2_group,def_table_sub$L3_group,sep="-"))
    for (c in L4_clusters){
      subdef <- def_table[(def_table$L0_group == str_sub(c,1,1)) & (def_table$L1_group == str_sub(c,3,3)) & (def_table$L2_group == str_sub(c,5,5)) & (def_table$L3_group == str_sub(c,7,7)) & (!is.na(def_table$L4_group)),"def_decomp"]
      if(" " %in% subdef | "+" %in% subdef){
        subdef2 <- subdef[(subdef != " ") & (subdef != "+")]
        subdef2_code <- subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)]
        subdef2_number <- as.numeric(subdef2[!subdef2 %in% subdef2_code])
        subdef2_expression <- expression_table[rownames(expression_table) %in% subdef2_code,sample]
        if(length(subdef2_expression) > 0){
          value <- mean(subdef2_expression)
        }else{
          value <- 0
        }
      } else if("," %in% subdef){
        subdef2 <- subdef[subdef != ","]
        subdef2_code <- subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)]
        subdef2_number <- as.numeric(subdef2[!subdef2 %in% subdef2_code])
        subdef2_expression <- expression_table[rownames(expression_table) %in% subdef2_code,sample]
        if(length(subdef2_expression) > 0){
          value <- max(subdef2_expression)
        }else{
          value <- 0
        }
      } else {
        value <- 0
      }
      #Replace definition by value
      definition_expression <- gsub(paste(c("(",subdef,")"),collapse=""),value,definition_expression, fixed = TRUE)
    }
  }


    if (level == "L3_group"){
      #L3
      def_table_sub <- def_table[,c("L0_group","L1_group","L2_group","L3_group")][complete.cases(def_table[,c("L0_group","L1_group","L2_group","L3_group")]),]
      L3_clusters <- unique(paste(def_table_sub$L0_group,def_table_sub$L1_group,def_table_sub$L2_group,sep="-"))
      for (c in L3_clusters){
        subdef <- def_table[(def_table$L0_group == str_sub(c,1,1)) & (def_table$L1_group == str_sub(c,3,3)) & (def_table$L2_group == str_sub(c,5,5)) & (!is.na(def_table$L3_group)),"def_decomp"]
        if(" " %in% subdef | "+" %in% subdef){
          subdef2 <- subdef[(subdef != " ") & (subdef != "+")]
          subdef2_code <- subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)]
          subdef2_number <- as.numeric(subdef2[!subdef2 %in% subdef2_code])
          subdef2_expression <- expression_table[rownames(expression_table) %in% subdef2_code,sample]
          if(length(subdef2_expression) > 0){
            value <- mean(subdef2_expression)
          }else{
            value <- 0
          }
        } else if("," %in% subdef){
          subdef2 <- subdef[subdef != ","]
          subdef2_code <- subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)]
          subdef2_number <- as.numeric(subdef2[!subdef2 %in% subdef2_code])
          subdef2_expression <- expression_table[rownames(expression_table) %in% subdef2_code,sample]
          if(length(subdef2_expression) > 0){
            value <- max(subdef2_expression)
          }else{
            value <- 0
          }
        } else {
          value <- 0
        }
        #Replace definition by value
        definition_expression <- gsub(paste(c("(",subdef,")"),collapse=""),value,definition_expression, fixed = TRUE)
      }
    }

    if (level == "L2_group"){
      #L2
      def_table_sub <- def_table[,c("L0_group","L1_group","L2_group")][complete.cases(def_table[,c("L0_group","L1_group","L2_group")]),]
      L2_clusters <- unique(paste(def_table_sub$L0_group,def_table_sub$L1_group,sep="-"))
      for (c in L2_clusters){
        subdef <- def_table[(def_table$L0_group == str_sub(c,1,1)) &(def_table$L1_group == str_sub(c,3,3)) & (!is.na(def_table$L2_group)),"def_decomp"]
        if(" " %in% subdef | "+" %in% subdef){
          subdef2 <- subdef[(subdef != " ") & (subdef != "+")]
          subdef2_code <- subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)]
          subdef2_number <- as.numeric(subdef2[!subdef2 %in% subdef2_code])
          subdef2_expression <- expression_table[rownames(expression_table) %in% subdef2_code,sample]
          if(length(subdef2_expression) > 0){
            value <- mean(subdef2_expression)
          }else{
            value <- 0
          }
        } else if("," %in% subdef){
          subdef2 <- subdef[subdef != ","]
          subdef2_code <- subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)]
          subdef2_number <- as.numeric(subdef2[!subdef2 %in% subdef2_code])
          subdef2_expression <- expression_table[rownames(expression_table) %in% subdef2_code,sample]
          if(length(subdef2_expression) > 0){
            value <- max(subdef2_expression)
          }else{
            value <- 0
          }
        } else {
          value <- 0
        }
        #Replace definition by value
        definition_expression <- gsub(paste(c("(",subdef,")"),collapse=""),value,definition_expression, fixed = TRUE)
      }
    }

    if (level == "L1_group"){
      #L1
      def_table_sub <- def_table[,c("L0_group","L1_group")][complete.cases(def_table[,c("L0_group","L1_group")]),]
      L1_clusters <- unique(def_table_sub$L0_group)
      for (c in L1_clusters){
        subdef <- def_table[(def_table$L0_group == c) & (!is.na(def_table$L1_group)),"def_decomp"]
        if(" " %in% subdef | "+" %in% subdef){
          subdef2 <- subdef[(subdef != " ") & (subdef != "+")]
          subdef2_code <- subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)]
          subdef2_number <- as.numeric(subdef2[!subdef2 %in% subdef2_code])
          subdef2_expression <- expression_table[rownames(expression_table) %in% subdef2_code,sample]
          if(length(subdef2_expression) > 0){
            value <- mean(subdef2_expression)
          }else{
            value <- 0
          }
        } else if("," %in% subdef){
          subdef2 <- subdef[subdef != ","]
          subdef2_code <- subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)]
          subdef2_number <- as.numeric(subdef2[!subdef2 %in% subdef2_code])
          subdef2_expression <- expression_table[rownames(expression_table) %in% subdef2_code,sample]
          if(length(subdef2_expression) > 0){
            value <- max(subdef2_expression)
          }else{
            value <- 0
          }
        } else {
          value <- 0
        }
        #Replace definition by value
        definition_expression <- gsub(paste(c("(",subdef,")"),collapse=""),value,definition_expression, fixed = TRUE)

      }
    }

    if (level == "L0_group"){
      #L0
      subdef <- def_table[!is.na(def_table$L0_group),"def_decomp"]
      if(" " %in% subdef | "+" %in% subdef){
        subdef2 <- subdef[(subdef != " ") & (subdef != "+")]
        subdef2_code <- subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)]
        subdef2_number <- as.numeric(subdef2[!subdef2 %in% subdef2_code])
        subdef2_expression <- expression_table[rownames(expression_table) %in% subdef2_code,sample]
        if(length(subdef2_expression) > 0){
          value <- mean(c(subdef2_number,subdef2_expression))
        }else{
          value <- 0
        }
      } else if("," %in% subdef){
        subdef2 <- subdef[subdef != ","]
        subdef2_code <- subdef2[grepl("_", subdef2, fixed = TRUE) | grepl("[A-Z]", subdef2, fixed = FALSE)]
        subdef2_number <- as.numeric(subdef2[!subdef2 %in% subdef2_code])
        subdef2_expression <- expression_table[rownames(expression_table) %in% subdef2_code,sample]
        if(length(subdef2_expression) > 0){
          value <- max(c(subdef2_number,subdef2_expression))
        }else{
          value <- 0
        }
      } else {
        value <- 0
      }
      definition_expression <- gsub(paste(subdef,collapse=""),value,definition_expression, fixed = TRUE)
    }

    return(definition_expression)
}
