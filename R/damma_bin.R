#' Convert the fullness table of metabolic pathways/modules into a binary (present/absent) matrix
#'
#' @param fullness_table Table containing fullness scores of functions (X axis) per MAG (Y axis)
#' @param threshold Minimum fullness value to consider a function (pathway/module) to be present in a MAG
#' @return A binary function matrix
#' @examples
#' damma_binary(fullness_table)
#' damma_binary(fullness_table,0.9)
#' damma_binary(fullness_table,threshold=0.9)
#' @export

damma_binary <- function(fullness_table,threshold=0.9){

    fullness_table_binary <- fullness_table
    fullness_table_binary[fullness_table_binary >= threshold] <- 1
    fullness_table_binary[fullness_table_binary < threshold] <- 0

  return(fullness_table_binary)
}
