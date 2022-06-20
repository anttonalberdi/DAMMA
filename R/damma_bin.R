#' Convert the fullness table of metabolic pathways/modules into a binary (present/absent) matrix
#'
#' @param fullness_table Table containing fullness scores of functions (X axis) per MAG (Y axis)
#' @param threshold Minimum fullness value to consider a function (pathway/module) to be present in a MAG
#' @param completeness Two-column table indicating the completeness (eg. CheckM) percentage (column 2) of each MAG (column 1)
#' @return A binary function matrix
#' @examples
#' damma_bin(fullness_table)
#' damma_bin(fullness_table,0.9)
#' damma_bin(fullness_table,threshold=0.9,completeness)
#' @export

damma_bin <- function(fullness_table,threshold=0.9,completeness){
  if(missing(completeness)){

  }

  if(missing(completeness)){
    #IGNORING MAG COMPLETENESS
    #Convert table into binary
    fullness_table_binary <- fullness_table
    fullness_table_binary[fullness_table_binary >= threshold] <- 1
    fullness_table_binary[fullness_table_binary < threshold] <- 0

  }else{
    #ACCOUNTING FOR MAG COMPLETENESS

    if(nrow(fullness_table) != nrow(completeness)) stop("Dimensions of function fullness and MAG completeness data do not match")

    #Adjust thresholds to each MAG
    threshold_corrected <- completeness[,2]/100 * threshold

    #Convert table into binary
    fullness_table_binary <- fullness_table[completeness[,1],]
    for(r in c(1:nrow(fullness_table_binary))){
      row <- fullness_table_binary[r, ]
      row[row >= threshold_corrected[r]] <- 1
      row[row < threshold_corrected[r]] <- 0
      fullness_table_binary[r, ] <- row
    }
  }

  return(fullness_table_binary)
}
