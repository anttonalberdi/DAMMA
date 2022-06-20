#' Convert the fullness table of metabolic pathways/modules into a binary (present/absent) matrix
#'
#' @param fullness_table Table containing MAG identifiers and annotation codes
#' @param threshold Table containing definitions and metadata of metabolic functions (provided by GAMMA)
#' @param completeness Column index (number) of the annotations table containing the MAG identifiers
#' @return A binary function matrix
#' @examples
#' damma_bin(fullness_table)
#' damma_bin(fullness_table,0.9)
#' damma_bin(fullness_table,threshold=0.9,completeness)
#' @export

damma_bin <- function(fullness_table,threshold=0.9,completeness){
  if(missing(completeness)){
    if(nrow(fullness_table) != nrow(completeness)) stop("Dimensions of function fullness and MAG completeness data do not match")
  }

  if(missing(completeness)){
    #IGNORING MAG COMPLETENESS
    #Convert table into binary
    fullness_table_binary <- fullness_table
    fullness_table_binary[fullness_table_binary >= threshold] <- 1
    fullness_table_binary[fullness_table_binary < threshold] <- 0

  }else{
    #ACCOUNTING FOR MAG COMPLETENESS
    #Adjust thresholds to each MAG
    threshold_corrected <- completeness[,1]/100 * threshold

    #Convert table into binary
    fullness_table_binary <- fullness_table
    for(r in c(1:nrow(fullness_table_binary))){
      row <- fullness_table_binary[r, ]
      row[row >= threshold_corrected[r]] <- 1
      row[row < threshold_corrected[r]] <- 0
      fullness_table_binary[r, ] <- row
    }
  }

  return(fullness_table_binary)
}
