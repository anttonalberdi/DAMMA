#' Models the relationship between function fullness and MAG completeness and applies a correction factor
#'
#' @param MCI_table Table containing MCI scores outputed by damma()
#' @param genome_completeness A matrix containing MAG names (1st column) and completemess values (2nd column)
#' @param stats Whether to print correction statistics on screen or not. Default=TRUE
#' @return A corrected fullness table
#' @import tidyverse
#' @examples
#' damma_correction(MCI_table, genome_completeness)
#' @export

damma_correction <- function(MCI_table,genome_completeness,stats=TRUE){

  #Sort Genomes and test matching
  MCI_table <- MCI_table[order(rownames(MCI_table)),]
  genome_completeness <- genome_completeness[order(genome_completeness[,1]),]
  if(all(as.character(rownames(MCI_table)) != as.character(genome_completeness[,1]))) stop("Pathway table is missing")

  #Get completeness values
  Genome_completeness <- as.numeric(genome_completeness[,2])

  #Create corrected fullness matrix
  MCI_table_corrected <- matrix(0,nrow = nrow(MCI_table),ncol = ncol(MCI_table))
  colnames(MCI_table_corrected)=colnames(MCI_table)
  rownames(MCI_table_corrected)=rownames(MCI_table)

  #Iterate modelling and correction for each function
  suppressWarnings(
    for(i in 1:ncol(MCI_table)){
      Model <- glm(MCI_table[,i]~Genome_completeness,family = "binomial")
      slope_coef <- coef(Model)[2]
      if(slope_coef > 0){
        for(j in 1:nrow(MCI_table)){
          # Model prediction of fullness if completeness was 100%
          pred_100 <- round(predict(Model,newdata = data.frame(Genome_completeness=100),type = "response"),1)
          # Model prediction of fullness for actual completeness of the focal MAG
          pred_focal <- round(predict(Model,newdata = data.frame(Genome_completeness=Genome_completeness[j]),type = "response"),1)
          # The expected change in function fullness if focal MAG was 100% complete
          pred_diff <- pred_100-pred_focal
          MCI_table_corrected[j,i] = MCI_table[j,i]+pred_diff
        }
      }else if(slope_coef <= 0){
        MCI_table_corrected[,i] <- MCI_table[,i]
      }
    }
  )

  # If corrected fullness >1, convert it to 1.
  MCI_table_corrected[MCI_table_corrected>1] <- 1

  #Outout overall correction statistics on screen
  total <- nrow(MCI_table)*ncol(MCI_table)
  changes <- c(MCI_table == MCI_table_corrected)
  changes <- length(changes[!(changes)])
  percentage <- round(changes / total * 100,1)
  cat(paste0(changes," out of ",total," (",percentage,"%) fullness values were corrected\n"))

  if(stats == TRUE){
    #Outout overall correction statistics on screen
    total <- ncol(MCI_table_corrected)
    for(r in rownames(MCI_table_corrected)){
        completeness <- round(as.numeric(genome_completeness[genome_completeness[,1] == r,2]),1)
        changes <- MCI_table[r,] == MCI_table_corrected[r,]
        changes <- length(changes[!(changes)])
        percentage <- round(changes / total * 100,1)
        cat(paste0("\t",r," (",completeness,"%): ",changes,"/",total," (",percentage,"%) fullness values were corrected\n"))
    }
  }

  #Output corrected table
  return(MCI_table_corrected)

}
