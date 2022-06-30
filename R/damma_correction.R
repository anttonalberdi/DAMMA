#' Models the relationship between function fullness and MAG completeness and applies a correction factor
#'
#' @param fullness_table Table containing fullness scores of functions (X axis) per MAG (Y axis)
#' @param mag_completeness A matrix containing MAG names (1st column) and completemess values ()
#' @return A modified fullness table
#' @examples
#' damma_correction(fullness_table, mag_completeness)
#' @export

damma_correction <- function(fullness_table,mag_completeness){

  #### UNDER DEVELOPMENT ####
  MAG_completeness <- as.numeric(mag_completeness[,2])

  #Create corrected fullness matrix
  fullness_table_corrected <- matrix(0,nrow = nrow(fullness_table),ncol = ncol(fullness_table))
  colnames(fullness_table_corrected)=colnames(fullness_table)
  rownames(fullness_table_corrected)=rownames(fullness_table)

  #Iterate modelling and correction for each function
  for(i in 1:ncol(fullness_table)){
    # Arbitrary number of 10 steps is used as weights
    Model <- glm(fullness_table[,i]~MAG_completeness,family = "binomial",weights = rep(10,length(MAG_completeness)))
    slope_coef <- coef(Model)[2]
    if(slope_coef > 0){
      for(j in 1:nrow(fullness_table)){
        # Model prediction of fullness if completeness was 100%
        pred_100 <- round(predict(Model,newdata = data.frame(MAG_completeness=100,weights=10),type = "response"),1)
        # Model prediction of fullness for actual completeness of the focal MAG
        pred_focal <- round(predict(Model,newdata = data.frame(MAG_completeness=MAG_completeness[j],weights=10),type = "response"),1)
        # The expected change in function fullness if focal MAG was 100% complete
        pred_diff <- pred_100-pred_focal
        fullness_table_corrected[j,i] = fullness_table[j,i]+pred_diff
      }
    }else if(slope_coef <= 0){
      fullness_table_corrected[,i] <- fullness_table[,i]
    }
  }

  # If corrected fullness >1, convert it to 1.
  fullness_table_corrected[fullness_table_corrected>1] <- 1

  return(fullness_table_corrected)

}
