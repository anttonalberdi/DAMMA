#' Calculate correlations between compounds and core functions
#'
#' @param fullness_table Table containing fullness scores of functions (X axis) per MAG (Y axis)
#' @param method Correlation method to be used. Default="pearson"
#' @param output Whether to output a table or a plot. Default="plot"
#' @return A correlation matrix or plot
#' @import ggplot2
#' @import ggcorrplot
#' @examples
#' damma_cor(fullness_table)
#' damma_cor(fullness_table,output="table")
#' @export

damma_correlation <- function(fullness_table,method="spearman", output="plot"){

  #### UNDER DEVELOPMENT ####
  cor_matrix <-cor(data.frame(fullness_table),method = method)
  cor_matrix[is.na(cor_matrix)] <- 0
  ggcorrplot(cor_matrix,
             outline.color = "black",
             show.diag  = F,
             hc.order = FALSE,
             type = "upper",
             lab = F,
             digits = 1,
             insig = "blank",
             ggtheme = theme_bw())
}
