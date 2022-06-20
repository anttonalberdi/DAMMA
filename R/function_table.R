#' Reference table containing function descriptions and metadata
#'
#' Information compiled from databases KEGG, MetaSys and MEROPS.
#'
#' Data created from version 20220619
#' function_table <- read.table("data/DAMMA_functions_20220619.tsv", header=TRUE, sep="\t")
#' save(function_table, file="data/function_table.RData")
#'
#' @docType data
#' @usage data(function_table)
#' @format A dataframe.
#' @keywords datasets
#' @examples
#' data(function_table)
#' @export
"function_table"
