#' Aggregates compound values into indices at the function level
#'
#' @param compound_table Compound table outputted by aggregate_compounds()
#' @param functions Table containing definitions and metadata of metabolic functions (provided by GAMMA)
#' @param normalise Whether to normalise the scale to the maximum values in the dataset. Default=FALSE
#' @return A function table aggregated at the compound level
#' @examples
#' aggregate_functions(compound_table)
#' @export

aggregate_functions <- function(compounds_table,functions_table,normalise=FALSE){

  #Create empty table
  aggregated_table <- c()

  #Aggregate dietary carbohydrates
  compounds <- unique(functions_table[functions_table$Function == "Dietary carbohydrate degradation","Compound"])
  compounds_table_sub <- compounds_table[,compounds]
  compounds_table_sub <- rowSums(compounds_table_sub)
  if(normalise == FALSE){
    compounds_table_sub <- compounds_table_sub/length(compounds)
  }else{
    compounds_table_sub <- compounds_table_sub/max(compounds_table_sub)
  }
  aggregated_table <- cbind(aggregated_table,"Dietary carbohydrate degradation"=compounds_table_sub)

  #Aggregate dietary lipids
  compounds <- unique(functions_table[functions_table$Function == "Dietary lipid degradation","Compound"])
  compounds_table_sub <- compounds_table[,compounds]
  compounds_table_sub <- rowSums(compounds_table_sub)
  if(normalise == FALSE){
    compounds_table_sub <- compounds_table_sub/length(compounds)
  }else{
    compounds_table_sub <- compounds_table_sub/max(compounds_table_sub)
  }
  aggregated_table <- cbind(aggregated_table,"Dietary lipid degradation"=compounds_table_sub)

  #Aggregate dietary proteins
  compounds_table_sub <- t(t(compounds_table[,"Peptides"]))
  if(normalise == FALSE){
    compounds_table_sub <- compounds_table_sub/length(compounds)
  }else{
    compounds_table_sub <- compounds_table_sub/max(compounds_table_sub)
  }
  colnames(compounds_table_sub) <- "Dietary protein degradation"

  #Aggregate mucins
  compounds <- unique(functions_table[functions_table$Function == "Mucin degradation","Compound"])
  compounds_table_sub <- compounds_table[,compounds]
  compounds_table_sub <- rowSums(compounds_table_sub)
  if(normalise == FALSE){
    compounds_table_sub <- compounds_table_sub/length(compounds)
  }else{
    compounds_table_sub <- compounds_table_sub/max(compounds_table_sub)
  }
  aggregated_table <- cbind(aggregated_table,"Mucin degradation"=compounds_table_sub)

  #Aggregate SCFAs
  compounds <- unique(functions_table[functions_table$Function == "SCFA production","Compound"])
  compounds_table_sub <- compounds_table[,compounds]
  compounds_table_sub <- rowSums(compounds_table_sub)
  if(normalise == FALSE){
    compounds_table_sub <- compounds_table_sub/length(compounds)
  }else{
    compounds_table_sub <- compounds_table_sub/max(compounds_table_sub)
  }
  aggregated_table <- cbind(aggregated_table,"SCFA production"=compounds_table_sub)

  #Aggregate organic anions
  compounds <- unique(functions_table[functions_table$Function == "Organic anion production","Compound"])
  compounds_table_sub <- compounds_table[,compounds]
  compounds_table_sub <- rowSums(compounds_table_sub)
  if(normalise == FALSE){
    compounds_table_sub <- compounds_table_sub/length(compounds)
  }else{
    compounds_table_sub <- compounds_table_sub/max(compounds_table_sub)
  }
  aggregated_table <- cbind(aggregated_table,"Organic anion production"=compounds_table_sub)

  #Aggregate secondary bile acids
  compounds <- unique(functions_table[functions_table$Function == "Secondary bile acid production","Compound"])
  compounds_table_sub <- compounds_table[,compounds]
  compounds_table_sub <- rowSums(compounds_table_sub)
  if(normalise == FALSE){
    compounds_table_sub <- compounds_table_sub/length(compounds)
  }else{
    compounds_table_sub <- compounds_table_sub/max(compounds_table_sub)
  }
  aggregated_table <- cbind(aggregated_table,"Secondary bile acid production"=compounds_table_sub)

  #Aggregate amino acids
  compounds <- unique(functions_table[functions_table$Function == "Amino acid production","Compound"])
  compounds_table_sub <- compounds_table[,compounds]
  compounds_table_sub <- rowSums(compounds_table_sub)
  if(normalise == FALSE){
    compounds_table_sub <- compounds_table_sub/length(compounds)
  }else{
    compounds_table_sub <- compounds_table_sub/max(compounds_table_sub)
  }
  aggregated_table <- cbind(aggregated_table,"Amino acid production"=compounds_table_sub)

  #Aggregate amino acid derivatives
  compounds <- unique(functions_table[functions_table$Function == "Amino acid derivative production","Compound"])
  compounds_table_sub <- compounds_table[,compounds]
  compounds_table_sub <- rowSums(compounds_table_sub)
  if(normalise == FALSE){
    compounds_table_sub <- compounds_table_sub/length(compounds)
  }else{
    compounds_table_sub <- compounds_table_sub/max(compounds_table_sub)
  }
  aggregated_table <- cbind(aggregated_table,"Amino acid derivative production"=compounds_table_sub)

  #Aggregate vitamins
  compounds <- unique(functions_table[functions_table$Function == "Vitamin production","Compound"])
  compounds_table_sub <- compounds_table[,compounds]
  compounds_table_sub <- rowSums(compounds_table_sub)
  if(normalise == FALSE){
    compounds_table_sub <- compounds_table_sub/length(compounds)
  }else{
    compounds_table_sub <- compounds_table_sub/max(compounds_table_sub)
  }
  aggregated_table <- cbind(aggregated_table,"Vitamin production"=compounds_table_sub)
  
  return(aggregated_table)

}
