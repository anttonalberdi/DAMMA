#' Aggregates compound values into indices at the function level
#'
#' @param compound_table Compound table outputted by aggregate_compounds()
#' @param functions Table containing definitions and metadata of metabolic functions (provided by GAMMA)
#' @param normalise Whether to normalise the scale to the maximum values in the dataset. Default=FALSE
#' @return A function table aggregated at the compound level
#' @examples
#' dama_functions(compound_table)
#' @export

damma_functions <- function(compounds_table,functions_table,normalise=FALSE){

  functions_table <- as.data.frame(functions_table)

  #Create empty table
  aggregated_table <- c()

  #Aggregate polysaccharides
  compounds <- unique(functions_table[functions_table$Function == "Polysaccharide degradation","Compound"])
  compounds_table_sub <- compounds_table[,compounds]
  compounds_table_sub <- rowSums(compounds_table_sub)
  if(normalise == FALSE){
    compounds_table_sub <- compounds_table_sub/length(compounds)
  }else{
    compounds_table_sub <- compounds_table_sub/max(compounds_table_sub)
  }
  aggregated_table <- cbind(aggregated_table,"Polysaccharide degradation"=compounds_table_sub)

  #Aggregate sugars
  compounds <- unique(functions_table[functions_table$Function == "Sugar degradation","Compound"])
  compounds_table_sub <- compounds_table[,compounds]
  compounds_table_sub <- rowSums(compounds_table_sub)
  if(normalise == FALSE){
    compounds_table_sub <- compounds_table_sub/length(compounds)
  }else{
    compounds_table_sub <- compounds_table_sub/max(compounds_table_sub)
  }
  aggregated_table <- cbind(aggregated_table,"Sugar degradation"=compounds_table_sub)

  #Aggregate lipids
  compounds <- unique(functions_table[functions_table$Function == "Lipid degradation","Compound"])
  compounds_table_sub <- compounds_table[,compounds]
  compounds_table_sub <- rowSums(compounds_table_sub)
  if(normalise == FALSE){
    compounds_table_sub <- compounds_table_sub/length(compounds)
  }else{
    compounds_table_sub <- compounds_table_sub/max(compounds_table_sub)
  }
  aggregated_table <- cbind(aggregated_table,"Lipid degradation"=compounds_table_sub)

  #Aggregate proteins
  compounds_table_sub <- t(t(compounds_table[,"Peptides"]))
  if(normalise == FALSE){
    compounds_table_sub <- compounds_table_sub
  }else{
    compounds_table_sub <- compounds_table_sub/max(compounds_table_sub)
  }
  colnames(compounds_table_sub) <- "Protein degradation"

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
