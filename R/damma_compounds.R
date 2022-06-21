#' Aggregates raw fullness values into meaningful indices at the compound level
#'
#' @param distilled_table Function distillation table outputted by damma()
#' @return A function table aggregated at the compound level
#' @param functions Table containing definitions and metadata of metabolic functions (provided by GAMMA)
#' @examples
#' damma_compounds(compound_table)
#' @export

damma_compounds <- function(distilled_table,functions){

  #Create empty table
  aggregated_table <- c()

  #Aggregate dietary carbohydrates
  codes <- functions[functions$Function == "Dietary carbohydrate degradation","Code"]
  compounds <- functions[functions$Function == "Dietary carbohydrate degradation","Compound"]
  distilled_table_sub <- distilled_table[,codes]
  colnames(distilled_table_sub) <- compounds
  aggregated_table <- cbind(aggregated_table,distilled_table_sub)

  #Aggregate dietary lipids
  codes <- functions[functions$Function == "Dietary lipid degradation","Code"]
  compounds <- functions[functions$Function == "Dietary lipid degradation","Compound"]
  distilled_table_sub <- distilled_table[,codes]
  colnames(distilled_table_sub) <- compounds
  aggregated_table <- cbind(aggregated_table,distilled_table_sub)

  #Aggregate dietary proteins
  codes <- functions[functions$Function == "Protein degradation","Code"]
  enzyme_groups <- functions[functions$Function == "Protein degradation","Compound"]
  distilled_table_sub <- distilled_table[,codes]
  distilled_table_sub <- rowMeans(distilled_table_sub)
  distilled_table_sub <- t(t(distilled_table_sub))
  colnames(distilled_table_sub) <- "Peptides"
  aggregated_table <- cbind(aggregated_table,distilled_table_sub)

  #Aggregate mucins
  codes <- functions[functions$Function == "Mucin degradation","Code"]
  compounds <- functions[functions$Function == "Mucin degradation","Compound"]
  distilled_table_sub <- distilled_table[,codes]
  colnames(distilled_table_sub) <- compounds
  aggregated_table <- cbind(aggregated_table,distilled_table_sub)

  #Aggregate SCFAs
  codes <- functions[functions$Function == "SCFA production","Code"]
  compounds <- functions[functions$Function == "SCFA production","Compound"]
  distilled_table_sub <- distilled_table[,codes]
  distilled_table_mer <- merge(t(distilled_table_sub),functions[functions$Function == "SCFA production",c("Code","Compound")],by.x="row.names",by.y="Code")
  distilled_table_agg <- aggregate(subset(distilled_table_mer, select = -c(Row.names,Compound) ),by=list(distilled_table_mer$Compound),FUN=max)
  rownames(distilled_table_agg) <- distilled_table_agg[,1]
  distilled_table_agg <- t(distilled_table_agg[,-1])
  aggregated_table <- cbind(aggregated_table,distilled_table_agg)

  #Aggregate organic anions
  codes <- functions[functions$Function == "Organic anion production","Code"]
  distilled_table_sub <- distilled_table[,codes]
  distilled_table_mer <- merge(t(distilled_table_sub),functions[functions$Function == "Organic anion production",c("Code","Compound")],by.x="row.names",by.y="Code")
  distilled_table_agg <- aggregate(subset(distilled_table_mer, select = -c(Row.names,Compound) ),by=list(distilled_table_mer$Compound),FUN=max)
  rownames(distilled_table_agg) <- distilled_table_agg[,1]
  distilled_table_agg <- t(distilled_table_agg[,-1])
  aggregated_table <- cbind(aggregated_table,distilled_table_agg)

  #Aggregate secondary bile acids
  codes <- functions[functions$Function == "Secondary bile acid production","Code"]
  compounds <- functions[functions$Function == "Secondary bile acid production","Compound"]
  distilled_table_sub <- distilled_table[,codes]
  colnames(distilled_table_sub) <- compounds
  aggregated_table <- cbind(aggregated_table,distilled_table_sub)

  #Aggregate amino acids
  codes <- functions[functions$Function == "Amino acid production","Code"]
  compounds <- functions[functions$Function == "Amino acid production","Compound"]
  distilled_table_sub <- distilled_table[,codes]
  distilled_table_mer <- merge(t(distilled_table_sub),functions[functions$Function == "Amino acid production",c("Code","Compound")],by.x="row.names",by.y="Code")
  distilled_table_agg <- aggregate(subset(distilled_table_mer, select = -c(Row.names,Compound) ),by=list(distilled_table_mer$Compound),FUN=max)
  rownames(distilled_table_agg) <- distilled_table_agg[,1]
  distilled_table_agg <- t(distilled_table_agg[,-1])
  aggregated_table <- cbind(aggregated_table,distilled_table_agg)

  #Aggregate amino acid derivatives
  codes <- functions[functions$Function == "Amino acid derivative production","Code"]
  compounds <- functions[functions$Function == "Amino acid derivative production","Compound"]
  distilled_table_sub <- distilled_table[,codes]
  colnames(distilled_table_sub) <- compounds
  aggregated_table <- cbind(aggregated_table,distilled_table_sub)

  #Aggregate vitamins
  codes <- functions[functions$Function == "Vitamin production","Code"]
  compounds <- functions[functions$Function == "Vitamin production","Compound"]
  distilled_table_sub <- distilled_table[,codes]
  distilled_table_mer <- merge(t(distilled_table_sub),functions[functions$Function == "Vitamin production",c("Code","Compound")],by.x="row.names",by.y="Code")
  distilled_table_agg <- aggregate(subset(distilled_table_mer, select = -c(Row.names,Compound) ),by=list(distilled_table_mer$Compound),FUN=max)
  rownames(distilled_table_agg) <- distilled_table_agg[,1]
  distilled_table_agg <- t(distilled_table_agg[,-1])
  aggregated_table <- cbind(aggregated_table,distilled_table_agg)

  return(aggregated_table)
}
