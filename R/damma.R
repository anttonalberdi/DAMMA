#' Generates a gene presence-based Metabolic Capacity Index (MCI) table from a bacterial genome annotation table
#'
#' @param annotation_table Table containing Genome identifiers and gene annotations
#' @param pathway_table Table containing definitions and metadata of metabolic functions (provided by DAMMA)
#' @param genomecol Column index (number) of the annotation_table containing the genome identifiers
#' @param keggcol Column index(es) of the annotation_table in which to search for KEGG KO annotations
#' @param eccol Column index(es) of the annotation_table in which to search for Enzyme Commision (EC) annotations
#' @param pepcol Column index(es) of the annotation_table in which to search for Peptidase annotations
#' @importFrom stringr str_extract str_match_all
#' @return A pathway-level MCI table
#' @examples
#' damma(annotation_table,pathway_table,genomecol,keggcol,eccol,pepcol)
#' damma(annotation_table,pathway_table,genomecol=2,keggcol=9,eccol=c(10,19),pepcol=12)
#' @export

damma <- function(annotation_table,pathway_table,genomecol,keggcol,eccol,pepcol){

  #Sanity check
  if(missing(annotation_table)) stop("Genome annotation table is missing")
  if(missing(pathway_table)) stop("Pathway table is missing")
  if(missing(genomecol)) stop("Specify a column containing Genome identifiers")
  if(length(genomecol)!=1) stop("The argument genomecol must be an integer indicating the number of the column containing the Genome identifiers in the annotations table")
  if(missing(keggcol) & missing(eccol) & missing(pepcol)) stop("Specify at least one column containing functional annotations")

  #Convert annotation table to data frame
  annotation_table <- as.data.frame(annotation_table)

  #Convert pathway table to data frame
  pathway_table <- as.data.frame(pathway_table)

  #List Genomes
  Genomes <- unique(annotation_table[,genomecol])

  #Calculate MCI's for each Genome
  MCI_table <- c()
  cat("Calculating MCIs for Genome:\n")
  m=0
  for(Genome in Genomes){
    m=m+1
    cat("\t",Genome," (",m,"/",length(Genomes),")\n", sep = "")
    #Fetch Genome annotations
    annotations_Genome <- annotation_table[annotation_table[,genomecol] == Genome,]

    #Declare vector of identifiers
    Identifier_vector <- c()

    #KEGG identifiers
    #K00000
    if(!missing(keggcol)){
    kegg <- str_extract(c(unlist(c(annotations_Genome[,keggcol]))), "K[0-9]+")
    kegg <- unique(kegg[!is.na(kegg)])
    }else{
    kegg <- c()
    }
    Identifier_vector <- c(Identifier_vector,kegg)

    #Enzyme Commission codes
    #[EC:0.0.0.0]
    if(!missing(eccol)){
    EC <- unlist(str_match_all(c(unlist(c(annotations_Genome[,eccol]))), "(?<=\\[EC:).+?(?=\\])")) #Extract ECs
    EC <- unique(unlist(strsplit(EC, " "))) #Dereplicate
    EC <- EC[!grepl("-", EC, fixed = TRUE)] #Remove ambiguous codes
    EC <- EC[grepl(".", EC, fixed = TRUE)] #Remove NAs and inproperly formatted codes
    }else{
    EC <- c()
    }
    Identifier_vector <- c(Identifier_vector,EC)

    #Peptidases
    if(!missing(pepcol)){
    pep <- unique(c(unlist(c(annotations_Genome[,pepcol]))))
    pep <- pep[pep != ""]
    }else{
    pep <- c()
    }
    Identifier_vector <- c(Identifier_vector,pep)

    #Calculate MCI's for each Pathway and append to vector
    MCI_vector <- c()
    suppressWarnings(
      for(f in c(1:nrow(pathway_table))){
        definition=pathway_table[f,"Definition"]
        MCI <- compute_MCI(definition,Identifier_vector)
        MCI_vector <- c(MCI_vector,MCI)
      }
    )
    #Append MCI vector of the Genome to the MCI table containing MCI values of all Genomes
    MCI_table <- rbind(MCI_table,MCI_vector)
  }

  #Format output MCI table
  rownames(MCI_table) <- Genomes
  colnames(MCI_table) <- pathway_table$Code_pathway
  MCI_table[is.na(MCI_table)] <- 0

  #Output MCI table
  return(MCI_table)
}
