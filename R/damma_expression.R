#' Generates a gene abundance/expression-based Metabolic Capacity Index (MCI) table from a bacterial genome annotation table and a related gene expression table
#'
#' @param expression_table Table containing normalised gene expression data with gene identifiers in rows and sample identifiers in columns
#' @param annotation_table Table containing gene and genome identifiers, and gene annotations
#' @param pathway_table Table containing definitions and metadata of metabolic functions (provided by DAMMA)
#' @param genecol Column index (number) of the annotations table containing the gene identifiers
#' @param genomecol Column index (number) of annotation_table containing the genome identifiers
#' @param keggcol Column index(es) of annotation_table in which to search for KEGG KO annotations
#' @param eccol Column index(es) of annotation_table in which to search for Enzyme Commision (EC) annotations
#' @param pepcol Column index(es) of annotation_table in which to search for Peptidase annotations
#' @importFrom stringr str_extract str_match_all
#' @return A list of pathway-expression matrices (one table per genome)
#' @examples
#' damma_expression(expression_table,annotation_table,pathway_table,genecol,genomecol,keggcol,eccol,pepcol)
#' damma_expression(expression_table,annotation_table,pathway_table,genecol=1,genomecol=2,keggcol=9,eccol=c(10,19),pepcol=12)
#' @export

#UNDER DEVELOPMENT
damma_expression <- function(expression_table,annotation_table,pathway_table,genecol,genomecol,keggcol,eccol,pepcol){

  #Sanity check
  if(missing(expression_table)) stop("Gene expression table is missing")
  if(missing(annotation_table)) stop("Genome annotation table is missing")
  if(missing(pathway_table)) stop("Pathway table is missing")
  if(missing(genecol)) stop("Specify a column containing Gene identifiers")
  if(missing(genomecol)) stop("Specify a column containing Genome identifiers")
  if(length(genecol)!=1) stop("The argument genecol must be an integer indicating the number of the column containing the Gene identifiers in the annotations table")
  if(length(genomecol)!=1) stop("The argument genomecol must be an integer indicating the number of the column containing the Genome identifiers in the annotations table")
  if(missing(keggcol) & missing(eccol) & missing(pepcol)) stop("Specify at least one column containing functional annotations")

  cat("Starting DAMMA expression analysis\n(Note this may take a while)...\n")

  #Convert annotation table to data frame
  annotation_table <- as.data.frame(annotation_table)

  #Validate expression table
  sharedgenes <- intersect(rownames(expression_table),annotation_table[,genecol])

  #Filter expression table
  expression_table <- as.data.frame(expression_table[sharedgenes,])

  #Filter annotations table
  annotation_table <- annotation_table[annotation_table[,genecol] %in% sharedgenes,]

  #List Genomes
  Genomes <- unique(annotation_table[,genomecol])

  #Calculate expression values for each Genome
  cat("Calculating gene expression-based MCIs for Genome:\n")
  m=0
  expression_MCI_table_list <- list()
  for(Genome in Genomes){
    m=m+1
    cat("\t",Genome," (",m,"/",length(Genomes),")\n", sep = "")

    #Subset annotation data for the specific Genome
    annotations_Genome <- annotation_table[annotation_table[,genomecol] == Genome,]

    #Declare expression table
    expression_MCI_table <- data.frame()

    #KEGG identifiers
    if(!missing(keggcol)){
      cat("\t\tProcessing KEGG annotations...\n", sep = "")
      kegg <- str_extract(c(unlist(c(annotations_Genome[,keggcol]))), "K[0-9]+")
      kegg <- sort(unique(kegg[!is.na(kegg)]))
      for(k in kegg){
        genes <- annotations_Genome[grep(k, annotations_Genome[,c(keggcol)]),genecol]
        expression_kegg <- expression_table[genes,]
        if(dim(expression_kegg)[1]>1){
          expression_kegg <- colSums(expression_kegg,na.rm=TRUE)
          expression_kegg <- t(expression_kegg)
          rownames(expression_kegg) <- k
          expression_MCI_table <- rbind(expression_MCI_table,expression_kegg)
        }
      }
    }

    #Enzyme Commission codes
    if(!missing(eccol)){
      cat("\t\tProcessing EC annotations...\n", sep = "")
      EC <- unlist(str_match_all(c(unlist(c(annotations_Genome[,eccol]))), "(?<=\\[EC:).+?(?=\\])")) #Extract ECs
      EC <- unique(unlist(strsplit(EC, " "))) #Dereplicate
      EC <- EC[!grepl("-", EC, fixed = TRUE)] #Remove ambiguous codes
      EC <- sort(EC[grepl(".", EC, fixed = TRUE)]) #Remove NAs and inproperly formatted codes and sort codes
      for(e in EC){
        genes <- annotations_Genome[grep(e, annotations_Genome[,c(eccol)]),genecol]
        expression_EC <- expression_table[genes,]
        if(dim(expression_EC)[1]>1){
          expression_EC <- colSums(expression_EC,na.rm=TRUE)
          expression_EC <- t(expression_EC)
          rownames(expression_EC) <- e
          expression_MCI_table <- rbind(expression_MCI_table,expression_EC)
        }
      }
    }

    #Peptidases
    if(!missing(pepcol)){
      cat("\t\tProcessing peptidase annotations...\n", sep = "")
      pep <- unique(c(unlist(c(annotations_Genome[,pepcol]))))
      pep <- pep[pep != ""]
      for(p in pep){
        genes <- annotations_Genome[grep(p, annotations_Genome[,c(pepcol)]),genecol]
        expression_pep <- expression_table[genes,]
        if(dim(expression_pep)[1]>1){
          expression_pep <- colSums(expression_pep,na.rm=TRUE)
          expression_pep <- t(expression_pep)
          rownames(expression_pep) <- p
          expression_MCI_table <- rbind(expression_MCI_table,expression_pep)
        }
      }
    }

    #Compute expression scores
    cat("\t\tCalculating gene expression-based MCIs for\n")
    cat("\t\t",nrow(pathway_table),"pathways in",ncol(expression_MCI_table),"samples...\n")
    suppressWarnings(
      for(f in c(1:nrow(pathway_table))){
        definition=pathway_table[f,"Definition"]
        expression_MCI <- compute_fullness_expression(definition,expression_MCI_table)
        if(f == 1){
          #Create list if it is the first function
          expression_MCI_list <- expression_MCI
        }else{
          #Append to list if it is not the first function
          expression_MCI_list <- Map(c, expression_MCI_list, expression_MCI)
        }
      }
    )
    #Convert sample list to matrix
    expression_MCI_list <- lapply(expression_MCI_list,function(x) as.numeric(x))
    expression_MCI_table <- do.call(rbind, expression_MCI_list)
    colnames(expression_MCI_table) <- pathway_table$Code

    #Append to Genome list
    expression_MCI_table_list[[Genome]] <- expression_MCI_table
  }

  return(expression_MCI_table_list)

}
