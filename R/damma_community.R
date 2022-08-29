#' Generate community-level MCIs from set of genomes, optionally accounting for relative abundance of each genome
#'
#' @param annotation_table Table containing Genome identifiers and annotation codes
#' @param pathway_table Table containing definitions and metadata of metabolic functions (included in DAMMA package)
#' @param abundance_table (Relative) abundance table with samples in columns and Genomes in rows. Required for computing sample-specific MCI values
#' @param MCI_table Gene presence-based MCI table provided by damma() function
#' @param completeness_table Genome completeness table to correct MCIs for varying completeness values
#' @param genomecol Column index (number) of annotation_table containing the MAG identifiers
#' @param keggcol Column index(es) of annotation_table in which to search for KEGG KO annotations
#' @param eccol Column index(es) of annotation_table in which to search for Enzyme Commision (EC) annotations
#' @param pepcol Column index(es) of annotation_table in which to search for Peptidase annotations
#' @importFrom stringr str_extract str_match_all
#' @return A pathway fullness vector (if no abundance data are provided) or matrix (if abundance data are provided)
#' @examples
#' damma_community(annotation_table,pathway_table,abundance_table,fullness_table,genomecol,keggcol,eccol,pepcol)
#' damma_community(annotation_table,pathway_table,genomecol,keggcol,eccol,pepcol)
#' @export

damma_community <- function(annotation_table,pathway_table,abundance_table,MCI_table,completeness_table,genomecol,keggcol,eccol,pepcol){

  #Sanity check
  if(missing(annotation_table)) stop("Genome annotation table is missing")
  if(missing(pathway_table)) stop("Pathway table is missing")
  if(missing(genomecol)) stop("Specify a column containing Genome identifiers")
  if(length(genomecol)!=1) stop("The argument genomecol must be an integer indicating the number of the column containing the Genome identifiers in the annotations table")
  if(missing(keggcol) & missing(eccol) & missing(pepcol)) stop("Specify at least one column containing functional annotations")

  #Declare TSS function
  tss <- function(abund){sweep(abund, 2, colSums(abund), FUN="/")}

  #Convert tables into data frames
  annotation_table <- as.data.frame(annotation_table)
  pathway_table <- as.data.frame(pathway_table)
  if(!missing(abundance_table)){abundance_table <- as.data.frame(abundance_table)}
  if(!missing(completeness_table)){completeness_table <- as.data.frame(completeness_table)}

  #####
  # Generate Genome- and Community-level MCIs
  #####

  #Create MCI table if it is not provided
  if(missing(MCI_table)){
    # If fullness table is not provided, run damma() function to create it
    cat("Running DAMMA script for creating the fullness table...\n")
    cat("\tNote: if you have already run DAMMA on this data set,\n")
    cat("\tyou can include the fullness table to avoid this step.\n")
    MCI_table <- damma(annotation_table,pathway_table,genomecol,keggcol,eccol,pepcol)
    #Apply correction factor if completeness_table is provided

    if(!missing(completeness_table)){
      MCI_table <- damma_correction(MCI_table,completeness_table)
    }
  }else{
    cat("Fullness table has been provided.\n")
    MCI_table <- as.data.frame(MCI_table)
  }

  #Compute overall community-level MCI table
  cat("Computing overall community-level MCIs...\n")
  cat("\tExtracting annotations...\n")

  #Declare vector of identifiers
  Identifier_vector <- c()

  #KEGG identifiers
  #K00000
  if(!missing(keggcol)){
  kegg <- str_extract(c(unlist(c(annotation_table[,keggcol]))), "K[0-9]+")
  kegg <- unique(kegg[!is.na(kegg)])
  }else{
  kegg <- c()
  }
  Identifier_vector <- c(Identifier_vector,kegg)

  #Enzyme Commission codes
  #[EC:0.0.0.0]
  if(!missing(eccol)){
  EC <- unlist(str_match_all(c(unlist(c(annotation_table[,eccol]))), "(?<=\\[EC:).+?(?=\\])")) #Extract ECs
  EC <- unique(unlist(strsplit(EC, " "))) #Dereplicate
  EC <- EC[!grepl("-", EC, fixed = TRUE)] #Remove ambiguous codes
  EC <- EC[grepl(".", EC, fixed = TRUE)] #Remove NAs and inproperly formatted codes
  }else{
  EC <- c()
  }
  Identifier_vector <- c(Identifier_vector,EC)

  #Peptidases
  if(!missing(pepcol)){
  pep <- unique(c(unlist(c(annotation_table[,pepcol]))))
  pep <- pep[pep != ""]
  }else{
  pep <- c()
  }
  Identifier_vector <- c(Identifier_vector,pep)

  #Compute MCIs
  cat("\tCalculatting MCIs...\n")
  community_MCI_vector <- c()
  suppressWarnings(
    for(f in c(1:nrow(pathway_table))){
      definition=pathway_table[f,"Definition"]
      MCI <- compute_fullness(definition,Identifier_vector)
      community_MCI_vector <- c(community_MCI_vector,MCI)
    }
  )
  names(community_MCI_vector) <- pathway_table$Code
  community_MCI_vector[is.na(community_MCI_vector)] <- 0

  #####
  # Apply community-level pathway chunckness penalties with or without accounting for relative abundance data
  #####

  if(!missing(abundance_table)){
    cat("Weighing MCIs by relative abundances...\n")
    abundance_table <- tss(as.data.frame(abundance_table))

    samples <- colnames(abundance_table)

    community_MCI_table <- c()
    for(s in samples){
      relabun <- abundance_table[,s]
      weighedMCI <- sweep(MCI_table, 1, relabun, FUN = "*")
      community_row <- mean(c(colSums(weighedMCI),community_MCI_vector))
      community_MCI_table <- rbind(community_MCI_table,community_row)
    }
    rownames(community_MCI_table) <- samples

  }else{
    cat("Only the overall MCI of the Genome catalogue has been calculated.\n")
    cat("\tInput a Genome abundance table to calculate sample-specific MCIs\n")
    cat("\tthat account for the relative abundances of each Genome in each sample.\n")

    community_MCI_table <- mean(c(colMeans(MCI_table),community_MCI_vector))

  }

  return(community_MCI_table)
}
