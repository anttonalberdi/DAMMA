#' Generate the fullness table of metabolic pathways/modules from a MAG annotation table
#'
#' @param annotations Table containing Genome identifiers and annotation codes
#' @param functions_table Table containing definitions and metadata of metabolic functions (included in DAMMA package)
#' @param abundance_table (Relative) abundance table with samples in columns and Genomes in rows
#' @param fullness_table Pathway fullness table provided by damma() function
#' @param genomecol Column index (number) of the annotations table containing the MAG identifiers
#' @param keggcol Column index(es) of the annotations table in which to search for KEGG KO annotations
#' @param eccol Column index(es) of the annotations table in which to search for Enzyme Commision (EC) annotations
#' @param pepcol Column index(es) of the annotations table in which to search for Peptidase annotations
#' @importFrom stringr str_extract str_match_all
#' @return A fullness matrix
#' @examples
#' damma_community(fullness_table,abundance_table,functions_table)
#' @export

damma_community <- function(annotations,functions_table,abundance_table,fullness_table,genomecol,keggcol,eccol,pepcol){

  tss <- function(abund){sweep(abund, 2, colSums(abund), FUN="/")}

  #Prepare tables
  annotations <- as.data.frame(annotations)
  annotations2 <- annotations[,c(genomecol,keggcol,eccol,pepcol)]
  colnames(annotations2) <- c("MAGs",paste0("K",c(1:length(keggcol))),paste0("E",c(1:length(eccol))),paste0("P",c(1:length(pepcol))))
  functions_table2 <- as.data.frame(functions_table)

  #####
  # Generate Genome- and Community-level fullness scores
  #####

  if(missing(fullness_table)){
    # If fullness table is not provided, run damma() function to create it
    cat("Running DAMMA script for creating the fullness table...\n")
    cat("\tNote: if you have already run DAMMA on this data set,\n")
    cat("\tyou can include the fullness table to avoid this step.\n")
    fullness_table2 <- damma(annotations,functions_table,genomecol,keggcol,eccol,pepcol)
  }else{
    cat("Fullness table has been provided.\n")
    fullness_table2 <- as.data.frame(fullness_table)
  }

  #Compute overall community-level fullness table
  cat("Computing overall community-level MCIs...\n")
  cat("\tExtracting annotations...\n")
  #K00000
  kegg <- str_extract(annotations2$K1, "K[0-9]+")
  kegg <- unique(kegg[!is.na(kegg)])
  #[EC:0.0.0.0]
  EC <- unlist(str_match_all(annotations2$E1, "(?<=\\[EC:).+?(?=\\])"))
  EC <- unique(unlist(strsplit(EC, " ")))
  EC <- EC[!grepl("-", EC, fixed = TRUE)]
  #(EC 0.0.0.0)
  EC2 <- unlist(str_match_all(annotations2$E2, "(?<=\\(EC ).+?(?=\\))"))
  EC2 <- unique(unlist(strsplit(EC2, " ")))
  EC2 <- EC2[!grepl("-", EC2, fixed = TRUE)]
  #Peptidases
  pep <- unique(annotations2$P1)
  pep <- pep[pep != ""]
  #Concatenate all annotations
  present <- unique(c(kegg,EC,EC2,pep))
  #Compute fullness scores
  cat("\tCalculatting fullness scores...\n")
  community_fullness_vector <- c()
  suppressWarnings(
    for(f in c(1:nrow(functions_table2))){
      definition=functions_table2[f,"Definition"]
      fullness <- compute_fullness(definition,present)
      community_fullness_vector <- c(community_fullness_vector,fullness)
    }
  )
  names(community_fullness_vector) <-  functions_table2$Code
  community_fullness_vector[is.na(community_fullness_vector)] <- 0

  #####
  # Generate Genome- and Community-level fullness scores
  #####

  if(!missing(abundance_table)){
    cat("Weighing fullness scores by relative abundances...\n")
    abundance_table2 <- tss(as.data.frame(abundance_table))

    samples <- colnames(abundance_table)

    community_table <- c()
    for(s in samples){
      relabun <- abundance_table2[,s]
      weigedfullness <- sweep(fullness_table2, 1, relabun, FUN = "*")
      community_row <- colSums(weigedfullness) * community_fullness_vector
      community_table <- rbind(community_table,community_row)
    }
    rownames(community_table) <- samples

  }else{
    cat("Only the overall MCI of the Genome catalogue has been calculated.\n")
    cat("\tInput a Genome abundance table to calculate sample-specific MCIs\n")
    cat("\that account for the relative abundances of each Genome in each sample.\n")

    community_table <- colMeans(fullness_table2) * community_fullness_vector

  }

  return(community_table)
}
