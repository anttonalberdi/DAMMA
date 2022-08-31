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
#' @importFrom stringr str_extract str_match_all str_detect
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

  cat("Starting DAMMA commmunity analysis\n(Note this may take a while)...\n")

  #Declare TSS function
  tss <- function(abund){sweep(abund, 2, colSums(abund), FUN="/")}

  #Convert tables into data frames
  annotation_table <- as.data.frame(annotation_table)
  pathway_table <- as.data.frame(pathway_table)
  if(!missing(abundance_table)){abundance_table <- as.data.frame(abundance_table)}
  if(!missing(completeness_table)){completeness_table <- as.data.frame(completeness_table)}

  if(missing(abundance_table)){
    #If abundance table does not exist, create a mock abundance table of a single even community
    abundance_table <- rep(1/length(unique(annotation_table[,genomecol])),length(unique(annotation_table[,genomecol])))
    names(abundance_table) <- unique(annotation_table[,genomecol])
    abundance_table <- t(t(abundance_table))
    colnames(abundance_table) <- "Community"

    #Declare single community
    communities <- "Community"
  }else{
    #Declare communities from abundance table
    communities <- colnames(abundance_table)
  }

  #Merge annotations and relative abundance information
  cat("\tMerging annotations and relative abundance data...\n")
  annotation_abundance_table <- merge(annotation_table,tss(abundance_table),by.x=genomecol,by.y="row.names")

  #Declare index (column numbers) of the relative abundance data
  relabun_index <- grep(paste(communities,collapse="|"), colnames(annotation_abundance_table))

  #Filter annotations of 0 abundance genomes
  annotation_abundance_table <- annotation_abundance_table[rowSums(annotation_abundance_table[,relabun_index]) != 0,]

####
# Prepare relative abundance table
####

  id_relabun_table <- c()
  identifier_vector <- c()

  #KEGG identifiers
  #K00000
  if(!missing(keggcol)){
    cat("\t\tExtracting relative abundance data for KEGG identifiers...\n")
    for(col in keggcol){
      column <- annotation_abundance_table[,col]
      kegg_detect <- str_detect(column, "K[0-9]+")
      kegg_detect[is.na(kegg_detect)] <- FALSE
      column_sub <- column[kegg_detect]
      kegg_codes <- unlist(str_match_all(column_sub, "K[0-9]+"))
      annotation_abundance_table_sub <- annotation_abundance_table[kegg_detect,c(col,1,relabun_index)]
      annotation_abundance_table_sub[,1] <- kegg_codes
      colnames(annotation_abundance_table_sub)[1] <- "ID"
      if(nrow(annotation_abundance_table_sub)>0){
        id_relabun_table <- rbind(id_relabun_table,annotation_abundance_table_sub)
      }
    }
  }

  #Enzyme Commission codes
  #[EC:0.0.0.0]
  if(!missing(eccol)){
    cat("\t\tExtracting relative abundance data for EC identifiers...\n")
    for(col in eccol){
      column <- annotation_abundance_table[,col]
      EC_detect <- str_detect(column, "(?<=\\[EC:).+?(?=\\])")
      EC_detect[is.na(EC_detect)] <- FALSE
      column_sub <- column[EC_detect]
      EC_codes <- unlist(str_match_all(column_sub, "(?<=\\[EC:).+?(?=\\])"))
      annotation_abundance_table_sub <- annotation_abundance_table[EC_detect,c(col,1,relabun_index)]
      annotation_abundance_table_sub[,1] <- EC_codes
      colnames(annotation_abundance_table_sub)[1] <- "ID"
      if(nrow(annotation_abundance_table_sub)>0){
        id_relabun_table <- rbind(id_relabun_table,annotation_abundance_table_sub)
      }
    }
  }

  #Peptidases
  #[EC:0.0.0.0]
  if(!missing(pepcol)){
    cat("\t\tExtracting relative abundance data for peptidase family identifiers...\n")
    for(col in pepcol){
      column <- annotation_abundance_table[,col]
      pep_codes <- unique(c(unlist(c(annotation_abundance_table[,col]))))
      pep_codes <- pep_codes[!is.na(pep_codes)]
      pep_codes <- pep_codes[pep_codes != ""]
      annotation_abundance_table_sub <- annotation_abundance_table[annotation_abundance_table[,col] %in% pep_codes, c(col,1,relabun_index)]
      colnames(annotation_abundance_table_sub)[1] <- "ID"
      if(nrow(annotation_abundance_table_sub)>0){
        id_relabun_table <- rbind(id_relabun_table,annotation_abundance_table_sub)
      }
    }
  }

####
# Resolve ambiguities and duplications
####

cat("\tResolving ambiguities and redundancy...\n")

#Remove redundancy
id_relabun_table <- unique(id_relabun_table)

#Split ambiguous rows
ambiguities <- id_relabun_table[grepl(" ", id_relabun_table$ID),"ID"]
for(a in ambiguities){
  elements <- unlist(strsplit(a, " "))
  rows <- id_relabun_table[id_relabun_table$ID == a,]
  if(nrow(rows)>0){
    rows1 <- rows
    rows1$ID <- elements[1]
    #Rename original rows
    id_relabun_table[id_relabun_table$ID == a,] <- rows1
    #Create new rows
    for(e in c(2:length(elements))){
      rows2 <- rows
      rows2$ID <- elements[e]
      id_relabun_table <- rbind(id_relabun_table,rows2)
    }
  }
}

#Remove redundancy (again)
id_relabun_table <- unique(id_relabun_table)

#Remove ambiguous ECs
id_relabun_table <- id_relabun_table[!grepl("-", id_relabun_table$ID),]

cat("\tCalculating community-weighed gene representation values...\n")
#Remove redundancy
id_relabun_table_agg <- aggregate(id_relabun_table[,c(3:ncol(id_relabun_table))],by=list(id_relabun_table[,"ID"]),FUN=sum)
rownames(id_relabun_table_agg) <- id_relabun_table_agg[,1]
id_relabun_table_agg <- id_relabun_table_agg[,-1]

####
# Generate community-specific MCIs
####
cat("\tCalculating community-level MCIs for community:\n")

MCI_table <- c()
m=0
for(community in communities){
  m=m+1
  cat("\t\t",community," (",m,"/",length(communities),")\n", sep = "")
  comm_abun <- id_relabun_table_agg[,community]
  names(comm_abun) <- rownames(id_relabun_table_agg)

  MCI_vector <- c()
  suppressWarnings(
    for(f in c(1:nrow(pathway_table))){
      definition=pathway_table[f,"Definition"]
      MCI <- compute_fullness_community(definition,comm_abun)
      MCI_vector <- c(MCI_vector,MCI)
    }
  )
  #Append MCI vector of the Genome to the MCI table containing MCI values of all Genomes
  MCI_table <- rbind(MCI_table,MCI_vector)
}

  #Format output MCI table
  rownames(MCI_table) <- communities
  colnames(MCI_table) <- pathway_table$Code
  MCI_table[is.na(MCI_table)] <- 0

  #Output MCI table
  return(MCI_table)

}
