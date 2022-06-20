#' Decomposes a definition string of a metabolic pathway/module into its information units
#'
#' @param Definition Definition string of a given metabolic pathway/module
#' @return A decomoposed definition string
#' @examples
#' decompose_definition(definition)
#' decompose_definition("K01580 (K13524,K07250,K00823,K16871) (K00135,K00139,K17761)")

decompose_definition <- function(definition){
    #Separate strings
    def_decomp <- unlist(strsplit(definition, split=''))
    #Merge units
    def_decomp2 <- c()
    unit<-c()
    for (i in c(1:length(def_decomp))){
      chr = def_decomp[i]
      if(grepl("[[:alpha:]]", chr) | grepl("[[:digit:]]", chr) | grepl("_", chr) | grepl("\\.", chr)){
        unit <- paste(unit,chr,collapse = '')
      }else{
        unit <- gsub(" ","",unit)
        def_decomp2 <- c(def_decomp2,unit)
        unit<-c()
        def_decomp2 <- c(def_decomp2,chr)
      }
      if(i == length(def_decomp)){
        unit <- gsub(" ","",unit)
        def_decomp2 <- c(def_decomp2,unit)
      }
    }
    return(def_decomp2)
  }
