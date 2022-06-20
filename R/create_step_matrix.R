#' Creates a matrix that specifies the hierarchical structure of the definition of a metabolic pathway/module
#'
#' @param def_decomp Vector of decomposed units of a definition string of a metabolic pathway/module (produced by decompose_definition())
#' @param def_level Vector of the hierarchical level of decomposed units of a definition string of a metabolic pathway/module (produced by set_levels())
#' @return A decomposed hierarchy matrix
#' @examples
#' create_step_matrix(definition_decomposed,definition_levels)

create_step_matrix <- function(def_decomp,def_level){
  def_table <- as.data.frame(cbind(def_decomp,def_level))
  #L0 units
  L0_splits <- rownames(def_table[((def_decomp == " ") | (def_decomp == "+") | (def_decomp == ",")) & (def_level == 0),])
  L0_group <- c()
  group <- 1
  for (i in c(1:length(def_decomp))){
    if(i %in% L0_splits){group=group+1}
    L0_group <- c(L0_group,group)
  }
  def_table[,"L0_group"] <- L0_group

  #L1 units
  L1_splits <- rownames(def_table[((def_decomp == " ") | (def_decomp == "+") | (def_decomp == ",")) & (def_level <= 1),])
  L1_group <- c()
  group <- 1
  for (i in c(1:length(def_decomp))){
    if(i %in% L1_splits){group=group+1}
    if(def_table[i,"def_level"]>=1){
      L1_group <- c(L1_group,group)
    }else{
      L1_group <- c(L1_group,NA)
      group=1
    }
  }
  def_table[,"L1_group"] <- L1_group

  #L2 units
  L2_splits <- rownames(def_table[((def_decomp == " ") | (def_decomp == "+") | (def_decomp == ",")) & (def_level <= 2),])
  L2_group <- c()
  group <- 1
  for (i in c(1:length(def_decomp))){
    if(i %in% L2_splits){group=group+1}
    if(def_table[i,"def_level"]>=2){
      L2_group <- c(L2_group,group)
    }else{
      L2_group <- c(L2_group,NA)
      group=1
    }
  }
  def_table[,"L2_group"] <- L2_group

  #L3 units
  L3_splits <- rownames(def_table[((def_decomp == " ") | (def_decomp == "+") | (def_decomp == ",")) & (def_level <= 3),])
  L3_group <- c()
  group <- 1
  for (i in c(1:length(def_decomp))){
    if(i %in% L3_splits){group=group+1}
    if(def_table[i,"def_level"]>=3){
      L3_group <- c(L3_group,group)
    }else{
      L3_group <- c(L3_group,NA)
      group=1
    }
  }
  def_table[,"L3_group"] <- L3_group

  #L4 units
  L4_splits <- rownames(def_table[((def_decomp == " ") | (def_decomp == "+") | (def_decomp == ",")) & (def_level <= 4),])
  L4_group <- c()
  group <- 1
  for (i in c(1:length(def_decomp))){
    if(i %in% L4_splits){group=group+1}
    if(def_table[i,"def_level"]>=4){
      L4_group <- c(L4_group,group)
    }else{
      L4_group <- c(L4_group,NA)
      group=1
    }
  }
  def_table[,"L4_group"] <- L4_group

  #L5 units
  L5_splits <- rownames(def_table[((def_decomp == " ") | (def_decomp == "+") | (def_decomp == ",")) & (def_level <= 5),])
  L5_group <- c()
  group <- 1
  for (i in c(1:length(def_decomp))){
    if(i %in% L5_splits){group=group+1}
    if(def_table[i,"def_level"]>=5){
      L5_group <- c(L5_group,group)
    }else{
      L5_group <- c(L5_group,NA)
      group=1
    }
  }
  def_table[,"L5_group"] <- L5_group

  #CLEAN TABLE
  def_table[(def_table$def_decomp == "(") | (def_table$def_decomp == ")"),] <- NA
  return(def_table)
}
