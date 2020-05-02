GMYC.results <- function(GMYC.empirical.tree, GMYC.result.empiric, empirical_tree){
  
  if(!grepl("n.s.", GMYC.result.empiric[7])){
    N.species <- as.numeric(sub("\t.+\t", "", GMYC.result.empiric[12]))
  } else {
    N.species <- 1
  }
  
  Species.result <- data.frame(matrix(NA, N.species, 2))
  colnames(Species.result) <- c("Species", "N.samples")
  
  if(N.species > 1){
    
    GMYC_res <- spec.list(GMYC.empirical.tree)
    
    for (u in 1:N.species){
      Species.result[u,1] <- u
      Species.result[u,2] <- length(which(GMYC_res[,1] == u))
    }} else {
      Species.result[1,1] <- N.species
      Species.result[1,2] <- length(empirical_tree$tip.label)
    }
  
  return (Species.result)
}