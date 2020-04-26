GMYC.pvalue <- function(empirical.tree, Result.n.species, Bootstrap_number_species, nboot){
  
  n.tips <- length(empirical.tree$tip.label)
  n_species <- Result.n.species[1,1]
  species.threshold <- round(n.tips * 0.5)
  
  if (n_species <= species.threshold){
    Penalty <- 1
  } else {
    Penalty <-  1/(abs(n_species - n.tips)^2)
  }
  
  p.value <- length(which(Bootstrap_number_species[,1] <= Result.n.species[1,1]*0.1*Penalty))/length(Bootstrap_number_species[,1])
  
  graphic.data <- data.frame(matrix(NA,nboot+1,1))
  colnames(graphic.data) <- "P2C2M_GMYC.results"
  
  graphic.data[1,1] <- Result.n.species[1,1]*0.1*Penalty
  graphic.data[2:(nboot+1),1] <- Bootstrap_number_species[,1]
  
  P2C2M_GMYC.result <- list(p.value,graphic.data)
  
  return(P2C2M_GMYC.result)
  
}
