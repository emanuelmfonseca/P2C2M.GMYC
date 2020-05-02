GMYC.pvalue <- function(Result.n.species, Bootstrap_number_species,perc.treshold){

  p.value <- length(which(Bootstrap_number_species[,1] <= Result.n.species[1,1]*perc.treshold))/length(Bootstrap_number_species[,1])

  P2C2M_GMYC.results <- list()

  P2C2M_GMYC.results[[1]] <- p.value
  P2C2M_GMYC.results[[2]] <- Result.n.species[1,1]*perc.treshold
  P2C2M_GMYC.results[[3]] <- Bootstrap_number_species[,1]

  names(P2C2M_GMYC.results) <- c("p.value","Threshold","Null.distribution")

  return(P2C2M_GMYC.results)

}
