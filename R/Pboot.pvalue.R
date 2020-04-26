GMYC.pvalue <- function(empirical.tree, Result.n.species, Bootstrap_number_species, nboot, nspecies.penalty,perc.treshold){

  n.tips <- length(empirical.tree$tip.label)
  n_species <- Result.n.species[1,1]
  species.threshold <- round(n.tips * nspecies.penalty)

  if (n_species <= species.threshold){
    Penalty <- 1
  } else {
    Penalty <-  1/(abs(n_species - n.tips)^2)
  }

  p.value <- length(which(Bootstrap_number_species[,1] <= Result.n.species[1,1]*perc.treshold*Penalty))/length(Bootstrap_number_species[,1])

  P2C2M_GMYC.results <- list()

  P2C2M_GMYC.results[[1]] <- p.value
  P2C2M_GMYC.results[[2]] <- Result.n.species[1,1]*perc.treshold*Penalty
  P2C2M_GMYC.results[[3]] <- Bootstrap_number_species[,1]

  names(P2C2M_GMYC.results) <- c("p.value","Threshold","Null.distribution")

  return(P2C2M_GMYC.results)

}
