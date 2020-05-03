Simulated.n.species.Pboot <- function(RAxML.Trees, i){

  Simulated_tree <- read.tree(text = RAxML.Trees[[i]])
  Simulated_tree <- force.ultrametric(Simulated_tree, method="nnls")
  Simulated_tree <- multi2di(Simulated_tree)
  edge <- which(Simulated_tree$edge.length < 0.00001)
  Simulated_tree$edge.length[edge] <- 0.001
  Simulated_tree <- force.ultrametric(Simulated_tree, method="nnls")
  Simulated_tree$edge.length <- Simulated_tree$edge.length/max(nodeHeights(Simulated_tree)[,2])

  GMYC.simulated.tree<-gmyc(Simulated_tree, method="single",quiet=T)
  GMYC.result.simulation <- capture.output(summary.gmyc(GMYC.simulated.tree))

  if(!grepl("n.s.", GMYC.result.simulation[7])){
    N.species <- as.numeric(sub("	.+	", "", GMYC.result.simulation[12]))
  } else {
    N.species <- 1
  }

  return(N.species)

}
