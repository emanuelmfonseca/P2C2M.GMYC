Simulated.n.species.PPS <- function(RAxML.Trees, i){

  Simulated_tree <- read.tree(text = RAxML.Trees[[i]])
  Simulated_tree <- force.ultrametric(Simulated_tree, method="nnls")
  Simulated_tree <- multi2di(Simulated_tree)
  edge <- which(Simulated_tree$edge.length < 0.00001)
  Simulated_tree$edge.length[edge] <- 0.001
  Simulated_tree <- force.ultrametric(Simulated_tree, method="nnls")
  Simulated_tree$edge.length <- Simulated_tree$edge.length/max(nodeHeights(Simulated_tree)[,2])

  n.tips <- length(Simulated_tree$tip.label)

  test <- bgmyc.singlephy(Simulated_tree, mcmc=100000, burnin=90000, thinning=100, py2=1.2, t1=1,t2=n.tips,start=c(1,1,30), scale=c(15,20,0.5))

  results.probmat <- spec.probmat(test)
  results.spec <- bgmyc.point(results.probmat, ppcutoff=0.5)

  N.species <- length(results.spec)

  return(N.species)

}
