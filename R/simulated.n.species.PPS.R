Simulated.n.species.PPS <- function(RAxML.Trees, i, mcmc, burnin, thinning, py2, t1, t2, start, scale){

  Simulated_tree <- read.tree(text = RAxML.Trees[[i]])
  Simulated_tree <- force.ultrametric(Simulated_tree, method="nnls")
  Simulated_tree <- multi2di(Simulated_tree)
  edge <- which(Simulated_tree$edge.length < 0.00001)
  Simulated_tree$edge.length[edge] <- 0.001
  Simulated_tree <- force.ultrametric(Simulated_tree, method="nnls")
  Simulated_tree$edge.length <- Simulated_tree$edge.length/max(nodeHeights(Simulated_tree)[,2])

  n.tips <- length(Simulated_tree$tip.label)

  test <- bgmyc.singlephy(Simulated_tree, mcmc=mcmc, burnin=burnin, thinning=thinning, py2=py2, t1=t1,t2=t2,start=start, scale=scale)

  results.probmat <- spec.probmat(test)
  results.spec <- bgmyc.point(results.probmat, ppcutoff=0.5)

  N.species <- length(results.spec)

  return(N.species)

}
