GMYC.simulated.tree <- function(Yule.tree, Species.result, n.tips){

  if(length(Species.result[,1]) != n.tips & length(Species.result[,1]) != 1){

  number = 0
  for (u in 1:length(Species.result[,1])){
    if (Species.result[u,2] > 1){
      if (number == 0){
        number = u
        sim <- sim.coaltree(nspecies=Species.result[u,2],theta=Species.result[u,3])
        sim <-read.tree(text=paste(sim,";",sep=""))
        sim$tip.label<-paste(rep(paste("Tip", Species.result[u,1], ".", sep=""), e=Species.result[u,2]), 1:Species.result[u,2], sep="")

        if (any(Species.result[,1] > 1)){
          Grafted.tree <- bind.tree(Yule.tree, sim, where=which(Yule.tree$tip.label == paste("Tip", u, ".1", sep="")))
        }} else {
          sim <- sim.coaltree(nspecies=Species.result[u,2],theta=Species.result[u,3])
          sim<-read.tree(text=paste(sim,";",sep=""))
          sim$tip.label<-paste(rep(paste("Tip", Species.result[u,1], ".", sep=""), e=Species.result[u,2]), 1:Species.result[u,2], sep="")
          Grafted.tree <- bind.tree(Grafted.tree, sim, where=which(Grafted.tree$tip.label == paste("Tip", u, ".1", sep="")))
          Grafted.tree <- force.ultrametric(Grafted.tree, method="nnls")
          Grafted.tree$edge.length<- Grafted.tree$edge.length/max(nodeHeights(Grafted.tree)[,2])
        }}}
  } else if (length(Species.result[,1]) == n.tips){
    Grafted.tree <- Yule.tree
  } else {
    sim <- sim.coaltree(nspecies=Species.result[1,2],theta=Species.result[1,3])
    sim <-read.tree(text=paste(sim,";",sep=""))
    sim$tip.label<-paste(rep(paste("Tip", Species.result[1,1], ".", sep=""), e=Species.result[1,2]), 1:Species.result[1,2], sep="")
    Grafted.tree <- sim
  }
  ci <- sum(dist(coalescent.intervals(Grafted.tree)$interval.length))
  cat(ci, file="/Users/emanuelfonseca/Desktop/P2C2M_low_popsize.txt", append=TRUE, sep = "\n")

  return (Grafted.tree)

}
