GMYC.simulated.tree <- function(Yule.tree, Species.result){
  
  if(length(Species.result[,1]) > 1){
  
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
        } else {
          Grafted.tree <- sim
        }} else {
          sim <- sim.coaltree(nspecies=Species.result[u,2],theta=Species.result[u,3])
          sim<-read.tree(text=paste(sim,";",sep=""))
          sim$tip.label<-paste(rep(paste("Tip", Species.result[u,1], ".", sep=""), e=Species.result[u,2]), 1:Species.result[u,2], sep="")
          Grafted.tree <- bind.tree(Grafted.tree, sim, where=which(Grafted.tree$tip.label == paste("Tip", u, ".1", sep="")))
          Grafted.tree <- force.ultrametric(Grafted.tree, method="nnls")
          Grafted.tree$edge.length<- Grafted.tree$edge.length/max(nodeHeights(Grafted.tree)[,2])
        }}}
  } else {
    Grafted.tree <- Yule.tree
  }
  
  return (Grafted.tree)
  
}
