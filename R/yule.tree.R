Yule.tree <- function(Species.result, empirical.tree){

  tips <- empirical.tree$tip.label
  t <- drop.tip(empirical.tree, tips,trim.internal=F)
  yule <- yule(t)
  yule_lambda <- rnorm(1,yule$lambda, yule$se)

  if (length(Species.result[,1]) > 1){
    yule.tree<-sim.bd.taxa(n=length(Species.result[,1]), numbsim=1, lambda=yule_lambda, mu=0, frac = 1, complete = TRUE, stochsampling = FALSE)[[1]]
    yule.tree$tip.label<-paste(rep(paste("Tip", 1:length(Species.result[,1]), ".1", sep="")))
    yule.tree$tip.label<-sample(yule.tree$tip.label)
    yule.tree$tip.label <- sample(yule.tree$tip.label)
    yule.tree$edge.length[yule.tree$edge[,2] <= Ntip(yule.tree)]<- yule.tree$edge.length[yule.tree$edge[,2] <= Ntip(yule.tree)] + 0.05
    yule.tree$edge.length<- yule.tree$edge.length/max(nodeHeights(yule.tree)[,2])

    return(yule.tree)

  } else {

    return(NULL)

  }}
