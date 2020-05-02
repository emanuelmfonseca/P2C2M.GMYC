Bootstrap.n.species <- function(Final_result, n.boot, n.sim){
  for (i in 1:n.boot){
    
    if (i == 1){
      Bootstrap <- data.frame(matrix(NA,n.boot,1))
      colnames(Bootstrap) <- "Bootstrap.n.species"
    }
    
    n <- 1:n.sim
    r <- sample(1:n.sim,n.sim*0.2,r=F)  
    n2 <- n[-r]
    n3 <- c(n2, sample(n2,n.sim*0.2,r=T))
    n4 <- Final_result[n3,1]
    Bootstrap[i,1] <- sd(n4)
  }
  
  return(Bootstrap)
  
}
