UPGMA.trees <- function(Simulated.tree, seq.length, i){
  
  Sequence <- capture.output(capture.output(seqgen(opts = paste0("-mGTR -l", seq.length, " -s1"), newick.tree = write.tree(Simulated.tree)), type="message"))
  
  n.tips <- length(Simulated.tree$tip.label)
  Sequence <- Sequence[1:(n.tips+1)]
  Sequence <- Sequence[-1]
  
  df <- data.frame(matrix(NA,length(Sequence),2))
  
  for (j in 1:length(Sequence)){
    df[j,1] <- sub("(\\d+)\\s+(\\w+)", "\\1", Sequence[j])
    df[j,2] <-  sub(".*\\s+(\\w+)", "\\1", Sequence[j])
  }
  
  df_seq <- t(sapply(strsplit(df[,2],""), tolower))
  rownames(df_seq) <- df[,1]
  
  seq <- as.DNAbin(df_seq)
  seq <- phyDat(seq, type = "DNA", levels = NULL)
  
  D <- dist.ml(seq, exclude = "pairwise")
  D[is.infinite(D)] <- 1
  D[is.nan(D)] <- 0
  
  tr <- phangorn::upgma(D)
  edge <- which(tr$edge.length < 0.00001)
  tr$edge.length[edge] <- 0.01
  tr$edge.length <- tr$edge.length/max(nodeHeights(tr)[,2])
  tr <- force.ultrametric(tr, method="nnls")
  tr$edge.length <- tr$edge.length/max(nodeHeights(tr)[,2])

  return(tr)
}
