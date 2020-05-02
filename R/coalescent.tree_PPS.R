coaslescent_tree.PPS <- function(GMYC.empirical.tree, N.species, seq){

  Coal.result <- data.frame(matrix(NA, N.species, 1))
  colnames(Coal.result) <- "N.Species"

  if(N.species == 1){
    fasta <- read.dna(seq, format="fasta")
    Coal.result[,1] <- nuc.div(fasta)
  } else {

    GMYC.result <- GMYC.empirical.tree

    accessions <- list()
    number = 0
    for (k in 1:length(unique(GMYC.result[,1]))){
      if (length(which(GMYC.result[,1] == k)) > 1){
        codes <- which(GMYC.result[,1] == k)
        list <- character()
        for (code in codes){
          list <- c(list, as.character(GMYC.result[code,2]))
        }
        number = number + 1
        accessions[[number]] <- list
      } else {
        number = number + 1
        code <- which(GMYC.result[,1] == k)
        accessions[[number]] <- as.character(GMYC.result[code,2])
      }}

    alignment <- readLines(seq)

    for (k in 1:length(accessions)){
      names <- paste(">", accessions[[k]], sep="")
      for (x in 1:length(names)){
        if (x == 1){
          new_sequences <- character()
        }
        reference <- which(names[x] == alignment)
        new_sequences <- c(new_sequences, alignment[reference])
        new_sequences <- c(new_sequences, alignment[reference + 1])
      }

      x <- new_sequences[seq(1, length(new_sequences),2)]
      y <- new_sequences[seq(2, length(new_sequences),2)]

      x1 <- c(x,y)
      x2 <- structure(c(x1), .Dim = c(as.numeric(length(x)), 2))

      z <- t(sapply(strsplit(x2[,2],""), tolower))
      rownames(z) <- x2[,1]

      fasta <- as.DNAbin(z)

      if (length(labels(fasta)) > 1){
        Coal.result[k,1] <- nuc.div(fasta)
      } else {
        Coal.result[k,1] <- NA
      }}

  }

  return(Coal.result)
}
