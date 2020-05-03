#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param tree.input PARAM_DESCRIPTION
#' @param tree.format PARAM_DESCRIPTION
#' @param seq PARAM_DESCRIPTION
#' @param nsim PARAM_DESCRIPTION, Default: NULL
#' @param nboot PARAM_DESCRIPTION, Default: NULL
#' @param perc.treshold PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[devtools]{remote-reexports}}
#' @rdname P2C2M_GMYC.pboot
#' @export
#' @importFrom devtools install_github

P2C2M_GMYC.pboot <- function(tree.input,
                             tree.format,
                             seq,
                             nsim=NULL,
                             nboot=NULL,
                             perc.treshold=NULL){

  list.of.packages <- c("ape",
                        "TreeSim",
                        "pegas",
                        "phyclust",
                        "phangorn",
                        "devtools",
                        "paran")

  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)

  if(!"splits" %in% installed.packages()){
    install.packages("splits", repos="http://R-Forge.R-project.org", type="source")
  }

  if(!"phytools" %in% installed.packages()){
    devtools::install_github("liamrevell/phytools")
  }

  suppressMessages(library(ape))
  suppressMessages(library(TreeSim))
  suppressMessages(library(pegas))
  #suppressMessages(library(phybase))
  suppressMessages(library(phyclust))
  suppressMessages(library(pegas))
  suppressMessages(library(phangorn))
  suppressMessages(library(splits))
  suppressMessages(library(phytools))

  quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }

  if(is.null(nsim)){
    nsim <- 100
  }

  if(is.null(nboot)){
    nboot <- 100
  }

  if(is.null(perc.treshold)){
    perc.treshold <- 0.1
  }

  Result.n.species <- data.frame(matrix(NA,nsim+1,1))
  colnames(Result.n.species) <- "N.species"

  Final.result <- data.frame(matrix(NA,nsim,1))
  colnames(Final.result) <- c("N.species")

  seq_length <- nchar(readLines(seq)[2])

  if (tree.format == "nexus"){
    tryCatch(empirical.tree <- read.nexus(tree.input),
             error = function(e) {
               print("Tree format not recognized. Tree must be in nexus or newick format.")
             })
  }else if (tree.format == "newick"){
    tryCatch(empirical.tree <- read.tree(tree.input),
             error = function(e) {
               print("Tree format not recognized. Tree must be in nexus or newick format.")
             })
  }

  if (substr(readLines(seq)[1], 1, 1) == ">" & substr(readLines(seq)[3], 1, 1) == ">"){
    invisible()
  } else {
    stop("Sequence input must be a sequencial fasta file")
  }

  edge <- which(empirical.tree$edge.length < 0.0000001)
  empirical.tree$edge.length[edge] <- 0.001
  empirical.tree <- force.ultrametric(empirical.tree, method="nnls")
  empirical.tree$edge.length<- empirical.tree$edge.length/max(nodeHeights(empirical.tree)[,2])

  n.tips <- length(empirical.tree$tip.label)

  GMYC.empirical.tree <- quiet(gmyc(empirical.tree, method="single", quiet = TRUE))
  GMYC.result.empirical <- capture.output(summary.gmyc(GMYC.empirical.tree))

  print(noquote("Analysis progress:"))

  print(noquote("1/5 - Running GMYC for the empirical dataset"))
  progress <- txtProgressBar(min = 0, max = 100, style = 3)
  setTxtProgressBar(progress,100)
  close(progress)

  Species_result <- GMYC.results(GMYC.empirical.tree, GMYC.result.empirical,empirical.tree)

  N.species <- length(Species_result[,1])

  Result.n.species[1,1] <- N.species

  Coal_result <- coaslescent_tree.Pboot(GMYC.empirical.tree, N.species, seq)

  Species_result <- cbind(Species_result, Coal_result)

  Simulated_trees <- character()

  for (i in 1:nsim){

    if (i == 1){
      print(noquote("2/5 - Simulating trees under the correct model"))
      progress <- txtProgressBar(min = 0, max = nsim, style = 3)
    }

    Yule_tree <- Yule.tree(Species_result, empirical.tree)
    GMYC_simulated.tree <- GMYC.simulated.tree(Yule_tree, Species_result, n.tips)
    Simulated_trees <- c(Simulated_trees, write.tree(GMYC_simulated.tree))

    setTxtProgressBar(progress, i)

    if (i ==nsim){
      close(progress)
    }
  }

  for (i in 1:nsim){

    if (i == 1){
      print(noquote("3/5 - Estimating UPGMA trees"))
      progress <- txtProgressBar(min = 0, max = nsim, style = 3)
      UPGMA_trees <- character()
    }

    Simulated_tree <- read.tree(text = Simulated_trees[[i]])
    tree_UPGMA <- UPGMA.trees(Simulated_tree, seq_length, i)
    UPGMA_trees <- c(UPGMA_trees, write.tree(tree_UPGMA))

    setTxtProgressBar(progress, i)

    if (i == nsim){
      close(progress)
    }}

  for (i in 1:nsim){

    if (i == 1){
      print(noquote("4/5 - Running GMYC for the simulated datasets"))
      progress <- txtProgressBar(min = 0, max = nsim, style = 3)
    }

    Result.n.species[i+1,1] <- quiet(Simulated.n.species.Pboot(UPGMA_trees, i))
    setTxtProgressBar(progress, i)

    if (i ==nsim){
      close(progress)
    }
  }

  for (i in 1:nsim){

    if (i == 1){
      print(noquote("5/5 - Calculating p-value"))
      progress <- txtProgressBar(min = 0, max = nsim, style = 3)
    }

    Final.result[i,1] <- abs(as.numeric(Result.n.species[1,1]) - as.numeric(Result.n.species[i+1,1]))

    setTxtProgressBar(progress,i)

    if (i ==nsim){
      close(progress)
    }
  }

  Bootstrap_n_species <- Bootstrap.n.species(Final.result, nboot, nsim)

  P2C2M_GMYC.res <- GMYC.pvalue(Result.n.species, Bootstrap_n_species, perc.treshold)

  p_value <- P2C2M_GMYC.res$p.value

  if (p_value <= 0.05){
    print(noquote("P2C2M_GMYC result:"))
    print(noquote(paste0("Your data does not violate GMYC model: p-value = ", p_value)))
    print(noquote("Significance: 0.05"))
  } else {
    print(noquote("P2C2M_GMYC result:"))
    print(noquote(paste0("Your data does not violate GMYC model: p-value = ", p_value)))
    print(noquote("Significance: 0.05"))
  }

  assign("P2C2M_GMYC.results", P2C2M_GMYC.res, envir = .GlobalEnv)

}
