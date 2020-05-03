P2C2M_GMYC.PPS <- function(tree.input,
                           seq,
                           nsim=NULL,
                           nboot=NULL,
                           ntree = NULL,
                           perc.treshold = NULL,
                           mcmc = NULL,
                           burnin = NULL,
                           thinning = NULL,
                           py1 = NULL,
                           py2 = NULL,
                           pc1 = NULL,
                           pc2 = NULL,
                           t1 = NULL,
                           t2 = NULL,
                           scale = NULL,
                           start = NULL,
                           ppcutoff= NULL){

  list.of.packages <- c("ape",
                        "TreeSim",
                        "pegas",
                        "phyclust",
                        "phangorn",
                        "devtools")

  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)


  if(!"phytools" %in% installed.packages()){
    devtools::install_github("liamrevell/phytools")
  }

  if(!"bGMYC" %in% installed.packages()){
    install.packages(system.file("bGMYC_1.0.2.tar", package="P2C2M.GMYC"), repos = NULL, type="source")
  }

  suppressMessages(library(ape))
  suppressMessages(library(TreeSim))
  suppressMessages(library(pegas))
  #suppressMessages(library(phybase))
  suppressMessages(library(phyclust))
  suppressMessages(library(phangorn))
  suppressMessages(library(phytools))
  suppressMessages(library(bGMYC))

  quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }

  path.files <- getwd()

  params <- c("nsim", "nboot", "ntree", "perc.treshold", "mcmc", "burnin", "thinning", "ppcutoff", "py1", "py2", "pc1", "pc2", "t1", "scale", "start")

  params.default <- c("100", "100", "100", "0.1", "100000", "90000", '100', "0.5", "0", '2', "0", "2", "2", "c(20, 10, 5)", "c(1, 0.5, 50)")


  for (x in 1:length(params)){
    if(is.null(get(params[x]))){
    assign(params[x], eval(parse(text=params.default[x])))
  }}

  seq_length <- nchar(readLines(seq)[2])

  Result.n.species <- data.frame(matrix(NA,nsim+1,1))
  colnames(Result.n.species) <- "N.species"

  Final.result <- data.frame(matrix(NA,nsim,1))
  colnames(Final.result) <- c("N.species")

  empirical.tree <- read.nexus(tree.input)
  empirical.tree <- empirical.tree[-c(1:round((length(empirical.tree)*0.1)))]
  empirical.tree <- empirical.tree[sample(1:9000,ntree)]

  n.tips <- length(empirical.tree[[1]]$tip.label)

  if(is.null(t2)){
    t2 <- n.tips
  }

  if (substr(readLines(seq)[1], 1, 1) == ">" & substr(readLines(seq)[3], 1, 1) == ">"){
    invisible()
  } else {
    stop("Sequence input must be a sequencial fasta file")
  }

  setwd(path.files)

  bGMYC.results <- bgmyc.multiphylo(empirical.tree, mcmc=mcmc, burnin=burnin, thinning=thinning, py2=py2, t1=t1,t2=t2,start=start, scale=scale)

  results.probmat <- spec.probmat(bGMYC.results)
  results.spec <- bgmyc.point(results.probmat, ppcutoff=0.5)

  results <- data.frame(matrix(NA,n.tips,2))

  number <- 0
  for (x in 1:length(results.spec)){
    sp <- results.spec[x]
    for (i in 1:length(sp[[1]])){
      number <- number + 1
      results[number,1] <- x
      results[number,2] <- sp[[1]][i]
    }
  }

  Species_result <- data.frame(matrix(NA,length(unique(results[,1])),2))
  Species_result[,1] <- 1:length(unique(results[,1]))
  Species_result[,2] <- table(results[,1])

  N.species <- length(Species_result[,1])

  Result.n.species[1,1] <- N.species

  Coal_result <- coaslescent_tree.PPS(results, N.species, seq)

  Species_result <- cbind(Species_result, Coal_result)
  Species_result[,3][is.na(Species_result[,3])] <- 0

  Simulated_trees <- character()

  for (i in 1:nsim){

    if (i == 1){
      print(noquote("2/5 - Simulating trees under the correct model"))
      progress <- txtProgressBar(min = 0, max = nsim, style = 3)
    }

    empirical.Tree <- empirical.tree[[i]]
    Yule_tree <- Yule.tree(Species_result, empirical.Tree)
    GMYC_simulated.tree <- GMYC.simulated.tree(Yule_tree, Species_result)
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
    Result.n.species[i+1,1] <- Simulated.n.species.PPS(UPGMA_trees, i, mcmc, burnin, thinning, py2, t1, t2, start, scale)
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

  Bootstrap_n.species <- Bootstrap.n.species(Final.result, nboot, nsim)

  P2C2M_GMYC.res <- GMYC.pvalue(Result.n.species, Bootstrap_n.species, perc.treshold)

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

