% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/P2C2M_GMYC_PPS.R
\name{P2C2M_GMYC.PPS}
\alias{P2C2M_GMYC.PPS}
\title{Run P2C2M.GMYC - posterior predictive simulations}
\usage{
P2C2M_GMYC.PPS(tree.input, seq, nsim = NULL, nboot = NULL,
     ntree = NULL, perc.treshold = NULL, mcmc = NULL,
     burnin = NULL, thinning = NULL, py1 = NULL, py2 = NULL,
     pc1 = NULL, pc2 = NULL, t1 = NULL, t2 = NULL,
     scale = NULL, start = NULL, ppcutoff = NULL)
}
\arguments{
\item{tree.input}{A file name which contains the posterior distribution of gene genealogies.}

\item{seq}{A file name which contains a sequential FASTA alignment.}

\item{nsim}{Number of simulated trees, Default: 100 trees.}

\item{nboot}{The size of the null distribution, Default: 100.}

\item{ntree}{Number of trees sampled from the posterior distribution of gene genealogies, Default: 100}

\item{perc.treshold}{Percentage of the number of species recovered by the empirical dataset used to set the threshold to calculate the p-value, Default: 0.1 (10\%).}

\item{mcmc}{bGMYC package: number of samples to take from the Markov Chain, Default: 100000.}

\item{burnin}{bGMYC package: the number of samples to discard as burn-in, Default: 90000.}

\item{thinning}{bGMYC package: the interval at which samples are retained from the Markov Chain, Default: 100}

\item{py1}{bGMYC package: governs the prior on the Yule (speciation) rate change parameter. using the default prior distribution, this is the lower bound of a uniform distribution. this can be the most influential prior of the three. rate change is parameterized as n^py where n is the number of lineages in a waiting interval (see Pons et al. 2006). if there are 50 sequences in an analysis and the Yule rate change parameter is 2, this allows for a potential 50-fold increase in speciation rate. this unrealistic parameter value can cause the threshold between Yule and Coalescent process to be difficult to distinguish. are more reasonable upper bound for the prior would probably be less than 1.5 (a potential 7-fold increase). Or you could modify the prior function to use a different distribution entirely, Default: 0.}

\item{py2}{bGMYC package: governs the prior on the Yule rate change parameter. using the default prior distribution, this is the upper bound of a uniform distribution., Default: 2.}

\item{pc1}{bGMYC package: governs the prior on the coalescent rate change parameter. using the default prior distribution, this is the lower bound of a uniform distribution. rate change is parameterized as (n(n-1))^pc where n is the number of lineages in a waiting interval (see Pons et al. 2006). In principle pc can be interpreted as change in effective population size (pc<1 decline, pc>1 growth) but because identical haplotypes must be excluded from this analysis an accurate biological interpretation is not possible, Default: 0.}

\item{pc2}{bGMYC package: governs the prior on the coalescent rate change parameter. using the default prior distribution, this is the upper bound of a uniform distribution, Default: 2.}

\item{t1}{bGMYC package: governs the prior on the threshold parameter. the lower bound of a uniform distribution. the bounds of this uniform distribution should not be below 1 or greater than the number of unique haplotypes in the analysis, Default: 1}

\item{t2}{bGMYC package: governs the prior on the threshold parameter. the upper bound of a uniform distribution, Default: number of tips in the phylogenetic tree}

\item{scale}{bGMYC package: a vector of scale parameters governing the proposal distributions for the markov chain. the first to are the Yule and coalescent rate change parameters. increasing them makes the proposals more conservative. the third is the threshold parameter. increasing it makes the proposals more liberal, Default: c(20, 10, 5)}

\item{start}{bGMYC package: a vector of starting parameters in the same order as the scale parameters, py, pc, t. t may need to be set so that it is not impossible given the dataset, Default: c(1, 0.5, 50)}

\item{ppcutoff}{Conspecificity probability threshold, Default: 0.5}
}
\value{
A list containing the p-value, the threshold used to calculate the p-value, and the null distribution.
}
\description{
Identifying model violations under the Generalized Mixed Yule Coalescent (GMYC) model using posterior predictive simulations.
}
\details{
This function take a ultrametric phylogenetic tree and a sequential FASTA alignment and test for model violation of the GMYC model using posterior predictive simulations.
}
\examples{
## Lygodactylus dataset

P2C2M_GMYC.PPS(tree.input=system.file("extdata", "Lygodactylus.trees", package="P2C2M.GMYC"),
               seq=system.file("extdata", "Lygodactylus.fas", package="P2C2M.GMYC"))
}

