# Coverage the variational bayes posterior to the mcmc posterior
#
library(survival)
library(survival.svb)
library(survival.ss)

source("00-functions.R")
source("../application/tcga_data.R", chdir=T)


d <- function(a, b) if(Sys.getenv(a) != "") as(Sys.getenv(a), class(b)) else b 
CORES <- d("CORES", 8)                 # number of CPU cores to use
SIMNUM <- d("SIMNUM", 1)               # which simulation parameters to use
DGP <- d("DGP", 1:5)                   # which data generating process to use


# ----------------------------------------
# Simulation settings
# ----------------------------------------
n <- c(200, 200)[SIMNUM]     		# number of samples
p <- c(1e3, 1e3)[SIMNUM]     		# number of features
s <- c(10, 10)[SIMNUM]       		# number of non-zero betas
cl <- c(0.25, 0.4)[SIMNUM]   		# censoring
r <- 100                     		# runs

S <- cor(tcga.rna[ , 1:p])
data_gen_settings <- list(
    # data generating processes
    dgp=c(data_gen_diag,               # Independent design
	  data_gen_diag,               # Diagonal design
	  data_gen_block,	       # Block design
	  data_gen_cov,                # Cov est from real data
	  data_gen_design),            # Psuedo response simulations 

    # parameters for data generating process
    par=list(0, 0.6, 0.6, S, tcga.rna)
)


# ----------------------------------------
# Run simulation
# ----------------------------------------
for (i in DGP) 
{
    # SVB (ours) 
    rname <- sprintf("c_%d_%d_%d", i, SIMNUM, 1)
    assign(rname,
	simulation_svb(n, p, s, cl,
	data_gen_settings$par[[i]], data_gen_settings$dgp[[i]], runs=r,
	inclusion_threshold=0.5, lambda=1, a0=1, b0=p, CORES=CORES))
    save(list=c(rname), file=sprintf("../../RData/comparison/%s.RData", rname))

    # MCMC
    rname <- sprintf("c_%d_%d_%d", i, SIMNUM, 2)
    assign(rname,
	simulation_mcmc(n, p, s, cl,
	data_gen_settings$par[[i]], data_gen_settings$dgp[[i]], runs=r,
	inclusion_threshold=0.5, lambda=1, a0=1, b0=p, CORES=CORES, 
	mcmc_samples=1e4))
    save(list=c(rname), file=sprintf("../../RData/comparison/%s.RData", rname))
}

