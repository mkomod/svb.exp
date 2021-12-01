# Run simulations for various settings
#
library(survival)
library(BhGLM)
library(survival.svb)
library(BVSNLP)
library(glmnet)

source("./00-functions.R")
source("../application/tcga_data.R", chdir=T)

d <- function(a, b) if(Sys.getenv(a) != "") as(Sys.getenv(a), class(b)) else b 
CORES <- d("CORES", 8)                 # number of CPU cores to use
SIMNUM <- d("SIMNUM", 1)               # which simulation parameters to use
DGP <- d("DGP", 1:5)                   # which data generating process to use


# ----------------------------------------
# Simulation settings
# ----------------------------------------
n <- c(5e2, 5e2, 1e3, 1e3, 5e2, 5e2)[SIMNUM]      # number of samples
p <- c(5e3, 5e3, 1e4, 1e4, 1e4, 1e4)[SIMNUM]      # number of features
s <- c(30, 30, 60, 60, 60, 60)[SIMNUM]            # number of non-zero betas
cl <- c(0.25, 0.4, 0.25, 0.4, 0.25, 0.4)[SIMNUM]  # censoring
r <- 100                                          # runs
alpha <- 1                                        # glmnet param, 1: LASSO

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
# Each method is ran for each setting in data_gen_settings
for (i in DGP) {
    
    # Sims 3 and 4 are too large for the TCGA data
    if (n > nrow(tcga.rna) & i == 5) next

    # BgGLM
    rname <- sprintf("r_%d_%d_%d", i, SIMNUM, 1)
    assign(rname,
	simulation_bhglm(n, p, s, cl, 
	data_gen_settings$par[[i]], data_gen_settings$dgp[[i]], runs=r,
	inclusion_threshold=0.5, ss=c(0.04, 0.5), CORES=CORES))
    save(list=c(rname), file=sprintf("../../RData/simulations/%s.RData", rname))

    # BVSNLP
    rname <- sprintf("r_%d_%d_%d", i, SIMNUM, 2)
    assign(rname,
	simulation_bvsnlp(n, p, s, cl, 
	data_gen_settings$par[[i]], data_gen_settings$dgp[[i]], runs=r,
	inclusion_threshold=0.5, CORES=CORES))
    save(list=c(rname), file=sprintf("../../RData/simulations/%s.RData", rname))

    # SVB (ours) 
    rname <- sprintf("r_%d_%d_%d", i, SIMNUM, 3)
    assign(rname,
	simulation_svb(n, p, s, cl,
	data_gen_settings$par[[i]], data_gen_settings$dgp[[i]], runs=r,
	inclusion_threshold=0.5, lambda=1, a0=1, b0=p, CORES=CORES,
	alpha=alpha))
    save(list=c(rname), file=sprintf("../../RData/simulations/%s.RData", rname))

    # LASSO
    # rname <- sprintf("r_%d_%d_%d", i, SIMNUM, 3)
    # assign(rname,
	# simulation_lasso(n, p, s, cl, 
	# data_gen_settings$par[[i]], data_gen_settings$dgp[[i]], runs=r,
	# inclusion_threshold=0.5, CORES=CORES))
    # save(list=c(rname), file=sprintf("../../RData/simulations/%s.RData", rname))

}

