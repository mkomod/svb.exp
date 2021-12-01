library(survival)
library(survival.svb)

source("./00-functions.R")
source("../application/tcga_data.R", chdir=T)

d <- function(a, b) if(Sys.getenv(a) != "") as(Sys.getenv(a), class(b)) else b 
CORES <- d("CORES", 8)                 # number of CPU cores to use
SIMNUM <- d("SIMNUM", 1)               # which simulation parameters to use
DGP <- d("DGP", 1:5)                   # which data generating process to use


# ----------------------------------------
# Simulation settings
# ----------------------------------------
n <- c(5e2, 5e2, 1e3, 1e3)[SIMNUM]     # number of samples
p <- c(5e3, 5e3, 1e4, 1e4)[SIMNUM]     # number of features
s <- c(30, 30, 60, 60)[SIMNUM]         # number of non-zero betas
cl <- c(0.25, 0.4, 0.25, 0.4)[SIMNUM]  # censoring
r <- 100                               # runs
alpha <- 1                             # glmnet param

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

lambdas <- c(0.25, 0.5, 1, 2, 4, 8, 10, 20)
a0s <- c(1, 10, 25, 50, 100, 250, 500)
vals <- expand.grid(lambdas=lambdas, a0s=a0s)

# ----------------------------------------
# Run simulation
# ----------------------------------------
# Each method is ran for each setting in data_gen_settings
for (i in DGP) {
    
    # Sims 3 and 4 are too large for the TCGA data
    if (SIMNUM > 2 & i == 5) next

    rname <- sprintf("m%d.%d", i, SIMNUM)
    res <- c()
    for (v in 1:nrow(vals)) {
	lambda <- vals[v, 1]
	a0 <-  vals[v, 2]

	x <- simulation_svb(n, p, s, cl,
	    data_gen_settings$par[[i]], data_gen_settings$dgp[[i]], runs=r,
	    inclusion_threshold=0.5, lambda=lambda, a0=a0, b0=p, CORES=CORES,
	    alpha=alpha, eval_model=TRUE)

	x <- cbind(lambda=lambda, a0=a0, x)
	res <- rbind(res, x)
    }
    assign(rname, res)

    save(list=c(rname), file=sprintf("../../RData/sensitivity/%s.RData", rname))
}

