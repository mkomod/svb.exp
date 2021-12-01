# Comparison of MCMC to Variational Bayes
library(survival)
library(survival.svb)
library(survival.ss)

source("00-functions.R")
source("../application/tcga_data.R", chdir=T)


# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------
LOAD <- TRUE
n <- 200
p <- 1000
s <- 10
M <- 1e4
censoring <- 0.25

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
# Fit models
# ----------------------------------------
# for (i in 1:5)
for (i in 1) 
{
    d <- data_gen_settings$dgp[[i]](n, p, s, censoring, data_gen_settings$par[[i]])
    fname <- sprintf("../../RData/comparison/m%d.RData", i)

    if (LOAD && file.exists(fname)) {
	load(fname)
    } else {
	m <- survival.ss::run_sampler(d$Y, d$d, d$X, 1, mcmc_samples=M)
	m$m <- apply(m$b[ , 1e3:M], 1, mean)
	m$g <- apply(m$z[ , 1e3:M], 1, mean)

	mname <- sprintf("m%d", i)
	assign(mname, m)
	save(list=c(mname), file=fname)
    }

    s <- survival.svb::svb.fit(d$Y, d$d, d$X)

    assign(sprintf("s%d", i), s)
    assign(sprintf("d%d", i), d)
}

