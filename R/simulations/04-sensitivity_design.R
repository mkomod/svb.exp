# Explore the effect of varaince in the design and the method's
# performance
#
library(survival)
library(survival.svb)
library(BhGLM)
library(BVSNLP)

source("00-functions.R")


# ----------------------------------------
# Settings
# ----------------------------------------
set.seed(1)

reps <- 10
n <- 200                   # number of samples
p <- 1000                  # number of parameters
s <- 10                    # number of non-zero betas
clvl <- 0.4                # censoring level

S <- sample(1:p, s)        # indices of non-zero beta
S <- sort(S)
Sigma <- rep(1, p)
b <- rep(0, p)
b.vals <- seq(-3, 3, by=0.25)[-13]  # get rid of 0


# ----------------------------------------
# Design varaince and method performance
# ----------------------------------------
d.1 <- d.2 <- d.3 <- matrix(0, ncol=0, nrow=length(b.vals))
for (i in 1:3) {
    Sigma[S] <- seq(0.05 + (0.5*i-0.5) , 0.5*i, by=0.05)
    g.1 <- g.2 <- g.3 <- matrix(0, ncol=10, nrow=length(b.vals))

    for (seed in 1:reps) {
	for (i in seq_along(b.vals)) {
	    b[S] <- b.vals[i]
	    
	    # generate the dataset
	    d <- data_gen_cov(n, p, s, clvl, diag(Sigma), b=b, seed=seed)
	    y <- Surv(as.matrix(d$Y), as.matrix(0+d$d))
	    resp <- matrix(c(d$Y, d$d), ncol=2)
	    
	    # fit the different models
	    f.1 <- BhGLM::bmlasso(d$X, y, family="cox")
	    f.2 <- BVSNLP::bvs(data.frame(d$X), resp, family="survival", ncpu=1)
	    f.3 <- survival.svb::svb.fit(d$Y, d$d, d$X, verbose=F, alpha=1)
	    
	    g.1[i, ] <- g.1[i, ] + (S %in% which(f.1$p > 0.5))
	    g.2[i, ] <- g.2[i, ] + (S %in% which(f.2$inc_probs > 0.5))
	    g.3[i, ] <- g.3[i, ] + (S %in% which(f.3$g > 0.5))
	}
    }
    d.1 <- cbind(d.1, g.1)
    d.2 <- cbind(d.2, g.2)
    d.3 <- cbind(d.3, g.3)
}

save(d.1, file="../../RData/sensitivity/d.1.RData")
save(d.2, file="../../RData/sensitivity/d.2.RData")
save(d.3, file="../../RData/sensitivity/d.3.RData")
