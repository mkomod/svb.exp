# Functions for setting up simulations

# ----------------------------------------
# Methods
# ----------------------------------------
#  - BhGLM
#  - SVB (ours)
#  - BVSNLP
#  - LASSO
#  - MCMC spike-and-slab
simulation_bhglm <- function(n, p, s, censoring, corr, data_gen, runs, 
    inclusion_threshold, ss=c(0.04, 0.5), b=NULL, CORES=CORES)
{

    cl <- parallel::makeCluster(getOption("cl.cores", CORES), outfile="")
    parallel::clusterExport(cl,
	varlist=c("n", "p", "s", "censoring", "corr", "data_gen", "ss", "b",
		  "inclusion_threshold", "method_summary", "method_roc",
		  "is.Surv"), 
	envir=environment())
    on.exit(parallel::stopCluster(cl))


    res <- parallel::parSapply(cl, 1:runs, function(run) {
	cat(run)
	d <- data_gen(n, p, s, censoring, corr, seed=run, b=b)
	y <- survival::Surv(as.matrix(d$Y), as.matrix(as.numeric(d$d)))

	fit.time <- system.time({
	    fit <- BhGLM::bmlasso(d$X, y, family="cox", ss=ss)
	})	

	r1 <- method_summary(d$beta, fit$beta, fit$p, inclusion_threshold)
	r2 <- method_roc(d$beta, fit$p)
	
	c(unlist(r1), unlist(fit.time[3]), auc=r2$auc, tpr=r2$TPR, fpr=r2$FPR)
    })

    return(t(res))
}


simulation_svb <- function(n, p, s, censoring, corr, data_gen, runs,
    inclusion_threshold, lambda, a0, b0, b=NULL, CORES=CORES, alpha=1,
    eval_model=FALSE)
{
    # create a cluster and export variables to each
    cl <- parallel::makeCluster(getOption("cl.cores", CORES), outfile="")
    parallel::clusterExport(cl,
	varlist=c("n", "p", "s", "censoring", "corr", "data_gen", "lambda",
		  "inclusion_threshold", "a0", "b0", "method_summary", 
		  "concordance.index", "method_coverage", "b", "method_roc", 
		  "alpha", "eval_model", "svb.credible_interval"), 
	envir=environment())
    on.exit(parallel::stopCluster(cl))


    res <- parallel::parSapply(cl, 1:runs, function(run) {
	# generate data, seed set to run number
	d <- data_gen(n, p, s, censoring, corr, seed=run, b=b)

	fit.time <- system.time({
	fit <- survival.svb::svb.fit(d$Y, d$d, d$X, 
	    lambda=lambda, a0=a0, b0=b0,
	    maxiter=1500, tol=1e-3, alpha=alpha)
	})

	r1 <- method_summary(d$beta, fit$mu*fit$gamma, fit$gamma, 
		inclusion_threshold)
	r2 <- method_roc(d$beta, fit$gamma)
	r3 <- method_coverage(fit, d)

	if (eval_model) {
	    m.elbo <- survival.svb::elbo(d$Y, d$d, d$X, fit)
	    c.index <- concordance.index(d$Y, d$d, d$X, fit$m * fit$g)

	    return(c(unlist(r1), unlist(fit.time[3]), auc=r2$auc, tpr=r2$TPR, fpr=r2$FPR,
	      unlist(r3), unlist(m.elbo), cindex=c.index))
	}
		
	return(c(unlist(r1), unlist(fit.time[3]), auc=r2$auc, tpr=r2$TPR, fpr=r2$FPR,
	  unlist(r3)))
    })

    return(t(res))
}


simulation_bvsnlp <- function(n, p, s, censoring, corr, data_gen, runs,
    inclusion_threshold, b=NULL, CORES=CORES)
{
    # create a cluster and export variables to each
    cl <- parallel::makeCluster(getOption("cl.cores", CORES), outfile="")
    parallel::clusterExport(cl,
	varlist=c("n", "p", "s", "censoring", "corr", "data_gen",
		  "inclusion_threshold", "method_summary", "b", "method_roc"), 
	envir=environment())
    on.exit(parallel::stopCluster(cl))
    
    res <- parallel::parSapply(cl, 1:runs, function(run) {
	# generate data
	d <- data_gen(n, p, s, censoring, corr, seed=run, b=b)
	resp <- matrix(c(d$Y, d$d), ncol=2)

	fit.time <- system.time({
	    fit <- BVSNLP::bvs(data.frame(d$X), resp, family="survival", ncpu=1)
	})

	beta.hat <- rep(0, p)
	if (length(fit$MPM) != 0) {
	    # prevents errors
	    beta.hat[fit$MPM] <- fit$beta_hat
	}

	r1 <- method_summary(d$beta, beta.hat, fit$inc_probs, 
		inclusion_threshold)
	r2 <- method_roc(d$beta, fit$inc_probs)
	
	c(unlist(r1), unlist(fit.time[3]), auc=r2$auc, tpr=r2$TPR, fpr=r2$FPR)
    })

    return(t(res))
}


simulation_lasso <- function(n, p, s, censoring, corr, data_gen, runs,
    inclusion_threshold, b=NULL, CORES=CORES)
{
    # create a cluster and export variables to each
    cl <- parallel::makeCluster(getOption("cl.cores", CORES), outfile="")
    parallel::clusterExport(cl,
	varlist=c("n", "p", "s", "censoring", "corr", "data_gen",
		  "inclusion_threshold", "method_summary", "b", "method_roc"), 
	envir=environment())
    on.exit(parallel::stopCluster(cl))
    
    res <- parallel::parSapply(cl, 1:runs, function(run) {
	# generate data
	d <- data_gen(n, p, s, censoring, corr, seed=run, b=b)

	y <- survival::Surv(matrix(d$Y), matrix(d$d + 0)) 

	fit.time <- system.time({
	    fit <- glmnet::cv.glmnet(d$X, y, family="cox")
	    beta.hat <- fit$glmnet.fit$beta[ , fit$index[1]]
	})

	r1 <- method_summary(d$beta, beta.hat, beta.hat != 0, inclusion_threshold)
	r2 <- method_roc(d$beta, beta.hat != 0)
	
	c(unlist(r1), unlist(fit.time[3]), auc=r2$auc, tpr=r2$TPR, fpr=r2$FPR)
    })

    return(t(res))
}


simulation_mcmc <- function(n, p, s, censoring, corr, data_gen, runs,
    inclusion_threshold, lambda, a0, b0, b=NULL, CORES=CORES, mcmc_samples=1e4)
{
    # create a cluster and export variables to each
    cl <- parallel::makeCluster(getOption("cl.cores", CORES), outfile="")
    parallel::clusterExport(cl,
	varlist=c("n", "p", "s", "censoring", "corr", "data_gen", "lambda",
		  "inclusion_threshold", "a0", "b0", "method_summary", 
		  "mcmc_samples", "method_roc", "method_coverage", 
		  "mcmc.credible_interval"), 
	envir=environment())
    on.exit(parallel::stopCluster(cl))

    res <- parallel::parSapply(cl, 1:runs, function(run) {
	# generate data, seed set to run number
	d <- data_gen(n, p, s, censoring, corr, seed=run, b=b)

	fit.time <- system.time({
	fit <- survival.ss::run_sampler(d$Y, d$d, d$X,
	    lambda=lambda, a_0=a0, b_0=b0, mcmc_samples=mcmc_samples, 
	    verbose=F)
	})

	fit$m <- apply(fit$b[ , 1e3:mcmc_samples], 1, mean)
	fit$g <- apply(fit$z[ , 1e3:mcmc_samples], 1, mean)

	r1 <- method_summary(d$beta, fit$m*fit$g, fit$g, inclusion_threshold)
	r2 <- method_roc(d$beta, fit$g)
	r3 <- method_coverage(fit, d, mcmc=T)

	return(c(unlist(r1), unlist(fit.time[3]), auc=r2$auc, tpr=r2$TPR, 
	    fpr=r2$FPR, unlist(r3)))
    })

    return(t(res))
}


# ----------------------------------------
# Simulation data generation
# ----------------------------------------
# Covariance matrices
# - block: where blocks of observations are correlated
# - diag: where there is correlation along the diagonal of the covariance matrix
# - real: a real design is used


# Generate data where observations within blocks are correlated
# and observations between blocks are independent
data_gen_block <- function(n, p, s, censoring, corr, 
    omega=1, b=NULL, seed=1, block_size=50) 
{
    set.seed(seed)

    if ((p %% block_size) != 0) stop("'block_size' must be multiple of p")

    # generate true coefs is not given
    if (is.null(b))
	b <- sample(c(sample(c(-1, 1), s, replace=T) * runif(s, 0.5, 2.0),
		rep(0, p - s)))

    if (corr == 0) {
	X <- matrix(rnorm(n * p), nrow=n)
    } else {
	X <- matrix(nrow=n, ncol=0)
	for (block in 1:(p / block_size)) {
	    S <- matrix(corr, nrow=block_size, ncol=block_size)
	    diag(S) <- 1
	    X <- cbind(X, mvtnorm::rmvnorm(n, mean=rep(0, block_size), S))
	}
    }

    Y <- log(1 - runif(n)) / - (exp(X %*% b) * omega)
    d <- runif(n) > censoring
    Y[!d] <- Y[!d] * runif(sum(!d)) 

    return(list(Y=Y, d=d, X=X, beta=b))
}


# Generate data where the diagonal of the design matrix is correlated
# setting corr=0 gives uncorrelated design
data_gen_diag <- function(n, p, s, censoring, corr, 
    omega=1, b=NULL, seed=1)
{
    set.seed(seed)
    
    # generate a random beta if one is not supplied
    if (is.null(b))
	b <- sample(c(sample(c(-1, 1), s, replace=T) * runif(s, 0.5, 2.0),
		rep(0, p - s)))

    if (corr == 0) {
	X <- matrix(rnorm(n * p), nrow=n)
    } else {
	S <- outer(1:p, 1:p, function(i, j) corr^(abs(i - j)))
	X <- mvtnorm::rmvnorm(n, rep(0, p), sigma=S)
    }

    Y <- log(1 - runif(n)) / - (exp(X %*% b) * omega)
    d  <- runif(n) > censoring
    Y[!d] <- Y[!d] * runif(sum(!d)) 

    return(list(Y=Y, d=d, X=X, n=n, beta=b))
}


# Generate data for a given covariance matrix
data_gen_cov <- function(n, p, s, censoring, corr, 
    omega=1, b=NULL, seed=1)
{
    set.seed(seed)
    
    # generate a random beta if one is not supplied
    if (is.null(b))
	b <- sample(c(sample(c(-1, 1), s, replace=T) * runif(s, 0.5, 2.0),
		rep(0, p - s)))

    X <- mvtnorm::rmvnorm(n, rep(0, p), sigma=corr)
    Y <- log(1 - runif(n)) / - (exp(X %*% b) * omega)
    d  <- runif(n) > censoring
    Y[!d] <- Y[!d] * runif(sum(!d)) 

    return(list(Y=Y, d=d, X=X, n=n, beta=b))
}


# Generate data from a design matrix
data_gen_design <- function(n, p, s, censoring, X, omega=1, 
    b=NULL, seed=1)
{
    set.seed(seed)

    if (is.null(b))
	b <- sample(c(sample(c(-1, 1), s, replace=T) * runif(s, 0.5, 2.0),
		rep(0, p - s)))
    
    X <- X[sample(1:nrow(X), n), sample(1:ncol(X), p)]
    Y <- log(1 - runif(n)) / - (exp(X %*% b) * omega)
    d <- runif(n) > censoring
    Y[!d] <- Y[!d] * runif(sum(!d)) 

    return(list(Y=Y, d=d, X=X, beta=b))
}


# ----------------------------------------
# Method evaluation
# ----------------------------------------

# summarise the method output
method_summary <- function(beta.true, beta.hat, inclusion.prob, thresh) 
{
    g <- beta.true != 0
    tab <- table(inclusion.prob > thresh, g)

    if (nrow(tab) == 2) {
	TP <- tab[2,2]; FP <- tab[2,1]
	FN <- tab[1,2]; TN <- tab[1,1]
    } else if(all(inclusion.prob > thresh)) {
	TP <- sum(g);   FP <- sum(1 - g)
	FN <- 0;        TN <- 0
    } else {
	TP <- 0;        FP <- 0
	FN <- sum(g);   TN <- sum(1 - g)
    }

    return(list(
	 acc = (TN + TP) / (TN + TP + FN + FP), 
	 err = (FP + FN) / (TN + TP + FN + FP), 
	 tpr = TP / (TP + FN), 
	 tnr = TN / (TN + FP), 
	 fpr = FP / (TN + FP),
	 ppv = ifelse(TP + FP == 0, NA, TP / (TP + FP)),
	 fdr = ifelse(TP + FP == 0, NA, FP / (TP + FP)),
	 l1 = sum(abs(beta.true - beta.hat)),
	 l2 = sum((beta.true - beta.hat)^2)
    ))
}


# generate the ROC curve and AUC statistic
method_roc <- function(beta.true, inclusion.prob)
{
    g <- beta.true != 0
    tpr <- numeric(0)
    fpr <- numeric(0)

    for (thresh in seq(0, 1, by=0.01)) {
	tab <- table(inclusion.prob >= thresh, g)
	if (nrow(tab) == 2) {
	    TP <- tab[2,2]; FP <- tab[2,1]
	    FN <- tab[1,2]; TN <- tab[1,1]
	} else if(all(inclusion.prob >= thresh)) {
	    TP <- sum(g);   FP <- sum(1 - g)
	    FN <- 0; 	    TN <- 0
	} else {
	    TP <- 0;        FP <- 0
	    FN <- sum(g);   TN <- sum(1 - g)
	}
	tpr <- c(tpr, TP / (TP + FN))
	fpr <- c(fpr, FP / (TN + FP))
    }
    
    # compute the AUC using trapizium int
    auc <- -sum((tpr[1:100] + tpr[2:101]) / 2 * diff(fpr))

    return(list(auc=auc, FPR=fpr, TPR=tpr))
}


# Compute the coverage
method_coverage <- function(fit, d, a=0.05, threshold=0.5, conditional=FALSE, 
	mcmc=FALSE)
{
    if (mcmc) {
	ci <- mcmc.credible_interval(fit, a, conditional=conditional)
    } else {
	ci <- svb.credible_interval(fit, a, conditional=conditional)
    }
    l <- ci[ , 1]
    u <- ci[ , 2]
    dirac <- ci[ , 3]

    coverage <- (((l <= d$beta) & (d$beta <= u)) | ((d$beta == 0) & (dirac == T)))
    bs <- d$beta != 0

    return(list(
	coverage.n0 = mean(coverage[bs]),
	length.n0 = mean(u[bs] - l[bs]),
	coverage.0 = mean(coverage[!bs]),
	length.0 =  mean(u[!bs] - l[!bs])
    ))
}


# Compute the c-index (slow)
concordance.index <- function(y, d, X, b)
{
    b <- as.matrix(b)
    num <- 0
    den <- 0
    
    for (i in 2:nrow(X)) {
	x.i <- X[i, ]
	y.i <- y[i]
	d.i <- d[i]
	for (j in 1:(i-1)) {
	    x.j <- X[j, ]
	    y.j <- y[j]
	    d.j <- d[j]
	    x.ij <- t(b) %*% x.i - t(b) %*% x.j
	    num <- num +
		(y.i < y.j) * (x.ij > 0) * d.i + (y.j < y.i) * (x.ij < 0) * d.j
	    den <- den +
		(y.i < y.j) * d.i + (y.j < y.i) * d.j
	}
    }
    return(num / den)
}


# compute the credible intervals
svb.credible_interval <- function(fit, a=0.05, conditional=FALSE)
{
    credible.interval <- sapply(1:length(fit$g), function(i) {
	g <- fit$g[i]
	m <- fit$m[i]
	s <- fit$s[i]
	
	if (conditional) {
	    return(qnorm(c(a/2, 1-a/2), m, s))
	}

	if (!conditional) {
	    if (g > 1 - a) {
		# compute the interval that contains 1 - a.g of the mass
		# i.e. if g = 0.97 then the interval needs to be wider to
		# contain 95% of the total mass
		a.g <- 1 - (1-a)/g
		interval <- qnorm(c(a.g/2, 1-a.g/2), m, s)
		contains.dirac <- FALSE
		
		if (interval[1] <= 0 && interval[2] >= 0) {
		    # if interval contains Dirac mass it needs to be smaller
		    interval <- qnorm(c(a.g/2 + (1-g)/2, 1 - a.g/2 - (1-g)/2), m, s)
		    contins.dirac <- TRUE
		}

		return(c(lower=interval[1], upper=interval[2], 
			 contains.dirac=contains.dirac))
	    } else if (g < a) {
		# if the spike contains (1-a)% of the mass then we take 
		# the Dirac mass at 0
		return(c(lower=0, upper=0, contains.dirac=T))
	    } else {
		# will always contain the Dirac
		# so we remove the density accounted for by the
		# Dirac mass from the interval
		interval <- qnorm(c(a/2 + (1-g)/2, 1-a/2 -(1-g)/2), m, s)

		return(c(lower=interval[1], upper=interval[2], 
			contains.dirac=TRUE))
	    }
	}
    })

    return(t(credible.interval))
}


# compute the credible intervals
mcmc.credible_interval <- function(fit, a=0.05, conditional=FALSE, burnin=1e3)
{
    p <- ncol(fit$b)

    credible.interval <- sapply(1:length(fit$g), function(i) {
	g <- fit$g[i]	

	if (conditional) {
	    m <- fit$b[i, burnin:p]
	    z <- !!fit$z[i, burnin:p]
	    m <- m[z]
	    f <- density(m)
	    cdf <- ecdf(m)
	    return(quantile(cdf, c(a/2, 1 - a/2)))
	}

	if (!conditional) {
	    if (g > 1 - a) {
		# Slab contains 1 - a of mass
		a.g <- 1 - (1-a)/g
		m <- fit$b[i, burnin:p]
		z <- !!fit$z[i, burnin:p]
		m <- m[z]
		f <- density(m)
		cdf <- ecdf(m)
		interval <- quantile(cdf, c(a.g/2, 1 - a.g/2))
		contains.dirac <- FALSE
		
		if (interval[1] <= 0 && interval[2] >= 0) {
		    # if interval contains Dirac mass it needs to be smaller
		    interval <- quantile(cdf, c(a.g/2+(1-g)/2, 1-a.g/2-(1-g)/2))
		    contins.dirac <- TRUE
		}

		return(c(lower=interval[1], upper=interval[2], 
			 contains.dirac=contains.dirac))
	    } else if (g < a) {
		# Dirac contains 1 - a of mass
		return(c(lower=0, upper=0, contains.dirac=T))
	    } else {
		# will always contain the Dirac
		m <- fit$b[i, burnin:p]
		z <- !!fit$z[i, burnin:p]
		m <- m[z]
		f <- density(m)
		cdf <- ecdf(m)
		interval <- quantile(cdf, c(a/2 +(1-g)/2, 1 - a/2 - (1-g)/2))

		return(c(lower=interval[1], upper=interval[2], 
			 contains.dirac=TRUE))
	    }
	}
    })

    return(t(credible.interval))
}

