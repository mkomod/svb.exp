#
# Functions for fitting models to real world datasets
#
library(Rcpp)


# Fit models for a grid of lambda values
fit_svb <- function(Y, d, X, lambda, a0, b0, folds, 
    CORES=CORES, dname="", remove_ties=FALSE, maxiter=1500, seed=1)
{
    if (remove_ties) {
	set.seed(seed)
	Y <- jitter(Y)
    }

    # cross validation set
    set.seed(seed)
    ts.samples <- list(); tr.samples <- list()
    n <- length(Y)

    if (folds > 1) {
	k <- floor(n / folds)
	index <- sample(1:n)
	
	for (f in 1:folds) {
	    test <- index[((f-1)*k + 1):(f * k)]
	    train <- setdiff(index, test)

	    ts.samples[[f]] <- test
	    tr.samples[[f]] <- train
	}
    } else {
	ts.samples[[1]] <- 0
	tr.samples[[1]] <- 1:n
    }


    # setup cluster
    cl <- parallel::makeCluster(getOption("cl.cores", CORES), outfile="")
    parallel::clusterExport(cl, varlist=c("Y", "d", "X", "a0", "b0", 
	    "dname", "folds", "tr.samples", "ts.samples"), 
	    envir=environment())
    on.exit(parallel::stopCluster(cl))
  

    # start fitting
    parallel::parSapplyLB(cl, lambda, function(l) 
    {
	fname <- sprintf("../../RData/models/%s/%s_l_%s.RData", 
	    dname, dname, l)
	fit <- list()
	
	# initialise the fit
	for (i in 1:folds) 
	{
	    tr <- tr.samples[[i]]
	    ts <- ts.samples[[i]]
	    
	    # init model
	    f <- survival.svb::svb.fit(Y[tr], d[tr], X[tr, ], 
		lambda=l, a0=a0, b0=b0, maxiter=0)

	    # track the ELBO
	    tr.elbo.mean <- c(); tr.elbo.sd <- c(); tr.elbo.kl <- c()
	    ts.elbo.mean <- c(); ts.elbo.sd <- c(); ts.elbo.kl <- c()
	    iter <- 0
	    
	    while (!f$converged && (iter < maxiter)) {
		# fit 5 iterations of model
		f <- survival.svb::svb.fit(Y[tr], d[tr], X[tr, ], 
		    lambda=f$lambda, a0=f$a0, b0=f$b0, 
		    mu.init=f$m, s.init=f$s, g.init=f$g, maxiter=5)
		
		# compute the ELBO
		tr.elbo <- survival.svb::elbo(Y[tr], d[tr], X[tr, ], f, nrep=5e2)
		tr.elbo.mean <- c(tr.elbo.mean, tr.elbo$mean)
		tr.elbo.sd <-   c(tr.elbo.sd,   tr.elbo$sd)
		tr.elbo.kl <-   c(tr.elbo.kl,   tr.elbo$kl)

		ts.elbo <- survival.svb::elbo(Y[ts], d[ts], X[ts, ], f, nrep=5e2)
		ts.elbo.mean <- c(ts.elbo.mean, ts.elbo$mean)
		ts.elbo.sd <-   c(ts.elbo.sd,   ts.elbo$sd)
		ts.elbo.kl <-   c(ts.elbo.kl,   ts.elbo$kl)
		
		cat("lambda: ", l, "\n", "ELBO", tr.elbo.mean, "\n")
		iter <- iter + 5
	    }
	    
	    # save ELBO to f
	    f$tr.elbo.mean <- tr.elbo.mean
	    f$tr.elbo.kl <-   tr.elbo.kl
	    f$tr.elbo.ll <-   tr.elbo.mean + tr.elbo.kl
	    f$tr.elbo.sd <-   tr.elbo.sd

	    f$ts.elbo.mean <- ts.elbo.mean
	    f$ts.elbo.kl <-   ts.elbo.kl
	    f$ts.elbo.ll <-   ts.elbo.mean + ts.elbo.kl
	    f$ts.elbo.sd <-   ts.elbo.sd
	    
	    fit[[i]] <- f
	}

	if (folds == 1) {
	    fit <- f
	} else {
	    fit$tr.samples <- tr.samples
	    fit$ts.samples <- ts.samples
	    fit$folds <- folds
	}

	save(fit, file=fname)
	return(0)
    })

    invisible(NULL)
}


ci <- function(l, m, s, width=2, ...) {
    arrows(l, m + width * s, l, m - width * s, code=0, ...)
}


model_summary <- function(Y, delta, X, fit) 
{
    # get x-val samples
    tr.samples <- fit$tr.samples
    ts.samples <- fit$ts.samples

    res <- sapply(1:fit$folds, function(i) {
	f <- fit[[i]]
	f$beta_hat <- f$m * f$g

	tr <- tr.samples[[i]]
	ts <- ts.samples[[i]]

	iter <- length(f$tr.elbo.mean)

	# training
	tr.elbo.mean <- f$tr.elbo.mean[iter]
	tr.elbo.kl <-   f$tr.elbo.kl[iter]
	tr.elbo.ll <-   f$tr.elbo.ll[iter]
	tr.elbo.sd <-   f$tr.elbo.sd[iter]  	
	tr.cin <- cindex(Y[tr], delta[tr], X[tr, ] %*% matrix(f$beta_hat))$cindex
	
	nvars <- sum(f$g > 0.5)

	# validation
	ts.elbo.mean <- f$ts.elbo.mean[iter]
	ts.elbo.kl <- 	f$ts.elbo.kl[iter]
	ts.elbo.ll <- 	f$ts.elbo.ll[iter]
	ts.elbo.sd <- 	f$ts.elbo.sd[iter]  	
	ts.cin <- cindex(Y[ts], delta[ts], X[ts, ] %*% matrix(f$beta_hat))$cindex
	
	return(c(
	    tr.elbo.mean=tr.elbo.mean, tr.elbo.kl=tr.elbo.kl, 
	    tr.elbo.ll=tr.elbo.ll, tr.elbo.sd=tr.elbo.sd, tr.cin=tr.cin,
	    nvars=nvars,
	    ts.elbo.mean=ts.elbo.mean, ts.elbo.kl=ts.elbo.kl, 
	    ts.elbo.ll=ts.elbo.ll, ts.elbo.sd=ts.elbo.sd, ts.cin=ts.cin
	))
    })

    return(res)
}


# compute the concordance
Rcpp::sourceCpp("concordance.cpp")


plot_fit <- function(fit, val, col="lightgrey", lwd=2, 
    add_mean=TRUE, m.lwd=2, m.col=2, negate=FALSE, ...)
{
    l <- list()
    m <- c()
    lwr <- c()
    upr <- c()

    for (i in 1:fit$folds) {
	f <- fit[[i]]
	x <- f[[val]] * ifelse(negate, -1, 1)
	l[[i]] <- x

	m <- c(m, length(x))
	lwr <- c(lwr, min(x))
	upr <- c(upr, max(x))
    }

    k <- which.max(m)
    xs <- l[[k]]
    index <- 1:max(m)*5

	
    plot(index, l[[k]], ylim=range(lwr, upr),
	 type="l", lwd=lwd, col=col, ...)

    for (i in setdiff(1:fit$folds, k)) {
	x <- l[[i]]
	newx <- setdiff(1:max(m), 1:length(x))
	x[newx] <- x[length(x)]

	lines(index, x, col=col, lwd=lwd)

	if (add_mean) 
	    xs <- rbind(xs, x)
    }
    
    if (add_mean) {
	lines(index, apply(xs, 2, mean), col=m.col, lwd=m.lwd)
    }
}


ci.vi <- function(fit, b0, a=0.95, size=3, offset=0.5, add_line=T,
	lty.col=1, lty.wd=2, side="both", ...) 
{
    m <- fit$m
    s <- fit$s
    a <- 1 - a
    width <- qnorm(1 - a/2)
    for (i in which(b0 != 0)) {
	xs <- seq(m[i] - s[i] * width, m[i] + s[i] * width, by=0.005)
	ys <- dnorm(xs, m[i], s[i])
	if (side == "right" || side == "both") {
	    polygon(i + ys * size - min(ys) * size + offset, xs, ...) 
	    if (add_line)
		lines(c(i, i) + offset, range(xs), col=lty.col, lwd=lty.wd)
	} 
	if (side == "left" || side == "both") {
	    polygon(i - ys * size + min(ys) * size - offset, xs, ...) 
	    if (add_line)
		lines(c(i, i) - offset, range(xs), col=lty.col, lwd=lty.wd)
	}
    }
}

