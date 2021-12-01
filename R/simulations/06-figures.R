# Create plots of simulations
#
library(lattice)
library(latex2exp)

source("03-mcmc_models.R", chdir=T)
source("../application/tcga_data.R", chdir=T)


# ----------------------------------------
# Functions
# ----------------------------------------
ci.vi <- function(fit, b0, a=0.95, size=3, offset=0.5, add_line=T,
	lty.col=1, lty.wd=2, ...) 
{
    m <- fit$m
    s <- fit$s
    g <- fit$g

    for (i in which(b0 != 0)) 
    {
	if (g[i] <= a) next

	a.g <- 1 - a/g[i]
	width <- qnorm(1 - a.g/2)

	xs <- seq(m[i] - s[i] * width, m[i] + s[i] * width, by=0.005)
	ys <- dnorm(xs, m[i], s[i])
	polygon(i + ys * size - min(ys) * size + offset, xs, ...) 
	if (add_line) {
	    lines(c(i, i) + offset, range(xs), col=lty.col, lwd=lty.wd)
	}
    }
}


ci.mcmc <- function(fit, b0, a=0.95, size=3, offset=0.5, add_line=T,
	lty.col=1, lty.wd=2, ...) 
{
    a <- 1 - a
    for (i in which(b0 != 0)) {
	m <- fit$b[i, !!fit$z[i, ]]
	f <- density(m)
	cdf <- ecdf(m)
	qs <- quantile(cdf, c(a/2, 1 - a/2))

	xs <- f$x[f$x > qs[1] & f$x < qs[2]]
	ys <- f$y[f$x > qs[1] & f$x < qs[2]]

	polygon(i - ys * size + min(ys) * size - offset, xs, ...) 

	if (add_line) {
	    lines(c(i, i) - offset, range(xs), col=lty.col, lwd=lty.wd)
	}
    }
}


# ----------------------------------------
# 01-mcmc_comp
# ----------------------------------------
# for (comp in c(4:5)) 
for (comp in c(1)) 
{
    d <- get(sprintf("d%d", comp))
    s <- get(sprintf("s%d", comp))
    m <- get(sprintf("m%d", comp))
    m.m <- get(sprintf("m%d", comp))$m
    m.g <- get(sprintf("m%d", comp))$g

    bs <- which(d$b != 0)

    pdf(file=sprintf("../../figures/comparison_0%d_beta.pdf", comp), 
	width=6, height=4)
	par(family="Times", mar=c(4, 2, 0.2, 0.2))
	plot(d$b, pch=8, ylim=c(-2.5, max(m.m * m.g) + 0.5), 
	     ylab=expression(beta), las=1)
	points(m.g * m.m, col="darkorange")
	points(s$m * s$g, col="darkblue", pch=20)
	ci.vi(s, s$g > 0.95, size=15, col=rgb(1, 0.549, 0, 0.2), border=NA, 
	      lty.col=2, offset=3, lty.wd=3)
	ci.mcmc(m, m.g > 0.95, size=15, col=rgb(0, 0, 0.545, 0.2), border=NA, 
		lty.col="blue", offset=3)
	legend("topleft", legend=c(TeX("$\\hat{\\beta}$ MCMC$"), 
	       TeX("$\\hat{\\beta}$ VA"), TeX("$\\beta_0$")), bty="n",
	       pch=c(1, 20, 8), col=c("darkblue", "darkorange", "black"))
    dev.off()

    plot(m.g, col="darkblue", lwd=2, ylab="Inclusion Probability")
    points(s$g, col="darkorange", pch=20, cex=0.8)
    points(bs, rep(1.035, 10), pch=25, bg="darkorchid1", col=0)
    legend("topleft", legend=c("MCMC", "VA", "True locations"), pch=c(1, 20, 25),
	   col=c("darkblue", "darkorange", "darkorchid1"), pt.bg=c(NA, NA, "darkorchid1"),
	   bty="n")
}


# ----------------------------------------
# 04-design-variance
# ----------------------------------------
cols <- colorRampPalette(c(rgb(.95,.95,.95,.1), c(rgb(0.6,0,0,0.9))))

for (i in 1:3) {
    load(sprintf("../../RData/sensitivity/d.%d.RData", i))
    z <- get(sprintf("d.%d", i))
    z <- rbind(z[1:12, ], 0, z[13:24, ])
    z <- 1 - z / 10

    pdf(file=sprintf("../../figures/design_var-%d.pdf", i), width=9, height=6)
    par(family="Times", mar=c(4, 4, 1, 1))
    p1 <- levelplot(t(z),
	    xlab=expression(sigma^2), ylab=expression(beta), 
	    col.regions = cols, useRaster=TRUE,
	    scales=list(
		x=list(at=seq(0, ncol(z), length.out=11),
		labels=round(seq(0.0, 1.5, length.out=11), 2)),
		y=list(at=seq(1, nrow(z), length.out=13),
		labels=round(seq(-3, 3, length.out=13), 2))
		))
    update(p1, aspect=0.65)
    print(p1)
    dev.off()
}


# Differences
x <- d.3
x <- rbind(x[1:12, ], 0, x[13:24, ])
x <- 1 - x / 10

for (i in 1:2) {
    load(sprintf("../../RData/sensitivity/d.%d.RData", i))
    y <- get(sprintf("d.%d", i))
    y <- rbind(y[1:12, ], 0, y[13:24, ])
    y <- 1 - y / 10
    z <- x - y

    pdf(file=sprintf("../../figures/design_var_diff-3_%d.pdf", i), width=9, height=6)
    par(family="Times", mar=c(4, 4, 1, 1))
    p1 <- levelplot(t(z),
	    xlab=expression(sigma^2), ylab=expression(beta), 
	    at = seq(-1, 1, by=0.1),
	    col.regions = hcl.colors(41, palette="Blue-Red"),
	    useRaster=TRUE,
	    scales=list(
		x=list(at=seq(0, ncol(z), length.out=11),
		labels=round(seq(0.0, 1.5, length.out=11), 2)),
		y=list(at=seq(1, nrow(z), length.out=13),
		labels=round(seq(-3, 3, length.out=13), 2))
		))
    update(p1, aspect=0.65)
    print(p1)
    dev.off()
}


# ----------------------------------------
# 05-sensitivity :: lambda
# ----------------------------------------
pdf(file="../../figures/sensitivity_lambda.pdf", width=9, height=9)
par(family="Times", mar=c(4, 4, 3, 1))

metrics <- c("l2", "l1", "tpr", "fdr", "auc", "elapsed", 
	     "mean", "expected.likelihood", "cindex")
metrics.lab <- c(TeX("$l_2^2$-error"), TeX("$l_1$-error"), "TPR", "FDR", 
		 "AUC", "Time (s)", "ELBO", "ELL", "c-index")
layout(matrix(1:9, nrow=3, byrow=T))
for (s in 1) {
    for (m in seq_along(metrics)) {
	xs <- c()
	for (i in 1:5) {
	    load(file=sprintf("../../RData/sensitivity/m%d.%d.RData", s, i))
	    x <- get(sprintf("m%d.%d", s, i))
	    a0 <- x[ , "a0"] == 1
	    lambdas <- x[ , "lambda"][a0]
	    x.mean <- as.numeric(by(x[a0, metrics[m]], lambdas, mean, simplify=T))
	    xs <- rbind(xs, x.mean)
	}
	lambdas <- unique(lambdas)
	matplot(lambdas, t(xs), log="x", type="l", lwd=2, col=2:6,
	    xlab=expression(lambda), ylab=metrics.lab[m], 
	    main=metrics.lab[m], las=1)
	mtext(sprintf("(%s)", letters[m]), 3, line=0.5, adj=-0.15)

	if (any(m == c(3, 5, 7:9))) {
	    for (r in 1:nrow(xs)) {
		maxs <- which(xs[r, ] == max(xs[r, ]))
		points(lambdas[maxs], xs[r, maxs], pch=20+r, col=1,
		       bg=r+1, cex=1.3)
	    }
	}
	if (any(m == c(1:2, 4))) {
	    for (r in 1:nrow(xs)) {
		mins <- which(xs[r, ] == min(xs[r, ]))
		points(lambdas[mins], xs[r, mins], pch=20+r, 
		       col=1, bg=r+1, cex=1.3)
	    }
	}

	if (m == 1)
	    legend("topleft", legend=sprintf("Setting %d", 1:5), 
		   col=2:6, lwd=2, lty=1:5, bty="n")
    }
}
dev.off()


# ----------------------------------------
# 05-sensitivity :: alpha
# ----------------------------------------
pdf(file="../../figures/sensitivity_alpha.pdf", width=9, height=9)
par(family="Times", mar=c(4, 4, 3, 1))

metrics <- c("l2", "l1", "tpr", "fdr", "auc", "elapsed", 
	     "mean", "expected.likelihood", "cindex")
metrics.lab <- c(TeX("$l_2^2$-error"), TeX("$l_1$-error"), "TPR", "FDR", 
		 "AUC", "Time (s)", "ELBO", "ELL", "c-index")
layout(matrix(1:9, nrow=3, byrow=T))
for (s in 1) {
    for (m in seq_along(metrics)) {
	xs <- c()
	for (i in 1:5) {
	    load(file=sprintf("../../RData/sensitivity/m%d.%d.RData", s, i))
	    x <- get(sprintf("m%d.%d", s, i))
	    lambda <- x[ , "lambda"] == 1
	    a0 <- x[ , "a0"][lambda]
	    x.mean <- as.numeric(by(x[lambda, metrics[m]], a0, mean, simplify=T))
	    xs <- rbind(xs, x.mean)
	}
	a0 <- unique(a0)
	matplot(a0, t(xs), log="x", type="l", lwd=2, col=2:6,
	    xlab=expression(alpha[0]), ylab=metrics.lab[m], 
	    main=metrics.lab[m], las=1)
	mtext(sprintf("(%s)", letters[m]), 3, line=0.5, adj=-0.15)

	if (any(m == c(3, 5, 7:9))) {
	    for (r in 1:nrow(xs)) {
		maxs <- which(xs[r, ] == max(xs[r, ]))
		points(a0[maxs], xs[r, maxs], pch=20+r, col=1,
		       bg=r+1, cex=1.3)
	    }
	}
	if (any(m == c(1:2, 4))) {
	    for (r in 1:nrow(xs)) {
		mins <- which(xs[r, ] == min(xs[r, ]))
		points(a0[mins], xs[r, mins], pch=20+r, 
		       col=1, bg=r+1, cex=1.3)
	    }
	}

	if (m == 1)
	    legend("topleft", legend=sprintf("Setting %d", 1:5), 
		   col=2:6, lwd=2, lty=1:5, bty="n")
    }
}
dev.off()

# ----------------------------------------
# 05-sensitivity :: lambda, alpha
# ----------------------------------------
metrics <- c("l2", "l1", "tpr", "fdr", "auc", "elapsed", 
	     "mean", "expected.likelihood", "cindex")
metrics.lab <- c(TeX("$l_2^2$-error"), TeX("$l_1$-error"), "TPR", "FDR", 
		 "AUC", "Time (s)", "ELBO", "ELL", "c-index")
for (i in 1:5) {
for (m in seq_along(metrics)) {
    pdf(file=sprintf("../../figures/sensitivity_alpha_lambda_%d_%d.pdf", i, m), 
	width=6, height=4.5)
    par(family="Times", mar=c(4, 4, 3, 1))
    # plot.new()

    load(file=sprintf("../../RData/sensitivity/m%d.%d.RData", s, i))
    x <- get(sprintf("m%d.%d", 1, i))

    lambda <- x[ , "lambda"]
    a0 <- x[ , "a0"]

    x.means <- outer(unique(a0), unique(lambda), Vectorize(function(a, l) {
	mean(x[ , metrics[m]][a0 == a & lambda == l])
    }))

    x <- unique(a0); 	y <- unique(lambda)
    xtick <- x;	 	ytick <- y
    x <- log(x);	y <- log(y)


    if (any(m == c(3, 5, 7:9))) {
	rev <- F
	mm <- which(x.means == max(x.means), arr.ind=T)
    }
    if (any(m == c(1:2, 4, 6))) {
	rev <- T
	mm <- which(x.means == min(x.means), arr.ind=T)
    }

    f <- filled.contour(x, y, x.means,
	color.palette = function(n) hcl.colors(n, rev=rev),
	plot.title = title(main = metrics.lab[m],
	    xlab = expression(alpha[0]), ylab=expression(lambda)),
	plot.axes = { 
	    axis(1, at=log(xtick), label=xtick);
	    axis(2, at=log(ytick), label=ytick);
	    points(x[mm[, 1]], y[mm[, 2]], pch=18, col=2, cex=1.2)
	}, useRaster=T,
	key.title = title(main=metrics.lab[m]))

    mtext(sprintf("(%s)", letters[m]), 3, line=1, adj=-0.15)
    
    dev.off()
}
}

