#
# Miscellaneous figures
#

# ----------------------------------------
# Spike-and-Slab distribution
# ----------------------------------------
laplace <- function(x, lambda=0.5) 0.5 * lambda * exp(-abs(x) * lambda)

xs <- seq(-3, 3, 0.05)
ys <- laplace(xs, lambda=2)

pdf(file="../figures/prior_spike_slab.pdf", width=4, height=2.5)
par(family="Times", mar=c(3, 3, 0, 0))
    plot(xs, ys, lwd=3, type="l", ylim=c(0, 7), axes=FALSE,
	  ylab="", xlab="", frame.plot=F)
    axis(side=1, at=c(-3, 0, 3), labels=FALSE)
    axis(side=1, at=c(0), labels=TRUE)
    arrows(0, 0, 0, 5, code=0, lwd=3, col="blue")
    points(0, 5, pch=20, cex=1.2, col="blue")
    legend("topright", legend=c("spike", "slab"), 
	   lty=c(1,1), col=c("blue", "black"), lwd=c(3,3), bty="n", cex=1)
dev.off()


# ----------------------------------------
# Spike-and-Slab LASSO
# ----------------------------------------
xs <- seq(-3, 3, 0.05)
y1 <- laplace(xs, lambda=2)
y2 <- laplace(xs, lambda=12)

pdf(file="../figures/prior_spike_slab_lasso.pdf", width=4, height=2.5)
par(family="Times", mar=c(3, 3, 0, 0))
    plot(xs, y1, lwd=3, type="l", ylim=c(0, 7), axes=FALSE,
	  ylab="", xlab="", frame.plot=F)
    lines(xs, y2, lwd=2, col="blue")
    axis(side=1, at=c(-3, 0, 3), labels=FALSE)
    axis(side=1, at=c(0), labels=TRUE)
    # axis(side=2, at=c(0, 3), labels=FALSE)
    legend("topright", legend=c("spike", "slab"), 
	   lty=c(1,1), col=c("blue", "black"), lwd=c(3,3), bty="n", cex=1)
dev.off()


# ----------------------------------------
# piMOM prior
# ----------------------------------------
pimom <- function(x, r=1, tau=0.5) ifelse(x!=0, tau^(r/2)/gamma(r/2)*abs(x)^-(r+1)*exp(-tau/x^2), 0)

xs <- seq(-3, 3, 0.05)
ys <- pimom(xs)

pdf(file="../figures/prior_pimom.pdf", width=4, height=2.5)
par(family="Times", mar=c(3, 3, 0, 0))
    plot(xs, ys, lwd=3, type="l", ylim=c(0, 0.6), axes=FALSE,
	  ylab="", xlab="", frame.plot=F)
    # lines(xs, dnorm(xs, 0, 0.5), lwd=2, col="red", lty=2)
    axis(side=1, at=c(-3, 0, 3), labels=FALSE)
    axis(side=1, at=c(0), labels=TRUE)
    # axis(side=2, at=c(0, 0.4), labels=FALSE)
dev.off()


# ----------------------------------------
# Laplace prior
# ----------------------------------------
xs <- seq(-3, 3, 0.05)
ys <- laplace(xs, 2)

pdf(file="../figures/prior_laplace.pdf", width=4, height=2.5)
par(family="Times", mar=c(3, 3, 0, 0))
    plot(xs, ys, lwd=3, type="l", ylim=c(0, 1.2), axes=FALSE,
	  ylab="", xlab="", frame.plot=F)
    # lines(xs, dnorm(xs, 0, 0.5), lwd=2, col="red", lty=2)
    axis(side=1, at=c(-3, 0, 3), labels=FALSE)
    axis(side=1, at=c(0), labels=TRUE)
    # axis(side=2, at=c(0, 0.4), labels=FALSE)
dev.off()


# ----------------------------------------
# Comparison of Normal to Laplace
# ----------------------------------------
xs <- seq(-5, 5, 0.05)
plot(xs, dnorm(xs), type="l", col="red", lwd=2, ylim=c(0, 0.8), lty=2,
     ylab="density")
lines(xs, laplace(xs, lambda=1), col="green", lwd=2)
lines(xs, laplace(xs, lambda=0.5), col="orange", lwd=2)
legend("topright", 
       legend=c(expression(N(0,1)), expression(L(lambda==1.0)), 
		expression(L(lambda==0.5))),
       lty=c(2, 1, 1), lwd=c(2, 2, 2),  col=c("red", "orange", "green"))

