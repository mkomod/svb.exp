# Tables and figures for TCGA OV data
#
library(survival)
library(latex2exp)

source("tcga_data.R", chdir=T)
source("00-functions.R")

# ----------------------------------------
# Setup
# ----------------------------------------
lambda  <- c(0.05, 0.1, 0.25, 0.5, 0.75, 
	     1.0, 1.25, 1.5, 1.75, 2.0, 
	     2.50, 3, 4, 5)            # lambda grid


# ----------------------------------------
# Model selection
# ----------------------------------------
m0 <- c(); s0 <- c(); mMax <- c(); genes <- c()

for (l in as.character(lambda)) {
    load(sprintf("../../RData/models/tcga/tcga_l_%s.RData", l))
    fit.summary <- model_summary(tcga.OS$OS.time, tcga.OS$OS, tcga.rna, fit)

    m <- apply(fit.summary, 1, mean); m0 <- rbind(m0, m)
    s <- apply(fit.summary, 1, sd);   s0 <- rbind(s0, s)
    m <- apply(fit.summary, 1, max);  mMax <- rbind(mMax, m)

    for (i in 1:fit$folds)
	genes <- c(genes, colnames(tcga.rna)[fit[[i]]$g > 0.5])
}

rownames(m0) <- lambda
round(m0, 2)

rownames(mMax) <- lambda
round(mMax, 3)

gt <- sort(table(genes), T) / (length(lambda)*fit$folds)
round(gt[gt > 0.05], 3)


load("../../RData/models/tcga/tcga_l_1.RData")
model_summary(tcga.OS$OS.time, tcga.OS$OS, tcga.rna, fit)
f <- fit[[9]]

plot(f$m)
plot(f$m * f$g, ylim=c(-1, 1))
plot(f$s)
plot(f$g)


# ----------------------------------------
# Figures
# ----------------------------------------
# fit lambda=1
pdf("../../figures/tcga_convergence_diagnostics.pdf", width=6, height=9)
    par(family="Times", mar=c(4, 4, 2, 1)) 
    layout(matrix(1:6, nrow=3, byrow=T))
    plot_fit(fit, "tr.elbo.mean",  main="ELBO - Training", cex.main=1.1,
	 ylab=expression(hat(L)[Q]), xlab="Iteration", las=1)
    mtext("A", 3, line=0.5, adj=-0.2, font=1)
    plot_fit(fit, "ts.elbo.mean",  main="ELBO - Validation", cex.main=1.1,
	 ylab=expression(hat(L)[Q]), xlab="Iteration", las=1)
    mtext("B", 3, line=0.5, adj=-0.2, font=1)
    plot_fit(fit, "tr.elbo.ll", main="Exp. log-lik - Training", cex.main=1.1, 
	 ylab=TeX("$E(l(\\beta; D)$"), xlab="Iteration", las=1)
    mtext("C", 3, line=0.5, adj=-0.2, font=1)
    plot_fit(fit, "ts.elbo.ll", main="Exp. log-lik - Validation", cex.main=1.1, 
	 ylab=TeX("$E(l(\\beta; D)$"), xlab="Iteration", las=1)
    mtext("D", 3, line=0.5, adj=-0.2, font=1)
    plot_fit(fit, "tr.elbo.kl", negate=T, main="Neg. KL - Training", cex.main=1.1,
	 ylab=TeX("$KL(Q ||\\Pi)"), xlab="Iteration", las=1)
    mtext("E", 3, line=0.5, adj=-0.2, font=1)
    plot_fit(fit, "ts.elbo.kl", negate=T, main="Neg. KL - Validation", cex.main=1.1,
	 ylab=TeX("$KL(Q ||\\Pi)"), xlab="Iteration", las=1)
    mtext("F", 3, line=0.5, adj=-0.2, font=1)
dev.off()


pdf("../../figures/tcga_model_selection.pdf", width=6, height=9)
par(mar = c(4, 4, .5, 4) + 0.3)  # Leave space for z axis
    plot(lambda, m0[ , "ts.elbo.mean"], ylim=c(-350, 0), type="b", pch=20, col="darkblue", 
	 ylab="ELBO", xlab=expression(lamdba))
    par(new = TRUE)
    plot(lambda, m0[ , "ts.cin"], ylim=c(0.5, 1), type="b", bty="n", axes=F,
	 pch=20, col="darkorange", xlab="", ylab="")
    axis(side=4, at = pretty(range(0.5, 1)))
    mtext("c-index", side=4, line=3)
dev.off()



# ----------------------------------------
# Tables
# ----------------------------------------
# Goodness of fit across different models
for (i in 1:length(lambda)) 
{
    m <- m0[i, c("tr.elbo.mean", "tr.elbo.ll", "tr.elbo.kl", "tr.cin", 
		 "nvars", "ts.elbo.mean", "ts.elbo.ll", "ts.elbo.kl", "ts.cin")]
    s <- s0[i, c("tr.elbo.mean", "tr.elbo.ll", "tr.elbo.kl", "tr.cin",
		 "nvars", "ts.elbo.mean", "ts.elbo.ll", "ts.elbo.kl", "ts.cin")]
    x <- sapply(1:length(m), function(i) c(m[i], s[i]))
    
    cat(
	sprintf("%.2f", lambda[i]), 
	sprintf("& %.1f (%.1f)", x[1, 1:3], x[2, 1:3]),
	sprintf("& %.3f (%.3f)", x[1, 4], x[2, 4]),
	sprintf("& %.1f (%.1f)", x[1, 5:8], x[2, 5:8]),
	sprintf("& %.3f (%.3f)", x[1, 9], x[2, 9]),
	"\\\\ \n"
    )
}

# Genes 
genes <- gt[gt >= 0.05]
per_line <- 10

for (i in 1:ceiling(length(genes) / per_line)) {
    g <- genes[(1 + (i - 1)*per_line):(min(i*per_line, length(genes)))]
    cat(paste0(names(g), collapse=" & "), "\\\\\n",
	paste0(round(g, 3), collapse=" & "), "\\\\ \\hline \n")
}


# --------------------
# Prognositc index
# --------------------
y <- survival::Surv(matrix(tcga.OS$OS.time), matrix(0+tcga.OS$OS))

tr.samples <- fit$tr.samples[[9]]
ts.samples <- fit$ts.samples[[9]]
eta <- tcga.rna[tr.samples, ] %*% as.matrix(f$m * f$g)
eta.0 <- tcga.rna[tr.samples, ] %*% as.matrix(f$m * f$g)
risk <- eta > median(eta)
s <- survfit(y[tr.samples] ~ risk, data=as.data.frame(risk))
s0 <- survfit(y[tr.samples] ~ -1)
lr <- survdiff(y[tr.samples] ~ risk, data=as.data.frame(risk), rho=0)
pval <- pchisq(lr$chisq, df=1, lower.tail=FALSE)
txt <- ifelse(pval > 1e-3, "log-rank: %.3f", "log-rank: %.2e")

pdf(file="../../figures/tcga_km_training.pdf", width=4, height=4)
# par(family="Times", mar=c(4, 4, 3, 1))
    plot(s, main="Training", ylab="Surv. prob", xlab="Time (years)",
	col=c("blue", "red"), lwd=1.2, conf.int=TRUE, pch=4,
	cex=0.7, conf.times=0:12*500, mark.time=T, xscale=365.25)
    lines(s0, conf.int=FALSE, lwd=1)
    legend("topright",
	legend=c("Percentile 0-50", "Percentile 50-100", "All", 
		 sprintf(txt, pval)),
	col=c("blue", "red", "black", 0), lty=c(1,1,1, 0), cex=.8, bty="n")
dev.off()


eta <- tcga.rna[ts.samples, ] %*% as.matrix(f$m * f$g)
risk <- eta > median(eta.0)
s <- survfit(y[ts.samples] ~ risk, data=as.data.frame(risk))
s0 <- survfit(y[ts.samples] ~ -1)
lr <- survdiff(y[ts.samples] ~ risk, data=as.data.frame(risk), rho=0)
pval <- pchisq(lr$chisq, df=1, lower.tail=FALSE)
txt <- ifelse(pval > 1e-3, "log-rank: %.3f", "log-rank: %.2e")

pdf(file="../../figures/tcga_km_validation.pdf", width=4, height=4)
par(family="Times", mar=c(4, 4, 1, 1))
    plot(s, main="Validation", ylab="Surv. prob", xlab="Time (years)",
	col=c("blue", "red"), lwd=1.2, conf.int=TRUE, pch=4,
	cex=0.7, conf.times=0:12*500, mark.time=T, xscale=365.25)
    lines(s0, conf.int=FALSE, lwd=1)
    legend("topright",
	legend=c("Percentile 0-50", "Percentile 50-100", "All", 
		 sprintf(txt, pval)),
	col=c("blue", "red", "black", 0), lty=c(1,1,1, 0), cex=.8, bty="n")
dev.off()


