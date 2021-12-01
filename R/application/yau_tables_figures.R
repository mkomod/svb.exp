# Tables and figures for Yau BC dataset
library(survival)
library(latex2exp)

source("yau_data.R", chdir=T)
source("00-functions.R")


# ----------------------------------------
# Setup 
# ----------------------------------------
lambda <- c(0.05, 0.1, 0.25, 0.5, 0.75, 
	     1.0, 1.25, 1.5, 1.75, 2.0, 
	     2.50, 3, 4, 5)            # lambda grid


# ----------------------------------------
# Model selection
# ----------------------------------------
m0 <- c(); s0 <- c(); mMax <- c(); genes <- c()

for (l in as.character(lambda)) {
    load(sprintf("../../RData/models/yau/yau_l_%s.RData", l))
    fit.summary <- model_summary(yau.OS$OS.time, yau.OS$OS, yau.rna, fit)

    m <- apply(fit.summary, 1, mean);   m0 <- rbind(m0, m)
    m <- apply(fit.summary, 1, max);  mMax <- rbind(mMax, m)
    s <- apply(fit.summary, 1, sd);     s0 <- rbind(s0, s)

    for (i in 1:fit$folds)
	genes <- c(genes, colnames(yau.rna)[fit[[i]]$g > 0.5])
}

# metrics
rownames(m0) <- lambda
round(m0, 2)

rownames(mMax) <- lambda
round(mMax, 3)

gt <- sort(table(genes), T) / (length(lambda)*fit$folds)
round(gt[gt >= 0.05], 3)

# pick a model, largest validation c-index
load("../../RData/models/yau/yau_l_1.25.RData")

# ----------------------------------------
# Figures
# ----------------------------------------
pdf("../../figures/yau_convergence_diagnostics.pdf", width=6, height=9)
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



# ----------------------------------------
# Tables
# ----------------------------------------
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
per_line <- 8

for (i in 1:ceiling(length(genes) / per_line)) {
    g <- genes[(1 + (i - 1)*per_line):(min(i*per_line, length(genes)))]
    cat(paste0(names(g), collapse=" & "), "\\\\\n",
	paste0(round(g, 3), collapse=" & "), "\\\\ \\hline \n")
}


# --------------------
# Prognositc index
# --------------------
y <- survival::Surv(matrix(yau.OS$OS.time), matrix(0+yau.OS$OS))

model_summary(yau.OS$OS.time, yau.OS$OS, yau.rna, fit)
tr.samples <- fit$tr.samples[[2]]
ts.samples <- fit$ts.samples[[2]]
f <- fit[[2]]

eta <- yau.rna[tr.samples, ] %*% as.matrix(f$m * f$g)
eta.0 <- yau.rna[tr.samples, ] %*% as.matrix(f$m * f$g)
risk <- eta > median(eta)
s <-  survfit(y[tr.samples] ~ risk, data=as.data.frame(risk))
s0 <- survfit(y[tr.samples] ~ -1)
lr <- survdiff(y[tr.samples] ~ risk, data=as.data.frame(risk), rho=0)
pval <- pchisq(lr$chisq, df=1, lower.tail=FALSE)
txt <- ifelse(pval > 1e-3, "log-rank: %.3f", "log-rank: %.2e")

pdf(file="../../figures/yau_km_training.pdf", width=4, height=4)
par(family="Times", mar=c(4, 4, 1, 1))
    plot(s, main="Training", ylab="Surv. prob", xlab="Time (years)",
	col=c("blue", "red"), lwd=1.2, conf.int=TRUE, pch=4,
	cex=0.7, conf.times=0:12*2, mark.time=T)
    lines(s0, conf.int=FALSE, lwd=1)
    legend("bottomleft",
	legend=c("Percentile 0-50", "Percentile 50-100", "All", 
		 sprintf(txt, pval)),
	col=c("blue", "red", "black", 0), lty=c(1,1,1, 0), cex=.8, bty="n")
dev.off()

eta <- yau.rna[ts.samples, ] %*% as.matrix(f$m * f$g)
risk <- eta > median(eta.0)
s <-  survfit(y[ts.samples] ~ risk, data=as.data.frame(risk))
s0 <- survfit(y[ts.samples] ~ -1)
lr <- survdiff(y[ts.samples] ~ risk, data=as.data.frame(risk), rho=0)
pval <- pchisq(lr$chisq, df=1, lower.tail=FALSE)
txt <- ifelse(pval > 1e-3, "log-rank: %.3f", "log-rank: %.2e")

pdf(file="../../figures/yau_km_validation.pdf", width=4, height=4)
par(family="Times", mar=c(4, 4, 1, 1))
    plot(s, main="Validation", ylab="Surv. prob", xlab="Time (years)",
	col=c("blue", "red"), lwd=1.2, conf.int=TRUE, pch=4,
	cex=0.7, conf.times=0:12*2, mark.time=T)
    lines(s0, conf.int=FALSE, lwd=1)
    legend("bottomleft",
	legend=c("Percentile 0-50", "Percentile 50-100", "All", 
		 sprintf(txt, pval)),
	col=c("blue", "red", "black", 0), lty=c(1,1,1, 0), cex=.8, bty="n")
dev.off()


