#
# Application of SVB to Yau et al 2010 dataset
#
library(survival.svb)

source("00-functions.R")
source("yau_data.R", chdir=T)


# ----------------------------------------
# Setup
# ----------------------------------------
CORES <- 15                            # number of CPU cores
lambda  <- c(0.05, 0.1, 0.25, 0.5, 0.75, 
	     1.0, 1.25, 1.5, 1.75, 2.0, 
	     2.50, 3, 4, 5)            # lambda grid
n <- nrow(yau.rna)


# ----------------------------------------
# Fit SVB
# ----------------------------------------
fit_svb(yau.OS$OS.time, yau.OS$OS, yau.rna, 
	lambda, a0=ncol(yau.rna)/100, b0=ncol(yau.rna),
	folds=10, CORES=CORES, dname="yau", remove_ties=T)


