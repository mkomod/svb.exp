#
# Application of SVB to TCGA Ovarian cancer data
#
library(survival)
library(survival.svb)
library(BhGLM)
library(parallel)

source("00-functions.R")
source("tcga_data.R", chdir=T)


# ----------------------------------------
# Setup
# ----------------------------------------
CORES <- 15                            # Number of CPU cores
n <- nrow(tcga.rna)
lambda  <- c(0.05, 0.1, 0.25, 0.5, 0.75, 
	     1.0, 1.25, 1.5, 1.75, 2.0, 
	     2.50, 3, 4, 5)            # lambda grid


# ----------------------------------------
# Fit SVB
# ----------------------------------------
fit_svb(tcga.OS$OS.time, tcga.OS$OS, tcga.rna, 
	lambda, a0=ncol(tcga.rna)/100, b0=ncol(tcga.rna), 
	folds=10, CORES=CORES, dname="tcga", remove_ties=T)


