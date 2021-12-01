# Processing of TCGA Ovarian cancer dataset
#
# The dataset used is publicly available and can be downloaded from
# https://xenabrowser.net/datapages/?dataset=TCGA.OV.sampleMap%2FHT_HG-U133A&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#
library(data.table)
SAVE <- FALSE                          # save data
LOAD <- TRUE                           # load processed data


# ----------------------------------------
# Process data
# ----------------------------------------
if (!LOAD) {

tcga.rna <- data.table::fread(file="../../data/TCGA_RNA_seq.csv")
tcga.surv <- data.table::fread(file="../../data/TCGA_OV_survival.tsv")

# RNA seq data
tcga.rna.samples <- colnames(tcga.rna)[-1]
tcga.rna.gene_names <- tcga.rna$sample
tcga.rna <- t(tcga.rna[, -1])
colnames(tcga.rna) <- tcga.rna.gene_names
tcga.rna <- cbind(data.table::data.table(sampleID=tcga.rna.samples), tcga.rna)

# Clinical data
colnames(tcga.surv)[1:2] <- c("sampleID", "patient")
tcga.surv <- tcga.surv[, c("sampleID", "OS", "OS.time")]

# Merge, aligns times and rna expression data
tcga <- merge(tcga.surv, tcga.rna, by="sampleID")
tcga <- tcga[complete.cases(tcga), ]

} else {
    load("../../RData/data/tcga.RData")
}

if (SAVE)
    save(tcga, file="../../RData/data/tcga.RData")


# ----------------------------------------
# Create dataframes
# ----------------------------------------
tcga.OS <- tcga[ , 2:3]
tcga.rna <- as.matrix(tcga[ , 4:ncol(tcga)])

# cleanup
rm(list=c("LOAD", "SAVE", "tcga"))

