# Application of SVB to Yau et al. Breast cancer data
#
# The dataset is publicly available and was downloaded from
# https://xenabrowser.net/datapages/?cohort=Breast%20Cancer%20(Yau%202010)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#
library(data.table)
SAVE <- FALSE                          # save the processed dataset
LOAD <- TRUE                           # load saved dataset


# ----------------------------------------
# Process datasets
# ----------------------------------------
if (!LOAD) {

    yau.clinical <- fread("../../data/YAU_clinical.tsv")
    yau.rna <- fread("../../data/YAU_RNA_seq.tsv")

    # RNA seq data
    yau.rna.gene_names <- yau.rna$probe
    yau.rna <- t(yau.rna[, -1])
    colnames(yau.rna) <- yau.rna.gene_names
    yau.rna <- cbind(data.table(sampleID=rownames(yau.rna)), yau.rna)

    # Clinical data
    yau.samples <- colnames(yau.rna)[-1]
    yau.clinical <- yau.clinical[ , c("sampleID", "DMFS", "DMFS.time", "DMFS.unit")]
    colnames(yau.clinical)[2:4] <- c("OS", "OS.time", "OS.unit")

    # Merge, aligns times and rna expression data
    yau <- merge(yau.clinical, yau.rna, by="sampleID")
    yau <- yau[complete.cases(yau) , ] # rm missing values

} else {
    load("../../RData/data/yau.RData")
}

if (SAVE)
    save(yau, file="../../RData/data/yau.RData")


# ----------------------------------------
# Create dataframes
# ----------------------------------------
yau.OS <- yau[ , 2:3]
yau.rna <- as.matrix(yau[ , 5:ncol(yau)])


# cleanup
rm(list=c("SAVE", "LOAD", "yau"))
