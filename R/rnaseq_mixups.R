# study mix-ups with RNA-seq data

library(lineup2)
library(qtl2)

# cross object
file <- "Data/attieDO_v1.rds"
if(file.exists(file)) {
    do <- readRDS(file)
} else {
    do <- readRDS("Data/attieDO_v0.rds")
    tab <- table(do$covar$mouse)
    dup_mice <- names(tab[tab>1])
    dup_mice <- sprintf("DO%03d", as.numeric(sub("DO", "", dup_mice)))
    nm <- n_missing(do, summary="prop")
    nm_before <- nm[dup_mice]
    nm_after <- nm[paste0(dup_mice, "_b")]
    keep <- names(nm_before)[nm_before < nm_after]
    omit <- dup_mice
    omit[omit %in% keep] <- paste0(omit[omit %in% keep], "_b")
    ids <- ind_ids(do)
    do <- do[ids[!(ids %in% omit)],]
    for(i in seq_along(do$geno)) {
        rownames(do$geno[[i]]) <- sub("_b$", "", rownames(do$geno[[i]]))
    }
    rownames(do$covar) <- sub("_b$", "", rownames(do$covar))
    rownames(do$cross_info) <- sub("_b$", "", rownames(do$cross_info))
    names(do$is_female) <- sub("_b$", "", names(do$is_female))
    stopifnot(check_cross2(do))
    saveRDS(do, file)
}

gmap <- insert_pseudomarkers(do$gmap, step=0.2, stepwidth="max")
file <- "Data/attieDO_v1_apr.rds"
if(file.exists(file)) {
    apr <- readRDS(file)
} else {
    for(i in c(1:19,"X")) {
        cat("Chr ", i, "\n")
        pr <- calc_genoprob(do[,i], gmap[i], error_prob=0.002, map_function="c-f", cores=0)
        saveRDS(pr, file=paste0("Data/Genoprobs/attieDO_v1_pr_", i, ".rds"))
        this_apr <- genoprob_to_alleleprob(pr)
        if(i==1) {
            apr <- this_apr
        } else {
            apr <- cbind(apr, this_apr)
        }
    }
    saveRDS(apr, file)
}

# grab expression data
file <- "Data/rnaseq_pheno.rds"
if(file.exists(file)) {
    expr.mrna <- readRDS(file)
} else {
    load("~/Projects/AttieDOv2/RawData/RNAseq/DO378_islet_v2.RData")
    rownames(expr.mrna) <- sub("\\-", "", rownames(expr.mrna))
    saveRDS(expr.mrna, file)
}

# peak QTL
file <- "Data/rnaseq_peaks.rds"
if(file.exists(file)) {
    peaks <- readRDS(file)
} else {
    out <- scan1(apr, expr.mrna, cores=0)
    peaks <- find_peaks(apr, gmap, threshold=25)
    saveRDS(peaks, file)
}
