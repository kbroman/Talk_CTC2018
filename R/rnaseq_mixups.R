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
saveRDS(gmap, "Data/attieDO_gmap.rds")
pmap <- interp_map(gmap, do$gmap, do$pmap)
saveRDS(pmap, "Data/attieDO_pmap.rds")

file <- "Data/attieDO_v1_apr.rds"
if(file.exists(file)) {
    apr <- readRDS(file)
} else {
    for(i in c(1:19,"X")) {
        cat("chr", i, "\n")
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
    batches <- batch_cols(expr.mrna, 50)
    peaks <- vector("list", length(batches))
    for(i in seq_along(batches)) {
        cat("batch", i, "\n")
        out <- scan1(apr, expr.mrna[,batches[[i]]$cols], cores=0)
        peaks[[i]] <- find_peaks(out, gmap, threshold=25)
        saveRDS(peaks, file)
    }
    saveRDS(peaks, file)
}

file <- "Data/rnaseq_mixups_results.RData"
if(file.exists(file)) {
    load(file)
} else {
    peaks <- do.call("rbind", peaks)
    peaks <- peaks[peaks$lod > 100,]
    expr_obs <- expr.mrna[,peaks$lodcolumn]
    library(parallel)
    expr_exp <- mclapply(1:nrow(peaks), function(i) {
        p <- pull_genoprobpos(apr, find_marker(gmap, peaks$chr[i], peaks$pos[i]))
        out <- fit1(p, expr_obs[,i])
        p %*% out$coef }, mc.cores=detectCores())
    expr_exp <- matrix(unlist(expr_exp), ncol=ncol(expr_obs))
    dimnames(expr_exp) <- list(rownames(apr[[1]]), colnames(expr_obs))

    d <- lineup2::dist_betw_matrices(expr_obs, expr_exp)
    d <- d[order(as.numeric(sub("DO", "", rownames(d)))), order(as.numeric(sub("DO", "", colnames(d))))]

    save(expr_obs, expr_exp, d, file=file)
}
