# convert genotypes to fst format

library(qtl2)
library(qtl2fst)

# load genoprobs
for(i in c(1:19,"X")) {
    if(i==1) {
        pr <- readRDS("Data/Genoprobs/attieDO_v1_pr_1.rds")
    } else {
        pr <- cbind(pr, readRDS(paste0("Data/Genoprobs/", "attieDO_v1_pr_", i, ".rds")))
    }
    cat(i, "\n")
}

# save in fst format
pr_fst <- fst_genoprob(pr, "attieDO_v1_pr", "Data/GenoprobsFST")

# save the key object
saveRDS(pr_fst, "Data/GenoprobsFST/attieDO_v1_pr_fst.rds")

# do the same with the allele probs
apr <- readRDS("Data/attieDO_v1_apr.rds")
apr_fst <- fst_genoprob(apr, "attieDO_v1_apr", "Data/AprobsFST")
saveRDS(apr_fst, "Data/AprobsFST/attieDO_v1_apr_fst.rds")
