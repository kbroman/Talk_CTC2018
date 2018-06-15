library(qtl2fst)
library(qtl2)
library(broman)

do <- readRDS("Data/attieDO_v1.rds")
gmap <- readRDS("Data/attieDO_gmap.rds")
pr_fst <- readRDS("Data/GenoprobsFST/attieDO_v1_pr_fst.rds")
errlod <- do.call("cbind", readRDS("R/_cache/errlod.rds"))
snpint <- fst::read.fst("Data/attieDO_allint.fst")
omit <- c("DO126", "DO133", "DO134", "DO135", "DO136", "DO138", "DO139",
          "DO140", "DO141", "DO144", "DO145", "DO149", "DO150", "DO165",
          "DO171", "DO202", "DO204", "DO021", "DO212", "DO213", "DO214",
          "DO215", "DO340", "DO357", "DO397", "DO402_b", "DO403_b")
snpint <- snpint[,!(names(snpint) %in% omit)]
names(snpint) <- sub("_b$", "", names(snpint))
snpint <- snpint[,c(1:2, 2+order(as.numeric(sub("DO", "", names(snpint)[-(1:2)]))))]

g <- do.call("cbind", do$geno)
g[g==0] <- 4
fg <- do.call("cbind", do$founder_geno)

mar <- "UNCHS021230"
mar <- "UNCHS027234"
mar <- "UNCHS006327"

mar <- "UNCrs232949530" # (F founder missing; looks like it should be 1)
                        # ... tried both and this fit the allele intensities clusters better
pos <- find_markerpos(gmap, mar)
this_fg <- t(do$founder_geno[[pos$chr]][,mar,drop=FALSE])
this_fg[,"F"] <- 1
sdp <- calc_sdp(this_fg)
this_fg[,"F"] <- 3
alt_sdp <- calc_sdp(this_fg)
snpinfo <- data.frame(chr=pos$chr, pos=pos$pos, snp=mar,
                      sdp=sdp)
snpinfo <- index_snps(gmap, snpinfo)
snp_pr <- genoprob_to_snpprob(pr_fst, snpinfo)
snp_inf <- maxmarg(snp_pr)[[1]][,1]
alt_snpinfo <- snpinfo;alt_snpinfo$sdp <- alt_sdp
alt_snp_pr <- genoprob_to_snpprob(pr_fst, alt_snpinfo)
alt_snpinf <- maxmarg(alt_snp_pr)[[1]][,1]

a <- unlist(snpint[snpint[,1]==mar & snpint[,2]=="X",-(1:2)])
b <- unlist(snpint[snpint[,1]==mar & snpint[,2]=="Y",-(1:2)])
grayplot(a, b, bg=brocolors("f2")[g[,mar]])
grayplot(a, b, bg=brocolors("f2")[snp_inf])
grayplot(a, b, bg=brocolors("f2")[alt_snpinf])
grayplot(a, b, bg=brocolors("sex")[2-(errlod[,mar]>4)])


mar <- "Mit0045" # genotypes have no relationship to the allele intensities

pos <- find_markerpos(gmap, mar)
this_fg <- t(do$founder_geno[[pos$chr]][,mar,drop=FALSE])
sdp <- calc_sdp(this_fg)
snpinfo <- data.frame(chr=pos$chr, pos=pos$pos, snp=mar,
                      sdp=sdp)
snpinfo <- index_snps(gmap, snpinfo)
snp_pr <- genoprob_to_snpprob(pr_fst, snpinfo)
snp_inf <- maxmarg(snp_pr)[[1]][,1]
snp_inf[is.na(snp_inf)] <- 4

a <- unlist(snpint[snpint[,1]==mar & snpint[,2]=="X",-(1:2)])
b <- unlist(snpint[snpint[,1]==mar & snpint[,2]=="Y",-(1:2)])
grayplot(a, b, bg=brocolors("f2")[g[,mar]])
grayplot(a, b, bg=brocolors("f2")[snp_inf])
grayplot(a, b, bg=brocolors("sex")[2-(errlod[,mar]>4)])

mar <- "UNCHS028710"

pos <- find_markerpos(gmap, mar)
this_fg <- t(do$founder_geno[[pos$chr]][,mar,drop=FALSE])
sdp <- calc_sdp(this_fg)
snpinfo <- data.frame(chr=pos$chr, pos=pos$pos, snp=mar,
                      sdp=sdp)
snpinfo <- index_snps(gmap, snpinfo)
snp_pr <- genoprob_to_snpprob(pr_fst, snpinfo)
snp_inf <- maxmarg(snp_pr)[[1]][,1]
snp_inf[is.na(snp_inf)] <- 4

a <- unlist(snpint[snpint[,1]==mar & snpint[,2]=="X",-(1:2)])
b <- unlist(snpint[snpint[,1]==mar & snpint[,2]=="Y",-(1:2)])
grayplot(a, b, bg=brocolors("f2")[g[,mar]])
grayplot(a, b, bg=brocolors("f2")[snp_inf])
grayplot(a, b, bg=brocolors("sex")[2-(errlod[,mar]>4)])

mar <- "UNCHS000060"

pos <- find_markerpos(gmap, mar)
this_fg <- t(do$founder_geno[[pos$chr]][,mar,drop=FALSE])
sdp <- calc_sdp(this_fg)
snpinfo <- data.frame(chr=pos$chr, pos=pos$pos, snp=mar,
                      sdp=sdp)
snpinfo <- index_snps(gmap, snpinfo)
snp_pr <- genoprob_to_snpprob(pr_fst, snpinfo)
snp_inf <- maxmarg(snp_pr)[[1]][,1]
snp_inf[is.na(snp_inf)] <- 4

a <- unlist(snpint[snpint[,1]==mar & snpint[,2]=="X",-(1:2)])
b <- unlist(snpint[snpint[,1]==mar & snpint[,2]=="Y",-(1:2)])
grayplot(a, b, bg=brocolors("f2")[g[,mar]], xlim=c(0, 2))
grayplot(a, b, bg=brocolors("f2")[snp_inf], xlim=c(0, 2))
grayplot(a, b, bg=brocolors("sex")[2-(errlod[,mar]>4)], xlim=c(0,2))

mar <- "UNCJPD002501" # inferred snp genotypes have no relationship to snp intensity plot
                      # ...because the snp doesn't actually belong at the cited position
pos <- find_markerpos(gmap, mar)
this_fg <- t(do$founder_geno[[pos$chr]][,mar,drop=FALSE])
sdp <- calc_sdp(this_fg)
snpinfo <- data.frame(chr=pos$chr, pos=pos$pos, snp=mar,
                      sdp=sdp)
snpinfo <- index_snps(gmap, snpinfo)
snp_pr <- genoprob_to_snpprob(pr_fst, snpinfo)
snp_inf <- maxmarg(snp_pr)[[1]][,1]
snp_inf[is.na(snp_inf)] <- 4

a <- unlist(snpint[snpint[,1]==mar & snpint[,2]=="X",-(1:2)])
b <- unlist(snpint[snpint[,1]==mar & snpint[,2]=="Y",-(1:2)])
grayplot(a, b, bg=brocolors("f2")[g[,mar]])
grayplot(a, b, bg=brocolors("f2")[snp_inf])
grayplot(a, b, bg=brocolors("sex")[2-(errlod[,mar]>4)])

mar <- "UNCrs50250732" # observed genotype calls really odd, but in any way this is a hard one
pos <- find_markerpos(gmap, mar)
this_fg <- t(do$founder_geno[[pos$chr]][,mar,drop=FALSE])
sdp <- calc_sdp(this_fg)
snpinfo <- data.frame(chr=pos$chr, pos=pos$pos, snp=mar,
                      sdp=sdp)
snpinfo <- index_snps(gmap, snpinfo)
snp_pr <- genoprob_to_snpprob(pr_fst, snpinfo)
snp_inf <- maxmarg(snp_pr)[[1]][,1]
snp_inf[is.na(snp_inf)] <- 4

a <- unlist(snpint[snpint[,1]==mar & snpint[,2]=="X",-(1:2)])
b <- unlist(snpint[snpint[,1]==mar & snpint[,2]=="Y",-(1:2)])
grayplot(a, b, bg=brocolors("f2")[g[,mar]])
grayplot(a, b, bg=brocolors("f2")[snp_inf])
grayplot(a, b, bg=brocolors("sex")[2-(errlod[,mar]>4)])

calc_snpprob <-
    function(marker, cross=do, genoprobs=pr_fst, map=gmap)
{
}


plot_snpint <-
    function(marker, use=c("observed", "inferred", "error"))
{


}
