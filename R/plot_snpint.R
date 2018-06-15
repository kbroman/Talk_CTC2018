## Functions to get inferred SNP genotypes from genotype probabilities
## and then to plot SNP allele intensities, colored by observed or inferred genotypes,
## or to indicate genotyping errors


# library(qtl2fst)
# library(qtl2)
# library(broman)

# do <- readRDS("Data/attieDO_v1.rds")
# gmap <- readRDS("Data/attieDO_gmap.rds")
# pr_fst <- readRDS("Data/GenoprobsFST/attieDO_v1_pr_fst.rds")
# errlod <- do.call("cbind", readRDS("R/_cache/errlod.rds"))
# snpint <- fst::read.fst("Data/attieDO_allint.fst")
# omit <- c("DO126", "DO133", "DO134", "DO135", "DO136", "DO138", "DO139",
#           "DO140", "DO141", "DO144", "DO145", "DO149", "DO150", "DO165",
#           "DO171", "DO202", "DO204", "DO021", "DO212", "DO213", "DO214",
#           "DO215", "DO340", "DO357", "DO397", "DO402_b", "DO403_b")
# snpint <- snpint[,!(names(snpint) %in% omit)]
# names(snpint) <- sub("_b$", "", names(snpint))
# snpint <- snpint[,c(1:2, 2+order(as.numeric(sub("DO", "", names(snpint)[-(1:2)]))))]
#
# g <- do.call("cbind", do$geno)
# g[g==0] <- 4
# fg <- do.call("cbind", do$founder_geno)
#
# mar <- "UNCHS021230"
# mar <- "UNCHS027234"
# mar <- "UNCHS006327"
# mar <- "UNCrs232949530" # (F founder missing; looks like it should be 1)
#                         # ... tried both and this fit the allele intensities clusters better
# mar <- "Mit0045" # genotypes have no relationship to the allele intensities (my guess this is mitochondrial?)
# mar <- "UNCHS028710"
# mar <- "UNCHS000060"
# mar <- "UNCJPD002501" # inferred snp genotypes have no relationship to snp intensity plot
#                       # ...because the snp doesn't actually belong at the cited position
# mar <- "UNCrs50250732" # observed genotype calls really odd, but in any way this is a hard one
# mar <- "ICR4840"  # founder D is missing; clearly should be 3
# mar <- "UNCHS032160" # extra blob; should be clear but terrible calls ... last 2 waves all NA
# mar <- "UNCHS006705" # hard to call but mostly right; well-placed blobs but they bleed together
# mar <- "UNC9180222"  # ugly partly due to low MAF
# mar <- "UNCHS022987" # puzzling no calls; mostly in waves 3 and 4
# mar <- "UNCHS029321" # similar...no calls mostly in wave 2

infer_snp <-
    function(marker, cross=do, genoprobs=pr_fst, map=gmap, sdp=NULL, quiet=FALSE, ...)
{
    require(qtl2)
    require(qtl2fst)
    pos <- find_markerpos(map, mar)

    if(is.null(sdp)) {
        this_fg <- t(cross$founder_geno[[pos$chr]][,mar,drop=FALSE])

        if(any(this_fg==0)) {
            # loop over all possible ways of filling in the data
            wh <- which(this_fg==0)
            bin <- broman:::binary.v(length(wh))*2+1
            result <- vector("list", ncol(bin))
            for(i in 1:ncol(bin)) {
                tmp <- this_fg
                tmp[tmp==0] <- bin[,i]
                sdp <- calc_sdp(tmp)
                if(!quiet) message("sdp = ", sdp, " (", i, " of ", ncol(bin), ")")
                result[[i]] <- infer_snp(marker, cross, genoprobs, gmap, sdp)
            }
            names(result) <- apply(bin, 2, paste, collapse="")
            return(result)
        }
        else {
            sdp <- calc_sdp(this_fg)
        }
    }

    snpinfo <- data.frame(chr=pos$chr, pos=pos$pos, snp=mar,
                          sdp=sdp)
    snpinfo <- index_snps(map, snpinfo)
    snp_pr <- genoprob_to_snpprob(genoprobs, snpinfo)
    snp_inf <- maxmarg(snp_pr)[[1]][,1]

    snp_inf
}


plot_snpint <-
    function(marker, snp_inf=NULL, use=c("observed", "inferred", "error", "geno_and_error", "genoinf_and_error"),
             cross=do, snp_int=snpint, genoprobs=pr_fst, map=gmap,
             error_lod=errlod, pointcolors=NULL, error_threshold=4, ...)
{
    use <- match.arg(use)
    if(use=="inferred" && is.null(snp_inf)) {
        snp_inf <- infer_snp(marker, cross=cross, genoprobs=genoprobs, map=map, quiet=FALSE)
    }

    intA <- unlist(snp_int[snp_int[,1]==marker & snp_int[,2]=="X", -(1:2)])
    intB <- unlist(snp_int[snp_int[,1]==marker & snp_int[,2]=="Y", -(1:2)])

    if(use=="inferred" || use=="genoinf_and_error") geno <- snp_inf
    else if(use=="observed" || use=="geno_and_error") geno <- do.call("cbind", cross$geno)[,marker]
    else if(use=="error") geno <- 1 + (error_lod[,mar] > error_threshold)

    if(use != "error") geno[geno==0 | is.na(geno)] <- 4

    if(is.null(pointcolors)) {
        if(use != "error") pointcolors <- broman::brocolors("f2")
        else pointcolors <- c(broman::brocolors("f2")[4], broman::brocolors("sex")[1])
    }

    if(use=="geno_and_error" || use=="genoinf_and_error") {
        err <- 1 + (errlod[,mar] > error_threshold)
        linewidth <- c(1,2)[err]
        linecolor <- c("black", broman::brocolors("sex")[1])[err]
    } else {
        linewidth <- 1
        linecolor <- "black"
    }


    id <- qtl2::get_common_ids(intA, intB, geno)

    this_plot <- function(x, y, pch=21, bg=pointcolors[geno[id]], col=linecolor, lwd=linewidth,
                          xlab="intensity A", ylab="intensity B", ...)
        grayplot(x, y, pch=pch, bg=bg, xlab=xlab, ylab=ylab, col=col, lwd=lwd, ...)

    this_plot(intA, intB, ...)
}
