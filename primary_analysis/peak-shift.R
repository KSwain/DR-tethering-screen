options(stringsAsFactors=FALSE)
library("stats4")

datadir <- "./work"

if (!exists("sgd")) {
    sgdfile <- sprintf("%s/SGD_features.tab", datadir)
    if (!file.exists(sgdfile)) {
        sgd <- download.file('https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab', destfile=sgdfile)
    }
    sgd <- read.delim(sgdfile, header=FALSE, quote="",
                      col.names=c("sgdid", "type", "qual", "name", "gene", "alias",
                                  "parent", "sgdid2", "chrom", "start", "end",
                                  "strand", "genpos", "cver", "sver", "desc"))
}

peakRun <- function(run) {
    cts <- read.csv(sprintf("%s/%s-assigned-counts.csv", datadir, run))

    cts$gene <- sgd[match(cts$yorf, sgd$name),"gene"]

    pseudo <- data.frame(FL = cts$FL / sum(cts$FL),
                         L = cts$L / sum(cts$L),
                         R = cts$R / sum(cts$R),
                         FR = cts$FR / sum(cts$FR))
    
    ## bounds are c(LOWER, UPPER)
    loglikeBin <- function(bounds, peak) {
        log(pnorm(bounds[[2]], mean=peak) - pnorm(bounds[[1]], mean=peak))
    }
    
    ## binEdges <- qnorm(p=c(0.25, 0.50, 0.75))
    binBounds = list(FR=qnorm(c(0.001, 0.225)),
                     R=qnorm(c(0.275,0.475)),
                     L=qnorm(c(0.525,0.725)),
                     FL=qnorm(c(0.775,0.999)))

    mlePeak <- function(pseudos, bounds) {
                                        #  peak0 <- weighted.mean(
                                        #    x=c(mean(binBounds$FL), mean(binBounds$L), mean(binBounds$R), mean(binBounds$FR)),
                                        #    w=c(sum(pseudos$FL), sum(pseudos$L), sum(pseudos$R), sum(pseudos$FR))
                                        #  )
        
        peak0 <- 0
        
        nll <- function(peak) {
            ll <- c(sum(pseudos$FL) * loglikeBin(bounds$FL, peak),
                    sum(pseudos$L)  * loglikeBin(bounds$L,  peak),
                    sum(pseudos$R)  * loglikeBin(bounds$R,  peak),
                    sum(pseudos$FR) * loglikeBin(bounds$FR, peak))
            ll[!is.finite(ll)] <- -36
            res <- -sum(1e8*ll, na.rm=TRUE)
            res
        }
        
        mle(minuslogl = nll,
            start=list(peak=peak0),
            method="L-BFGS-B",
            lower=-3, upper=3,
            control=(REPORT=1))
    }
    
    foo <- lapply(seq(1,nrow(pseudo)), function(i) { mlePeak(pseudo[i,], binBounds)})
    
    cts$mlePeak <- sapply(foo, coef)

    write.csv(x=cts, file=sprintf("%s/%s-barcode-mle-peak.csv", datadir, run))
    
    multiCts2 <- cts[duplicated(cts$frag),]
    multiCts2 <- multiCts2[!duplicated(multiCts2$frag),]
    multiCts1 <- cts[!duplicated(cts$frag),]
    multiCts1 <- multiCts1[match(multiCts2$frag, multiCts1$frag),]
    
    multi <- data.frame(barcode1 = multiCts1$barcode,
                        barcode2 = multiCts2$barcode,
                        frag = multiCts1$frag,
                        yorf = multiCts1$yorf,
                        gene = multiCts1$gene,
                        FL1 = multiCts1$FL,
                        L1 = multiCts1$L,
                        R1 = multiCts1$R,
                        FR1 = multiCts1$FR,
                        unsort1 = multiCts1$unsort,
                        FL2 = multiCts2$FL,
                        L2 = multiCts2$L,
                        R2 = multiCts2$R,
                        FR2 = multiCts2$FR,
                        unsort2 = multiCts2$unsort,
                        mlePeak1 = multiCts1$mlePeak,
                        mlePeak2 = multiCts2$mlePeak)
    
    ## Verify correct matching!
    table(cts[match(multi$barcode1, cts$barcode),"frag"] == cts[match(multi$barcode2, cts$barcode),"frag"])
    
    medianPeak <- aggregate(cts$mlePeak, by=list(frag=cts$frag), median)
    nbc <- aggregate(cts$mlePeak, by=list(frag=cts$frag), length)
    nread <- aggregate(cts$FL + cts$L + cts$R + cts$FR,
                       by=list(frag=cts$frag), sum)
    
    byfrag <- data.frame(frag=medianPeak$frag,
                         mlePeak = medianPeak$x,
                         nbc = nbc$x,
                         nread = nread$x)
    byfrag$yorf <- cts[match(byfrag$frag, cts$frag),"yorf"]
    byfrag$gene <- sgd[match(byfrag$yorf, sgd$name),"gene"]
    byfrag$desc <- sgd[match(byfrag$yorf, sgd$name),"desc"]
    
    write.csv(x=byfrag, file=sprintf("%s/%s-frag-mle-peak.csv", datadir, run))
    write.csv(x=byfrag[byfrag$nbc>1,], file=sprintf("%s/%s-fragmulti-mle-peak.csv", datadir, run))
}

peakRun("niks015")
peakRun("niks018")
