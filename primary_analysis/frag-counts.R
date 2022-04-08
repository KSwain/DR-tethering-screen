options(stringsAsFactors=FALSE)

datadir <- "./work"
assigndir <- "./work"

seen <- read.delim(pipe(sprintf("cut -f1-4 %s/pacbio-190731-barcode-assign-all.txt", assigndir)), header=FALSE)
assign <- read.csv(sprintf("%s/pacbio-190731-facs-assign-yorf.txt", assigndir))

assignRun <- function(name) {
    ctall <- read.delim(sprintf("%s/%s.txt", datadir, name))
    colnames(ctall) <- sub(".count.txt$", "", colnames(ctall))
    colnames(ctall) <- sub("P[0-9]", "", colnames(ctall))
    colnames(ctall) <- sub(".*niks01[58].", "", colnames(ctall))

    ctall$assign <- seen[match(ctall$barcode, seen$V1), "V4"]
    ctall$assign <- ifelse(is.na(ctall$assign), "NotSeen", ctall$assign)
    
    write.csv(x=table(2**floor(log2(rowSums(ctall[,c("FL", "L", "R", "FR")]))),
                      ctall$assign),
              file=sprintf("%s/%s-assignment-stats.csv", datadir, name),
              quote=FALSE, row.names=TRUE)

    unassigned <- ctall[!(ctall$barcode %in% assign$barcode),]
    unassigned <- unassigned[order(rowSums(unassigned[,c("FL", "L", "R", "FR")]), decreasing=TRUE),]
    write.csv(x=unassigned, file=sprintf("%s/%s-unassigned.csv", datadir, name),
              quote=FALSE, row.names=FALSE)
    
    cts <- ctall[ctall$barcode %in% assign$barcode,]
    
    cts$frag <- assign[match(cts$barcode, assign$barcode),"frag"]
    cts$yorf <- assign[match(cts$barcode, assign$barcode),"yorf"]
    
    write.csv(cts, file=sprintf("%s/%s-assigned-counts-all.csv", datadir, name), quote=FALSE, row.names=FALSE)

    ## ctshi <- cts[cts$unsort >= 32,]
    ## unsort is not very good in NIKS018
    ctshi <- cts
    ctshi <- ctshi[rowSums(ctshi[,c("FL", "L", "R", "FR")]) >= 128,]
    write.csv(ctshi, file=sprintf("%s/%s-assigned-counts.csv", datadir, name), quote=FALSE, row.names=FALSE)
}

assignRun("niks015")
assignRun("niks018")
