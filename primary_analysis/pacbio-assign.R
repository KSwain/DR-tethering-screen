options(stringsAsFactors=FALSE)

datadir <- "~/Tethering/PacBio190731"

for (sample in c("pacbio-190731-facs")) {
    bed <- read.delim(sprintf("%s/%s-assign-gene.bed", datadir, sample), header=FALSE)
    bed$frag <- sprintf("%s:%d-%d(%s)", bed$V1, bed$V2, bed$V3, bed$V6)
    fragcts <- table(bed$frag)
    fragdist <- table(2**floor(log2(fragcts)))
    
    write.table(fragdist, file=sprintf("%s/%s-bc-per-frag.txt", datadir, sample),
                quote=FALSE, sep="\t", row.names=FALSE,
                col.names=FALSE)
    
    bedFragUniq <- bed[!duplicated(bed$frag),]
    yorfcts <- table(bedFragUniq$V10)
    yorfdist <- table(2**floor(log2(yorfcts)))

    write.table(yorfdist, file=sprintf("%s/%s-frag-per-yorf.txt", datadir, sample),
                quote=FALSE, sep="\t", row.names=FALSE,
                col.names=FALSE)
    
    assign <- data.frame(barcode=bed$V4,
                         ref=bed$V1,
                         start=bed$V2,
                         end=bed$V3,
                         str=bed$V6,
                         frag=bed$frag,
                         yorf=bed$V10,
                         cts=bed$V5)
    
    write.csv(assign, file=sprintf("%s/%s-assign-yorf.txt", datadir, sample),
              quote=FALSE, row.names=FALSE)
}

