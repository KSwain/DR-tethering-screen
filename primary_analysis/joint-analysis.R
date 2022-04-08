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

bc15 <- read.csv(sprintf("%s/niks015-barcode-mle-peak.csv", datadir), row.names=1)
bc18 <- read.csv(sprintf("%s/niks018-barcode-mle-peak.csv", datadir), row.names=1)

bc15merge <- bc15
bc15merge$frag <- NULL
bc15merge$yorf <- NULL
bc15merge$gene <- NULL
bcjoint <- merge(bc15merge, bc18, by=c("barcode"), suffixes=c(".15", ".18"))
bcall <- merge(bc15, bc18, by=c("barcode"), suffixes=c(".15", ".18"), all=TRUE)
bcall$frag <- ifelse(is.na(bcall$frag.15), bcall$frag.18, bcall$frag.15)
bcall$yorf <- ifelse(is.na(bcall$yorf.15), bcall$yorf.18, bcall$yorf.15)
bcall$gene <- ifelse(is.na(bcall$gene.15), bcall$gene.18, bcall$gene.15)
bcall$frag.15 <- NULL
bcall$yorf.15 <- NULL
bcall$gene.15 <- NULL
bcall$frag.18 <- NULL
bcall$yorf.18 <- NULL
bcall$gene.18 <- NULL

write.csv(bcall, sprintf("%s/joint-barcode-mle-peak.csv", datadir), row.names=FALSE)
write.csv(bcjoint, sprintf("%s/joint-both-barcode-mle-peak.csv", datadir), row.names=FALSE)

f15 <- read.csv(sprintf("%s/niks015-frag-mle-peak.csv", datadir), row.names=1)
f18 <- read.csv(sprintf("%s/niks018-frag-mle-peak.csv", datadir), row.names=1)

f15merge <- f15
f15merge$yorf <- NULL
f15merge$gene <- NULL
f15merge$desc <- NULL
fjoint <- merge(f15merge, f18, by=c("frag"), suffixes=c(".15", ".18"))
fall <- merge(f15, f18, by=c("frag"), suffixes=c(".15", ".18"), all=TRUE)
fall$yorf <- ifelse(is.na(fall$yorf.15), fall$yorf.18, fall$yorf.15)
fall$yorf.15 <- NULL
fall$yorf.18 <- NULL
fall$gene <- ifelse(is.na(fall$gene.15), fall$gene.18, fall$gene.15)
fall$gene.15 <- NULL
fall$gene.18 <- NULL
fall$desc <- ifelse(is.na(fall$desc.15), fall$desc.18, fall$desc.15)
fall$desc.15 <- NULL
fall$desc.18 <- NULL

fall$mlePeak <- ifelse(abs(fall$mlePeak.15) < abs(fall$mlePeak.18),
                       fall$mlePeak.15, fall$mlePeak.18)

write.csv(fall, sprintf("%s/joint-frag-mle-peak.csv", datadir), row.names=FALSE)
write.csv(fjoint, sprintf("%s/joint-both-frag-mle-peak.csv", datadir), row.names=FALSE)


