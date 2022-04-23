options(stringsAsFactors=FALSE)
library(dplyr)

##read in barcode count and activity data from two screen reps
joint_frag <- read.csv("../primary_analysis/work-cached/joint-frag-mle-peak.csv", header = TRUE, sep=",")

head(joint_frag)

joint_frag$nbc.15[is.na(joint_frag$nbc.15)] <- 0
joint_frag$nbc.18[is.na(joint_frag$nbc.18)] <- 0
joint_frag$nread.15[is.na(joint_frag$nread.15)] <- 0
joint_frag$nread.18[is.na(joint_frag$nread.18)] <- 0
head(joint_frag)

plot(log10(joint_frag$nread.15), log10(joint_frag$nread.18))
cor(joint_frag$mlePeak.15, joint_frag$mlePeak.18)
dim(joint_frag)

## create df with second screen data only
mle18 <- joint_frag[,c("frag", "mlePeak.18", "nbc.18", "nread.18", "yorf", "gene")]
head(mle18)
dim(mle18)
library(tidyr)

mle18 <- mle18 %>% drop_na(mlePeak.18)
mle15 <- joint_frag[,c("frag", "mlePeak.15", "nbc.15", "nread.15", "yorf", "gene")]
head(mle15)
dim(mle15)
mle15 <- mle15 %>% drop_na(mlePeak.15)

mle18_kde <- density(mle18$mlePeak.18)
plot(mle18_kde)
head(mle18)

mle15_kde <- density(mle15$mlePeak.15)
plot(mle15_kde)

## set read cutoff to 500 or more                                               
present_18 <- subset(mle18, mle18$nread.18 > 499)
# #switch mle peak values so repressors are negative and activators are positive 
# present_18$mlePeak.18 <- present_18$mlePeak.18 * -1
mle18$mlePeak.18 <- mle18$mlePeak.18 * -1
mle15$mlePeak.15 <- mle15$mlePeak.15 * -1

dim(present_18)
present18_kde <- density(present_18$mlePeak.18)
plot(present18_kde)

Ded1_18 <- subset(mle18, mle18$gene=="DED1", select=c(1:6))
dim(Ded1_18)
Ded1_18

plot(Ded1_18$mlePeak.18, Ded1_18$nread.18)
plot(Ded1_18$mlePeak.18, log2(Ded1_18$nbc.18))
rug(Ded1_18$mlePeak.18, ticksize = 0.05, col = "blue", lwd=1.5)

library(ggplot2)
library(ggrepel) 

density_Ded1 <- density(Ded1_18$mlePeak.18)
plot(density_Ded1)
rug(Ded1_18$mlePeak.18, ticksize = 0.1, col = "blue",lwd=1.25)

ded_means <- Ded1_18 %>% summarize(mean=mean(mlePeak.18))
ded_means

png(file="output/density_ded1only.png",width=1200,height=1000,res=288)
plot_ded <- ggplot(Ded1_18, aes(x=mlePeak.18)) +
  geom_density(color="turquoise3", fill="lightcyan", size=1) +
  geom_vline(aes(xintercept=mean(mlePeak.18)),
             color="turquoise3", linetype="dashed", size=0.75) +
  geom_rug(size=1, alpha=0.6, color="turquoise3", length = unit(0.075, "npc")) + theme_classic() +
  xlim(-2.5,2.5) +
  theme(axis.title = element_blank(), axis.text = element_text(size=12),
        legend.position = c(0.15, 0.9), legend.text = element_text(size=12), legend.title = element_blank()) + 
  geom_text(data = ded_means, aes(x=mean, label = "1.07"), color="turquoise3", y = 0.65, vjust=-1, hjust=1.2, size=5)
plot_ded
dev.off()


## CCR4
Ccr4_18 <- subset(mle18, mle18$gene=="CCR4", select=c(1:6))
Ccr4_18
plot(Ccr4_18$mlePeak.18, Ccr4_18$nbc.18)
density_ccr4 <- density(Ccr4_18$mlePeak.18)
plot(density_ccr4)
rug(Ccr4_18$mlePeak.18, ticksize = 0.05, col = "red", lwd = 1.5)
## use ggplot for more controllable plot: 

ccr_means <- Ccr4_18 %>% summarize(mean=mean(mlePeak.18))
ccr_means
png(file="output/density_ccr4only.png",width=1200,height=1000,res=288)
plot_ccr <- ggplot(Ccr4_18, aes(x=mlePeak.18)) +
  geom_density(color="indianred1", fill="mistyrose1", size=1) +
  geom_vline(aes(xintercept=mean(mlePeak.18)),
             color="indianred1", linetype="dashed", size=0.75) +
  geom_rug(size=1, alpha=0.6, color="indianred1", length = unit(0.075, "npc")) + theme_classic() +
  xlim(-2.5,2.5) +
  theme(axis.title = element_blank(), axis.text = element_text(size=12),
        legend.position = c(0.15, 0.9), legend.text = element_text(size=12), legend.title = element_blank()) + 
  geom_text(data = ccr_means, aes(x=mean, label = "-0.46"), color="indianred2", y = 0.55, vjust=-1, hjust=1.1, size=5)
plot_ccr
dev.off()

histogram_18 <- ggplot(mle18, aes(x=mlePeak.18)) +
  geom_density(color="darkslategray3", fill = "darkslategray2") +
  geom_histogram(aes(y=..density..), colour="gray35", fill="white", bins = 60, alpha=0.4) +
  theme_classic()
histogram_18

#plot all fragments in screen according to activity score 
png(file="output/NIKS018_mlepeakdistribution.png",width=1200,height=1000,res=288)
all_18 <- ggplot(present_18, aes(x=mlePeak.18)) +
  geom_density(color="gray35", fill = "darkslategray2", size=0.75) +
  geom_vline(aes(xintercept=mean(mlePeak.18)),
             color="gray30", linetype="dashed", size=0.75) + theme_classic() +
  theme(axis.title = element_blank(), axis.text = element_text(size=12))
all_18
dev.off()


## kde plots for just Ccr4 v Ded1                                               
genes <- rbind(Ccr4_18, Ded1_18)
genes

# Find the mean of each group                                                   
library(plyr)
cdat <- ddply(genes, "gene", summarise, mle.mean=mean(mlePeak.18))
cdat

png(file="output/Ded1vCcr4.png",width=1200,height=1000,res=288)
plot_genes <- ggplot(genes, aes(x=mlePeak.18, color=gene)) +
  geom_density(size=0.75) +
  geom_vline(data=cdat, aes(xintercept=mle.mean, colour=gene),
             linetype="dashed", size=0.5, alpha=0.8) +
  geom_rug(size=1, alpha=0.6, length = unit(0.1, "npc")) + theme_classic() +
  xlim(-2.5,2.5) +
  theme(axis.title = element_blank(), axis.text = element_text(size=12),
        legend.position = c(0.15, 0.9), legend.text = element_text(size=12), legend.title = element_blank())
plot_genes
dev.off()


#### plot fragment size distribution

pacbio_frag <- read.csv("output/fragment_chrom_coord.csv", header = TRUE, sep= ",")
head(pacbio_frag)
pacbio_frag$X = NULL
pacbio_frag$length = pacbio_frag$start - pacbio_frag$end
pacbio_frag$length <- abs(pacbio_frag$length)

pacbio_means <- pacbio_frag %>% summarize(mean=mean(length))
pacbio_means$mean <- round(pacbio_means$mean)
pacbio_means


png(file="output/PacbiogDNA_size-distribution.png",width=600,height=500,res=144)
length_df <- ggplot(pacbio_frag, aes(x=length)) +
  geom_density(color="gray25", fill = "lightcoral", size=0.75) +
  geom_vline(aes(xintercept=mean(length)),
             color="gray30", linetype="dashed", size=0.75) + theme_classic() +
  geom_text(data = pacbio_means, aes(x=mean, label = mean), y = 0.004, vjust=-1, hjust=1.2, size=5) +
  xlim(-0,1000) + ylim(0,0.0045) +
  theme(axis.title = element_blank(), axis.text = element_text(size=12))
plot(length_df)
dev.off()

### working with identified domains:                                            

all_domains <- read.csv("output/all_domains.csv", header=TRUE, sep=",")
head(all_domains)
dim(all_domains)

valids <- all_domains[all_domains$pvalue <= 0.05,]
valids$padj <- p.adjust(valids$pvalue, method="BH")

## table of all domains with significant p values and p adjusted values         
valids
valids$fullname = c("Adenylosuccinate synthase", "Initiation factor 2 subunit family", "Endo/Exonuclease/phosphatase",
                    "Catalase", "FAD dependent oxidoreductase", "	Electron transfer flavoprotein", "Lysophospholipase",
                    "Amino-transferase class IV", "Glyceraldehyde 3-phosphate dehydrogenase", "BAR", "ANTH", "Proteasome", 
                    "Mannitol dehydrogenase C-terminal", "Betaglucosidase (SUN)", "Family of unknown function", 
                    "Hyphally regulated cell wall protein", "Aminotransferase class I and II", "DEAD/DEAH box helicase",
                    "RhoGAP", "Spb1 C-terminal", "tRNA synthetases class I catalytic","Glycolipid 2-α-mannosyltransferase",
                    "AICARFT/IMPCHase bienzyme", "RasGEF", "Formyl transferase", "Malic enzyme, N-terminal", 
                    "Eukaryotic rRNA processing protein EBP2", "DUF2722", "eIF3 subunit")
valids$fullname <- factor(valids$fullname, levels = valids$fullname[order(-valids$mlePeak)])
valids

png(file="output/domaincount_padj_hires.png", width=2400,height=2400,res=250)
valids %>%
  ggplot(aes(x=mlePeak,
             y=fullname,
             colour=padj,
             size=count)) + xlim(-1.0, 1.5) +
  geom_point() + labs(x="Mean activity score", y="Protein domain", colour="Adjusted P value", size="Count") +
  theme_minimal() + scale_color_gradient(low="olivedrab2", high="maroon1") +
  theme(axis.title= element_blank(), axis.text = element_text(size=14), 
        legend.text = element_text(size=12), legend.title = element_text(size=12)) +
  geom_vline(xintercept = 0, size=0.5, color="grey76", linetype="dashed")
dev.off()


## Load in gene annotation file
if (!file.exists("SGD_features.tab")) {
  sgd <- download.file('https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab', destfile="SGD_features.tab")
}

sgd <- read.delim("SGD_features.tab", header=FALSE, quote="",
                  col.names=c("sgdid", "type", "qual", "name", "gene", "alias",
                              "parent", "sgdid2", "chrom", "start", "end",
                              "strand", "genpos", "cver", "sver", "desc"))
head(sgd)
# 

joint_barcode <- read.csv("../primary_analysis/work-cached/joint-barcode-mle-peak.csv", header=TRUE, sep=",")
dim(joint_barcode)
head(joint_barcode)
# 
complete18 <- joint_barcode[!is.na(joint_barcode$mlePeak.18),]
dim(complete18)
add15 <- joint_barcode[is.na(joint_barcode$mlePeak.18),]
dim(add15)
# 
head(complete18)
complete18 <- complete18[,c("barcode", "mlePeak.18", "yorf", "gene")]
names(complete18)[2] <- "mlePeak"
complete18$origin <- "18"
# 
add15 <- add15[,c("barcode", "mlePeak.15", "yorf", "gene")]
head(add15)
names(add15)[2] <- "mlePeak"
add15$origin <- "15"
# 
mlescores<- rbind(complete18, add15)
dim(mlescores)
head(mlescores)

## load in non-collapsed by 90% identity fragments: 
frags_mle <- read.csv("output/frags_coordinates_mlep.csv", header = TRUE, sep = ",")
dim(frags_mle)
head(frags_mle)

minus <- subset(frags_mle, frags_mle$strand_x=="-", select=c(1:21))
plus <- subset(frags_mle, frags_mle$strand_x=="+", select=c(1:21))
head(plus)

dim(plus)
dim(minus)

head(minus)
names(minus)[3] <- 'fragEnd'
names(minus)[4] <- 'fragStart'
names(minus)[8] <- 'geneEnd'
names(minus)[9] <- 'geneStart'
## correct for proteins encoded from sequences that are oriented in the minus direction
minus$proteinStart <- (minus$geneStart - (minus$fragStart -1))/3
minus$proteinEnd <- (minus$geneStart - (minus$fragEnd - 1))/3
minus$proteinStart <- round(minus$proteinStart)
minus$proteinEnd <- round(minus$proteinEnd)

head(minus) 

dim(minus)   
dim(plus)
minus <- minus[,c(1:7,9,8,10:21)]


head(minus)
head(plus)

strandcorrect_db <- rbind(plus,minus)
dim(strandcorrect_db)
head(strandcorrect_db)

## save new database with corrected protein coordinates for - strand: 
strandcorrect_db <- strandcorrect_db[order(strandcorrect_db$chrom, strandcorrect_db$geneStart),]
head(strandcorrect_db)
strandcorrect_db$mlePeak.15 <- strandcorrect_db$mlePeak.15 * -1
strandcorrect_db$X = NULL

write.csv(strandcorrect_db, file= "output/fragment_protein_coordinates-minusstrandcorrect.csv")

## strandcorrect_db = coordinates of all fragments in our in-frame database
head(strandcorrect_db)
fragment_db <- strandcorrect_db[, c("barcode", "proteinStart", "proteinEnd", "mlePeak.15", "mlePeak.18", "yorf", "gene", "overlap")]
head(fragment_db)

fragment_db$UTRcorrect <- ifelse(fragment_db$proteinStart <0, 0, fragment_db$proteinStart)
fragment_db$protlength <- fragment_db$proteinEnd - fragment_db$proteinStart
fragment_db$UTRcorrect = NULL

fragment_db$desc <- sgd[match(fragment_db$yorf, sgd$name), "desc"]
fragment_db <- fragment_db[,c("barcode", "gene", "yorf","proteinStart", "proteinEnd", "protlength", "mlePeak.18", 
                              "mlePeak.15", "desc", "overlap")]
head(fragment_db)
dim(fragment_db)

cor(fragment_db$mlePeak.15, fragment_db$mlePeak.18)
write.csv(fragment_db, file= "output/fragment_db_all-minuscorrect.csv")

mlevFC <- read.csv("data/mlevFC.csv", header = TRUE, sep=",")
mlevFC
mlevFC$fragment = factor(mlevFC$fragment, levels=mlevFC[order(mlevFC$mlePeak18), "fragment"])
mlevFC$logFC = log2(mlevFC$FC)
cor.test(mlevFC$mlePeak18, mlevFC$FC)

#Pearson's product-moment correlation

#data:  mlevFC$mlePeak18 and mlevFC$FC
#t = 6.8474, df = 10, p-value = 4.473e-05
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.6972354 0.9741883
#sample estimates:
 #     cor 
#0.9078613 

mytheme <- theme(axis.line = element_line(size = 1.0, colour = "gray30"),
                 panel.background = element_rect(fill = "white"))


png(file="output/mlevFCscore.png", width=1800,height=1800,res=432)
mFCdots <- ggplot(mlevFC, aes(x=mlePeak18, y=logFC)) + 
  geom_point(aes(colour = factor(fragment)), size=2.5) + 
  geom_smooth(method=lm, se=FALSE, color="gray70", linetype="dashed", size=0.75)  + theme_minimal() +
  ylim(-1.5, 1.5) +
  geom_hline(yintercept=0, size=0.5, color = "grey70") + 
  geom_vline(xintercept = 0, size=0.5, color="grey70")+ 
  theme(axis.title=element_blank(), axis.text=element_text(size=14), 
        legend.text = element_text(size=12), legend.title = element_blank(), 
        legend.key.size = unit(1.5, 'lines')) + 
  scale_colour_manual(values = c("firebrick4", "firebrick3", "red", "darkgoldenrod1", "gold1", "khaki",
                                 "greenyellow","limegreen", "lightseagreen", "steelblue2",  "royalblue", "mediumblue"))
mFCdots
dev.off()


## looking at 90% sim 

head(fragment_db)
Ded1_all <- subset(fragment_db, fragment_db$gene=="DED1")
Ded1_all <- Ded1_all[,c(1,4,5,7,8)]
head(Ded1_all)


png(file="output/Ded1_dotplot.png", width=1200,height=1000,res=288)
ded1plot <- ggplot(Ded1_all, aes(x=proteinEnd, y=mlePeak.18)) + geom_point(size=1, color="seagreen3") + 
    geom_smooth(method=lm, se=TRUE, color="palevioletred2", linetype="dashed", size=1, fill="mistyrose1") +
  xlim(0,110) + theme_minimal() +
  theme(axis.title=element_blank(), axis.text=element_text(size=14))
ded1plot
dev.off()

head(fragment_db)
Ccr4_all <- subset(fragment_db, fragment_db$gene=="CCR4", select=c(1,4,5,6,7,8,10))
dim(Ccr4_all)
Ccr4_all <- Ccr4_all[order(Ccr4_all$proteinStart),]

Ccr4_all

C_nterm <- Ccr4_all[1:6,]
C_nterm
Cn_means <- C_nterm %>% summarize(mean=mean(mlePeak.18))
Cn_means

C_mid <- Ccr4_all[7:17,]
C_mid
Cmid_means <- C_mid %>% summarize(mean=mean(mlePeak.18))
Cmid_means

C_llr <- Ccr4_all[18:25,]
C_llr
Cllr_means <- C_llr %>% summarize(mean=mean(mlePeak.18))
Cllr_means

plot(density(C_nterm$mlePeak.18))
plot(density(C_mid$mlePeak.18))
plot(density(C_llr$mlePeak.18))


png(file="output/density_ccr4nterm.png",width=1200,height=1000,res=288)
plot_n <- ggplot(C_nterm, aes(x=mlePeak.18)) +
  geom_density(color="tomato2", size=1) +
  geom_vline(aes(xintercept=mean(mlePeak.18)),
             color="tomato2", linetype="dashed", size=0.75) +
  geom_rug(size=1, alpha=0.6, color="tomato2", length = unit(0.075, "npc")) + theme_classic() +
  xlim(-2.5,2.5) + ylim(0,0.65) +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        legend.position = c(0.15, 0.9), legend.text = element_text(size=12), legend.title = element_blank()) + 
  geom_text(data=Cn_means, aes(x=mean, label = "-0.9"), color="tomato2", y = 0.63, vjust=0, hjust=-.3, size=5)
plot_n
dev.off()

png(file="output/density_ccr4mid.png",width=1200,height=1000,res=288)
plot_mid <- ggplot(C_mid, aes(x=mlePeak.18)) +
  geom_density(color="sienna2", size=1) +
  geom_vline(aes(xintercept=mean(mlePeak.18)),
             color="sienna2", linetype="dashed", size=0.75) +
  geom_rug(size=1, alpha=0.6, color="sienna2", length = unit(0.075, "npc")) + theme_classic() +
  xlim(-2.5,2.5) + ylim(0,0.75) +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        legend.position = c(0.15, 0.9), legend.text = element_text(size=12), legend.title = element_blank()) + 
  geom_text(data = Cmid_means, aes(x=mean, label = "-0.6"), color="sienna2", y = 0.72, vjust=0, hjust=-.3, size=5)
plot_mid
dev.off()

png(file="output/density_ccr4llr.png",width=1200,height=1000,res=288)
plot_llr <- ggplot(C_llr, aes(x=mlePeak.18)) +
  geom_density(color="tan1", size=1) +
  geom_vline(aes(xintercept=mean(mlePeak.18)),
             color="tan1", linetype="dashed", size=0.75) +
  geom_rug(size=1, alpha=0.6, color="tan1", length = unit(0.075, "npc")) + theme_classic() +
  xlim(-1,1) +
  theme(axis.title = element_blank(), axis.text = element_blank()) + 
  geom_text(data = Cllr_means, aes(x=mean, label = "-0.1"), color="tan1", y = 7.35, vjust=0, hjust=1.1, size=5)
plot_llr
dev.off()

neg_motif <- read.csv("output/neg_mlep_motif_genes.csv", header=TRUE, sep=",")
head(neg_motif)

neg_motif$gene <- sgd[match(neg_motif$sequence_name, sgd$name), "gene"]


pos_motif <- read.csv("output/pos_mlep_motif_genes.csv", header = TRUE, sep=",")
head(pos_motif)
pos_motif <- subset(pos_motif,
                    pos_motif$q.value < 0.08)

## remove Ubiquitin and Y helicase motifs #2 and #6: 
pos_motif<- subset(pos_motif, 
                   pos_motif$motif_id!="6")
pos_motif<- subset(pos_motif, 
                   pos_motif$motif_id!="2")
## get rid of y helicase motifs (did multiple) and GAG retrotransposon motifs
pos_motif <- subset(pos_motif,
                    pos_motif$matched_sequence!="SQKPIIYKVHRDNNNLSPVQNEQKSWNKTQKKSNKVYNSK")
pos_motif <- subset(pos_motif,
                    pos_motif$matched_sequence!="SQKPIIYKVHRDNNHLSPVQNEQKSWNKTQKRSNKVYNSK")

act_motif <- pos_motif
act_motif$motif_id[which(act_motif$motif_id == "3")] <- "2"
act_motif$motif_id[which(act_motif$motif_id == "4")] <- "3"
act_motif$motif_id[which(act_motif$motif_id == "5")] <- "4"
act_motif$motif_id[which(act_motif$motif_id == "7")] <- "5"
act_motif$motif_id[which(act_motif$motif_id == "8")] <- "6"
act_motif$X = NULL
write.csv(act_motif, "~/screen_manuscript/screenII-only/activator_motifs.csv")

v = nrow(subset(act_motif, 
                act_motif$motif_id==6))
v
# motif 1 = 157
# motif 2 = 563
# motif 3 = 38
# motif 4 = 9
# motif 5 = 63
# motif 6 = 618

## repressors motifs:
neg_motif <- read.csv("output/neg_mlep_motif_genes.csv", header = TRUE, sep=",")
head(neg_motif)
neg_motif <- subset(neg_motif,
                    neg_motif$q.value < 0.08)

neg_motif <- subset(neg_motif,
                    neg_motif$matched_sequence!="YPQQRMMTPQQ")
rep_motif <- neg_motif
rep_motif$X = NULL
write.csv(rep_motif, "~/screen_manuscript/screenII-only/repressor_motifs.csv")


v = nrow(subset(rep_motif, 
                rep_motif$motif_id==6))
v
# motif 1 = 106
# motif 2 = 524
# motif 3 = 1
# motif 4 = 3
# motif 5 = 9
# motif 6 = 483

## interaction enrichment table
interaction_enrich <- read.csv("~/NIKS018/interaction-enrichment_new.csv", header = TRUE, sep=",")
head(interaction_enrich)

interact_table <- read.csv("~/NIKS018/interacting-activators-table.csv", header = TRUE, sep=",")
interact_list <- read.csv("~/NIKS018/interacting-activators-list.csv", header = TRUE, sep=",")
interact_list

head(interaction_enrich)
interaction_enrich$gene <- factor(interaction_enrich$gene, levels = interaction_enrich$gene[order(interaction_enrich$foldEnrich)])

png(file="~/NIKS018/Interactions_enriched.png", width=1600,height=1800,res=250)
interaction_enrich %>%
  ggplot(aes(x=foldEnrich,
             y=gene,
             colour=qEnrich,
             size=kActs)) + #xlim(-1.0, 1.5) +
  geom_point() + labs(x="Interactions with activator fragments", y="Gene", colour="Adjusted q-value", size="Interactors") +
  theme_minimal() + scale_color_gradient(low="maroon1", high="blue") +
  theme(axis.title= element_blank(), axis.text = element_text(size=18), 
        legend.text = element_text(size=14), legend.title = element_text(size=16))
dev.off()

barcode_counts <- read.csv("~/NIKS018/niks018-barcode-mle-peak.csv", header = TRUE, sep=",")
head(barcode_counts)

barcode_counts <- barcode_counts[,c("barcode", "yorf", "gene", "frag", "FL", "L", "R", "FR", "unsort", "mlePeak")]
head(barcode_counts)
barcode_counts$mlePeak <- barcode_counts$mlePeak * -1
## counts from 18
tif35 <- barcode_counts[barcode_counts$barcode == "GATAACACTTGCAGACCAATCTTTT",]
ebs1 <- barcode_counts[barcode_counts$barcode == "GGATACTAGCGCTCTGGTCATCTTA",]
ngr1 	<- barcode_counts[barcode_counts$barcode == "CGCCTTTATTGGGTTTATCTGTGGG",]
sro9 	<- barcode_counts[barcode_counts$barcode == "ACCGGGTTTTAAATTTTCTAAATCA",]
ded1 	<- barcode_counts[barcode_counts$barcode == "AAGTTCGGTCTGCACTAACCGTAAT",]
whi3 	<- barcode_counts[barcode_counts$barcode == "TTTCTGATTGACGCTAAGCTTACTG",]
yap1801	<- barcode_counts[barcode_counts$barcode == "GGTGCACCATCCCGGACTATCGTG",]
her1 	<- barcode_counts[barcode_counts$barcode == "GCGGGATCCACTGGGCGGCGAGATA",]
jsn1 	<- barcode_counts[barcode_counts$barcode == "CGGTGCGTATTTCGGCTCGGTGTGC",]
smy2 <- barcode_counts[barcode_counts$barcode == "GGGGTAGGTACGAATAGACCCGTGG",]

# counts from 15
barcode_counts15 <- read.csv("~/NIKS018/niks015-barcode-mle-peak.csv", header = TRUE, sep=",")
head(barcode_counts15)
dim(barcode_counts15)

barcode_counts15 <- barcode_counts15[,c("barcode", "yorf", "gene", "frag", "FL", "L", "R", "FR", "unsort", "mlePeak")]
head(barcode_counts15)
barcode_counts15$mlePeak <- barcode_counts15$mlePeak * -1
hsp26 <- barcode_counts[barcode_counts$barcode == "TGTGTTGTGGGTAGTTATTGAATAG",]
sbp1 	<- barcode_counts15[barcode_counts15$barcode == "TCTTCACTTACGGTGGAGTGGAATT",]

validated <- rbind(tif35, ebs1, ngr1, sro9, ded1, hsp26, jsn1, whi3, sbp1, yap1801, smy2, her1)
colnames(validated) <- c("barcode", "Yorf", "gene", "fragment start", "far left", "left", "right", "far right", "unsorted", "activity score")
validated 
write.csv(validated, "~/NIKS018/validated-fragments_AS.csv")



