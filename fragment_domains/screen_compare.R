options(stringsAsFactors=FALSE)

## naming fragment activity score - accidentally used barcode but it's only the fragments here
screens <- read.csv("output/frags_coordinates_mlep.csv", header = TRUE, sep=",")

## adjust for left sort = repressors  
screens$mlePeak.15 <- screens$mlePeak.15 * -1

head(screens)
dim(screens)
mean(screens$nread.15)
mean(screens$nread.18)

plot(screens$mlePeak.15, screens$mlePeak.18)
cor(screens$mlePeak.15, screens$mlePeak.18)
cor.test(screens$mlePeak.15, screens$mlePeak.18)


library(ggplot2)
library(ggrepel)

filter_sc <- subset(screens, screens$nread.15 > 1100 & screens$nread.18 > 1100)
plot(filter_sc$mlePeak.15, filter_sc$mlePeak.18)
cor(filter_sc$mlePeak.15, filter_sc$mlePeak.18)
cor.test(filter_sc$mlePeak.15, filter_sc$mlePeak.18)


png(file="output/libraryactscores_screencompare.png", width=1600,height=1600,res=250)
filter_plot <- ggplot(filter_sc, aes(x=mlePeak.15, y=mlePeak.18)) + geom_point(size=0.75, color="deepskyblue3") + 
  geom_smooth(method=lm, se=TRUE, color="gray20", linetype="dashed", size=0.75, fill="gray90") +
  theme_minimal() + theme(axis.title = element_text(size=14), axis.text= element_text(size=12)) + 
  labs(x="Library activity scores: screen rep I", y="Library activity scores: screen rep II")
filter_plot
dev.off()

#data rep I

barcode_counts15 <- read.csv("../primary_analysis/work-cached/niks015-barcode-mle-peak.csv", header = TRUE, sep=",", row.names=1)
head(barcode_counts15)
dim(barcode_counts15)

#data rep II

barcode_counts <- read.csv("../primary_analysis/work-cached/niks018-barcode-mle-peak.csv", header = TRUE, sep=",", row.names=1)
head(barcode_counts)

## counts from 18
tif35 <- barcode_counts[barcode_counts$barcode == "GATAACACTTGCAGACCAATCTTTT",]
ebs1 <- barcode_counts[barcode_counts$barcode == "GGGGTGATAATGGCCGTGAGGCCCA",]
ngr1 	<- barcode_counts[barcode_counts$barcode == "ACTGTATAAAGTAGGGAGTTAACAT",]
sro9 	<- barcode_counts[barcode_counts$barcode == "ACCGGGTTTTAAATTTTCTAAATCA",]
ded1 	<- barcode_counts[barcode_counts$barcode == "AAGTTCGGTCTGCACTAACCGTAAT",]
whi3 	<- barcode_counts[barcode_counts$barcode == "TTTCTGATTGACGCTAAGCTTACTG",]
yap1801	<- barcode_counts[barcode_counts$barcode == "GGTGCACCATCCCGGACTATCGTG",]
her1 	<- barcode_counts[barcode_counts$barcode == "GCGGGATCCACTGGGCGGCGAGATA",]
jsn1 	<- barcode_counts[barcode_counts$barcode == "CGGTGCGTATTTCGGCTCGGTGTGC",]
smy2 <- barcode_counts[barcode_counts$barcode == "GGGGTAGGTACGAATAGACCCGTGG",]
hsp26 <- barcode_counts[barcode_counts$barcode == "TGTGTTGTGGGTAGTTATTGAATAG",]
sbp1 	<- barcode_counts[barcode_counts$barcode == "TCTTCACTTACGGTGGAGTGGAATT",]

validated18 <- rbind(tif35, ebs1, ngr1, sro9, ded1, hsp26, jsn1, whi3, sbp1, yap1801, smy2, her1)
colnames(validated18) <- c("barcode", "R", "FL", "L", "FR", "unsort", "status", "frag", "yorf", "gene", "mlePeak18")
validated18
dim(validated18)

## counts from 15
tif35b <- barcode_counts15[barcode_counts15$barcode == "GATAACACTTGCAGACCAATCTTTT",]
ebs1b <- barcode_counts15[barcode_counts15$barcode == "GGGGTGATAATGGCCGTGAGGCCCA",]
ngr1b 	<- barcode_counts15[barcode_counts15$barcode == "ACTGTATAAAGTAGGGAGTTAACAT",]
sro9b 	<- barcode_counts15[barcode_counts15$barcode == "ACCGGGTTTTAAATTTTCTAAATCA",]
ded1b 	<- barcode_counts15[barcode_counts15$barcode == "AAGTTCGGTCTGCACTAACCGTAAT",]
whi3b 	<- barcode_counts15[barcode_counts15$barcode == "TTTCTGATTGACGCTAAGCTTACTG",]
yap1801b	<- barcode_counts15[barcode_counts15$barcode == "GGTGCACCATCCCGGACTATCGTG",]
her1b 	<- barcode_counts15[barcode_counts15$barcode == "GCGGGATCCACTGGGCGGCGAGATA",]
jsn1b 	<- barcode_counts15[barcode_counts15$barcode == "CGGTGCGTATTTCGGCTCGGTGTGC",]
smy2b <- barcode_counts15[barcode_counts15$barcode == "GGGGTAGGTACGAATAGACCCGTGG",]
hsp26b <- barcode_counts15[barcode_counts15$barcode == "TGTGTTGTGGGTAGTTATTGAATAG",]
sbp1b 	<- barcode_counts15[barcode_counts15$barcode == "TCTTCACTTACGGTGGAGTGGAATT",]


validated15 <- rbind(tif35b, ebs1b, ngr1b, sro9b, ded1b, hsp26b, jsn1b, whi3b, sbp1b, yap1801b, smy2b, her1b)
colnames(validated15) <- c("barcode", "R", "FL", "L", "FR", "unsort", "status", "frag", "yorf", "gene2", "mlePeak15")

dim(validated15)

validated18
validated15
validated <- merge(validated15, validated18[,c("barcode", "mlePeak18")], by="barcode")

head(validated)
plot(validated$mlePeak15, validated$mlePeak18)
cor(validated$mlePeak15, validated$mlePeak18)

cor.test(validated$mlePeak15, validated$mlePeak18)

# actscore <- validated[,c(2,8,9,10,11,22)]
actscore <- validated[,c("barcode", "frag", "yorf", "gene2", "mlePeak15", "mlePeak18")]

actscore

## Fragment names
fragnames <- data.frame(gene=c("TIF35", "EBS1", "NGR1", "SRO9", "DED1", "HSP26",
                               "JSN1", "WHI3", "SBP1", "YAP1801", "SMY2", "HER1"),
                        fragment=c("Tif35(115-239)", "Ebs1(691-876)", "Ngr1(474-626)", "Sro9(14-151)", "Ded1(1-108)", "Hsp26(1-182)",
                       "Jsn1(144-295)", "Whi3(518-662)", "Sbp1(14-178)", "Yap1801(374-527)", "Smy2(65-232)", "Her1(1025-1142)"))
actscore <- merge(actscore, fragnames, by.x="gene2", by.y="gene")

actscore$fragment = factor(actscore$fragment, levels=actscore[order(actscore$mlePeak18), "fragment"])

png(file="output/validated_screencompare.png", width=1800,height=1600,res=250)
valdots <- ggplot(actscore, aes(x=mlePeak15, y=mlePeak18)) + 
  geom_point(aes(colour = fragment), size=2.5) + 
  geom_smooth(method=lm, se=FALSE, color="gray70", linetype="dashed", size=0.75)  + theme_minimal() +
  geom_hline(yintercept=0, size=0.5, color = "grey70") + 
  geom_vline(xintercept = 0, size=0.5, color="grey70")+ 
  theme(axis.title = element_text(size=14), axis.text= element_text(size=12)) + 
  labs(x="Screen rep I activity score", y="Screen rep II activity score") +
  theme(legend.text = element_text(size=12), legend.title = element_blank(), 
        legend.key.size = unit(1.5, 'lines')) + 
  scale_colour_manual(values = c("firebrick4", "firebrick3", "red", "darkgoldenrod1", "gold1", "khaki",
                                 "greenyellow","limegreen", "lightseagreen", "steelblue2",  "royalblue", "mediumblue"))
valdots
dev.off()



