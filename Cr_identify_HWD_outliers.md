## Hardy-Weinberg Disequilibrium
```bash
# bash script
module load PLINK/2.00a6.9
plink2 --vcf Cr_filt_uniqID_noRel_exSNP_S194.vcf.gz --hardy --out Cr_filt_uniqID_noRel_exSNP_S194
plink2 --vcf Cr_filt_uniqID_noRel_exSNP_S194.vcf.gz --freq --out Cr_filt_uniqID_noRel_exSNP_S194

# R script
frq <- read.table("Cr_filt_uniqID_noRel_exSNP_S194.frq", skip = 1, header = F)
frq$MAF <- sapply(frq[,6], function(af) min(af, 1-af))
library(ggplot2)
ggplot(frq, aes(x = MAF)) +
  geom_histogram(binwidth = 0.005, fill = "steelblue", color = "black") +
  labs(title = "Frequency Distribution", x = "Frequency", y = "Count") +
  theme_minimal()
setwd("/path/to/Cr_raw_lcWGS_data/Cr_pop_analysis/stats")
pfrq <- read.table("Cr_filt_uniqID_noRel_exSNP_S194.afreq", skip = 1, header = F)
pfrq$MAF <- sapply(pfrq[,5], function(af) min(af, 1-af))
head(pfrq)

hardy <- read.table("Cr_filt_uniqID_noRel_exSNP_S194.hardy", skip = 1, header = F)
# Calculate the HWD ratio
hardy$HWD_ratio <- (hardy[,8] - hardy[,9]) 
#$/ hardy[,9]
head(hardy$HWD_ratio)

hwd_maf_df <- data.frame(x=pfrq$MAF , y = hardy$HWD_ratio)

library(ggplot2)
ggplot(hwd_maf_df, aes(x = x, y = y)) +
  geom_point(size = 1) +  theme_minimal() +
  geom_hline(yintercept = -0.05, color = "red", linetype = "dashed") + 
  labs(x = "MAF", y = "HWD")

writeLines(hardy[hwd_maf_df$y<=-0.05,2], "Cr_filt_uniqID_noRel_exSNP_S194_HWD_loci.txt")

# remove HWD
vcftools --gzvcf Cr_filt_uniqID_noRel_exSNP_S194.vcf.gz --exclude Cr_filt_uniqID_noRel_exSNP_S194_HWD_05_loci.txt --recode --recode-INFO-all --out Cr_filt_uniqID_noRel_exSNP_S194_HWD
```

## Outlier Selection
```bash

# PCAdapt

# LD pruned dataset

## Rangitāhua
genotypes <- "Cr_uniqID_noRel_filt_HWD_LD_RA.bed"
pca_matrix <- read.pcadapt(genotypes, type = "bed")
# choose the final K value for the rest of the analyses
x <- pcadapt(input = pca_matrix, K = 5, min.maf = 0.01)
plot(x)
plot(x, option = "screeplot")
pop_R <- c("21", "21", "22", "21", "22", "21", "21", "21", "21", "21", "22", "22", "22",
           "22", "22", "22", "22", "21", "21", "21", "22", "21", "21", "21", "21", "22", "22", "22")
pop_R <- c("21S109", "21S131", "22S132", "21S133", "22S142", "21S143", "21S16", "21S176", "21S177", "21S180",
  "22S191", "22S202", "22S206", "22S221", "22S22", "22S244", "22S270", "21S272", "21S30", "21S44", 
  "22S49", "21S50", "21S58", "21S62", "21S74", "22S7", "22S8", "22S92")
R <- plot(x, option = "scores")
R$layers <- R$layers[-1]  # remove default points
R+ geom_point(aes(shape = pop_R, colour = pop_R), size = 3) +
  scale_shape_manual(breaks = c("21", "22"), values = c("21" = 2, "22" = 3 )) + 
  scale_colour_manual( breaks = c("9", "10", "11", "12"), values = c("21" = "purple4","22" = "slateblue3")) +
theme(panel.grid = element_blank())



R = plot(x, option = "scores", pop = pop_R)
R + geom_text(aes(label = pop_R), size = 3, vjust = -1) 
plot(x, option = "scores", i = 3, j = 4, pop = pop_R)
summary(x)
plot(x , option = "manhattan", outliers = TRUE)
plot(x, option = "qqplot")
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 100, col = "red3")
plot(x, option = "stat.distribution")
x$gif
library(qvalue)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers <- which(qval < alpha)
length(outliers)
bim <- read.table("Cr_uniqID_noRel_filt_HWD_LD_RA.bim") 
outlier_snps_noR <- bim[outliers, 2]
write.table(outlier_snps_noR, file = "Cr_uniqID_noRel_filt_HWD_LD_RA_outliers_k5_a1.csv", sep = ",", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

## Australia
genotypes <- "Cr_uniqID_noRel_filt_HWD_LD_OZ.bed"
pca_matrix <- read.pcadapt(genotypes, type = "bed")
# choose the final K value for the rest of the analyses
x <- pcadapt(input = pca_matrix, K = 4, min.maf = 0.01)
plot(x)
plot(x, option = "screeplot")
pop_OZ <- c("3", "5", "4", "12", "11", "8", "11", "1", "2", "12", "8", "4", "10", "11", "9", "5", "7", "3", "10", "12", "4", "10", "7", 
            "13", "3", "4", "1", "3", "11", "5", "3", "1", "10", "12", "12", "3", "10", "9", "1", "10", "5", "7", "5", "6", "3", "3", 
            "1", "4", "12", "8", "5", "13", "11", "1", "8", "2", "1", "1", "2", "2", "13", "13", "1", "12", "7", "7", "7", "8", "13", 
            "9", "1", "1", "5", "6", "1", "9", "8", "12", "2", "8", "5", "11", "9", "7", "2", "6", "4", "9", "11", "4", "9", "4", "3", 
            "12", "7", "8", "4", "5", "5", "9", "13", "6", "10", "6", "2", "9", "7", "6", "7", "2", "11", "6", "4", "12", "12", "12", "12", 
            "12", "14", "14", "14", "14", "11", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "9", "9", "13", "3", "1", 
            "3", "3", "11", "9", "11", "8", "2", "7", "12", "3", "11", "3", "2", "11", "2", "9", "7", "2", "11", "10", "1", "2", "8", "10", "5", 
            "13", "13", "7", "9", "3", "2", "6", "6", "1", "5", "2", "6", "9")

pop_OZ <- c("3S100", "5S103", "4S104", "12S105", "11S106", "8S107", "11S108", "1S112", "2S115", "12S117", "8S11", "4S121", "10S122", "11S123", 
  "9S124", "5S126", "7S127", "3S12", "10S130", "12S138", "4S139", "10S13", "7S140", "13S141", "3S144", "4S146", "1S147", "3S149", 
  "11S14", "5S150", "3S151", "1S152", "10S153", "12S154", "12S155", "3S156", "10S158", "9S15", "1S161", "10S162", "5S164", "7S166", 
  "5S169", "6S170", "3S171", "3S172", "1S174", "4S175", "12S179", "8S182", "5S184", "13S185", "11S186", "1S187", "8S188", "2S190", 
  "1S195", "1S196", "2S197", "2S199", "13S19", "13S200", "1S205", "12S209", "7S20", "7S210", "7S211", "8S212", "13S213", "9S218", 
  "1S21", "1S220", "5S230", "6S231", "1S234", "9S236", "8S237", "12S238", "2S239", "8S240", "5S241", "11S242", "9S245", "7S246", 
  "2S247", "6S248", "4S24", "9S252", "11S253", "4S255", "9S257", "4S258", "3S25", "12S260", "7S262", "8S263", "4S265", "5S267", 
  "5S26", "9S271", "13S273", "6S274", "10S275", "6S276", "2S281", "9S282", "7S285", "6S286", "7S287", "2S288", "11S289", "6S28", 
  "4S290", "12S291", "12S292", "12S293", "12S294", "12S295", "14S296", "14S297", "14S298", "14S299", "11S2", "14S300", "14S301", 
  "14S302", "14S303", "14S304", "14S305", "14S306", "14S307", "14S308", "14S309", "14S310", "9S33", "9S34", "13S35", "3S37", "1S3", 
  "3S42", "3S43", "11S45", "9S47", "11S48", "8S4", "2S53", "7S5", "12S60", "3S63", "11S64", "3S65", "2S66", "11S68", "2S69", "9S6", 
  "7S70", "2S71", "11S72", "10S73", "1S75", "2S76", "8S78", "10S79", "5S81", "13S82", "13S84", "7S85", "9S86", "3S88", "2S90", "6S91", 
  "6S93", "1S94", "5S97", "2S98", "6S99", "9S9")
O <- plot(x, option = "scores")
O$layers <- O$layers[-1]
highlight <- c("4", "8", "10", "13")
fill_cols <- c("1"="#721E17", "2"="#95211B", "3"="#B8221E", "4"="#DA2222", "5"="#DF4828", "6"="#E4632D", "7"="#E67932", "8"="#E78C35",
               "9"="orange", "10"="goldenrod", "11"="goldenrod1", "12"="gold2", "13"="gold", "14"="yellow2")
outline_cols <- setNames(ifelse(names(fill_cols) %in% c("4", "8", "10", "13"), "black",fill_cols), names(fill_cols))
O + geom_point( aes(fill = pop_OZ, colour = pop_OZ),shape = 21, size = 3, stroke = 0.8) +
  scale_fill_manual(breaks = as.character(1:14), values = fill_cols, guide = guide_legend(override.aes = list(shape = 21,
  colour = ifelse(names(fill_cols) %in% highlight,"black", fill_cols), stroke = 0.8)))+
  scale_colour_manual(breaks = as.character(1:14), values = outline_cols) + theme(panel.grid = element_blank())

#to add labels to figure 
O <- plot(x, option = "scores", pop = pop_OZ)
O + geom_text(aes(label = pop_OZ), size = 3, vjust = -1)
plot(x, option = "scores", i = 5, j = 6, pop = pop_OZ)
summary(x)
x$gif
plot(x , option = "manhattan", qval = TRUE )
plot(x, option = "qqplot")
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 100, col = "orange")
plot(x, option = "stat.distribution")
library(qvalue)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers_OZ <- which(qval < alpha)
length(outliers_OZ)
bim <- read.table("Cr_uniqID_noRel_filt_HWD_LD_OZ.bim") 
outlier_snps_OZ <- bim[outliers_OZ, 2]
write.table(outlier_snps_OZ, file = "Cr_uniqID_noRel_filt_HWD_LD_OZ_outliers_k4_a1.csv", sep = ",", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

## Aotearoa
genotypes <- "Cr_uniqID_noRel_filt_HWD_LD_NZ.bed"
pca_matrix <- read.pcadapt(genotypes, type = "bed")
# choose the final K value for the rest of the analyses
x <- pcadapt(input = pca_matrix, K = 2, min.maf = 0.01)
plot(x)
plot(x, option = "screeplot")
pop_NZ <- c("17", "15", "18", "15", "17", "17", "15", "18", "16", "16", "19", "18", "20", "18", "18", "19", "20", "17", "20", "15", "19", "17", 
            "19", "19", "16", "17", "19", "15", "18", "20", "16", "19", "18", "20", "16", "19", "20", "18", "19", "20", "18", "15", "15", "19", 
            "18", "15", "20", "20", "16", "17", "16", "20", "15", "16", "16", "17", "16", "18", "17", "18", "17", "17", "19", "20", "15", "20", 
            "17", "19", "16", "20", "19", "20", "17", "18", "18", "15", "16", "15", "15", "16")

pop_NZ <- c("17S101", "15S102", "18S111", "15S113", "17S114", "17S118", "15S119", "18S120", "16S125", "16S128", "19S129", "18S134", "20S135",
  "18S136", "18S137", "19S145", "20S157", "17S160", "20S163", "15S167", "19S173", "17S178", "19S17", "19S181", "16S189", "17S18", 
  "19S193", "15S198", "18S1", "20S201", "16S203", "19S204", "18S208", "20S214", "16S216", "19S219", "20S222", "18S223", "19S224", 
  "20S226", "18S227", "15S228", "15S229", "19S23", "18S243", "15S249", "20S250", "20S251", "16S254", "17S256", "16S259", "20S264",
  "15S266", "16S268", "16S269", "17S277", "16S278", "18S279", "17S27", "18S280", "17S283", "17S29", "19S31", "20S32", "15S36",
  "20S38", "17S39", "19S40", "16S41", "20S46", "19S52", "20S54", "17S56", "18S57", "18S61", "15S67", "16S80", "15S83", "15S87", "16S96")

N <- plot(x, option = "scores")
N$layers <- N$layers[-1]  # remove default points
N+ geom_point(aes(shape = pop_NZ, colour = pop_NZ), size = 3) +
  scale_shape_manual(values = c(
    "15" = 0, "16" = 1, "17" = 2, "18" = 3 , "19" = 4, "20" = 5)) + 
  scale_colour_manual(values = c(
    "15" = "#549EB3", "16" = "#59A5A9", 
    "17" = "#60AB9E", "18" = "#69B190", "19" = "skyblue", "20"= "dodgerblue")) + 
  theme(panel.grid = element_blank())


N = plot(x, option = "scores", pop = pop_NZ)
N + geom_text(aes(label = pop_NZ), size = 3, vjust = -1)
plot(x, option = "scores", i = 3, j = 4, pop = pop_NZ)
summary(x)
x$gif
plot(x , option = "manhattan")
plot(x, option = "qqplot")
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 100, col = "green4", xaxt = "n")
ticks <- seq(0, 0.005, by = 0.001)
axis(side = 1, at = ticks)
plot(x, option = "stat.distribution")
library(qvalue)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers_NZ <- which(qval < alpha)
length(outliers_NZ)
bim <- read.table("Cr_uniqID_noRel_filt_HWD_LD_NZ.bim") 
outlier_snps_NZ <- bim[outliers_NZ, 2]
write.table(outlier_snps_NZ, file = "Cr_uniqID_noRel_filt_HWD_LD_NZ_outliers_k2_a1.csv", sep = ",", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# CSV files were combined and duplicates were removed

cd /path/to/Cr_raw_lcWGS_data/Cr_pop_analysis/LDprune
bcftools view -e 'ID=@Cr_combines_outliers.txt' Cr_uniqID_noRel_filt_HWD_LD_NZ.recode.vcf --output Cr_HWD_LD_noOuts_NZ.vcf.gz --output-type z
bcftools view -e 'ID=@Cr_combines_outliers.txt' Cr_uniqID_noRel_filt_HWD_LD_OZ.recode.vcf --output Cr_HWD_LD_noOuts_OZ.vcf.gz --output-type z
bcftools view -e 'ID=@Cr_combines_outliers.txt' Cr_uniqID_noRel_filt_HWD_LD_RA.recode.vcf --output Cr_HWD_LD_noOuts_RA.vcf.gz --output-type z
bcftools view -e 'ID=@Cr_combines_outliers.txt' Cr_uniqID_noRel_filt_HWD_LD.vcf.gz --output Cr_HWD_LD_noOuts_all.vcf.gz --output-type z

cd /path/to/Cr_raw_lcWGS_data/Cr_pop_analysis/thin
bcftools view -e 'ID=@Cr_thin_outliers.txt' Cr_uniqID_noRel_filt_HWD_th.recode.vcf --output Cr_uniqID_noRel_filt_HWD_th_noOuts_all.vcf.gz --output-type z
bcftools view -e 'ID=@Cr_thin_outliers.txt' Cr_uniqID_noRel_filt_HWD_th_NZ.recode.vcf --output Cr_uniqID_noRel_filt_HWD_th_noOuts_NZ.vcf.gz --output-type z
bcftools view -e 'ID=@Cr_thin_outliers.txt' Cr_uniqID_noRel_filt_HWD_th_OZ.recode.vcf --output Cr_uniqID_noRel_filt_HWD_th_noOuts_OZ.vcf.gz --output-type z
bcftools view -e 'ID=@Cr_thin_outliers.txt' Cr_uniqID_noRel_filt_HWD_th_RA.recode.vcf --output Cr_uniqID_noRel_filt_HWD_th_noOuts_RA.vcf.gz --output-type z
```
