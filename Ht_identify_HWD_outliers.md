## Hardy-Weinberg Disequilibrium
```bash
# bash script
module load PLINK/2.00a6.9
plink2 --vcf Ht_refilt_uniqID_noRel.recode.vcf --hardy --out Ht_refilt_uniqID_noRel
plink2 --vcf Ht_refilt_uniqID_noRel.recode.vcf --freq --out Ht_refilt_uniqID_noRel

# R script
#read in file
frq <- read.table("Ht_refilt_uniqID_noRel_2.frq", skip = 1, header = F)
frq$MAF <- sapply(frq[,6], function(af) min(af, 1-af))
library(ggplot2)
ggplot(frq, aes(x = MAF)) +
  geom_histogram(binwidth = 0.005, fill = "steelblue", color = "black") +
  labs(title = "Frequency Distribution", x = "Frequency", y = "Count") +
  theme_minimal()
setwd("/nesi/nobackup/ga03714/Ht_raw_lcWGS_data/Ht_pop_analysis/stats")
pfrq <- read.table("Ht_refilt_uniqID_noRel.afreq", skip = 1, header = F)
pfrq$MAF <- sapply(pfrq[,5], function(af) min(af, 1-af))
head(pfrq)

hardy <- read.table("Ht_refilt_uniqID_noRel.hardy", skip = 1, header = F)
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

writeLines(hardy[hwd_maf_df$y<=-0.05,2], "Ht_refilt_uniqID_noRel_HWD_loci.txt")
```

## Outlier Selection
````bash

# PCAdapt

# LD pruned dataset

## Rangitāhua
genotypes <- "Ht_filt_uniqID_noRel_HWD_LD_RA.bed"
pca_matrix <- read.pcadapt(genotypes, type = "bed")
# choose the final K value for the rest of the analyses
x <- pcadapt(input = pca_matrix, K = 5, min.maf = 0.01)
plot(x)
plot(x, option = "screeplot")
pop_R <- c("10", "12", "9", "10", "11", "10", "10", "9", "12", "10", "9", "9", "12", "11", "9", "9", "11", "11", 
           "9", "9", "10", "12", "10", "10", "12", "10", "10", "11", "11", "11", "9", "10", "11", "9", "11", "11",
           "12", "10", "12", "12", "9", "10", "9", "11", "11", "9", "9", "12", "12", "12", "12", "12", "10", "11", 
           "12", "9", "11", "10", "12")
R <- plot(x, option = "scores")
R$layers <- R$layers[-1]  # remove default points
R+ geom_point(aes(shape = pop_R, colour = pop_R), size = 3) +
  scale_shape_manual(breaks = c("9", "10", "11", "12"), values = c(
    "9" = 0, "10" = 1, "11" = 2, "12" = 3 )) + 
  scale_colour_manual( breaks = c("9", "10", "11", "12"), values = c(
    "9" = "#549EB3", "10" = "#59A5A9", "11" = "#60AB9E", "12" = "#69B190")) +
  theme(panel.grid = element_blank())

R = plot(x, option = "scores", pop = pop_R)
R + geom_text(aes(label = pop_R), size = 3, vjust = -1) 
plot(x, option = "scores", i = 6, j = 7, pop = pop_R)
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
bim <- read.table("Ht_filt_uniqID_noRel_HWD_LD_RA.bim") 
outlier_snps_noR <- bim[outliers, 2]
write.table(outlier_snps_noR, file = "Ht_refilt_uniqID_noRel_HWD_LD_RA_outliers_k5_a1.csv", sep = ",", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

## Australia
genotypes <- "Ht_filt_uniqID_noRel_HWD_LD_OZ.bed"
pca_matrix <- read.pcadapt(genotypes, type = "bed")
# choose the final K value for the rest of the analyses
x <- pcadapt(input = pca_matrix, K = 4, min.maf = 0.01)
plot(x)
plot(x, option = "screeplot")
pop_OZ <- c("5", "2", "1", "2", "4", "6", "5", "5", "3", "1", "4", "1", "1", "4", "1", "7", "3", "8", "8", "7", "8", 
            "4", "5", "6", "2", "3", "2", "6", "7", "6", "1", "2", "6", "7", "3", "6", "1", "7", "2", "5", "1", "5",
            "7", "6", "5", "8", "6", "8", "6", "1", "1", "2", "3", "3", "2", "3", "4", "2", "6", "6", "7", "1", "1", 
            "6", "2", "8", "2", "2", "6", "4", "4", "7", "8", "2", "8", "4", "3", "8", "1", "2", "8", "7", "7", "4", 
            "7", "4", "5", "5", "6", "3", "3", "5", "7", "3", "1", "7", "8", "2", "6", "5", "3", "7", "8", "8", "5", 
            "5", "7", "1", "3", "3", "4")



O <- plot(x, option = "scores")
O$layers <- O$layers[-1]  # remove default points
O+ geom_point(aes(shape = pop_OZ, colour = pop_OZ), size = 3) +
  scale_shape_manual(values = c(
    "1" = 0, "2" = 1, "3" = 2, "4" = 3 , "5" = 4, "6" = 5, "7" = 6 , "8" = 8 )) + 
  scale_colour_manual(values = c(
    "1" = "#721E17", "2" = "#95211B", "3" = "#B8221E", "4" = "#DA2222", 
    "5" = "#DF4828", "6" = "#E4632D", "7" = "#E67932", "8" = "#E78C35"))+
  theme(panel.grid = element_blank())

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
bim <- read.table("Ht_filt_uniqID_noRel_HWD_LD_OZ.bim") 
outlier_snps_OZ <- bim[outliers_OZ, 2]
write.table(outlier_snps_OZ, file = "Ht_filt_uniqID_noRel_HWD_LD_OZ_outliers_k4_a1.csv", sep = ",", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

## Aotearoa
genotypes <- "Ht_filt_uniqID_noRel_HWD_LD_NZ.bed"
pca_matrix <- read.pcadapt(genotypes, type = "bed")
# choose the final K value for the rest of the analyses
x <- pcadapt(input = pca_matrix, K = 2, min.maf = 0.01)
plot(x)
plot(x, option = "screeplot")
pop_NZ <- c("17", "16", "17", "13", "15", "18", "19", "18", "14", "13", "14", "19", "15", "17", "15", "19", "18", "13", 
            "14", "17", "18", "15", "16", "13", "18", "18", "17", "19", "13", "18", "13", "18", "19", "19", "14", "13", 
            "14", "13", "14", "18", "18", "17", "14", "16", "13", "17", "15", "17", "17", "13", "17", "14", "18", "17", 
            "14", "18", "13", "19", "14", "16", "15", "19", "14", "18", "17", "18", "13", "14", "18", "19", "14", "17", 
            "14", "13")

N <- plot(x, option = "scores")
N$layers <- N$layers[-1]  # remove default points
N+ geom_point(aes(shape = pop_NZ, colour = pop_NZ), size = 3) +
  scale_shape_manual(values = c(
    "13" = 0, "14" = 1, "15" = 2, "16" = 3 , "17" = 4, "18" = 5, "19" = 6)) + 
  scale_colour_manual(values = c(
    "13" = "midnightblue", "14" = "royalblue", "15" = "blue1", "16" = "dodgerblue", 
    "17" = "mediumblue", "18" = "darkslateblue", "19" = "skyblue"))+
  theme(panel.grid = element_blank())


N = plot(x, option = "scores", pop = pop_NZ)
N + geom_text(aes(label = pop_NZ), size = 3, vjust = -1)
plot(x, option = "scores", i = 5, j = 6, pop = pop_NZ)
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
bim <- read.table("Ht_filt_uniqID_noRel_HWD_LD_NZ.bim") 
outlier_snps_NZ <- bim[outliers_NZ, 2]
write.table(outlier_snps_NZ, file = "Ht_filt_uniqID_noRel_HWD_LD_NZ_outliers_k2_a1.csv", sep = ",", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# Thinned dataset 

## Rangitāhua
genotypes <- "Ht_filt_uniqID_noRel_HWD_th_RA.bed"
pca_matrix <- read.pcadapt(genotypes, type = "bed")
# choose the final K value for the rest of the analyses
x <- pcadapt(input = pca_matrix, K = 5, min.maf = 0.01)
plot(x)
plot(x, option = "screeplot")
pop_R <- c("10", "12", "9", "10", "11", "10", "10", "9", "12", "10", "9", "9", "12", "11", "9", "9", "11", "11", 
           "9", "9", "10", "12", "10", "10", "12", "10", "10", "11", "11", "11", "9", "10", "11", "9", "11", "11",
           "12", "10", "12", "12", "9", "10", "9", "11", "11", "9", "9", "12", "12", "12", "12", "12", "10", "11", 
           "12", "9", "11", "10", "12")
R <- plot(x, option = "scores")
R$layers <- R$layers[-1]  # remove default points
R+ geom_point(aes(shape = pop_R, colour = pop_R), size = 3) +
  scale_shape_manual(breaks = c("9", "10", "11", "12"), values = c(
    "9" = 0, "10" = 1, "11" = 2, "12" = 3 )) + 
  scale_colour_manual( breaks = c("9", "10", "11", "12"), values = c(
    "9" = "#549EB3", "10" = "#59A5A9", "11" = "#60AB9E", "12" = "#69B190")) 



R = plot(x, option = "scores", pop = pop_R)
R + geom_text(aes(label = pop_R), size = 3, vjust = -1) 
plot(x, option = "scores", i = 6, j = 7, pop = pop_R)
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
bim <- read.table("Ht_filt_uniqID_noRel_HWD_th_RA.bim") 
outlier_snps_noR <- bim[outliers, 2]
write.table(outlier_snps_noR, file = "Ht_filt_uniqID_noRel_HWD_th_RA_outliers_k5_a1.csv", sep = ",", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

## Australia
genotypes <- "Ht_filt_uniqID_noRel_HWD_th_OZ.bed"
pca_matrix <- read.pcadapt(genotypes, type = "bed")
# choose the final K value for the rest of the analyses
x <- pcadapt(input = pca_matrix, K = 4, min.maf = 0.01)
plot(x)
plot(x, option = "screeplot")
pop_OZ <- c("5", "2", "1", "2", "4", "6", "5", "5", "3", "1", "4", "1", "1", "4", "1", "7", "3", "8", "8", "7", "8", 
            "4", "5", "6", "2", "3", "2", "6", "7", "6", "1", "2", "6", "7", "3", "6", "1", "7", "2", "5", "1", "5",
            "7", "6", "5", "8", "6", "8", "6", "1", "1", "2", "3", "3", "2", "3", "4", "2", "6", "6", "7", "1", "1", 
            "6", "2", "8", "2", "2", "6", "4", "4", "7", "8", "2", "8", "4", "3", "8", "1", "2", "8", "7", "7", "4", 
            "7", "4", "5", "5", "6", "3", "3", "5", "7", "3", "1", "7", "8", "2", "6", "5", "3", "7", "8", "8", "5", 
            "5", "7", "1", "3", "3", "4")



O <- plot(x, option = "scores")
O$layers <- O$layers[-1]  # remove default points
O+ geom_point(aes(shape = pop_OZ, colour = pop_OZ), size = 3) +
  scale_shape_manual(values = c(
    "1" = 0, "2" = 1, "3" = 2, "4" = 3 , "5" = 4, "6" = 5, "7" = 6 , "8" = 8 )) + 
  scale_colour_manual(values = c(
    "1" = "#721E17", "2" = "#95211B", "3" = "#B8221E", "4" = "#DA2222", 
    "5" = "#DF4828", "6" = "#E4632D", "7" = "#E67932", "8" = "#E78C35")) 

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
bim <- read.table("Ht_filt_uniqID_noRel_HWD_th_OZ.bim") 
outlier_snps_OZ <- bim[outliers_OZ, 2]
write.table(outlier_snps_OZ, file = "Ht_filt_uniqID_noRel_HWD_th_OZ_outliers_k4_a1.csv", sep = ",", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

## Aotearoa
genotypes <- "Ht_filt_uniqID_noRel_HWD_th_NZ.bed"
pca_matrix <- read.pcadapt(genotypes, type = "bed")
# choose the final K value for the rest of the analyses
x <- pcadapt(input = pca_matrix, K = 2, min.maf = 0.01)
plot(x)
plot(x, option = "screeplot")
pop_NZ <- c("17", "16", "17", "13", "15", "18", "19", "18", "14", "13", "14", "19", "15", "17", "15", "19", "18", "13", 
            "14", "17", "18", "15", "16", "13", "18", "18", "17", "19", "13", "18", "13", "18", "19", "19", "14", "13", 
            "14", "13", "14", "18", "18", "17", "14", "16", "13", "17", "15", "17", "17", "13", "17", "14", "18", "17", 
            "14", "18", "13", "19", "14", "16", "15", "19", "14", "18", "17", "18", "13", "14", "18", "19", "14", "17", 
            "14", "13")

N <- plot(x, option = "scores")
N$layers <- N$layers[-1]  # remove default points
N+ geom_point(aes(shape = pop_NZ, colour = pop_NZ), size = 3) +
  scale_shape_manual(values = c(
    "13" = 0, "14" = 1, "15" = 2, "16" = 3 , "17" = 4, "18" = 5, "19" = 6)) + 
  scale_colour_manual(values = c(
    "13" = "midnightblue", "14" = "royalblue", "15" = "blue1", "16" = "dodgerblue", 
    "17" = "mediumblue", "18" = "darkslateblue", "19" = "skyblue")) 


N = plot(x, option = "scores", pop = pop_NZ)
N + geom_text(aes(label = pop_NZ), size = 3, vjust = -1)
plot(x, option = "scores", i = 4, j = 5, pop = pop_NZ)
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
bim <- read.table("Ht_filt_uniqID_noRel_HWD_th_NZ.bim") 
outlier_snps_NZ <- bim[outliers_NZ, 2]
write.table(outlier_snps_NZ, file = "Ht_filt_uniqID_noRel_HWD_th_NZ_outliers_k2_a1.csv", sep = ",", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# CSV files were combined and duplicates were removed

module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1 PLINK/2.00a6.9 BCFtools/1.22-GCC-12.3.0 Python/3.11.6-foss-2023a

cd /path/to/Ht_raw_lcWGS_data/Ht_pop_analysis/LDprune_data
bcftools view -e 'ID=@outlier_combined.txt' Ht_filt_uniqID_noRel_HWD_LD_NZ.recode.vcf --output Ht_HWD_LD_noOuts_NZ.vcf.gz --output-type z
bcftools view -e 'ID=@outlier_combined.txt' Ht_filt_uniqID_noRel_HWD_LD_OZ.recode.vcf --output Ht_HWD_LD_noOuts_OZ.vcf.gz --output-type z
bcftools view -e 'ID=@outlier_combined.txt' Ht_filt_uniqID_noRel_HWD_LD_RA.recode.vcf --output Ht_HWD_LD_noOuts_RA.vcf.gz --output-type z
bcftools view -e 'ID=@outlier_combined.txt' Ht_filt_uniqID_noRel_HWD_LD.vcf.gz --output Ht_HWD_LD_noOuts_all.vcf.gz --output-type z

bcftools view -e 'ID=@thin_outliers.txt' Ht_filt_uniqID_noRel_HWD_th.recode.vcf --output Ht_filt_uniqID_noRel_HWD_th_noOuts.vcf.gz --output-type z
bcftools view -e 'ID=@thin_outliers.txt' Ht_filt_uniqID_noRel_HWD_th_NZ.recode.vcf --output Ht_filt_uniqID_noRel_HWD_th_noOuts_NZ.vcf.gz --output-type z
bcftools view -e 'ID=@thin_outliers.txt' Ht_filt_uniqID_noRel_HWD_th_OZ.recode.vcf --output Ht_filt_uniqID_noRel_HWD_th_noOuts_OZ.vcf.gz --output-type z
bcftools view -e 'ID=@thin_outliers.txt' Ht_filt_uniqID_noRel_HWD_th_RA.recode.vcf --output Ht_filt_uniqID_noRel_HWD_th_noOuts_RA.vcf.gz --output-type z
```
