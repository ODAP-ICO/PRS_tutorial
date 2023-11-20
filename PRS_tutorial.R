# chage to cloned folder

getwd()

system("git clone https://github.com/ODAP-ICO/PRS_tutorial")

setwd("PRS_tutorial")

library(readxl)
library(data.table)
library(ggplot2)

# load data
CRC_SNP <- read_excel("data/resources/1_CRC_SNPs_list_205_Nature_2023.xlsx")

# The first 5 rows and first 7 columns are shown
CRC_SNP[c(1:5),c(1:7)]

# .txt file generation
SNP_file_CRC <- CRC_SNP[,c(1,2)]
head(SNP_file_CRC)

# two avoid scientific annotation
SNP_file_CRC$`POS (GRCh37)` <- format(SNP_file_CRC$`POS (GRCh37)`, scientific=FALSE)

# create folder data/1_snps_extraction/ and save as .txt
dir.create("data/1_snps_extraction/", recursive=TRUE)
write.table(SNP_file_CRC, "data/1_snps_extraction/1_SNP_file_CRC.txt", quote = F, col.names = F, row.names = F, qmethod = "double", sep = "\t")

# bash

system(" plink2 --vcf data/resources/GenRisk_all_CRC_SNPs_perm.vcf.gz --new-id-max-allele-len 50 --set-all-var-ids @:#_\\$r_\\$a --freq --make-pgen --out data/1_snps_extraction/raw_data" )

# re-load original CRC SNPs list.
CRC_SNP <- read_excel("data/resources/1_CRC_SNPs_list_205_Nature_2023.xlsx")

# Load the data of the SNPs extracted from the target data:
pvar<-fread("data/1_snps_extraction/raw_data.pvar")

colnames(pvar)<- c("CHR",   "POS",  "ID"    ,"REF", "ALT",  "FILTER",   "INFO")
head(pvar)

# Match variant format
CRC_SNP$id_pos <- paste0(CRC_SNP$CHR,":",CRC_SNP$`POS (GRCh37)`)
pvar$id_pos <- paste0(pvar$CHR,":",pvar$POS)


# SNPs not found
snpsNOTpvar<-CRC_SNP[!(CRC_SNP$id_pos %in% pvar$id_pos), ]
snpsNOTpvar


table(duplicated(pvar$id_pos))

multi_all <- pvar[duplicated(pvar$POS) | duplicated(pvar$POS, fromLast = TRUE), ]
multi_all


# the format of the variable ID between tables is equalized. From chr:position_other_allele/effect_allele to chr:position_other_allele_effect_allele
CRC_SNP$ID <- gsub("/","_",CRC_SNP$VARIANT)

#Which is the correct?
id_multi <- multi_all[!multi_all$ID %in% CRC_SNP$ID,] 
id_multi


#Save the incorrect SNP to delete in the QC step:
## Structure: CHROM  ID      REF     ALT
## ID = chr:pos_REF_ALT

id_multi_to_remove <- id_multi[,c(1,3:5)]

dir.create("data/2_SNPs_alignament_check")

write.table(id_multi_to_remove,"data/2_SNPs_alignament_check/multi_to_remove.txt", col.names = F, row.names = F, quote = F)

## SNPs that have been extracted: 
dim(CRC_SNP[CRC_SNP$ID %in% pvar$ID,])

## SNPs that have not been extracted: 
no_snps <- CRC_SNP[!CRC_SNP$ID %in% pvar$ID,]
no_snps

# load SNPs frequencyes  
freq<-read.table("data/1_snps_extraction/raw_data.afreq",header=F,sep="\t")
colnames(freq)<- c("CHROM", "ID",   "REF",  "ALT",  "ALT_FREQS",    "OBS_CT")
head(freq)

ambiguous_snps <- freq[freq$ALT_FREQS > 0.4 & freq$ALT_FREQS< 0.6 & 
                        ((freq$REF == "A" & freq$ALT == "T") | 
                         (freq$ALT == "T" & freq$REF == "A") | 
                         (freq$REF == "C" & freq$ALT == "G") | 
                         (freq$ALT == "G" & freq$REF == "C")), ]

ambiguous_snps

# Save
write.table(ambiguous_snps, "data/2_SNPs_alignament_check/exclrsIDs_ambiguous.txt", sep = "\t", quote = FALSE, row.names = FALSE)


head(pvar)

# Create column R2
pvar$INFO2 <- sub("ER2=.*", "", pvar$INFO)
pvar$R2 <- as.numeric(sub('.*R2=([0-9.]+);.*', '\\1', pvar$INFO2))

r2_04_remove <- pvar[pvar$R2 < 0.4,]
r2_04_remove


dir.create("data/3_QC")

# remove redundant SNPs 
system("plink2 --pfile data/1_snps_extraction/raw_data --indep-pairwise 200 1 0.3 --out data/3_QC/LD")

ld_snps <- read.table("data/3_QC/LD.prune.out",header=F,sep="\t")
ld_snps <- ld_snps$V1
length(ld_snps) 
ld_snps

system("plink2 --pfile data/1_snps_extraction/raw_data --hwe 1e-6 --write-snplist --out data/3_QC/HWE_QC")

multi_remove <-read.table("data/2_SNPs_alignament_check/multi_to_remove.txt", colClasses = c("character"))
colnames(multi_remove) <- c("CHROM", "ID", "REF", "ALT")
multi_remove

# The ID columns are selected
ambiguous_snps_ID <- ambiguous_snps$ID
r2_04_remove_ID <- r2_04_remove$ID
multi_remove_ID  <- multi_remove$ID


table(ambiguous_snps_ID == multi_remove_ID)
table(ambiguous_snps_ID %in% r2_04_remove_ID)
table(multi_remove_ID %in% r2_04_remove_ID)

table(ld_snps %in% multi_remove_ID)
table(ld_snps %in% r2_04_remove_ID)
table(ld_snps %in% ambiguous_snps_ID)

to_remove <- unique(c(ambiguous_snps_ID,multi_remove_ID,r2_04_remove_ID,ld_snps ))
to_remove

# Save
write.table(to_remove, "data/3_QC/to_remove.txt", quote = FALSE, row.names = FALSE, col.names = F)


system("plink2 --pfile data/1_snps_extraction/raw_data --exclude data/3_QC/to_remove.txt --maf 0.005 --write-snplist --make-pgen --out data/3_QC/dataQC")


# Load Info SNPs resultants del QC
snps_qc<-read.table("data/3_QC/dataQC.pvar")
colnames(snps_qc) <- c("CHROM", "POS",  "ID", "REF", "ALT", "FILTER", "INFO")

#Transform ID to: chr:position_all_1:all_2 format
#snps_qc$ID <- gsub("_",":",snps_qc$ID)

#Keep ID and ALT / RISK_ALLELE
names(snps_qc)[names(snps_qc) == "ALT"] <- "RISK_ALLELE"
names(snps_qc)[names(snps_qc) == "REF"] <- "OTHER_ALLELE"
snps_qc <- snps_qc[,c("ID", "RISK_ALLELE","OTHER_ALLELE" )]

head(snps_qc) 

# Load beta values.
betas <- read_excel("data/resources/1_CRC_SNPs_list_205_Nature_2023.xlsx")

#Transform ID to: chr:position_all_1:all_2 format -> 1:110365045_A_G
names(betas)[names(betas) == "VARIANT"] <- "ID"
names(betas)[names(betas) == "BETA"] <- "RISK_SCORE"
#betas$ID <- gsub("_",":",betas$ID)
betas$ID <- gsub("/","_",betas$ID)

betas <- betas[,c(3,7)]

head(betas)


table(snps_qc$ID %in% betas$ID)
score <- merge(snps_qc, betas, by = "ID", all.x = TRUE)
head(score)

score$RISK_ALLELE_2 <- score$RISK_ALLELE
head(score)

#Change Risk_allele <-> Other_allele in the negative beta values cases:
for (i in 1:nrow(score)){
      
    # iterate over the column nº4 ("RISK_SCORE")  and check if <0
    if(score[i,4]<0){
           
        # replace the value with OTHER_ALLELE (column 3) in the new column RISK_ALLEL_2 (5)
        score[i,5]<- score[i,3]
        }
}

#Change beta value sign - to -> +
score$RISK_SCORE_2 <- as.numeric(gsub("-", "",score$RISK_SCORE ))


#Delete some columns
score[,2:4] <- NULL

#Save the new score doc as score_final.txt 

dir.create("data/4_PRS_calculation",recursive = TRUE)

write.table(score,"data/4_PRS_calculation/score_final.txt", col.names = F, quote = F, row.names = F)

system("plink2 --pfile data/1_snps_extraction/raw_data --extract data/3_QC/dataQC.snplist --score data/4_PRS_calculation/score_final.txt no-mean-imputation --out data/4_PRS_calculation/PRS")

# Load results:
scores <- read.table("data/4_PRS_calculation/PRS.sscore")
colnames(scores)<- c("IID", "NMISS_ALLELE_CT",  "ALLELE_DOSAGE_SUM",    "SCORE1_AVG")

#Rescale the scores
scores$SCORE1_AVG_scale <- scale(scores$SCORE1_AVG, center = T, scale = T)
head(scores)

# Samples to extract
samples<-read.table("data/resources/GenRisk_ImputedData_CRC_Samples.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t")

head(samples)

table(samples$tipo_cancer)

#Unir columnes samples i PRS.score per les ID de pacients
scores_samples<-merge(scores,samples,by.x="IID",by.y="ID",all=FALSE)

head(scores_samples)

ttest<-t.test(scores_samples$SCORE1_AVG_scale~scores_samples$tipo_cancer)
tpval<- t.test(scores_samples$SCORE1_AVG_scale~scores_samples$tipo_cancer)$p.value
tpval

# Histogram:
histo <- ggplot(scores_samples, aes(x = SCORE1_AVG_scale))+ 
  geom_histogram(aes(y=..density.. ,fill = tipo_cancer, color = tipo_cancer), alpha = 0.5, position = "identity") +
  scale_fill_manual( values = c("steelblue4", "gray87")) +
   scale_color_manual(values = c("steelblue4", "gray50"))+
   theme_classic()+
  xlab("Polygenic risk score (PRS)") + 
  ylab("Density")+
 
  #labs(fill="xyz")+
   labs(title = "Colorectal cancer", subtitle = "p <0.0001 ")+
   theme(text = element_text(size = 16),
         legend.position = "bottom") +
theme(legend.title=element_blank())

histo

# Density plot
density <- ggplot(scores_samples, aes(x = SCORE1_AVG_scale, color = tipo_cancer))+ 
  geom_density()+
  scale_color_manual(values=c("#99B1C4", "gray50"))+
   labs(title = "Risk score distribution",
              #subtitle = "Plot of length by dose",
              caption = paste0("t-test pval=",tpval))+
   theme_minimal()+
  theme(text = element_text(size = 20)) 

density

# Boxplot
boxplot <- ggplot(scores_samples, aes(x=tipo_cancer, y=SCORE1_AVG_scale, fill=tipo_cancer)) +
  geom_boxplot()+
  labs(x="", y = "Risk Score")+
  scale_x_discrete(limits=c("Control", "Colorectal"))+
  scale_fill_manual(values=c("#99B1C4", "gray87"))+
  theme_classic()+
  theme(text = element_text(size = 20), legend.position = "none") 

boxplot

# qqplot
qwe1<-qqnorm(scores_samples$SCORE1_AVG_scale[scores_samples$tipo_cancer =="Control"],plot.it=FALSE)
qwe2<-qqnorm(scores_samples$SCORE1_AVG_scale[scores_samples$tipo_cancer=="Colorectal"],plot.it=FALSE)

plot(qwe1$x,qwe1$y,pch=19,las=1,main="",col="darkgrey",xlab="Theoretical Quantiles",ylab="Sample Quantiles")
qqline(scores_samples$SCORE1_AVG_scale[scores_samples$tipo_cancer== "Control"],col="darkgrey",lwd=2)
points(qwe2$x,qwe2$y,pch=19,col="#99B1C4")
qqline(scores_samples$SCORE1_AVG_scale[scores_samples$tipo_cancer== "Colorectal"],col="#99B1C4",lwd=2)
title("Risk score Q-Q Plot",sub=paste0("t-test pval=",tpval))
