require(data.table)
require(R.utils)
require(MungeSumstats)
require(readr)
require(dplyr)

######################################
##############PART 1##################
##European Continuous Traits Example###
#####################################

################################
###########Yengo 2022: Height###
################################
#downloaded file from: https://www.joelhirschhornlab.org/giant-consortium-results
Height<-fread("GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR",data.table=FALSE)

#convert to MAF
Height$MAF<-ifelse(Height$EFFECT_ALLELE_FREQ > 0.5, 1-Height$EFFECT_ALLELE_FREQ, Height$EFFECT_ALLELE_FREQ)

#subset to MAF > 1%
Height<-subset(Height, Height$MAF > 0.01)

#make SNPID null due to duplicate column that will be interested as rsid
Height$SNPID<-NULL

#output new file
write.table(Height, file = "Yengo_Height_EUR.txt", col.names=TRUE,quote = FALSE, row.names = FALSE)


################################
###########Pan UKB: BMI########
################################
#downloaded the file form terminal using wget
#wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-21001-both_sexes-irnt.tsv.bgz

#unzip the file
BMI <- gunzip("continuous-21001-both_sexes-irnt.tsv.bgz", "BMI_PanUKB.tsv")

#read in the tsv file
BMI <- read_tsv("BMI_PanUKB.tsv")

#convert to dataframe
BMI<-data.frame(BMI)

#convert allele frequency to MAF 
BMI$MAF<-ifelse(BMI$af_EUR > 0.5, 1-BMI$af_EUR, BMI$af_EUR)

#subset fo MAF > 1%
BMI<-subset(BMI, BMI$MAF > .01)

#subset to only european data
attach(BMI)
BMI<-data.frame(chr,pos,ref,alt,beta_EUR,se_EUR,pval_EUR,MAF)

#unzip and read in the rsID file from panUKB
rsid<-gunzip("full_variant_qc_metrics.txt.bgz", "panuKBrsID.txt")
RSid<-fread("panuKBrsID.txt",data.table=FALSE)

#save just the chromosome, base pair position, and rsID for merging with BMI
attach(RSid)
RSid<-data.frame(chrom,pos,rsid)

#rename column to match naming in BMI and merge on chromosome and base pair
colnames(RSid)[1]<-"chr"
BMI<-inner_join(RSid,BMI,by=c("chr","pos"))

#now we have a file in the format that standard LDSC expects
head(BMI)

#except that the _EUR flag will not be recognized, so let's also rename those columns
colnames(BMI)[6:8]<-c("beta","SE","pval")

#ready to go right? not so fast: pvalue column is not in format that is expected
#they report ln(p) 
#can convert back to regular pvalue by either calculating p on the Z-stat
#or backing out from the natural log, either works just as well
#BMI$pval<-2*pnorm(abs(BMI$beta)/BMI$SE,lower.tail=FALSE)
BMI$pval<-exp(1)^BMI$pval

#now let's output this file to be munged
write.table(BMI, file = "panUKB_BMI_EUR.txt", col.names=TRUE,quote = FALSE, row.names = FALSE)


#now for the practical just going to subset to chromosome 1
#NOTE THAT THIS PART YOU WOULD NOT DO IN PRACTICE
Height<-fread("Yengo_Height_EUR.txt",data.table=FALSE)
Height<-subset(Height, Height$CHR == 1)
write.table(Height, file = "Yengo_Height_EUR_chr1.txt", col.names=TRUE,quote = FALSE, row.names = FALSE)

BMI<-fread("panUKB_BMI_EUR.txt",data.table=FALSE)
BMI<-subset(BMI, BMI$chr == 1)
write.table(BMI, file = "panUKB_BMI_EUR_chr1.txt", col.names=TRUE,quote = FALSE, row.names = FALSE)



#########################
######PART 2##############
##East Asian Example###
##########################

################################
###########biobank Japan: BMI###
################################
#downloaded file from: https://pheweb.jp/
BMI_Japan<-fread("GWASsummary_BMI_Japanese_SakaueKanai2020.auto.txt.gz",data.table=FALSE)

#convert to MAF
BMI_Japan$MAF<-ifelse(BMI_Japan$A1FREQ > 0.5, 1-BMI_Japan$A1FREQ, BMI_Japan$A1FREQ)

#subset to MAF > 1%
BMI_Japan<-subset(BMI_Japan, BMI_Japan$MAF > 0.01)

attach(BMI_Japan)
BMI_Japan<-data.frame(SNP,CHR,BP,ALLELE1,ALLELE0,INFO,MAF,BETA,SE,P_BOLT_LMM_INF)

colnames(BMI_Japan)[10]<-"PVAL"
#output new file
write.table(BMI_Japan, file = "BMI_Japan.txt", col.names=TRUE,quote = FALSE, row.names = FALSE)
rm(BMI_Japan)

################################
###########Yengo 2022: Height########
################################
#downloaded file from: https://www.joelhirschhornlab.org/giant-consortium-results
Height<-fread("GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EAS.gz",data.table=FALSE)

#convert to MAF
Height$MAF<-ifelse(Height$EFFECT_ALLELE_FREQ > 0.5, 1-Height$EFFECT_ALLELE_FREQ, Height$EFFECT_ALLELE_FREQ)

#subset to MAF > 1%
Height<-subset(Height, Height$MAF > 0.01)

#make SNPID null due to duplicate column that will be interested as rsid
Height$SNPID<-NULL

#output new file
write.table(Height, file = "Yengo_Height_EAS.txt", col.names=TRUE,quote = FALSE, row.names = FALSE)
rm(Height)

#now for the practical just going to subset to chromosome 1
#NOTE THAT THIS PART YOU WOULD NOT DO IN PRACTICE
Height<-fread("Yengo_Height_EAS.txt",data.table=FALSE)
Height<-subset(Height, Height$CHR == 1)
write.table(Height, file = "Yengo_Height_EAS_chr1.txt", col.names=TRUE,quote = FALSE, row.names = FALSE)

BMI<-fread("BMI_Japan.txt",data.table=FALSE)
head(BMI)
BMI<-subset(BMI, BMI$CHR == 1)
write.table(BMI, file = "BMI_EAS_chr1.txt", col.names=TRUE,quote = FALSE, row.names = FALSE)




#prep the psychiatric traits for practical

setwd("~/Desktop/Workshops/ISGWorkshop_2023/Practical/EUR/")
BIP<-fread("pgc-bip2021-all.vcf.tsv",data.table=FALSE)
head(BIP)
colnames(BIP)[1]<-"CHR"
BIP<-subset(BIP, BIP$CHR == 1)
write.table(BIP, file = "BIP_chr1.txt")

SCZ<-fread("PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv",data.table=FALSE)
SCZ$NEFF<-SCZ$NEFF*2
SCZ<-subset(SCZ, SCZ$CHROM == 1)
write.table(SCZ, file = "SCZ_chr1.txt")
