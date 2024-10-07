library(GenomicSEM)












######################################

files <- c("/home/shadi/Desktop/S4_Project/GWAS_snps/ASD_Grove2019.txt",
           "/home/shadi/Desktop/S4_Project/GWAS_snps/ADHD_Demontis2023.txt",
           "/home/shadi/Desktop/S4_Project/GWAS_snps/SCZ_trubetskoy2022.tsv",
           "/home/shadi/Desktop/S4_Project/GWAS_snps/BD_Mullins2021.tsv",
           "/home/shadi/Desktop/S4_Project/GWAS_snps/MDD_meng2024.csv")
hm3 <- "/home/shadi/Desktop/S4_Project/GWAS_snps/w_hm3.snplist"
trait.names <- c("ASD","ADHD","SCZ","BD","MDD")

# list the sample sizes
# define the MAF filter
# define the imputation quality filter

N <- c(18381,38691,67390,41917,2074)
info.filter=0.9
maf.filter=0.01

## Run munge
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)

## Run LDSC

# traits = the name of the .sumstats.gz traits produced by munge
# ld = folder of LD scores 
# wld = folder of LD scores
# sample.prev = the proportion of cases to total sample size. For quantitative traits list NA
# population.prev = the population lifetime prevalence of the traits. For quantitative traits list NA
# trait.names = optional fifth argument to list trait names so they can be named in your model

traits <- c("ASD.sumstats.gz")
ld <- "/home/shadi/Desktop/S4_Project/GWAS_snps/1000G_Phase3_ldscores/LDscore"
wld <- "/home/shadi/Desktop/S4_Project/GWAS_snps/1000G_Phase3_ldscores/LDscore"

ld <- "/home/shadi/Desktop/S4_Project/GWAS_snps/eur_w_ld_chr"
wld <- "/home/shadi/Desktop/S4_Project/GWAS_snps/eur_w_ld_chr"


sample.prev <- c(NA,NA)
population.prev <- c(NA,NA)
trait.names <- c("ASD")

LDSC_1000P3 <- ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev, ld=ld, wld=wld,trait.names=trait.names)

LDSC_1000P3$m

save(LDSC_1000P3, file = "LDSC_EUR.RData")

# genetic covariance matrix
LDSC_1000P3$S

# genetic correlation matrix
cov2cor(LDSC_1000P3$S)

# lets grab teh standard errors from the V matrix
k <- nrow(LDSC_1000P3$S)
SE <- matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <- sqrt(diag(LDSC_1000P3$V))

# matrix of Z-stats
Z <- LDSC_1000P3$S/SE

# matrix of p-values
P <- 2*pnorm(Z,lower.tail=FALSE)

# create a genetic heatmap 
library(corrplot)
rownames(LDSC_1000P3$S) <- colnames(LDSC_1000P3$S)

corrplot(corr = cov2cor(LDSC_1000P3$S),
         method = "color",
         addCoef.col = "dark grey",
         add = F,
         bg = "white",
         diag = T,
         outline = T,
         mar = c(0,0,2,0),
         number.cex=2,
         cl.pos = "b",
         cl.ratio = 0.125,
         cl.align.text = "l",
         cl.offset = 0.2,
         tl.srt=45,
         tl.pos = "lt",
         tl.offset=0.2,
         tl.col = "black",
         pch.col = "white",
         addgrid.col = "black",
         xpd = T,
         tl.cex=1.3,
         is.corr=TRUE,
         title = "European")

rm(list = ls())

# East Asian Example########################
setwd("../EAS/")

files<-c("BMI_EAS_chr1.txt","Yengo_Height_EAS_chr1.txt")

#define the reference file being used to allign alleles across summary stats
#here we are using hapmap3
hm3<-"eas_ldscores/w_hm3.snplist"

#name the traits 
trait.names<-c("BMI","Height")

#list the sample sizes. 
#Yengo has sample size in the file 
#but pan UKB does not so we enter it manually
N<-c(163835, NA)

#definte the imputation quality filter
info.filter=0.9

#define the MAF filter
maf.filter=0.01

#run munge
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)

###Step 2: run LDSC

#traits = the name of the .sumstats.gz traits produced by munge
traits<-c("BMI.sumstats.gz","Height.sumstats.gz")

##ld = folder of LD scores 
ld <- "eas_ldscores/"

#wld = folder of LD scores
wld <- "eas_ldscores/"

#sample.prev = the proportion of cases to total sample size. For quantitative traits list NA
sample.prev<-c(NA,NA)

#population.prev = the population lifetime prevalence of the traits. For quantitative traits list NA
population.prev<-c(NA,NA)

#trait.names = optional fifth argument to list trait names so they can be named in your model
trait.names<-c("BMI", "Height")

LDSC_EAS <- ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev, ld=ld, wld=wld,trait.names=trait.names)
save(LDSC_EAS, file = "LDSC_EAS.RData")

#genetic covariance matrix
LDSC_EAS$S

#genetic correlation matrix
cov2cor(LDSC_EAS$S)

#standard errors
k<-nrow(LDSC_EAS$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSC_EAS$V))

#matrix of Z-stats
Z<-LDSC_EAS$S/SE

#matrix of p-values
P<-2*pnorm(Z,lower.tail=FALSE)

#create corrplot again
rownames(LDSC_EAS$S)<-colnames(LDSC_EAS$S)

corrplot(corr = cov2cor(LDSC_EAS$S),
         method = "color",
         addCoef.col = "dark grey",
         add = F,
         bg = "white",
         diag = T,
         outline = T,
         mar = c(0,0,2,0),
         number.cex=2,
         cl.pos = "b",
         cl.ratio = 0.125,
         cl.align.text = "l",
         cl.offset = 0.2,
         tl.srt=45,
         tl.pos = "lt",
         tl.offset=0.2,
         tl.col = "black",
         pch.col = "white",
         addgrid.col = "black",
         xpd = T,
         tl.cex=1.3,
         is.corr=TRUE,
         title = "East Asian")

rm(list = ls())



##European Binary (case/control) Traits Example###
#####################################


files <- c("/home/shadi/Desktop/S4_Project/GWAS_snps/SCZ_trubetskoy2022.tsv",
           "/home/shadi/Desktop/S4_Project/GWAS_snps/BD_Mullins2021.tsv")
hm3 <- "/home/shadi/Desktop/S4_Project/GWAS_snps/w_hm3.snplist"
trait.names <- c("SCZ","BD")

N <- c(67390,41917)
info.filter=0.9
maf.filter=0.01

#list the sample sizes. 
#we write NA here as these files already have Neff as columns
N<-c(NA, NA)


#run munge
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)

###Step 2: run LDSC
#traits = the name of the .sumstats.gz traits produced by munge
traits<-c("SCZ.sumstats.gz","BP.sumstats.gz")

##ld = folder of LD scores 
ld <- "eur_w_ld_chr/"

#wld = folder of LD scores
wld <- "eur_w_ld_chr/"

#sample.prev: enter as 0.5 because used sum of effecitve N, 
#that accounts for sample ascertainment
sample.prev<-c(.5,.5)

#population.prev = the population lifetime prevalence of the traits. 
population.prev<-c(.01,.02)

#trait.names = optional fifth argument to list trait names so they can be named in your model
trait.names<-c("SCZ", "BIP")

LDSC_SCZ_BIP <- ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev, ld=ld, wld=wld,trait.names=trait.names)
save(LDSC_SCZ_BIP , file = "LDSC_SCZ_BIP .RData")

#genetic covariance matrix
LDSC_SCZ_BIP$S

#genetic correlation matrix
cov2cor(LDSC_SCZ_BIP$S)

#lets grab teh standard errors from the V matrix
k<-nrow(LDSC_SCZ_BIP$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSC_SCZ_BIP$V))

#matrix of Z-stats
Z<-LDSC_SCZ_BIP$S/SE

#matrix of p-values
P<-2*pnorm(Z,lower.tail=FALSE)

#create a genetic heatmap 
rownames(LDSC_SCZ_BIP$S)<-colnames(LDSC_SCZ_BIP$S)

corrplot(corr = cov2cor(LDSC_SCZ_BIP$S),
         method = "color",
         addCoef.col = "dark grey",
         add = F,
         bg = "white",
         diag = T,
         outline = T,
         mar = c(0,0,2,0),
         number.cex=2,
         cl.pos = "b",
         cl.ratio = 0.125,
         cl.align.text = "l",
         cl.offset = 0.2,
         tl.srt=45,
         tl.pos = "lt",
         tl.offset=0.2,
         tl.col = "black",
         pch.col = "white",
         addgrid.col = "black",
         xpd = T,
         tl.cex=1.3,
         is.corr=TRUE,
        title = "SCZ / BIP")

rm(list = ls())

