#### The Roles of GRIK2 and HARs in the Evolution and Development of Synapses
#### Shadi Shahatit, Master's Thesis, IBENS, Paris 2024
# libraries ---------------------------------------------------------------

library("BSgenome.Hsapiens.UCSC.hg38")
library("topGO")
library("DESeq2")
library("motifbreakR")
library("atSNP")
library("MotifDb")
library("enrichR")
library("memes")
library("universalmotif")
library("JASPAR2024")
library("JASPAR2022")
library("Biostrings")
library("TFBSTools")
library("ggplot2")
library("forcats")
library("wesanderson")
library("stringr")
library("tidyr")
library("dplyr")
library("readr")
library("splitstackshape")
library("ggridges")
library("ggthemes")
library("karyoploteR")
library("BiocManager")
library("GenomicRanges")
library("cowplot")
library("forcats")
library("scales")
library("palmerpenguins")
library("ggbeeswarm")
library("systemfonts")
library("ggforce")
library("ggpubr")
library("rtracklayer")
library("ggcorrplot")
library("scMethrix")
library("BRGenomics")
library("minpack.lm")
library("sjPlot")
library("DHARMa")
library("RRPP")
library("devtools")
library("flexplot")
library("lme4")
library("betareg")
library("nnet")
library("MASS")
library("ordinal")
library("brant")
library("ggeffects")
library("lvplot")
library("caret")
library("patchwork")
library("ggrepel")
library("MuMIn")
library("lmtest")
library("fitdistrplus")
library("VennDiagram")
library("readxl")
library("maftools")
library("ggplot2")
library("forcats")
library("wesanderson")
library("stringr")
library("tidyr")
library("dplyr")
library("readr")
library("splitstackshape")
library("ggridges")
library("ggthemes")
library("karyoploteR")
library("BiocManager")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("GenometriCorr")
library("rGREAT")
library("regioneR")
library("plyranges")
library("ggmotif")

## note: replace your system's directory in sys_dir
sys_dir <- "/home/shadi/Desktop/S4_Project"

# how many GWAS snps are in HARs? -----------------------------------------

## get your HARs

## HAR_3100_grik2

HAR_3100 <- read.table(file = file.path(sys_dir,"HAR/GSE180714_HARs.bed"), header = T)
search_dataframe <- function(df, pattern) {
  matches <- apply(df, 1, function(row) any(grepl(pattern, row, ignore.case = TRUE)))
  return(df[matches, ])
}
HAR_3100_grik2 <- search_dataframe(HAR_3100,"GRIK2")

## hg19_HAR

hg19_HAR <- read.table(file.path(sys_dir,"HAR/hg19liftover_Girskis_3100_grik2.bed"))
colnames(hg19_HAR) <- c("chr","start","end","id")
hg19_HAR$size <- hg19_HAR$end-hg19_HAR$start+1

## HAR_coordinates

HAR_coordinates <- data.frame(chr = "chr6",
                              start = 101296282,
                              end = 103544094,
                              id = "All_HAR_coordinates")

## AD

AD_Jansen2019_sheet8 <- read_excel(file.path(sys_dir,"GWAS_snps/AD_Jansen2019_sheet8.xlsx"), sheet=8, skip=2) %>% 
  filter(Locus == 6) %>%
  filter(bp > HAR_coordinates$start, bp < HAR_coordinates$end) %>%
  filter(bp > grik2$start-1000, bp < grik2$end)

## output: 0

## ADHD

ADHD_Demontis2023 <- read.table(file=file.path(sys_dir,"GWAS_snps/ADHD_Demontis2023.txt"),header=T) %>%
  filter(CHR == 6) %>%
  filter(P < 0.05)
ADHD_Demontis2023_filtered <- ADHD_Demontis2023 %>%
  filter(BP > hg19_HAR$start[1]+10, BP < hg19_HAR$end[7]+10)
  # filter(BP > grik2$start-1000, BP < grik2$end)
ADHD_Demontis2023_filtered_2 <- ADHD_Demontis2023_filtered %>%
  filter(
    !((BP > hg19_HAR$start[1] & BP < hg19_HAR$end[1]) |
      (BP > hg19_HAR$start[2] & BP < hg19_HAR$end[2]) |
      (BP > hg19_HAR$start[3] & BP < hg19_HAR$end[3]) |
      (BP > hg19_HAR$start[4] & BP < hg19_HAR$end[4]) |
      (BP > hg19_HAR$start[5] & BP < hg19_HAR$end[5]) |
      (BP > hg19_HAR$start[6] & BP < hg19_HAR$end[6]) |
      (BP > hg19_HAR$start[7] & BP < hg19_HAR$end[7])))
nrow(ADHD_Demontis2023_filtered)
nrow(ADHD_Demontis2023_filtered_2)

## output: 0

## SCZ

SCZ_trubetskoy2022 <- read.table(file=file.path(sys_dir,"GWAS_snps/SCZ_trubetskoy2022_chr6_R.tsv"), header=F)
SCZ_trubetskoy2022_filtered <- SCZ_trubetskoy2022 %>% filter(V11 < 0.05)

SCZ_trubetskoy2022_filtered_2 <- SCZ_trubetskoy2022_filtered %>%
  filter(
    ((V3 > HAR_3100_grik2$start[1] & V3 < HAR_3100_grik2$end[1]) |
        (V3 > HAR_3100_grik2$start[2] & V3 < HAR_3100_grik2$end[2]) |
        (V3 > HAR_3100_grik2$start[3] & V3 < HAR_3100_grik2$end[3]) |
        (V3 > HAR_3100_grik2$start[4] & V3 < HAR_3100_grik2$end[4]) |
        (V3 > HAR_3100_grik2$start[5] & V3 < HAR_3100_grik2$end[5]) |
        (V3 > HAR_3100_grik2$start[6] & V3 < HAR_3100_grik2$end[6]) |
        (V3 > HAR_3100_grik2$start[7] & V3 < HAR_3100_grik2$end[7])))

nrow(SCZ_trubetskoy2022_filtered_2)

## output: 9

## ASD

ASD_Grove2019_chr6 <- read.table(file=file.path(sys_dir,"GWAS_snps/ASD_Grove2019_chr6.txt"), header=T) %>%
  filter(P < 0.05)

ASD_Grove2019_chr6_filtered <- ASD_Grove2019_chr6 %>%
  filter(
    ((BP > HAR_3100_grik2$start[1] & BP < HAR_3100_grik2$end[1]) |
       (BP > HAR_3100_grik2$start[2] & BP < HAR_3100_grik2$end[2]) |
       (BP > HAR_3100_grik2$start[3] & BP < HAR_3100_grik2$end[3]) |
       (BP > HAR_3100_grik2$start[4] & BP < HAR_3100_grik2$end[4]) |
       (BP > HAR_3100_grik2$start[5] & BP < HAR_3100_grik2$end[5]) |
       (BP > HAR_3100_grik2$start[6] & BP < HAR_3100_grik2$end[6]) |
       (BP > HAR_3100_grik2$start[7] & BP < HAR_3100_grik2$end[7])))

nrow(ASD_Grove2019_chr6_filtered)

## output: 4

## BD

BD_Mullins2021_chr6 <- read.table(file=file.path(sys_dir,"GWAS_snps/BD_Mullins2021_chr6.tsv"), header=F) %>%
  filter(V8 < 0.05)

BD_Mullins2021_chr6_filtered <- BD_Mullins2021_chr6 %>%
  filter(
    ((V2 > HAR_3100_grik2$start[1] & V2 < HAR_3100_grik2$end[1]) |
       (V2 > HAR_3100_grik2$start[2] & V2 < HAR_3100_grik2$end[2]) |
       (V2 > HAR_3100_grik2$start[3] & V2 < HAR_3100_grik2$end[3]) |
       (V2 > HAR_3100_grik2$start[4] & V2 < HAR_3100_grik2$end[4]) |
       (V2 > HAR_3100_grik2$start[5] & V2 < HAR_3100_grik2$end[5]) |
       (V2 > HAR_3100_grik2$start[6] & V2 < HAR_3100_grik2$end[6]) |
       (V2 > HAR_3100_grik2$start[7] & V2 < HAR_3100_grik2$end[7])))

nrow(BD_Mullins2021_chr6_filtered)

## output: 9

write.xlsx2(BD_Mullins2021_chr6_filtered, 'BD_Mullins2021_chr6_filtered.xlsx',sheetName = "Sheet1")
write.xlsx2(SCZ_trubetskoy2022_filtered_2, 'SCZ_trubetskoy2022_filtered_2.xlsx',sheetName = "Sheet1")
write.xlsx2(ASD_Grove2019_chr6_filtered, 'ASD_Grove2019_chr6_filtered.xlsx',sheetName = "Sheet1")

ASD_Grove2019_chr6_filtered$BP <- as.character(ASD_Grove2019_chr6_filtered$BP)

ASD_Grove2019_chr6_filtered_2 <- ASD_Grove2019_chr6_filtered %>% 
  select(c(X6,SNP,BP,OR,SE,P))
colnames(ASD_Grove2019_chr6_filtered_2) <- c("CHR","SNP","POS","OR","SE","P-value")
ASD_Grove2019_chr6_filtered_2$condition <- "ASD"

SCZ_trubetskoy2022_filtered_3 <- SCZ_trubetskoy2022_filtered_2 %>% 
  select(c(V1,V2,V3,V9,V10,V11))
colnames(SCZ_trubetskoy2022_filtered_3) <- c("CHR","SNP","POS","OR","SE","P-value")
SCZ_trubetskoy2022_filtered_3$condition <- "SCZ"

BD_Mullins2021_chr6_filtered_2 <- BD_Mullins2021_chr6_filtered %>% 
  select(c(V1,V3,V2,V6,V7,V8))
colnames(BD_Mullins2021_chr6_filtered_2) <- c("CHR","SNP","POS","OR","SE","P-value")
BD_Mullins2021_chr6_filtered_2$condition <- "BD"

GWAS_HARs <- rbind(ASD_Grove2019_chr6_filtered_2,SCZ_trubetskoy2022_filtered_3,BD_Mullins2021_chr6_filtered_2)

ggplot()+
  geom_point(GWAS_HARs, mapping=aes(x=GWAS_HARs$SNP,y=GWAS_HARs$OR,color=GWAS_HARs$`P-value`))+
  geom_errorbar(GWAS_HARs, mapping=aes(x=SNP,ymin=OR-SE, ymax=OR+SE,color=GWAS_HARs$`P-value`), width=0.3)+
  labs(y="SCZ enrichment", x="Functional categories")+ # Enrichment = SNPs h2 proportion / SNPs proportion
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## MDD

memory.limit()

MDD_meng2024_chr6 <- read.table(file=file.path(sys_dir,"GWAS_snps/MDD_meng2024_chr6_p0.05.csv"))

MDD_meng2024_chr6_filtered <- MDD_meng2024_chr6 %>%
  filter(
    ((V9 > HAR_3100_grik2$start[1] & V9 < HAR_3100_grik2$end[1]) |
       (V9 > HAR_3100_grik2$start[2] & V9 < HAR_3100_grik2$end[2]) |
       (V9 > HAR_3100_grik2$start[3] & V9 < HAR_3100_grik2$end[3]) |
       (V9 > HAR_3100_grik2$start[4] & V9 < HAR_3100_grik2$end[4]) |
       (V9 > HAR_3100_grik2$start[5] & V9 < HAR_3100_grik2$end[5]) |
       (V9 > HAR_3100_grik2$start[6] & V9 < HAR_3100_grik2$end[6]) |
       (V9 > HAR_3100_grik2$start[7] & V9 < HAR_3100_grik2$end[7])))

nrow(MDD_meng2024_chr6_filtered)

## output: 0

## MDD

intelligence_Savage2018_chr6 <- read_excel(file.path(sys_dir,"GWAS_snps/intelligence_Savage2018_sheet6.xlsx"), sheet=6, skip=4) %>% 
  filter(CHR == 6)

intelligence_Savage2018_chr6_filtered <- intelligence_Savage2018_chr6 %>%
  filter(
    ((BP > HAR_3100_grik2$start[1] & BP < HAR_3100_grik2$end[1]) |
       (BP > HAR_3100_grik2$start[2] & BP < HAR_3100_grik2$end[2]) |
       (BP > HAR_3100_grik2$start[3] & BP < HAR_3100_grik2$end[3]) |
       (BP > HAR_3100_grik2$start[4] & BP < HAR_3100_grik2$end[4]) |
       (BP > HAR_3100_grik2$start[5] & BP < HAR_3100_grik2$end[5]) |
       (BP > HAR_3100_grik2$start[6] & BP < HAR_3100_grik2$end[6]) |
       (BP > HAR_3100_grik2$start[7] & BP < HAR_3100_grik2$end[7])))

intelligence_Savage2018_chr6_filtered <- intelligence_Savage2018_chr6 %>%
  filter(BP >= HAR_coordinates$start, BP <= HAR_coordinates$end)

nrow(intelligence_Savage2018_chr6_filtered)

## output: 0

# ldsc - partitioned h2 ---------------------------------------------------

## ASD

ASD_h2_UCSC <- read.table(file = file.path(sys_dir,"GWAS_snps/ASD_h2perUCSC.txt"), header = T)
ASD_h2_HARall <- read.table(file = file.path(sys_dir,"GWAS_snps/ASD_HAR_baselineLD.results"), header = T)
ASD_h2_HARall[1, 1] <- "HAR_all"
ASD_h2_HARchr6 <- read.table(file = file.path(sys_dir,"GWAS_snps/ASD_HAR_baselineLD_chr6_changesscipt.results"), header = T)
ASD_h2_HARchr6[1, 1] <- "HAR_chr6"
ASD_h2_HARgrik2 <- read.table(file = file.path(sys_dir,"GWAS_snps/ASD_grik2_baselineLD_chr6_changesscipt.results"), header = T)
ASD_h2_HARgrik2[1, 1] <- "HAR_grik2"
ASD_h2_HARgrik2flanking <- read.table(file = file.path(sys_dir,"GWAS_snps/ASD_grik2_flanking_baselineLD_chr6_changesscipt.results"), header = T)
ASD_h2_HARgrik2flanking[1, 1] <- "HAR_grik2_flanking"
ASD_h2 <- rbind(ASD_h2_UCSC,ASD_h2_HARall,ASD_h2_HARchr6,ASD_h2_HARgrik2,ASD_h2_HARgrik2flanking)
ASD_h2$Phenotype <- "ASD"

ASD_h2$p_value <- ifelse(ASD_h2$Enrichment_p < 0.05, "p < 0.05","p > 0.05")
Category_order <- fct_relevel(ASD_h2$Category,"Coding_UCSCL2_0","Coding_UCSC.flanking.500L2_0",
                                              "Promoter_UCSCL2_0","Promoter_UCSC.flanking.500L2_0",
                                              "Intron_UCSCL2_0","Intron_UCSC.flanking.500L2_0",
                                              "UTR_3_UCSCL2_0","UTR_3_UCSC.flanking.500L2_0",
                                              "UTR_5_UCSCL2_0","UTR_5_UCSC.flanking.500L2_0",
                                              "HAR_all","HAR_chr6","HAR_grik2","HAR_grik2_flanking")

ggplot()+
  geom_point(ASD_h2, mapping=aes(x=Category_order,y=Enrichment,color=p_value))+
  geom_errorbar(ASD_h2, mapping=aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error,x=Category_order,color=p_value), width=0.3)+
  scale_color_manual(values = c("#E94B3CFF","black"))+
  labs(y="ASD enrichment", x="Functional categories")+ # Enrichment = SNPs h2 proportion / SNPs proportion
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## SCZ

SCZ_h2_UCSC <- read.table(file = file.path(sys_dir,"GWAS_snps/SCZ_h2perUCSC.txt"), header = T)
SCZ_h2_HARall <- read.table(file = file.path(sys_dir,"GWAS_snps/SCZ_HAR_baselineLD.results"), header = T)
SCZ_h2_HARall[1, 1] <- "HAR_all"
SCZ_h2_HARchr6 <- read.table(file = file.path(sys_dir,"GWAS_snps/SCZ_HAR_baselineLD_chr6_changesscipt.results"), header = T)
SCZ_h2_HARchr6[1, 1] <- "HAR_chr6"
SCZ_h2_HARgrik2 <- read.table(file = file.path(sys_dir,"GWAS_snps/SCZ_grik2_baselineLD_chr6_changesscipt.results"), header = T)
SCZ_h2_HARgrik2[1, 1] <- "HAR_grik2"
SCZ_h2_HARgrik2flanking <- read.table(file = file.path(sys_dir,"GWAS_snps/SCZ_grik2_flanking_baselineLD_chr6_changesscipt.results"), header = T)
SCZ_h2_HARgrik2flanking[1, 1] <- "HAR_grik2_flanking"
SCZ_h2 <- rbind(SCZ_h2_UCSC,SCZ_h2_HARall,SCZ_h2_HARchr6,SCZ_h2_HARgrik2,SCZ_h2_HARgrik2flanking)
SCZ_h2$Phenotype <- "SCZ"

SCZ_h2$p_value <- ifelse(SCZ_h2$Enrichment_p < 0.05, "p < 0.05","p > 0.05")
Category_order <- fct_relevel(SCZ_h2$Category,"Coding_UCSCL2_0","Coding_UCSC.flanking.500L2_0",
                              "Promoter_UCSCL2_0","Promoter_UCSC.flanking.500L2_0",
                              "Intron_UCSCL2_0","Intron_UCSC.flanking.500L2_0",
                              "UTR_3_UCSCL2_0","UTR_3_UCSC.flanking.500L2_0",
                              "UTR_5_UCSCL2_0","UTR_5_UCSC.flanking.500L2_0",
                              "HAR_all","HAR_chr6","HAR_grik2","HAR_grik2_flanking")

ggplot()+
  geom_point(SCZ_h2, mapping=aes(x=Category_order,y=Enrichment,color=p_value))+
  geom_errorbar(SCZ_h2, mapping=aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error,x=Category_order,color=p_value), width=0.3)+
  scale_color_manual(values = c("#E94B3CFF","black"))+
  labs(y="SCZ enrichment", x="Functional categories")+ # Enrichment = SNPs h2 proportion / SNPs proportion
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## BD

BD_h2_UCSC <- read.table(file = file.path(sys_dir,"GWAS_snps/BD_h2perUCSC.txt"), header = T)
BD_h2_HARall <- read.table(file = file.path(sys_dir,"GWAS_snps/BD_HAR_baselineLD.results"), header = T)
BD_h2_HARall[1, 1] <- "HAR_all"
BD_h2_HARchr6 <- read.table(file = file.path(sys_dir,"GWAS_snps/BD_HAR_baselineLD_chr6_changesscipt.results"), header = T)
BD_h2_HARchr6[1, 1] <- "HAR_chr6"
BD_h2_HARgrik2 <- read.table(file = file.path(sys_dir,"GWAS_snps/BD_grik2_baselineLD_chr6_changesscipt.results"), header = T)
BD_h2_HARgrik2[1, 1] <- "HAR_grik2"
BD_h2_HARgrik2flanking <- read.table(file = file.path(sys_dir,"GWAS_snps/BD_grik2_flanking_baselineLD_chr6_changesscipt.results"), header = T)
BD_h2_HARgrik2flanking[1, 1] <- "HAR_grik2_flanking"
BD_h2 <- rbind(BD_h2_UCSC,BD_h2_HARall,BD_h2_HARchr6,BD_h2_HARgrik2,BD_h2_HARgrik2flanking)
BD_h2$Phenotype <- "BD"

BD_h2$p_value <- ifelse(BD_h2$Enrichment_p < 0.05, "p < 0.05","p > 0.05")
Category_order <- fct_relevel(BD_h2$Category,"Coding_UCSCL2_0","Coding_UCSC.flanking.500L2_0",
                              "Promoter_UCSCL2_0","Promoter_UCSC.flanking.500L2_0",
                              "Intron_UCSCL2_0","Intron_UCSC.flanking.500L2_0",
                              "UTR_3_UCSCL2_0","UTR_3_UCSC.flanking.500L2_0",
                              "UTR_5_UCSCL2_0","UTR_5_UCSC.flanking.500L2_0",
                              "HAR_all","HAR_chr6","HAR_grik2","HAR_grik2_flanking")

ggplot()+
  geom_point(BD_h2, mapping=aes(x=Category_order,y=Enrichment,color=p_value))+
  geom_errorbar(BD_h2, mapping=aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error,x=Category_order,color=p_value), width=0.3)+
  scale_color_manual(values = c("#E94B3CFF","black"))+
  labs(y="BD enrichment", x="Functional categories")+ # Enrichment = SNPs h2 proportion / SNPs proportion
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## combine ASD_SCZ_BD

ASD_SCZ_BD_hs <- rbind(BD_h2,SCZ_h2,ASD_h2)
ASD_SCZ_BD_hs$p_value <- ifelse(ASD_SCZ_BD_hs$Enrichment_p < 0.05, "p < 0.05","p > 0.05")
Category_order <- fct_relevel(ASD_SCZ_BD_hs$Category,"Coding_UCSCL2_0","Coding_UCSC.flanking.500L2_0",
                              "Promoter_UCSCL2_0","Promoter_UCSC.flanking.500L2_0",
                              "Intron_UCSCL2_0","Intron_UCSC.flanking.500L2_0",
                              "UTR_3_UCSCL2_0","UTR_3_UCSC.flanking.500L2_0",
                              "UTR_5_UCSCL2_0","UTR_5_UCSC.flanking.500L2_0",
                              "HAR_all","HAR_chr6","HAR_grik2","HAR_grik2_flanking")

ASD_SCZ_BD_hs_subsetregion <- ASD_SCZ_BD_hs %>% filter(Category != "Coding_UCSC.flanking.500L2_0" &
                                                        Category != "Promoter_UCSC.flanking.500L2_0" &
                                                        Category != "Intron_UCSC.flanking.500L2_0" &
                                                        Category != "UTR_3_UCSC.flanking.500L2_0" &
                                                        Category != "UTR_5_UCSC.flanking.500L2_0" &
                                                        Category != "Intron_UCSCL2_0" &
                                                        Category != "UTR_3_UCSCL2_0" &
                                                        Category != "UTR_5_UCSCL2_0")
Category_subsetregion_order <- fct_relevel(ASD_SCZ_BD_hs_subsetregion$Category,"Coding_UCSCL2_0","Promoter_UCSCL2_0",
                              "HAR_all","HAR_chr6","HAR_grik2","HAR_grik2_flanking")

ggplot()+
  geom_point(ASD_SCZ_BD_hs_subsetregion, mapping=aes(x=Category_subsetregion_order,y=Enrichment,color=p_value))+
  geom_errorbar(ASD_SCZ_BD_hs_subsetregion, mapping=aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error,x=Category_subsetregion_order,color=p_value), width=0.3)+
  scale_color_manual(values = c("#E94B3CFF","black"))+
  labs(y="Enrichment", x="Functional categories")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~Phenotype, scales = "free")

# Supp_tabels -------------------------------------------------------------

Rdataframes_titles <- data.frame(
  
  Sheet = c("Sheet1", "Sheet2", "Sheet3",
            "Sheet4", "Sheet5", "Sheet6",
            "Sheet7", "Sheet8", "Sheet9",
            "Sheet10", "Sheet11",
            "Sheet12",
            
            "Sheet13", "Sheet14", "Sheet15", "Sheet16",
            
            "Sheet17", "Sheet18"
  ),
  
  DataFrameName = c("HAR_3100", "HAR_3100_counts", "HAR_3100_grik2",
                    "great_HAR_flank_table_MF_df", "great_HAR_flank_table_BP_df" , "great_HAR_flank_table_CC_df",
                    "TF_allHAR_brainRPKML_highexp_single_unique", "HAR_motifbreakrate_df_brain", "TF_motifbreakrate_overall_df",
                    "neanderthal_snps", "denisovan_snps",
                    "GWAS_HARs",
                    
                    "grik2_RPKML_age", "grik2_expcom", "test_results", "encode_RPKML_avgexp_human",
                    
                    "promoter_PhyloP", "promoter_PhyloP_motif_labeled"
  ),
  
  DataFrameNote = c("3100 HARs defined by Girskis et al., 2021",
                    "HARs counts per closest gene for Girskis 3100 HARs",
                    "grik2HARs", 
                    "rGREAT results - grik2HARs genomic regions enrichment table for GO Molecular Function",
                    "rGREAT results - grik2HARs genomic regions enrichment table for GO Biological Process",
                    "rGREAT results - grik2HARs genomic regions enrichment table for GO Cellular Component",
                    "motifbreakR results - TF expressed during synaptogenesis in cortical tissue with significantly (top 5%) altered binding motif between human and chimp alleles",
                    "motifbreakR results - add motif change rates and divergence per bp", 
                    "motifbreakR results - raw data that include TF NOT filtered for effect size nor brain expression",
                    "neanderthal variations info in grik2HARs",
                    "denisovan variations info in grik2HARs",
                    "GWAS snps that intersected with HARs, including snps for Autism Spectrum Disorder, Bipolar Disorder, and Schizophrenia",
                    
                    "PsychENCODE raw data for grik2 expression in human and macaque with sample metadata",
                    "PsychENCODE grik2 expression in cortical tissues with developmental timepoints labels",
                    "stats test results per tissue (t or wilcox test) for grik2 expression diffrerence between human and macaque at each developmental timepoint",
                    "PsychENCODE grik2 expression in cortical tissues in human with genes higher than the average brain expression levels (used in motifbreakR analysis)",
                    
                    "human grik2 promoter with conservation scores (PhyloP) from USCS",
                    "human grik2 promoter with conservation scores (PhyloP) binned across 10 bp and labeled with the highest scoring cis-regulatory module identified by cbust (inRange)"
  ))

write_xlsx(list(
  ContentTable = Rdataframes_titles,
  
  Sheet1 = HAR_3100,
  Sheet2 = HAR_3100_counts,
  Sheet3 = HAR_3100_grik2,
  Sheet4 = great_HAR_flank_table_MF_df,
  Sheet5 = great_HAR_flank_table_BP_df,
  Sheet6 = great_HAR_flank_table_CC_df,
  Sheet7 = TF_allHAR_brainRPKML_highexp_single_unique,
  Sheet8 = HAR_motifbreakrate_df_brain,
  Sheet9 = TF_motifbreakrate_overall_df,
  
  Sheet10 = neanderthal_snps,
  Sheet11 = denisovan_snps,
  Sheet12 = GWAS_HARs,
  
  Sheet13 = grik2_RPKML_age,
  Sheet14 = grik2_expcom,
  Sheet15 = test_results,
  Sheet16 = encode_RPKML_avgexp_human,
  
  Sheet17 = promoter_PhyloP,
  Sheet18 = promoter_PhyloP_motif_labeled
),
file.path(sys_dir,"supp_data/supp_data_Rdataframes.xlsx"))


