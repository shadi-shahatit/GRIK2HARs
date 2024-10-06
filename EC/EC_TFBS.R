#### Transcription factor binding site analysis of FRMPD2 and FRMPD2B promoters
#### September 2024, for Eirini

#### Content:
#### libraries
#### 1. FIMO
####    1.3 extract results
####        - top scoring (1% & 5%)
####        - uniqueness (complete & partial)
####    1.4 combine to make life easier
#### 2. ClusterBuster
####	  2.3 extract results
####        - top scoring cluster or cis regulatory module
####	  2.4 gene set enrichment
#### 3. karyoplot F2B and HAQERs

# libraries ---------------------------------------------------------------

library(writexl)
library(readxl)
library(dplyr)
library(tidyr)
library(GenomicDistributions)
library(GenomicDistributionsData)
library("org.Hs.eg.db")
library(topGO)
library("memes")
library("universalmotif")
library("JASPAR2024")
library("JASPAR2022")
library(Biostrings)
library(TFBSTools)
library("ggplot2")
library("forcats")
library(stringr)
library(tidyr)
library("dplyr")
library("readr")
library("splitstackshape")
library(ggridges)
library(ggthemes)
library(karyoploteR)
library(BiocManager)
library(GenomicRanges)
library(forcats)
library("scales")
library(palmerpenguins)
library(ggbeeswarm)
library(systemfonts)
library("ggforce")
library("ggpubr")
library(rtracklayer)
library("ggcorrplot")
library("scMethrix")
library("minpack.lm")
library("sjPlot")
library("RRPP")
library("devtools")
library("lme4")
library(betareg)
library("nnet")
library("MASS")
library("ordinal")
library(brant)
library(ggeffects)
library("lvplot")
library(caret)
library(patchwork)
library(ggrepel)
library(MuMIn)
library("lmtest")
library(readxl)

## note: replace your system's directory in sys_dir
sys_dir <- "/home/shadi/Desktop/S4_Project/EC"

# 1. FIMO --------------------------------------------------------------------

#### 1.3 extract results (top scoring and unique)

## filter based on top scoring (95 and 99 percentiles of q-value)

## For FRMPD2Bpromoter.txt

# define fimo output and remove duplicates TF motifs
Fimo_FRMPD2Bpromoter <- importFimo(file.path(sys_dir,"fimo_out_FRMPD2Bpromoter/fimo.tsv"))
Fimo_FRMPD2Bpromoter_df <- Fimo_FRMPD2Bpromoter %>% as.data.frame()
Fimo_FRMPD2Bpromoter_unique <- Fimo_FRMPD2Bpromoter_df %>%
  arrange(motif_alt_id, desc(score)) %>%
  distinct(motif_alt_id, .keep_all = TRUE)
# 95 and 99 percentile filtration
Fimo_FRMPD2Bpromoter_df_top5 <- Fimo_FRMPD2Bpromoter_unique %>% filter(qvalue <= quantile(qvalue, 0.05))
Fimo_FRMPD2Bpromoter_df_top1 <- Fimo_FRMPD2Bpromoter_unique %>% filter(qvalue <= quantile(qvalue, 0.01))
# count and collect
nrow(Fimo_FRMPD2Bpromoter_unique)
cat(Fimo_FRMPD2Bpromoter_unique$motif_alt_id, sep = "\n")
nrow(Fimo_FRMPD2Bpromoter_df_top5)
cat(Fimo_FRMPD2Bpromoter_df_top5$motif_alt_id, sep = "\n")
nrow(Fimo_FRMPD2Bpromoter_df_top1)
cat(Fimo_FRMPD2Bpromoter_df_top1$motif_alt_id, sep = "\n")

## For 6exFRMPD2humanpromoter.txt

# define fimo output and remove duplicates TF motifs
Fimo_6exFRMPD2humanpromoter <- importFimo(file.path(sys_dir,"fimo_out_6exFRMPD2humanpromoter/fimo.tsv"))
Fimo_6exFRMPD2humanpromoter_df <- Fimo_6exFRMPD2humanpromoter %>% as.data.frame()
Fimo_6exFRMPD2humanpromoter_unique <- Fimo_6exFRMPD2humanpromoter_df %>%
  arrange(motif_alt_id, desc(score)) %>%
  distinct(motif_alt_id, .keep_all = TRUE)
# 95 and 99 percentile filtration
Fimo_6exFRMPD2humanpromoter_df_top5 <- Fimo_6exFRMPD2humanpromoter_unique %>% filter(qvalue <= quantile(qvalue, 0.05))
Fimo_6exFRMPD2humanpromoter_df_top1 <- Fimo_6exFRMPD2humanpromoter_unique %>% filter(qvalue <= quantile(qvalue, 0.01))
# count and collect
nrow(Fimo_6exFRMPD2humanpromoter_unique)
cat(Fimo_6exFRMPD2humanpromoter_unique$motif_alt_id, sep = "\n")
nrow(Fimo_6exFRMPD2humanpromoter_df_top5)
cat(Fimo_6exFRMPD2humanpromoter_df_top5$motif_alt_id, sep = "\n")
nrow(Fimo_6exFRMPD2humanpromoter_df_top1)
cat(Fimo_6exFRMPD2humanpromoter_df_top1$motif_alt_id, sep = "\n")

## For 6exFRMPD2macaquepromoter.txt

# define fimo output and remove duplicates TF motifs
Fimo_6exFRMPD2macaquepromoter <- importFimo(file.path(sys_dir,"fimo_out_6exFRMPD2macaquepromoter/fimo.tsv"))
Fimo_6exFRMPD2macaquepromoter_df <- Fimo_6exFRMPD2macaquepromoter %>% as.data.frame()
Fimo_6exFRMPD2macaquepromoter_unique <- Fimo_6exFRMPD2macaquepromoter_df %>%
  arrange(motif_alt_id, desc(score)) %>%
  distinct(motif_alt_id, .keep_all = TRUE)
# 95 and 99 percentile filtration
Fimo_6exFRMPD2macaquepromoter_df_top5 <- Fimo_6exFRMPD2macaquepromoter_unique %>% filter(qvalue <= quantile(qvalue, 0.05))
Fimo_6exFRMPD2macaquepromoter_df_top1 <- Fimo_6exFRMPD2macaquepromoter_unique %>% filter(qvalue <= quantile(qvalue, 0.01))
# count and collect
nrow(Fimo_6exFRMPD2macaquepromoter_unique)
cat(Fimo_6exFRMPD2macaquepromoter_unique$motif_alt_id, sep = "\n")
nrow(Fimo_6exFRMPD2macaquepromoter_df_top5)
cat(Fimo_6exFRMPD2macaquepromoter_df_top5$motif_alt_id, sep = "\n")
nrow(Fimo_6exFRMPD2macaquepromoter_df_top1)
cat(Fimo_6exFRMPD2macaquepromoter_df_top1$motif_alt_id, sep = "\n")

## filter based on uniqueness (complete and partial)

# combine the results from the three promoters into one table with a sequence label (source)
Fimo_FRMPD2Bpromoter_unique <- Fimo_FRMPD2Bpromoter_unique %>% mutate(source = "FRMPD2Bpromoter")
Fimo_6exFRMPD2humanpromoter_unique <- Fimo_6exFRMPD2humanpromoter_unique %>% mutate(source = "FRMPD2humanpromoter")
Fimo_6exFRMPD2macaquepromoter_unique <- Fimo_6exFRMPD2macaquepromoter_unique %>% mutate(source = "FRMPD2macaquepromoter")
combined_fimo_out <- bind_rows(Fimo_FRMPD2Bpromoter_unique, Fimo_6exFRMPD2humanpromoter_unique, Fimo_6exFRMPD2macaquepromoter_unique)

# completely unique motifs (appear in only one promoter but not the others)

completely_unique <- combined_fimo_out %>%
  group_by(motif_alt_id) %>%
  filter(n_distinct(source) == 1) %>%
  ungroup()

unique_to_FRMPD2B <- completely_unique %>% filter(source == "FRMPD2Bpromoter")
unique_to_FRMPD2human <- completely_unique %>% filter(source == "FRMPD2humanpromoter")
unique_to_FRMPD2macaque <- completely_unique %>% filter(source == "FRMPD2macaquepromoter")

# count and collect
nrow(unique_to_FRMPD2B)
cat(unique_to_FRMPD2B$motif_alt_id, sep = "\n")
nrow(unique_to_FRMPD2human)
cat(unique_to_FRMPD2human$motif_alt_id, sep = "\n")
nrow(unique_to_FRMPD2macaque)
cat(unique_to_FRMPD2macaque$motif_alt_id, sep = "\n")

# partially unique motifs (shared between only two promoters, but absent in the third)

partially_unique <- combined_fimo_out %>%
  group_by(motif_alt_id) %>%
  filter(n_distinct(source) == 2) %>%
  ungroup()

shared_FRMPD2B_FRMPD2human <- partially_unique %>% filter(source %in% c("FRMPD2Bpromoter", "FRMPD2humanpromoter"))
shared_FRMPD2B_FRMPD2macaque <- partially_unique %>% filter(source %in% c("FRMPD2Bpromoter", "FRMPD2macaquepromoter"))
shared_FRMPD2human_FRMPD2macaque <- partially_unique %>% filter(source %in% c("FRMPD2humanpromoter", "FRMPD2macaquepromoter"))

# count and collect
nrow(shared_FRMPD2B_FRMPD2human)
cat(shared_FRMPD2B_FRMPD2human$motif_alt_id, sep = "\n")
nrow(shared_FRMPD2B_FRMPD2macaque)
cat(shared_FRMPD2B_FRMPD2macaque$motif_alt_id, sep = "\n")
nrow(shared_FRMPD2human_FRMPD2macaque)
cat(shared_FRMPD2human_FRMPD2macaque$motif_alt_id, sep = "\n")

#### 1.4 combine to make life easier

## top 5 and 1
Fimo_FRMPD2Bpromoter_df_top1 <- Fimo_FRMPD2Bpromoter_df_top1 %>% mutate(status = "Top1_FRMPD2Bpromoter")
Fimo_6exFRMPD2humanpromoter_df_top1 <- Fimo_6exFRMPD2humanpromoter_df_top1 %>% mutate(status = "Top1_FRMPD2humanpromoter")
Fimo_6exFRMPD2macaquepromoter_df_top1 <- Fimo_6exFRMPD2macaquepromoter_df_top1 %>% mutate(status = "Top1_FRMPD2macaquepromoter")
Fimo_FRMPD2Bpromoter_df_top5 <- Fimo_FRMPD2Bpromoter_df_top5 %>% mutate(status = "Top5_FRMPD2Bpromoter")
Fimo_6exFRMPD2humanpromoter_df_top5 <- Fimo_6exFRMPD2humanpromoter_df_top5 %>% mutate(status = "Top5_FRMPD2humanpromoter")
Fimo_6exFRMPD2macaquepromoter_df_top5 <- Fimo_6exFRMPD2macaquepromoter_df_top5 %>% mutate(status = "Top5_FRMPD2macaquepromoter")
## unique
unique_to_FRMPD2B <- unique_to_FRMPD2B %>% mutate(status = "unique_to_FRMPD2B")
unique_to_FRMPD2human <- unique_to_FRMPD2human %>% mutate(status = "unique_to_FRMPD2human")
unique_to_FRMPD2macaque <- unique_to_FRMPD2macaque %>% mutate(status = "unique_to_FRMPD2macaque")
shared_FRMPD2B_FRMPD2human <- shared_FRMPD2B_FRMPD2human %>% mutate(status = "PartiallyUnique_FRMPD2B_FRMPD2human")
shared_FRMPD2B_FRMPD2macaque <- shared_FRMPD2B_FRMPD2macaque %>% mutate(status = "PartiallyUnique_FRMPD2B_FRMPD2macaque")
shared_FRMPD2human_FRMPD2macaque <- shared_FRMPD2human_FRMPD2macaque %>% mutate(status = "PartiallyUnique_FRMPD2human_FRMPD2macaque")

fimo_filtered_full_info <- bind_rows(
  Fimo_FRMPD2Bpromoter_df_top1,
  Fimo_6exFRMPD2humanpromoter_df_top1,
  Fimo_6exFRMPD2macaquepromoter_df_top1,
  Fimo_FRMPD2Bpromoter_df_top5,
  Fimo_6exFRMPD2humanpromoter_df_top5,
  Fimo_6exFRMPD2macaquepromoter_df_top5,
  unique_to_FRMPD2B,
  unique_to_FRMPD2human,
  unique_to_FRMPD2macaque,
  shared_FRMPD2B_FRMPD2human,
  shared_FRMPD2B_FRMPD2macaque,
  shared_FRMPD2human_FRMPD2macaque)

write_xlsx(fimo_filtered_full_info, file.path(sys_dir,"fimo_filtered_full_info.xlsx"))

# 2. ClusterBuster -----------------------------------------------------------

## cbustout_FRMPD2Bpromoter_1_symbols.txt

# define cbust output and transpose
cbust_FRMPD2Bpro <- read.table(file.path(sys_dir,"cbustout_FRMPD2Bpromoter_1_symbols.txt"), header=T)
cbust_FRMPD2Bpro <- cbust_FRMPD2Bpro %>% dplyr::select(-c(Start,End))
cbust_FRMPD2Bpro_t <- as.data.frame(t(cbust_FRMPD2Bpro[-1]))
cbust_FRMPD2Bpro_t$gene_sym <- rownames(cbust_FRMPD2Bpro_t)
rownames(cbust_FRMPD2Bpro_t) <- seq(1:nrow(cbust_FRMPD2Bpro_t))
# 95 and 99 percentile filtration for the highest scoring cluster V1
cbust_FRMPD2Bpro_t_top5 <- cbust_FRMPD2Bpro_t %>% filter(V1 > quantile(V1,0.95))
cbust_FRMPD2Bpro_t_top1 <- cbust_FRMPD2Bpro_t %>% filter(V1 > quantile(V1,0.99))
# count and collect 
nrow(cbust_FRMPD2Bpro_t)
cat(cbust_FRMPD2Bpro_t$gene_sym, sep = "\n")
nrow(cbust_FRMPD2Bpro_t_top5)
cat(cbust_FRMPD2Bpro_t_top5$gene_sym, sep = "\n")
nrow(cbust_FRMPD2Bpro_t_top1)
cat(cbust_FRMPD2Bpro_t_top1$gene_sym, sep = "\n")

# top 1 in FRMPD2B promoter
# BCL11A
# EGR3
# EGR4
# NEUROG2.1
# ZNF140
# Sox5

## cbustout_6exFRMPD2humanpromoter_1_symbols.txt

# define cbust output and transpose
cbust_FRMPD2humpro <- read.table(file.path(sys_dir,"cbustout_6exFRMPD2humanpromoter_1_symbols.txt"), header=T)
cbust_FRMPD2humpro <- cbust_FRMPD2humpro %>% dplyr::select(-c(Start,End))
cbust_FRMPD2humpro_t <- as.data.frame(t(cbust_FRMPD2humpro[-1]))
cbust_FRMPD2humpro_t$gene_sym <- rownames(cbust_FRMPD2humpro_t)
rownames(cbust_FRMPD2humpro_t) <- seq(1:nrow(cbust_FRMPD2humpro_t))
# 95 and 99 percentile filtration for the highest scoring cluster V1
cbust_FRMPD2humpro_t_top5 <- cbust_FRMPD2humpro_t %>% filter(V1 > quantile(V1,0.95))
cbust_FRMPD2humpro_t_top1 <- cbust_FRMPD2humpro_t %>% filter(V1 > quantile(V1,0.99))
# count and collect 
nrow(cbust_FRMPD2humpro_t)
cat(cbust_FRMPD2humpro_t$gene_sym, sep = "\n")
nrow(cbust_FRMPD2humpro_t_top5)
cat(cbust_FRMPD2humpro_t_top5$gene_sym, sep = "\n")
nrow(cbust_FRMPD2humpro_t_top1)
cat(cbust_FRMPD2humpro_t_top1$gene_sym, sep = "\n")

# top 1 in FRMPD2hum promoter
# NEUROG2
# OLIG2
# OLIG1
# OLIG3
# BHLHA15
# TAL1..TCF3

## cbustout_6exFRMPD2macaquepromoter_1_symbols.txt

# define cbust output and transpose
cbust_FRMPD2macpro <- read.table(file.path(sys_dir,"cbustout_6exFRMPD2macaquepromoter_1_symbols.txt"), header=T)
cbust_FRMPD2macpro <- cbust_FRMPD2macpro %>% dplyr::select(-c(Start,End))
cbust_FRMPD2macpro_t <- as.data.frame(t(cbust_FRMPD2macpro[-1]))
cbust_FRMPD2macpro_t$gene_sym <- rownames(cbust_FRMPD2macpro_t)
rownames(cbust_FRMPD2macpro_t) <- seq(1:nrow(cbust_FRMPD2macpro_t))
# 95 and 99 percentile filtration for the highest scoring cluster V1
cbust_FRMPD2macpro_t_top5 <- cbust_FRMPD2macpro_t %>% filter(V1 > quantile(V1,0.95))
cbust_FRMPD2macpro_t_top1 <- cbust_FRMPD2macpro_t %>% filter(V1 > quantile(V1,0.99))
# count and collect 
nrow(cbust_FRMPD2macpro_t)
cat(cbust_FRMPD2macpro_t$gene_sym, sep = "\n")
nrow(cbust_FRMPD2macpro_t_top5)
cat(cbust_FRMPD2macpro_t_top5$gene_sym, sep = "\n")
nrow(cbust_FRMPD2macpro_t_top1)
cat(cbust_FRMPD2macpro_t_top1$gene_sym, sep = "\n")

# top 1 in FRMPD2mac promoter
# POU3F3
# ZNF766
# BHLHE40
# XBP1

####	  2.4 gene set enrichment with https://maayanlab.cloud/Enrichr/ aganist all other synaptogenesis TF or random background genes

# 3. karyoplot F2b and HAQERs -----------------------------

# define f2b coordinates in hg38
F2B <- data.frame(chr = "chr10",
                    start = 46870858,
                    end = 46894562,
                    hgnc_symbol = "FRMPD2B") 

# create a genomic range object
F2B_gr <- makeGRangesFromDataFrame(df=F2B,
                                     keep.extra.columns=T,
                                     ignore.strand=T,
                                     seqnames.field="chr",
                                     start.field="start",
                                     end.field="end",
                                     starts.in.df.are.0based=FALSE)

HAQER <- read_excel(file.path(sys_dir,"HAQERs_1-s2.0-S0092867422013587-mmc1.xlsx"), sheet = 1)
colnames(HAQER) <- c("NAME","Chr","Start","End",
                     "SCORE","GSM595920_DHS","GSM595922_DHS" ,                   
                     "GSM595926_DHS","HAR","DHS",
                     "KANTON_DA","Gained After Rhesus Split","NumDatasets",
                     "DSB Meiotic Recombination Hotspot" ,"StarrSeqConstructId","STARR-seqStatus",                  
                     "STARR-seq coord (hg238)","Nearby Gene","Notes")
HAQER <- HAQER %>% dplyr::select(c("Chr","Start","End")) %>% drop_na()
HAQER_gr <- makeGRangesFromDataFrame(df=HAQER,
                                     keep.extra.columns=T,
                                     ignore.strand=T,
                                     seqnames.field="Chr",
                                     start.field="Start",
                                     end.field="End",
                                     starts.in.df.are.0based=FALSE)

HAQER_kp <- plotKaryotype(genome="hg38", plot.type = 2, chromosomes = c("chr10"))
kpPlotMarkers(HAQER_kp, data=F2B_gr, labels=F2B_gr$hgnc_symbol, text.orientation = "horizontal",
              r1=0.3, cex=0.8, adjust.label.position = FALSE)
y1 = c(0)
kpPoints(HAQER_kp, data=HAQER_gr, y=y1, cex=2, pch=6 ,col="#EE353E")
kpAddMainTitle(HAQER_kp,"F2B + HAQER")

HAQER_chr10 <- HAQER %>% filter(Chr == "chr10")
write.table(HAQER_chr10, file = file.path("HAQER_chr10.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
HAQER_chr10$F2B_distance <- (((HAQER_chr10$Start+HAQER_chr10$End)/2)-(F2B$start))
sorted_HAQER_chr10 <- HAQER_chr10[order(abs(HAQER_chr10$F2B_distance)), ]


