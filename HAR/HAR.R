#### The Roles of GRIK2 and HARs in the Evolution and Development of Synapses
#### Shadi Shahatit, Master's Thesis, IBENS, Paris 2024
# libraries ---------------------------------------------------------------

library("ggsci")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("writexl")
library("org.Hs.eg.db")  
library("AnnotationDbi")
library("biomaRt")
library("BSgenome.Hsapiens.UCSC.hg38")
library("seqinr")
library("Biostrings")
library("viridis")
library("topGO")
library("motifbreakR")
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
library("minpack.lm")
library("sjPlot")
library("DHARMa")
library("RRPP")
library("devtools")
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
library("regioneR")
library("plyranges")
library("atSNP")
library("scMethrix")
library("DESeq2")
library("GenometriCorr")
library("rGREAT")
library("maftools")
library("flexplot")
library("BRGenomics")
library("ggmotif")

## note: replace your system's directory in sys_dir
sys_dir <- "/home/shadi/Desktop/S4_Project"

# karyoplot grik2 and zooHAR (Keough et al., 2023) -----------------------------

grik2 <- data.frame(chr = "chr6",
                    start = 100962701,
                    end = 102081622,
                    hgnc_symbol = "GRIK2") 
## change to ensemble_canonical with these coordinates chr6:101393708-102070083
grik2_gr <- makeGRangesFromDataFrame(df=grik2,
                                     keep.extra.columns=T,
                                     ignore.strand=T,
                                     seqnames.field="chr",
                                     start.field="start",
                                     end.field="end",
                                     starts.in.df.are.0based=FALSE)

zooHAR <- read_excel(file.path(sys_dir,"HAR/science.abm1696_tables_s1_to_s9/science.abm1696_table_s1.xlsx"),
                     sheet = 2)
zooHAR_chr6 <- zooHAR %>% filter(chrom == "chr6")
zooHAR_gr <- makeGRangesFromDataFrame(df=zooHAR,
                                     keep.extra.columns=T,
                                     ignore.strand=T,
                                     seqnames.field=c("chr","chrom"),
                                     start.field="start",
                                     end.field="end",
                                     starts.in.df.are.0based=FALSE)

grik2_zooHAR_kp <- plotKaryotype(genome="hg38", plot.type = 2, chromosomes = c("chr6"))
kpPlotMarkers(grik2_zooHAR_kp, data=grik2_gr, labels=grik2_gr$hgnc_symbol, text.orientation = "horizontal",
              r1=0.3, cex=0.8, adjust.label.position = FALSE)
y1 = c(0)
kpPoints(grik2_zooHAR_kp, data=zooHAR_gr, y=y1, cex=2, pch=6 ,col="#EE353E")
kpAddMainTitle(grik2_zooHAR_kp,"grik2 + zooHAR")
kpPlotDensity(grik2_zooHAR_kp, data=zooHAR_gr, col="#EE353E")
kpPlotCoverage(grik2_zooHAR_kp, data=zooHAR_gr, col="#EE353E")
# y0 <- rnorm(15, mean=0.3, sd=0)
# y1 <- c(0, 0.3)
# kpRect(grik2_zooHAR_kp, chr = c("chr6"), x0=zooHAR_chr6_gr@ranges@start+1000000,
#                                         x1=zooHAR_chr6_gr@ranges@start+zooHAR_chr6_gr@ranges@width,
#                                         # y0=zooHAR_chr6_gr$Genoa_score, y1=zooHAR_chr6_gr$Genoa_score+0.05,
#                                         y0=y0,y1=y1)
# kpPlotRegions(grik2_zooHAR_kp, data=zooHAR_gr, data.panel=2)

zooHAR_chr6$grik2_distance <- (((zooHAR_chr6$start+zooHAR_chr6$end)/2)-(grik2$start))
sorted_zooHAR_chr6 <- zooHAR_chr6[order(abs(zooHAR_chr6$grik2_distance)), ]

# karyoplot grik2 and HAR (Girskis et al., 2021) -------------------------------

HAR_3100 <- read.table(file = file.path(sys_dir,"HAR/GSE180714_HARs.bed"), header = T)

HAR_3100_counts <- HAR_3100 %>%
  mutate(nearest_gene = gsub("c\\(|\\)", "", nearest_gene)) %>%
  separate_rows(nearest_gene, sep = ";") %>%
  as.data.frame() %>%
  dplyr::select(c(chr,start,end,HAR_ID,nearest_gene)) %>% 
  group_by(nearest_gene) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  filter(!grepl("^LINC0", nearest_gene) &
          !nearest_gene %in% c("NONE", "42795", "42797", "42985"))
HAR_3100_counts <- HAR_3100_counts %>% filter(!grepl("^MIR", nearest_gene))
HAR_3100_counts <- HAR_3100_counts %>% filter(!grepl("^LOC", nearest_gene)) %>% as.data.frame()

table(HAR_3100_counts$count)

grik2_count <- 7
HARs_all <- length(HAR_3100_counts$count)
HARs_belowgrik2count <- length(HAR_3100_counts$count[HAR_3100_counts$count < grik2_count])
HARs_abovegrik2count <- length(HAR_3100_counts$count[HAR_3100_counts$count >= grik2_count])
# HARs_below7 <- sum(HAR_3100_counts$count < 7)
# HARs_above7 <- sum(HAR_3100_counts$count >= 7)
(HARs_belowgrik2count/HARs_all)*100
(HARs_abovegrik2count/HARs_all)*100

ggplot(HAR_3100_counts, aes(x=(reorder(nearest_gene,-count)), y=count))+
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=ifelse(nearest_gene == "GRIK2","GRIK2","")), color="#EE353E", vjust=-0.5)+
  labs(x="Genes closest to HARs", y="HAR count")+
  theme_classic()+
  theme(axis.text.x = element_blank())

ggplot(HAR_3100_counts, aes(x=count)) +
  geom_histogram(binwidth=1, color="black", alpha=0.7) +
  geom_vline(xintercept=grik2_count, color="#EE353E", linetype="dashed", size=0.5) +
  labs(x="HAR counts", y="Frequency") +
  theme_classic()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

ggsave(file.path(sys_dir,"thesis_fig/fig2_test.png"), plot=ggplot2::last_plot())

write_xlsx(HAR_3100_counts, file.path(sys_dir,"Figure_4A.xlsx"))
write_xlsx(HAR_3100_counts, file.path(sys_dir,"Figure_4A_left.xlsx"))
write_xlsx(HAR_flank_grik2, file.path(sys_dir,"Figure_4A_right.xlsx"))

## include gene size

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
tx_by_gene <- transcriptsBy(txdb, by="gene")
gene_lens <- max(width(tx_by_gene))
gene_lens <- width(tx_by_gene)
max_gene_lens <- sapply(gene_lens, max)
entrez_ids <- names(max_gene_lens)
hgnc_symbols <- mapIds(org.Hs.eg.db,
                       keys = entrez_ids,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")
gene_size_hgnc <- data.frame(HGNC_Symbol = hgnc_symbols, Max_Transcript_Length = max_gene_lens)
rownames(gene_size_hgnc) <- seq(1:nrow(gene_size_hgnc))

HAR_counts_genesize <- merge(HAR_3100_counts, gene_size_hgnc, by.x="nearest_gene", by.y="HGNC_Symbol") %>% 
  mutate(norm_count=-log10(count/Max_Transcript_Length))

colnames(gene_size_hgnc) <- c("HGNC_Symbol","Max_Transcript_Length")
colnames(HAR_3100_counts) <- c("HGNC_Symbol","count")
unmerged_HARs <- anti_join(HAR_3100_counts, gene_size_hgnc)
HAR_3100_counts_filtered <- anti_join(HAR_3100_counts, unmerged_HARs, by = "HGNC_Symbol")

table(HAR_3100_counts_filtered$count)

HARs_all <- length(HAR_3100_counts_filtered$count)
HARs_belowgrik2count <- length(HAR_3100_counts_filtered$count[HAR_3100_counts_filtered$count < 7])
HARs_abovegrik2count <- length(HAR_3100_counts_filtered$count[HAR_3100_counts_filtered$count >= 7])
(HARs_belowgrik2count/HARs_all)*100
(HARs_abovegrik2count/HARs_all)*100

table(HAR_counts_genesize$norm_count)
summary(HAR_counts_genesize$norm_count)

grik2_normcount <- 5.1382
HARs_all <- length(HAR_counts_genesize$norm_count)
HARs_abovegrik2_normcount <- length(HAR_counts_genesize$norm_count[HAR_counts_genesize$norm_count > grik2_normcount])
HARs_belowgrik2_normcount <- length(HAR_counts_genesize$norm_count[HAR_counts_genesize$norm_count <= grik2_normcount])
(HARs_belowgrik2_normcount/HARs_all)*100
(HARs_abovegrik2_normcount/HARs_all)*100

ggplot(HAR_counts_genesize, aes(x=norm_count)) +
  geom_histogram(binwidth=1, color="black", alpha=0.7) +
  geom_vline(xintercept=grik2_normcount, color="#EE353E", linetype="dashed", size=0.5) +
  labs(x="normalized HAR count", y="Frequency") +
  theme_classic()

## include windows

HAR_forwindows <- HAR_3100 %>%
  mutate(HAR_precence = 1) %>%
  # mutate(nearest_gene = gsub("c\\(|\\)", "", nearest_gene)) %>%
  # separate_rows(nearest_gene, sep = ";") %>%
  as.data.frame() %>%
  dplyr::select(c(chr,start,end,HAR_precence,nearest_gene))

bin          <- 5000000

df2          <- HAR_forwindows %>% 
  mutate(window = start %/% bin)
df2$window_CHR <- paste(df2$window, df2$chr, sep = "_") 
df2          <- df2 %>%
  group_by(window_CHR) %>%
  summarise(sum_nearest_gene = list(unlist(nearest_gene)),across(c(HAR_precence), sum))
df2[,4:5] <- stringr::str_split_fixed(df2$window_CHR, "_", 2)
# df2 <- df2[-c(1)] 
colnames(df2) <- c("ID","sum_nearest_gene","sum_HAR_precence","window","CHR")
df2$window <- as.numeric(df2$window)
df2          <- df2 %>% as.data.frame() %>%
  mutate(start = (window*bin)+1, stop=(window+1)*bin) %>% 
  relocate(c("window","CHR","start","stop","sum_nearest_gene","sum_HAR_precence","ID")) %>% as.data.frame()

df_out       <- c()
df2_tmp       <- c()
missing_win       <- c()
df_tmp       <- c()

for(i in unique(df2$CHR)) {
  df2_tmp      <- df2 %>% filter(CHR == i)
  missing_win  <- setdiff(0:max(df2_tmp$window),unique(df2_tmp$window))
  df_tmp       <- data.frame(window = missing_win) %>% 
    mutate(CHR    = i,
           sum_nearest_gene = NA,
           sum_HAR_precence = NA,
           start  = (window*bin)+1, 
           stop   = (window+1)*bin) %>% 
    bind_rows(df2_tmp) %>% 
    arrange(CHR,window)
  df_out       <- bind_rows(df_out,df_tmp)
}
HAR_forwindows_binned <- df_out %>% arrange(CHR, window) %>% drop_na()

df2 <- df2 %>% as.data.frame() %>% dplyr::select(c(sum_nearest_gene,sum_HAR_precence,ID))
grik2_window <- 10
HARs_all <- length(df2$sum_HAR_precence)
HARs_belowgrik2window <- length(df2$sum_HAR_precence[df2$sum_HAR_precence < grik2_window])
HARs_abovegrik2window <- length(df2$sum_HAR_precence[df2$sum_HAR_precence >= grik2_window])
(HARs_belowgrik2window/HARs_all)*100
(HARs_abovegrik2window/HARs_all)*100

ggplot(df2, aes(x=sum_HAR_precence)) +
  geom_histogram(binwidth=1, color="black", alpha=0.7) +
  geom_vline(xintercept=grik2_window, color="#EE353E", linetype="dashed", size=0.5) +
  labs(x="mormalized HAR count", y="Frequency") +
  theme_classic()

zooHAR_filtered <- zooHAR %>% dplyr::select("chr","start","end")
colnames(zooHAR_filtered) <- c("chr","start","end")
zooHAR_filtered_HAR_3100_merged <- merge(HAR_3100,zooHARR)
nrow(zooHAR_filtered_HAR_3100_merged)

HAR_3100_gr <- makeGRangesFromDataFrame(df=HAR_3100,
                                      keep.extra.columns=T,
                                      ignore.strand=T,
                                      seqnames.field=c("chr","chrom"),
                                      start.field="start",
                                      end.field="end",
                                      starts.in.df.are.0based=FALSE)
HAR_3100_kp <- plotKaryotype(genome="hg38")
kpPlotDensity(HAR_3100_kp, data=HAR_3100_gr, col="steelblue")
kpAddMainTitle(HAR_3100_kp,"HAR (Girskis et al., 2021)")
HAR_3100_kp <- plotKaryotype(genome="hg38", plot.type = 2, chromosomes = c("chr6"))
kpPlotMarkers(HAR_3100_kp, data=grik2_gr, labels=grik2_gr$hgnc_symbol, text.orientation = "horizontal",
              r1=0.3, cex=0.8, adjust.label.position = FALSE)
y1 = c(0)
kpPoints(HAR_3100_kp, data=HAR_3100_gr, y=y1, cex=2, pch=6 ,col="#EE353E")
kpAddMainTitle(HAR_3100_kp,"grik2 + HAR")

search_dataframe <- function(df, pattern) {
  matches <- apply(df, 1, function(row) any(grepl(pattern, row, ignore.case = TRUE)))
  return(df[matches, ])
}
HAR_3100_grik2 <- search_dataframe(HAR_3100,"GRIK2")
HAR_3100_grik2
write_xlsx(HAR_3100_grik2, file.path(sys_dir,"HAR/HAR_3100_grik2.xlsx"))
HAR_3100_grik2_gr <- makeGRangesFromDataFrame(df=HAR_3100_grik2,
                                        keep.extra.columns=T,
                                        ignore.strand=T,
                                        seqnames.field=c("chr","chrom"),
                                        start.field="start",
                                        end.field="end",
                                        starts.in.df.are.0based=FALSE)
HAR_3100_grik2_kp <- plotKaryotype(genome="hg38", plot.type = 2, chromosomes = c("chr6"))
kpPlotMarkers(HAR_3100_grik2_kp, data=grik2_gr, labels=grik2_gr$hgnc_symbol, text.orientation = "horizontal",
              r1=0.3, cex=0.8, adjust.label.position = FALSE)
kpPoints(HAR_3100_grik2_kp, data=HAR_3100_gr, y=y1, cex=2, pch=6 ,col="#DBDBDB")
kpPoints(HAR_3100_grik2_kp, data=HAR_3100_grik2_gr, y=y1, cex=2, pch=6 ,col="#EE353E")

## density

HAR_3100_grik2_kp <- plotKaryotype(genome="hg38", plot.type = 2, chromosomes = c("chr6"))
kpPlotMarkers(HAR_3100_grik2_kp, data=grik2_gr, labels=grik2_gr$hgnc_symbol, text.orientation = "horizontal",
              r1=1.5, cex=0.8, adjust.label.position = FALSE)
kpPlotDensity(HAR_3100_grik2_kp, data=HAR_3100_gr, col="#DBDBDB")
kpPlotDensity(HAR_3100_grik2_kp, data=HAR_3100_grik2_gr, col="#EE353E")

## zoom-in

zoom.region_grik2 <- toGRanges(data.frame("chr6", grik2_start-5000000, grik2_end+5000000))
HAR_grik2_kp <- plotKaryotype(chromosomes="chr6", zoom=zoom.region_grik2)
kpDataBackground(HAR_grik2_kp)
kpAddBaseNumbers(HAR_grik2_kp)
kpAddCytobandLabels(HAR_grik2_kp)
kpPlotMarkers(HAR_grik2_kp, data=grik2_gr, labels=grik2_gr$hgnc_symbol, text.orientation = "horizontal",r1=1, cex=0.8, adjust.label.position = FALSE)
kpPoints(HAR_grik2_kp, data=HAR_flank_grik2_gr, y=y1,  r0=0.1, r1=1, cex=2, pch=6, col="#EE353E")

## include metadata from the paper (epigenetic markers and MPRA)

HAR_3100_grik2$size <- abs(HAR_3100_grik2$start-HAR_3100_grik2$end)+1

HAR_3100_grik2$grik2_distance <- (((HAR_3100_grik2$start+HAR_3100_grik2$end)/2)-(grik2$start))
summary(abs(HAR_3100_grik2$grik2_distance))

nrow(HAR_3100_grik2)
summary(HAR_3100_grik2$size)
sd(HAR_3100_grik2$size)
summary(HAR_3100_grik2$grik2_distance)
sd(HAR_3100_grik2$grik2_distance)

ggplot()+
  geom_col(HAR_3100_grik2,mapping=aes(x=HAR_ID,y=size),fill="steelblue")+
  geom_text(HAR_3100_grik2,mapping=aes(x=HAR_ID,y=size,label=size), vjust = -0.5, color = "black")+
  labs(title = "",x = "HAR ID",y = "HAR size",fill = "") +
  theme_classic()

HAR_3100_grik2 <- HAR_3100_grik2 %>% dplyr::select(c(1,2,3))
write.table(HAR_3100_grik2, file = '/home/shadi/Desktop/S4_Project/HAR/HAR_3100_grik2.bed', sep = '\t', row.names = F,quote=F)

HAR_3100_sig <- read_excel(file.path(sys_dir,"HAR/1-s2.0-S0896627321005808-mmc7.xlsx"))

HAR_3100_sig_gr <- makeGRangesFromDataFrame(df=HAR_3100_sig,
                                        keep.extra.columns=T,
                                        ignore.strand=T,
                                        seqnames.field=c("chr","chrom"),
                                        start.field="start",
                                        end.field="stop",
                                        starts.in.df.are.0based=FALSE)

HAR_3100_sig_grik2 <- search_dataframe(HAR_3100_sig,"GRIK2")

## include metadata from the paper (epigenetic markers and MPRA) - prioritized HARs

HAR_210_sig <- read_excel(file.path(sys_dir,"HAR/1-s2.0-S0896627321005808-mmc7.xlsx"),sheet=2)
HAR_210_sig_gr <- makeGRangesFromDataFrame(df=HAR_210_sig,
                                            keep.extra.columns=T,
                                            ignore.strand=T,
                                            seqnames.field=c("chr","chrom"),
                                            start.field="start",
                                            end.field="stop",
                                            starts.in.df.are.0based=FALSE)
HAR_210_kp <- plotKaryotype(genome="hg38", plot.type = 2, chromosomes = c("chr6"))
kpPlotMarkers(HAR_210_kp, data=grik2_gr, labels=grik2_gr$hgnc_symbol, text.orientation = "horizontal",
              r1=0.3, cex=0.8, adjust.label.position = FALSE)
kpPoints(HAR_210_kp, data=HAR_3100_gr, y=y1, cex=2, pch=6 ,col="#DBDBDB")
kpPoints(HAR_210_kp, data=HAR_3100_grik2_gr, y=y1, cex=2, pch=6 ,col="#EE353E")
kpPoints(HAR_210_kp, data=HAR_210_sig_gr, y=y1, cex=2, pch=6 ,col="#6667AB")

HAR_210_sig_chr6 <- HAR_210_sig %>% filter(chr == "chr6")
HAR_210_sig_chr6$grik2_distance <- (((HAR_210_sig_chr6$start+HAR_210_sig_chr6$stop)/2)-(grik2$start))
sorted_HAR_210_sig_chr6 <- HAR_210_sig_chr6[order(abs(HAR_210_sig_chr6$grik2_distance)), ]
summary(abs(sorted_HAR_210_sig_chr6$grik2_distance))

# karyoplot grik2 and HAR (Doan et al., 2016) ----------------------------------

HAR_2737 <- read_excel(file.path(sys_dir,"HAR/1-s2.0-S0092867416311692-mmc2.xlsx"))

HAR_2737_gr <- makeGRangesFromDataFrame(df=HAR_2737,
                                        keep.extra.columns=T,
                                        ignore.strand=T,
                                        seqnames.field=c("Chr","chrom"),
                                        start.field="Start",
                                        end.field="End",
                                        starts.in.df.are.0based=FALSE)

HAR_2737_kp <- plotKaryotype(genome="hg19", plot.type = 2)
kpPlotDensity(HAR_2737_kp, data=HAR_2737_gr, col="steelblue")
HAR_2737_kp <- plotKaryotype(genome="hg19", plot.type = 2, chromosomes = c("chr6"))
# kpPlotMarkers(HAR_3100_kp, data=grik2_gr, labels=grik2_gr$hgnc_symbol, text.orientation = "horizontal",
#               r1=0.3, cex=0.8, adjust.label.position = FALSE) ## diffrernt ref
y1 = c(0)
kpPoints(HAR_2737_kp, data=HAR_2737_gr, y=y1, cex=2, pch=6 ,col="#EE353E")
kpAddMainTitle(HAR_3100_kp,"grik2 + HAR")

# karyoplot grik2 and HAR (Mangan et al., 2022)  --------------------------------------------------

HAQER <- read_excel(file.path(sys_dir,"HAR/1-s2.0-S0092867422013587-mmc1.xlsx"), sheet = 1)
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

HAQER_kp <- plotKaryotype(genome="hg38", plot.type = 2)
kpPlotDensity(HAQER_kp, data=HAQER_gr, col="steelblue")
HAQER_kp <- plotKaryotype(genome="hg38", plot.type = 2, chromosomes = c("chr6"))
kpPlotMarkers(HAQER_kp, data=grik2_gr, labels=grik2_gr$hgnc_symbol, text.orientation = "horizontal",
               r1=0.3, cex=0.8, adjust.label.position = FALSE)
y1 = c(0)
kpPoints(HAQER_kp, data=HAQER_gr, y=y1, cex=2, pch=6 ,col="#EE353E")
kpAddMainTitle(HAQER_kp,"grik2 + HAR")

HAQER_chr6 <- HAQER %>% filter(Chr == "chr6")
HAQER_chr6$grik2_distance <- (((HAQER_chr6$Start+HAQER_chr6$End)/2)-(grik2$start))
sorted_HAQER_chr6 <- HAQER_chr6[order(abs(HAQER_chr6$grik2_distance)), ]

# karyoplot grik2 and HAR (Gittelman et al. 2015) ----------------------------------

HAR_DHS <- read_excel(file.path(sys_dir,"HAR/Supplemental_Table_2.xlsx"))

HAR_DHS_gr <- makeGRangesFromDataFrame(df=HAR_DHS,
                                        keep.extra.columns=T,
                                        ignore.strand=T,
                                        seqnames.field=c("chr","chrom"),
                                        start.field="start",
                                        end.field="stop",
                                        starts.in.df.are.0based=FALSE)

HAR_DHS_kp <- plotKaryotype(genome="hg19", plot.type = 2, chromosomes = c("chr6"))
kpPlotMarkers(HAR_DHS_kp, data=grik2_gr, labels=grik2_gr$hgnc_symbol, text.orientation = "horizontal",
              r1=0.3, cex=0.8, adjust.label.position = FALSE)
y1 = c(0)
kpPoints(HAR_DHS_kp, data=HAR_DHS_gr, y=y1, cex=2, pch=6 ,col="#EE353E")
kpAddMainTitle(HAR_DHS_kp,"grik2 + HAR")

HAR_DHS_chr6 <- HAR_DHS %>% filter(chr == "chr6")
HAR_DHS_chr6$grik2_distance <- (((HAR_DHS_chr6$start+HAR_DHS_chr6$stop)/2)-(grik2$start))
sorted_HAR_DHS_chr6 <- HAR_DHS_chr6[order(abs(HAR_DHS_chr6$grik2_distance)), ]

# all HAR karyoplot -----------------------------------------------------------

allHAR_grik2_kp <- plotKaryotype(genome="hg38", plot.type = 2, chromosomes = c("chr6"))
kpPlotMarkers(allHAR_grik2_kp, data=grik2_gr, labels=grik2_gr$hgnc_symbol, text.orientation = "horizontal",
              r1=0.3, cex=0.8, adjust.label.position = FALSE,data.panel=2)
y1 = c(0)
kpPoints(allHAR_grik2_kp, data=zooHAR_gr, y=y1, cex=2, pch=6 ,col="#481568FF",r0=0, r1=0.25,data.panel=2)
kpPoints(allHAR_grik2_kp, data=HAR_3100_gr, y=y1, cex=2, pch=6 ,col="#39558CFF",r0=0.25, r1=0.5,data.panel=2)
kpPoints(allHAR_grik2_kp, data=HAR_2737_gr, y=y1, cex=2, pch=6 ,col="#3CBC75FF", r0=0.5, r1=0.75,data.panel=2)
kpPoints(allHAR_grik2_kp, data=HAQER_gr, y=y1, cex=2, pch=6 ,col="#DCE318FF", r0=0.75, r1=1,data.panel=2)
kpPlotDensity(allHAR_grik2_kp, data=zooHAR_gr, col="#481568FF",r0=0, r1=0.25)
kpPlotDensity(allHAR_grik2_kp, data=HAR_3100_gr, col="#39558CFF",r0=0.25, r1=0.5)
kpPlotDensity(allHAR_grik2_kp, data=HAR_2737_gr, col="#3CBC75FF", r0=0.5, r1=0.75)
kpPlotDensity(allHAR_grik2_kp, data=HAQER_gr, col="#DCE318FF", r0=0.75, r1=1)
kpAxis(allHAR_grik2_kp, ymax=1, r0=0, r1=0.25, numticks = 2, col="#666666", cex=0.5,text.col="white")
kpAxis(allHAR_grik2_kp, ymax=1, r0=0.25, r1=0.5, numticks = 2, col="#666666", cex=0.5,text.col="white")
kpAxis(allHAR_grik2_kp, ymax=1, r0=0.5, r1=0.75, numticks = 2, col="#666666", cex=0.5,text.col="white")
kpAxis(allHAR_grik2_kp, ymax=1, r0=0.75, r1=1, numticks = 2, col="#666666", cex=0.5,text.col="white")
kpAddMainTitle(grik2_zooHAR_kp,"grik2 + Different HAR lists")

# karyoplot grik2 and hCONDEL (Xue et al., 2023) ----------------------------

hCONDEL <- read_excel(file.path(sys_dir,"HAR/science.abn2253_table_s1.xlsx"),
                     sheet = 2) %>% dplyr::select(c("hg38_del_chr","hg38_del_start_pos","hg38_del_end_pos"))

hCONDEL_gr <- makeGRangesFromDataFrame(df=hCONDEL,
                                      keep.extra.columns=T,
                                      ignore.strand=T,
                                      seqnames.field=c("chr","hg38_del_chr"),
                                      start.field="hg38_del_start_pos",
                                      end.field="hg38_del_end_pos",
                                      starts.in.df.are.0based=FALSE)

hCONDEL_kp <- plotKaryotype(genome="hg38", plot.type = 2, chromosomes = c("chr6"))
kpPlotMarkers(grik2_zooHAR_kp, data=grik2_gr, labels=grik2_gr$hgnc_symbol, text.orientation = "horizontal",
              r1=0.3, cex=0.8, adjust.label.position = FALSE)
y1 = c(0)
kpPoints(hCONDEL_kp, data=hCONDEL_gr, y=y1, cex=2, pch=6 ,col="black")
kpAddMainTitle(hCONDEL_kp,"grik2 + hCONDEL")

# grik2 CRE (Meuleman et al., 2020) --------------------------------------------

## use ref: Index and biological spectrum of human DNase I hypersensitive sites
## get DHS file from: https://www.encodeproject.org/annotations/ENCSR857UZV/
## awk 'BEGIN {OFS="\t"} {gsub(" ", "", $10); gsub(" ", "", $11); gsub(" ", "", $12); print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10$11$12}' ENCFF503GCK.tsv > DHS_ENCFF503GCK.tsv 

DHS_hg38 <- read.table(file = file.path(sys_dir,"HAR/DHS_ENCFF503GCK.tsv"), header = T)
DHS_hg38_chr6_Neural <- DHS_hg38 %>% filter(seqname == "chr6") %>% filter(component == "Neural")

DHS_hg38_chr6_Neural_gr <- makeGRangesFromDataFrame(df=DHS_hg38_chr6_Neural,
                                        keep.extra.columns=T,
                                        ignore.strand=T,
                                        seqnames.field=c("chr","seqname"),
                                        start.field="start",
                                        end.field="end",
                                        starts.in.df.are.0based=FALSE)

DHS_kp <- plotKaryotype(genome="hg38", plot.type = 2, chromosomes = c("chr6"))
kpPlotDensity(DHS_kp, data=DHS_hg38_chr6_Neural_gr, col="steelblue")
kpPlotMarkers(DHS_kp, data=grik2_gr, labels=grik2_gr$hgnc_symbol, text.orientation = "horizontal",
              r1=0.3, cex=0.8, adjust.label.position = FALSE)
kpAddMainTitle(DHS_kp,"grik2 + DHS")

## extract CREs for grik which lies in 25 kb upstream and 3 kb downstream

upstream_25kb <- grik2$start - 25000
downstream_3kb <- grik2$end + 3000

CRE_grik2 <- DHS_hg38_chr6_Neural[DHS_hg38_chr6_Neural$start >= upstream_25kb & DHS_hg38_chr6_Neural$start <= downstream_3kb, ]

CRE_grik2_gr <- makeGRangesFromDataFrame(df=CRE_grik2,
                                                    keep.extra.columns=T,
                                                    ignore.strand=T,
                                                    seqnames.field=c("chr","seqname"),
                                                    start.field="start",
                                                    end.field="end",
                                                    starts.in.df.are.0based=FALSE)

CRE_grik2_kp <- plotKaryotype(genome="hg38", plot.type = 2, chromosomes = c("chr6"))
kpPlotDensity(CRE_grik2_kp, data=DHS_hg38_chr6_Neural_gr, col="#DBDBDB")
kpPlotDensity(CRE_grik2_kp, data=CRE_grik2_gr, col="black")
# y1 = c(0)
# kpPoints(CRE_grik2_kp, data=CRE_grik2_gr, y=y1, cex=2, pch=4 ,col="steelblue")
kpPlotMarkers(CRE_grik2_kp, data=grik2_gr, labels=grik2_gr$hgnc_symbol, text.orientation = "horizontal",
              r1=1.5, cex=0.8, adjust.label.position = FALSE)
kpPoints(CRE_grik2_kp, data=HAR_3100_grik2_gr, y=y1, cex=2, pch=6 ,col="#EE353E")
kpPoints(CRE_grik2_kp, data=HAR_210_sig_gr, y=y1, cex=2, pch=6 ,col="#6667AB")
kpAddMainTitle(CRE_grik2_kp,"grik2 + HAR + DHS")

CRE_HAR_overlaps <- findOverlaps(CRE_grik2_gr,HAR_3100_grik2_gr)
CRE_grik2$CRE_distance <- ((CRE_grik2$start)-(HAR_3100_grik2$start))
sorted_CRE_grik2 <- CRE_grik2[order(abs(CRE_grik2$CRE_distance)), ]

# grik2 HAR/DHS/hCONDEL overlap -------------------------------------------------------------------------

## test the following data sets for overlaps

zooHAR_gr # NO
HAR_3100_gr # YES 3 
HAR_2737_kp # YES 2
HAQER_gr # YES 1
HAR_DHS_gr # NO 
hCONDEL_gr # YES 3

grik2_overlap <- findOverlaps(grik2_gr,hCONDEL_gr)
grik2_overlap

HAR_overlap <- findOverlaps(HAR_3100_grik2_gr,zooHAR_gr)
HAR_overlap

# Girskis HAR --------------------------------------------

hg19_HAR <- read.table(file = file.path(sys_dir,"HAR/hg19liftover_Girskis_3100_grik2.bed"))
colnames(hg19_HAR) <- c("chr","start","end","id")
hg19_HAR$size <- hg19_HAR$end-hg19_HAR$start+1
hg19_HAR_gr <- makeGRangesFromDataFrame(df=hg19_HAR,
                                      keep.extra.columns=T,
                                      ignore.strand=T,
                                      seqnames.field=c("chr","chrom"),
                                      start.field="start",
                                      end.field="end",
                                      starts.in.df.are.0based=FALSE)
hg19_HAR_kp <- plotKaryotype(genome="hg19", plot.type = 2, chromosomes = c("chr6"))
y1 = c(0)
kpPoints(hg19_HAR_kp, data=hg19_HAR_gr, y=y1, cex=2, pch=6 ,col="#EE353E")

## full info list for Girskis HAR_3100_grik2

HAR_3100_grik2

HAR_coordinates <- data.frame(chr = "chr6",
                    start = 101296282,
                    end = 103544094,
                    id = "All_HAR_coordinates") 

# Neanderthal and Denisovan in HARs ---------------------------------------

## Neanderthal

ntSssSnps_intersected <- read.table(file = file.path(sys_dir,"HAR/Neanderthal/ntSssSnps_intersected.bed"), header = F)
ntSssZScorePMVar_intersected <- read.table(file = file.path(sys_dir,"HAR/Neanderthal/ntSssZScorePMVar_intersected.bed"), header = F)
colnames(ntSssSnps_intersected) <- c("chr","start_var","stop_var","genotype","score_1","strand","start_var","stop_var","score_2","chr","start","stop","id")
colnames(ntSssZScorePMVar_intersected) <- c("chr","start_var","stop_var","info","Zscore+variance","chr","start","stop","id")

odd_rows <- ntSssZScorePMVar_intersected[seq(1, nrow(ntSssZScorePMVar_intersected), 2), ] 
even_rows <- ntSssZScorePMVar_intersected[seq(2, nrow(ntSssZScorePMVar_intersected), 2), ] %>% dplyr::select(c("Zscore+variance"))
ntSssZScorePMVar_comb <- cbind(odd_rows,even_rows)
colnames(ntSssZScorePMVar_comb) <- c("chr","start_var","stop_var","info","Zscore","chr","start","stop","id","variance")
ntSssZScorePMVar_comb$variance <- abs(ntSssZScorePMVar_comb$variance)
ntSssSnps_intersected_info <- ntSssSnps_intersected %>% dplyr::select(c("genotype"))
neanderthal_snps <- cbind(ntSssZScorePMVar_comb,ntSssSnps_intersected_info) %>%
  dplyr::select(c("id","chr","start","stop","start_var","stop_var","genotype","Zscore","variance"))

## Denisovan

Denisovan_Pinky_intersected <- read.table(file = file.path(sys_dir,"HAR/Denisovan/Denisovan_Pinky_intersected.bed"), header = F)
colnames(Denisovan_Pinky_intersected) <- c("chr","start_var","stop_var","snp_id","chr_2","start","stop","snp_id_2","sample_info","id")
denisovan_snps <- Denisovan_Pinky_intersected %>%
  dplyr::select(c("chr","start_var","stop_var","snp_id","sample_info","id"))
length(unique(denisovan_snps$snp_id))
length(unique(denisovan_snps$id))

# Transcription factor binding site analysis ------------------------------

## The meme suite - FIMO

## human

Fimo <- importFimo(file.path(sys_dir,"HAR/TFBS/fimo_out/fimo.tsv"))
Fimo_df <- Fimo %>% as.data.frame()
Fimo_df$matched_sequence_length <- nchar(Fimo_df$matched_sequence)
Fimo_lowestq <- Fimo_df %>% filter(Fimo_df$qvalue < 0.05)

mean(Fimo_lowestq$matched_sequence_length)
sd(Fimo_lowestq$matched_sequence_length)

table(Fimo_lowestq$matched_sequence_length)
summary(Fimo_lowestq$matched_sequence_length)

Fimo_lowestq_highestms <- Fimo_lowestq %>% filter(Fimo_lowestq$matched_sequence_length == 21)
# plot_sequence_heatmap(Fimo_lowestq_highestms$matched_sequence)
Fimo_lowestq_freqms <- Fimo_lowestq %>% filter(Fimo_lowestq$matched_sequence_length == 11)
plot_sequence_heatmap(Fimo_lowestq_freqms$matched_sequence)

Fimo_df$motif_alt_id_count <- Fimo_df

Fimo_motif <- Fimo
colnames(Fimo_motif$matched_sequence) <- c("seqnames","start","end","width","strand","motif_id","motif_alt_id","score","pvalue","qvalue","motif","matched_sequence_length")

Fimo_top1_perhar <- Fimo_df %>%
  group_by(seqnames) %>%
  filter(qvalue <= quantile(qvalue, 0.01))
length(unique(Fimo_top1_perhar$motif_alt_id))

Fimo_df$num <- seq(1:nrow(Fimo_df))

ggplot()+
  geom_line(Fimo_df,mapping=aes(x=num,y=score),color="purple")+
  labs(x="TF",y="score")+
  theme_bw()

## does the above code works well? Yes

# Fimo_HAR2520 <- Fimo_df %>% filter(seqnames == "HARsv2_2520")
# summary(Fimo_HAR2520$qvalue)
# quantile(Fimo_HAR2520$qvalue, 0.01)
# Fimo_HAR2514_top1 <- Fimo_HAR2520 %>% filter(qvalue <= quantile(qvalue, 0.01))

## chimp

Fimo_chimp <- importFimo(file.path(sys_dir,"HAR/TFBS/fimo_out_chimp/fimo.tsv"))
Fimo_chimp_df <- Fimo_chimp %>% as.data.frame()
Fimo_chimp_df$matched_sequence_length <- nchar(Fimo_chimp_df$matched_sequence)
Fimo_chimp_lowestq <- Fimo_chimp_df %>% filter(Fimo_chimp_df$qvalue < 0.05)

mean(Fimo_chimp_lowestq$matched_sequence_length)
sd(Fimo_chimp_lowestq$matched_sequence_length)

table(Fimo_chimp_lowestq$matched_sequence_length)
summary(Fimo_chimp_lowestq$matched_sequence_length)

Fimo_chimp_top1_perhar <- Fimo_chimp_df %>%
  group_by(seqnames) %>%
  filter(qvalue <= quantile(qvalue, 0.01))
length(unique(Fimo_chimp_top1_perhar$motif_alt_id))

## compare human and chimp

Fimo_top1_perhar
Fimo_chimp_top1_perhar 

unique(Fimo_top1_perhar$motif_alt_id)
unique(Fimo_chimp_top1_perhar$motif_alt_id)
unique(Fimo_top1_perhar$matched_sequence)
unique(Fimo_chimp_top1_perhar$matched_sequence)

# comparison output: we have 5 different TF in HAR_15, HAR_17, HAR_19: 
# ONECUT1 -> same seq
# NR2F1 -> same seq
# ZNF85 -> seq is too short, Fimo_chimps filtered it out
# NRL -> diff; we lost it in a HAR snp
# TCF7L1 -> diff; we lost it in a HAR snp

Fimo_chimp_df_fil <- Fimo_chimp_df %>% dplyr::select(c("motif_alt_id","pvalue"))
colnames(Fimo_chimp_df_fil) <- c("motif_alt_id","pvalue_chimp")
Fimo_HC <- merge(Fimo_df,Fimo_chimp_df_fil)
Fimo_HC$pvalue_delta <- (-log10(abs((Fimo_HC$pvalue)-(Fimo_HC$pvalue_chimp)))) 
Fimo_HC$names <- seq(1,nrow(Fimo_HC))

unique_to_Fimo_df <- anti_join(Fimo_df, Fimo_chimp_df, by = c("motif_alt_id"))
unique_to_Fimo_chimp_df <- anti_join(Fimo_chimp_df, Fimo_df, by = c("motif_alt_id"))

cat(unique_to_Fimo_df$motif_alt_id, sep = "\n")
cat(unique_to_Fimo_chimp_df$motif_alt_id, sep = "\n")

ggplot()+
  # geom_point(Fimo_HC, mapping=aes(x=names,y=-log10(pvalue)),col="steelblue")+
  geom_point(Fimo_HC, mapping=aes(x=names,y=-log10(pvalue_chimp)),col="orange")
  # geom_point(Fimo_HC, mapping=aes(x=names,y=pvalue_delta),col="gray")+
  theme_classic()

ggplot()+
  geom_point(Fimo_top1_perhar, mapping=aes(y=motif_alt_id,x=-log10(qvalue)),col="steelblue")+
  theme_linedraw()

plot(Fimo_HC$names,Fimo_HC$pvalue_chimp)
plot(Fimo_HC$names,Fimo_HC$pvalue_delta)

length(unique(Fimo_df$motif_id))
length(unique(Fimo_chimp_df$motif_id))

# TFBS analysis - TF gene-set enrichment -------------------------------------

## Enrichr

## take top 5% of fimo_df

Fimo_top10_perhar <- Fimo_df %>%
  group_by(seqnames) %>%
  filter(qvalue <= quantile(qvalue, 0.10))
cat(Fimo_top10_perhar$motif_alt_id, sep = "\n")

Fimo_pcut4_perhar <- Fimo_df %>%
  group_by(seqnames) %>%
  filter(qvalue <= quantile(pvalue, 0.05))

## use online version at https://maayanlab.cloud/Enrichr/ 

## topGO

library(topGO)
library(ALL)
data(ALL)
data(geneList)

affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)

Fimo_top10_perhar <- Fimo_df %>%
  group_by(seqnames) %>%
  filter(qvalue <= quantile(qvalue, 0.10)) %>% 
  as.data.frame() %>%
  dplyr::select(c("motif_alt_id","qvalue"))

topGOresultFisher <- runTest(geneList, algorithm = "elim", statistic = "fisher", ontology="BP", nodeSize = 20)

# TFBS analysis - SNPs affecting TFBS (atSNP & motifbreakR) ----------------------------------------------------------------

## atSNP

## get jaspar motif id for lost TF in humans (TCF7L1 <- MA1421.1, NRL <- MA0842.3)

pwms_jaspar <- LoadMotifLibrary(filename = file.path(sys_dir,"HAR/TFBS/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt"))
# pwms <- LoadMotifLibrary(urlname="https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt")

target_motifs <- c("MA1421.1", "MA0842.3")
pwms_jaspar_tcf_nrl <- pwms_jaspar[names(pwms_jaspar) %in% target_motifs]

## get seq info of human and chimp of HAR19

snp_info <- LoadFastaData(ref.filename=file.path(sys_dir,"HAR/TFBS/h_HAR19.fasta"), 
                          snp.filename=file.path(sys_dir,"HAR/TFBS/c_HAR19.fasta"),
                          default.par = T)
str(snp_info)

## run atSNP
                         
atsnp.scores <- ComputeMotifScore(pwms_jaspar_tcf_nrl,
                                  snp_info,
                                  ncores = 1)
atsnp.scores$snp.tbl
atsnp.scores$motif.scores
atsnp.result <- ComputePValues(motif.lib = pwms_jaspar_tcf_nrl,
                               snp.info = snp_info,
                               motif.scores = atsnp.scores$motif.scores,
                               ncores = 1, testing.mc=TRUE)
atsnp.result
head(atsnp.result[order(atsnp.result$pval_rank), c("snpid", "motif", "pval_ref", "pval_snp", "pval_rank")])

## motifbreakR

## get snps in HAR in the right format

library(BSgenome.Hsapiens.UCSC.hg38)
HAR_snp <- file.path(sys_dir,"HAR/HAR_maf_vcf/merged_HAR_snp.bed")
read.table(HAR_snp, header = FALSE)
HAR_snp_gr <- snps.from.file(file = HAR_snp,
                                  search.genome = BSgenome.Hsapiens.UCSC.hg38,
                                  format = "bed")
## output: 35 snps in HARs_grik2 (not 41, after the correct alignment of HAR_20)

HARsnp_relative_loc 
colnames(HARsnp_relative_loc) <- c("chr","start","end","id","zero","sign")
HARsnp_relative_loc$nearest_snp <- sapply(HARsnp_relative_loc$start, function(x) {
  diffs <- abs(HARsnp_relative_loc$start - x)
  min_diff <- min(diffs[diffs > 0])
  ifelse(is.infinite(min_diff), 0, min_diff)
})
summary(HARsnp_relative_loc$nearest_snp)

# all_pos <- read.table(HAR_snp, header = FALSE)
# all_pos_2 <- all_pos$V3
# HAR_snp_gr_id <- unlist(HAR_snp_gr$SNP_id)
# not_included <- all_pos_2[!sapply(all_pos_2, function(x) any(grepl(x, HAR_snp_gr$SNP_id)))]
# print(not_included)
# HAR_snp <- read.table(file.path(sys_dir,"HAR/HAR_maf_vcf/merged_HAR_snp.bed"), header = FALSE)
# colnames(HAR_snp) <- c("chr","start","end","snp","ref","alt")
# HAR_snp_gr <- makeGRangesFromDataFrame(df=HAR_snp,
#                                      keep.extra.columns=T,
#                                      ignore.strand=T,
#                                      seqinfo=Seqinfo(genome="hg38"),
#                                      seqnames.field="chr",
#                                      start.field="start",
#                                      end.field="end",
#                                      starts.in.df.are.0based=T)
# names(HAR_snp_gr) <- HAR_snp_gr$snp

## get JASPAR motif in the right format

motifs_jaspar2022 <- query(MotifDb,
                           andStrings=c("hsapiens"),
                           orStrings=c("jaspar2022"))

## run motifbreakR

## HARsv2_2519 and chr6:103468164:G:C

motifbreakR_HAR19_res <- motifbreakR(snpList = HAR_snp_gr[25],
                                   pwmList = motifs_jaspar2022,
                                   threshold = 1e-4,
                                   method = "ic",
                                   bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                   BPPARAM = BiocParallel::bpparam())
plotMB(results = motifbreakR_HAR19_res[motifbreakR_HAR19_res$geneSymbol == "NRL"], rsid = "chr6:103468165:G:C", effect = "strong")

calculatePvalue(motifbreakR_results[motifbreakR_results$geneSymbol == "NRL"])

## HARsv2_2517 and chr6:102249170:A:G

motifbreakR_HAR17_res <- motifbreakR(snpList = HAR_snp_gr[14],
                                   pwmList = motifs_jaspar2022,
                                   threshold = 1e-4,
                                   method = "ic",
                                   bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                   BPPARAM = BiocParallel::bpparam())
# calculatePvalue(motifbreakR_HAR17_res[motifbreakR_HAR17_res$effect == "strong"])
summary(motifbreakR_HAR17_res$pctAlt)
motifres_pctAlt_99pen <- quantile(motifbreakR_HAR17_res$pctAlt, 0.99)
plotMB(results = motifbreakR_HAR17_res[motifbreakR_HAR17_res$pctAlt > motifres_pctAlt_99pen], rsid = "chr6:102249171:A:G", effect = "strong")

## per HAR and extract significant

HAR_snps_order <- c(6, 3, 4, 5, 4, 3, 10)
L1 <- 1
L2 <- 7
L3 <- 10
L4 <- 14
L5 <- 19
L6 <- 23
L7 <- 25
HAR14_snp_gr <- HAR_snp_gr[L1:(L1+HAR_snps_order[1]-1)]
HAR15_snp_gr <- HAR_snp_gr[L2:(L2+HAR_snps_order[2]-1)]
HAR16_snp_gr <- HAR_snp_gr[L3:(L3+HAR_snps_order[3]-1)]
HAR17_snp_gr <- HAR_snp_gr[L4:(L4+HAR_snps_order[4]-1)]
HAR18_snp_gr <- HAR_snp_gr[L5:(L5+HAR_snps_order[5]-1)]
HAR19_snp_gr <- HAR_snp_gr[L6:(L6+HAR_snps_order[6]-1)]
HAR20_snp_gr <- HAR_snp_gr[L7:(L7+HAR_snps_order[7]-1)]

motifbreakR_reslist <- list()
TF_sig_HAR_list <- list()

for (i in 14:20) {
  HAR_snp_gr <- get(paste0("HAR", i, "_snp_gr"))
  HAR_res <- motifbreakR(snpList = HAR_snp_gr,
                         pwmList = motifs_jaspar2022,
                         threshold = 1e-4,
                         method = "ic",
                         bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                         BPPARAM = BiocParallel::bpparam())
  
  ## extract quantiles (0.01+0.99)
  HAR_res_sig <- HAR_res[HAR_res$alleleEffectSize > quantile(HAR_res$alleleEffectSize, 0.99) | HAR_res$alleleEffectSize <  quantile(HAR_res$alleleEffectSize, 0.01),]
  ## extract the top 99% pctAlt scores
  # HAR_res_sig <- HAR_res[HAR_res$pctAlt > quantile(HAR_res$pctAlt, 0.99),]
  ## extract the two highest pctAlt scores
  # highest_indices <- order(HAR_res$pctAlt, decreasing = TRUE)[1:2] 
  # HAR_res_sig <- HAR_res[highest_indices, ]
  ## extract if the scaled alleleEffectSize is > 0.4
  # HAR_res_sig <- HAR_res[scale(HAR_res$alleleEffectSize) > 0.4,]
  ## extract if the scaled alleleDiff is > 0.4
  # HAR_res_sig <- HAR_res[(scale(HAR_res$scoreRef)-scale(HAR_res$scoreAlt)) > 0.4,]
  
  motifbreakR_reslist[[paste0("motifbreakR_HAR", i, "_res_sig")]] <- HAR_res_sig
}

for (i in 14:20) {
  HAR_res_sig <- motifbreakR_reslist[[paste0("motifbreakR_HAR", i, "_res_sig")]]
  TF_sig_HAR_list[[paste0("TF_sig_HAR", i)]] <- data.frame(
    geneSymbol = HAR_res_sig$geneSymbol,
    alleleEffectSize = HAR_res_sig$alleleEffectSize
  )
}

## use encode_RPKML_avgexp and encode_RPKML_avgexp_human instead of encode_RPKML_age_filtered

## for macaque and human gene expression

encode_RPKML_avgexp_macaque <- encode_RPKML_avgexp %>% filter(Species == "Macaque")
encode_RPKML_avgexp_macaque$RPKM_Zscore <- ((encode_RPKML_avgexp_macaque$avg_RPKM - (mean(encode_RPKML_avgexp_macaque$avg_RPKM)))/(sd(encode_RPKML_avgexp_macaque$avg_RPKM)))
colnames(encode_RPKML_avgexp_macaque) <- c("Species","avg_RPKM","sd_RPKM","ensemblID","geneSymbol","RPKM_Zscore")

encode_RPKML_avgexp_human_highexp <- encode_RPKML_avgexp_human %>% filter(avg_RPKM > mean(encode_RPKML_avgexp_human$avg_RPKM))
synaptogenesis_genes_human <- encode_RPKML_avgexp_human_highexp$geneSymbol
writeLines(synaptogenesis_genes_human, "synaptogenesis_genes_human.txt")

encode_RPKML_avgexp_macaque_highexp <- encode_RPKML_avgexp_macaque %>% filter(avg_RPKM > mean(encode_RPKML_avgexp_macaque$avg_RPKM))
synaptogenesis_genes_macaque <- encode_RPKML_avgexp_macaque_highexp$geneSymbol
writeLines(synaptogenesis_genes_macaque, "synaptogenesis_genes_macaque.txt")



TF_allHAR <- bind_rows(TF_sig_HAR_list, .id = "HAR_ID")
TF_allHAR_brainRPKML <- merge(TF_allHAR, encode_RPKML_avgexp_human)
TF_allHAR_brainRPKML$Effect <- ifelse(TF_allHAR_brainRPKML$alleleEffectSize > 0, "improving",
                                      ifelse(TF_allHAR_brainRPKML$alleleEffectSize < 0, "disrupting", "none"))

TF_allHAR_brainRPKML_highexp <- TF_allHAR_brainRPKML %>% filter(avg_RPKM > mean(TF_allHAR_brainRPKML$avg_RPKM))
nrow(TF_allHAR_brainRPKML)
nrow(TF_allHAR_brainRPKML_highexp)

TF_allHAR_brainRPKML_highexp_unique <- TF_allHAR_brainRPKML_highexp[!duplicated(TF_allHAR_brainRPKML_highexp$geneSymbol) &
                                                                      !duplicated(TF_allHAR_brainRPKML_highexp$geneSymbol, fromLast = TRUE), ]
nrow(TF_allHAR_brainRPKML_highexp)
nrow(TF_allHAR_brainRPKML_highexp_unique)

ggplot()+
  geom_rect(aes(xmin = mean(TF_allHAR_brainRPKML_highexp$avg_RPKM) - sd(TF_allHAR_brainRPKML_highexp$avg_RPKM), 
                xmax = mean(TF_allHAR_brainRPKML_highexp$avg_RPKM) + sd(TF_allHAR_brainRPKML_highexp$avg_RPKM), 
                ymin = -Inf, ymax = Inf), 
            fill = "gray", alpha = 0.8) + 
  geom_vline(xintercept = mean(TF_allHAR_brainRPKML_highexp$avg_RPKM), linetype = "solid", color = "black", size = 1)+
  geom_point(data=subset(TF_allHAR_brainRPKML_highexp, Effect == "improving"), mapping=aes(x=avg_RPKM,y=reorder(geneSymbol,avg_RPKM),col=Effect))+
  geom_errorbar(data=subset(TF_allHAR_brainRPKML_highexp, Effect == "improving"), mapping=aes(xmin=avg_RPKM-sd_RPKM, xmax=avg_RPKM+sd_RPKM,y=geneSymbol,col=Effect), width=0.3)+
  geom_point(data=subset(encode_RPKML_avgexp_human, geneSymbol == "GRIK2"), mapping=aes(x=avg_RPKM,y=reorder(geneSymbol,avg_RPKM),color=geneSymbol))+
  geom_errorbar(data=subset(encode_RPKML_avgexp_human, geneSymbol == "GRIK2"), mapping=aes(xmin=avg_RPKM-sd_RPKM, xmax=avg_RPKM+sd_RPKM,y=geneSymbol,color=geneSymbol), width=0.3)+
  # scale_color_manual(values = c("#FC766AFF","black","#5B84B1FF"))+
  scale_color_manual(values = c("black","#5B84B1FF"))+
  labs(y="TF",x="RPKM")+
  theme_linedraw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

ggplot()+
  geom_rect(aes(xmin = mean(TF_allHAR_brainRPKML_highexp$avg_RPKM) - sd(TF_allHAR_brainRPKML_highexp$avg_RPKM), 
                xmax = mean(TF_allHAR_brainRPKML_highexp$avg_RPKM) + sd(TF_allHAR_brainRPKML_highexp$avg_RPKM), 
                ymin = -Inf, ymax = Inf), 
            fill = "gray", alpha = 0.8) + 
  geom_vline(xintercept = mean(TF_allHAR_brainRPKML_highexp$avg_RPKM), linetype = "solid", color = "black", size = 1)+
  geom_point(data=subset(TF_allHAR_brainRPKML_highexp, Effect == "disrupting"), mapping=aes(x=avg_RPKM,y=reorder(geneSymbol,avg_RPKM),col=Effect))+
  geom_errorbar(data=subset(TF_allHAR_brainRPKML_highexp, Effect == "disrupting"), mapping=aes(xmin=avg_RPKM-sd_RPKM, xmax=avg_RPKM+sd_RPKM,y=geneSymbol,col=Effect), width=0.3)+
  geom_point(data=subset(encode_RPKML_avgexp_human, geneSymbol == "GRIK2"), mapping=aes(x=avg_RPKM,y=reorder(geneSymbol,avg_RPKM),color=geneSymbol))+
  geom_errorbar(data=subset(encode_RPKML_avgexp_human, geneSymbol == "GRIK2"), mapping=aes(xmin=avg_RPKM-sd_RPKM, xmax=avg_RPKM+sd_RPKM,y=geneSymbol,color=geneSymbol), width=0.3)+
  # scale_color_manual(values = c("#FC766AFF","black","#5B84B1FF"))+
  scale_color_manual(values = c("#FC766AFF","black"))+ 
  labs(y="TF",x="RPKM")+
  theme_linedraw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

# ggplot()+
#   geom_point(data=subset(TF_allHAR_brainRPKML_highexp, Effect == "disrupting"), mapping=aes(x=avg_RPKM,y=reorder(geneSymbol,avg_RPKM),col=Effect),size=0.000001)+
#   geom_errorbar(data=subset(TF_allHAR_brainRPKML_highexp, Effect == "disrupting"), mapping=aes(xmin=avg_RPKM-sd_RPKM, xmax=avg_RPKM+sd_RPKM,y=geneSymbol,col=Effect),width=0.0000001,size=0.000001)+
#   geom_rect(aes(xmin = mean(TF_allHAR_brainRPKML_highexp$avg_RPKM) - sd(TF_allHAR_brainRPKML_highexp$avg_RPKM),
#                 xmax = mean(TF_allHAR_brainRPKML_highexp$avg_RPKM) + sd(TF_allHAR_brainRPKML_highexp$avg_RPKM),
#                 ymin = -Inf, ymax = Inf),
#             fill = "gray", alpha = 0.8) +
#   geom_vline(xintercept = mean(TF_allHAR_brainRPKML_highexp$avg_RPKM), linetype = "solid", color = "black", size = 1)+
#   geom_point(data=subset(encode_RPKML_avgexp_human, geneSymbol == "GRIK2"), mapping=aes(x=avg_RPKM,y=reorder(geneSymbol,avg_RPKM),color=geneSymbol))+
#   geom_errorbar(data=subset(encode_RPKML_avgexp_human, geneSymbol == "GRIK2"), mapping=aes(xmin=avg_RPKM-sd_RPKM, xmax=avg_RPKM+sd_RPKM,y=geneSymbol,color=geneSymbol),width=0.3)+
#   # scale_color_manual(values = c("#FC766AFF","black","#5B84B1FF"))+
#   scale_color_manual(values = c("gray","black"))+
#   labs(y="TF",x="RPKM")+
#   theme_linedraw()+
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank())

## subset per HAR

ggplot()+
  geom_rect(aes(xmin = mean(TF_allHAR_brainRPKML_highexp$avg_RPKM) - sd(TF_allHAR_brainRPKML_highexp$avg_RPKM), 
                xmax = mean(TF_allHAR_brainRPKML_highexp$avg_RPKM) + sd(TF_allHAR_brainRPKML_highexp$avg_RPKM), 
                ymin = -Inf, ymax = Inf), 
            fill = "gray", alpha = 0.8) + 
  geom_vline(xintercept = mean(TF_allHAR_brainRPKML_highexp$avg_RPKM), linetype = "solid", color = "black", size = 1)+
  geom_point(data=subset(TF_allHAR_brainRPKML_highexp,HAR_ID == "TF_sig_HAR20"),
             mapping=aes(x=avg_RPKM,y=reorder(geneSymbol,avg_RPKM),col=Effect))+
  geom_errorbar(data=subset(TF_allHAR_brainRPKML_highexp, HAR_ID == "TF_sig_HAR20"),
                mapping=aes(xmin=avg_RPKM-sd_RPKM, xmax=avg_RPKM+sd_RPKM,y=geneSymbol,col=Effect), width=0.1)+
  scale_color_brewer()+
  labs(y="TF",x="RPKM")+
  facet_wrap(~HAR_ID, scales = "free")+
  theme_linedraw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank(), 
        strip.background = element_rect(color = "black", fill = NA),
        strip.text = element_text(size = 12),
        legend.position = "none")

p14 <- plotMB(results = (motifbreakR_reslist$motifbreakR_HAR14_res_sig),
              rsid = unique(motifbreakR_reslist$motifbreakR_HAR14_res_sig$SNP_id),
              effect = "strong")
p15 <- plotMB(results = (motifbreakR_reslist$motifbreakR_HAR15_res_sig),
              rsid = unique(motifbreakR_reslist$motifbreakR_HAR15_res_sig$SNP_id),
              effect = "strong")
p16.a <- plotMB(results = (motifbreakR_reslist$motifbreakR_HAR16_res_sig),
              rsid = (motifbreakR_reslist$motifbreakR_HAR16_res_sig$SNP_id[1]),
              effect = "strong")
p16.b <- plotMB(results = (motifbreakR_reslist$motifbreakR_HAR16_res_sig),
                rsid = (motifbreakR_reslist$motifbreakR_HAR16_res_sig$SNP_id[2]),
                effect = "strong")
p17.a <- plotMB(results = (motifbreakR_reslist$motifbreakR_HAR17_res_sig),
              rsid = (motifbreakR_reslist$motifbreakR_HAR17_res_sig$SNP_id[1]),
              effect = "strong")
p17.b <- plotMB(results = (motifbreakR_reslist$motifbreakR_HAR17_res_sig),
              rsid = (motifbreakR_reslist$motifbreakR_HAR17_res_sig$SNP_id[2]),
              effect = "strong")
p18 <- plotMB(results = (motifbreakR_reslist$motifbreakR_HAR18_res_sig),
              rsid = unique(motifbreakR_reslist$motifbreakR_HAR18_res_sig$SNP_id),
              effect = "strong")
p19.a <- plotMB(results = (motifbreakR_reslist$motifbreakR_HAR19_res_sig),
              rsid = (motifbreakR_reslist$motifbreakR_HAR19_res_sig$SNP_id[1]),
              effect = "strong")
p19.b <- plotMB(results = (motifbreakR_reslist$motifbreakR_HAR19_res_sig),
                rsid = (motifbreakR_reslist$motifbreakR_HAR19_res_sig$SNP_id[2]),
                effect = "strong")
p20 <- plotMB(results = (motifbreakR_reslist$motifbreakR_HAR20_res_sig),
              rsid = unique(motifbreakR_reslist$motifbreakR_HAR20_res_sig$SNP_id),
              effect = "strong")

TF_sig_HAR <- c(motifbreakR_reslist$motifbreakR_HAR14_res_sig$geneSymbol,
                motifbreakR_reslist$motifbreakR_HAR15_res_sig$geneSymbol,
                motifbreakR_reslist$motifbreakR_HAR16_res_sig$geneSymbol,
                motifbreakR_reslist$motifbreakR_HAR17_res_sig$geneSymbol,
                motifbreakR_reslist$motifbreakR_HAR18_res_sig$geneSymbol,
                motifbreakR_reslist$motifbreakR_HAR19_res_sig$geneSymbol,
                motifbreakR_reslist$motifbreakR_HAR20_res_sig$geneSymbol)

## include encode transcription data for all genes from /home/shadi/Desktop/S4_Project/PsychENCODE/PsychENCODE.R

## is TF_sig_HAR[#] present in encode_RPKML_age_filtered?

ggplot()+
  geom_point(data=subset(encode_RPKML_age_filtered, grepl(TF_sig_HAR[11], gene_name)),
             mapping=aes(x=Predicted_Age_log,y=RPKM_log,fill=Species,col=Species),size=0.1,shape=23)+
  geom_smooth(data=subset(encode_RPKML_age_filtered, grepl(TF_sig_HAR[11], gene_name)),
              mapping=aes(x=Predicted_Age_log,y=RPKM_log,col=Species),method = "loess")+
  labs(x="Predicted_Age",y="log2(RPKM+1)",labs="Species")+
  theme_classic()
# + facet_wrap(~Tissue)

for (TF in TF_sig_HAR) {
  matches <- encode_RPKML_age_filtered[grepl(TF, encode_RPKML_age_filtered$gene_name), ]
  encodeRPKML_TFsig <- rbind(encodeRPKML_TFsig, matches)
}
unique(encodeRPKML_TFsig$gene_name)
split_gene_name <- strsplit(as.character(encodeRPKML_TFsig$gene_name), split = "\\|")
split_df <- do.call(rbind, split_gene_name)
encodeRPKML_TFsig$ensemblID <- split_df[, 1]
encodeRPKML_TFsig$gene_symbol <- split_df[, 2]
encodeRPKML_TFsig$gene_name <- NULL
encodeRPKML_TFsig <- encodeRPKML_TFsig %>% filter(gene_symbol == TF_sig_HAR)
unique(encodeRPKML_TFsig$gene_symbol)

mean_RPKM <- mean(encodeRPKML_TFsig$RPKM)
sd_RPKM <- sd(encodeRPKML_TFsig$RPKM)
encodeRPKML_TFsig$RPKM_Zscore <- ((encodeRPKML_TFsig$RPKM - mean_RPKM)/sd_RPKM)

ggplot()+
  geom_rect(aes(xmin = mean(encodeRPKML_TFsig$RPKM_Zscore) - sd(encodeRPKML_TFsig$RPKM_Zscore), 
                xmax = mean(encodeRPKML_TFsig$RPKM_Zscore) + sd(encodeRPKML_TFsig$RPKM_Zscore), 
                ymin = -Inf, ymax = Inf), 
            fill = "gray", alpha = 0.5) + 
  geom_vline(xintercept = mean(encodeRPKML_TFsig$RPKM_Zscore), linetype = "solid", color = "black", size = 1)+
  geom_boxplot(data=subset(encodeRPKML_TFsig,Species = "Human"),mapping=aes(x=RPKM_Zscore,y=gene_symbol,fill="tomato"))+
  labs(y="TF",x="scaled RPKM")+
  theme_linedraw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

## filter based on motifbreakR & encode transcription data results

## test best filtration method for one HAR

## fast and dirty way to get TF expression in the human brain during development
rm(encode_RPKM,encode_RPKM_table,encode_RPKML_age)
encode_RPKML_human <- encode_RPKML_age_filtered %>%
  filter(encode_RPKML_age_filtered$Species == "Human")
rm(encode_RPKML_age_filtered)
encode_RPKML_gene_symbol <- encode_RPKML_human %>%
  group_by(gene_name) %>%
  summarize(average_expression = mean(RPKM))
encode_RPKML_gene_symbol <- encode_RPKML_gene_symbol %>% 
  filter(average_expression > 0)
split_gene_name <- strsplit(as.character(encode_RPKML_gene_symbol$gene_name), split = "\\|")
split_df <- do.call(rbind, split_gene_name)
encode_RPKML_gene_symbol$ensemblID <- split_df[, 1]
encode_RPKML_gene_symbol$gene_symbol <- split_df[, 2]
encode_RPKML_gene_symbol$gene_name <- NULL
encode_RPKML_gene_symbol <- encode_RPKML_gene_symbol %>% 
  dplyr::select(c("gene_symbol","average_expression")) %>%
  arrange(-average_expression) %>% 
  as.data.frame() %>%
  filter(average_expression > mean(average_expression)) 

for (TF in encode_RPKML_gene_symbol$gene_symbol) {
  subset_motifs <- subset(motifs_jaspar2022, geneSymbol == (TF))
  if (length(subset_motifs) == 0) {
    assign(paste0("subset_motifs_", TF), subset_motifs)
  }
}

motifbreakR_HAR14_test <- motifbreakR(HAR14_snp_gr,
                       pwmList = motifs_jaspar2022,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::bpparam())
quantile_99_test <- quantile(motifbreakR_HAR14_test$alleleEffectSize, 0.99)
quantile_1_test <- quantile(motifbreakR_HAR14_test$alleleEffectSize, 0.01)
motifbreakR_HAR14_test_sig <- motifbreakR_HAR14_test[motifbreakR_HAR14_test$alleleEffectSize > quantile_95_test | motifbreakR_HAR14_test$alleleEffectSize < quantile_5_test,]
motifbreakR_HAR14_test_sig <- motifbreakR_HAR14_test_sig[scale(motifbreakR_HAR14_test_sig$alleleEffectSize) > 0.4,]
motifbreakR_HAR14_test_sig <- motifbreakR_HAR14_test_sig[motifbreakR_HAR14_test_sig$effect == "strong",]
length(motifbreakR_HAR14_test)
length(motifbreakR_HAR14_test_sig)
TF_sig_HAR14_test <- motifbreakR_HAR14_test_sig$geneSymbol

encode_TF_sig_HAR14_test <- data.frame()
for (TF in TF_sig_HAR14_test) {
  matches <- encode_RPKML_gene_symbol[encode_RPKML_gene_symbol$gene_symbol == TF,]
  encode_TF_sig_HAR14_test <- rbind(encode_TF_sig_HAR14_test, matches)
}

TF_sig_encode_TF_sig_HAR14_test <- encode_TF_sig_HAR14_test$gene_symbol
encode_RPKML_sig_test <- data.frame()
for (TF in TF_sig_encode_TF_sig_HAR14_test) {
  matches <- encode_RPKML_age_filtered[grepl(TF, encode_RPKML_age_filtered$gene_name), ]
  encode_RPKML_sig_test <- rbind(encode_RPKML_sig_test, matches)
}

split_gene_name <- strsplit(as.character(encode_RPKML_sig_test$gene_name), split = "\\|")
split_df <- do.call(rbind, split_gene_name)
encode_RPKML_sig_test$ensemblID <- split_df[, 1]
encode_RPKML_sig_test$gene_symbol <- split_df[, 2]
encode_RPKML_sig_test$gene_name <- NULL
encode_RPKML_sig_test$RPKM_Zscore <- ((encode_RPKML_sig_test$RPKM - (mean(encode_RPKML_sig_test$RPKM)))/(sd(encode_RPKML_sig_test$RPKM)))

ggplot()+
  # geom_rect(aes(xmin = mean(encode_RPKML_sig_test$RPKM_log) - sd(encode_RPKML_sig_test$RPKM_log),
  #               xmax = mean(encode_RPKML_sig_test$RPKM_log) + sd(encode_RPKML_sig_test$RPKM_log),
  #               ymin = -Inf, ymax = Inf),
  #           fill = "gray", alpha = 0.5) +
  geom_rect(aes(xmin = mean(log2(encode_RPKML_gene_symbol$average_expression)) - sd(log2(encode_RPKML_gene_symbol$average_expression)),
                xmax = mean(log2(encode_RPKML_gene_symbol$average_expression)) + sd(log2(encode_RPKML_gene_symbol$average_expression)),
                ymin = -Inf, ymax = Inf),
            fill = "gray", alpha = 0.5) +
  geom_vline(xintercept = mean(log2(encode_RPKML_gene_symbol$average_expression)), linetype = "solid", color = "black", size = 1)+
  geom_boxplot(data=subset(encode_RPKML_sig_test, Species="Human"),
               mapping=aes(x=RPKM_log,y=reorder(gene_symbol,RPKM_log)),outlier.shape = "-")+
  labs(y="TF",x="scaled RPKM")+
  theme_linedraw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

TF_HAR14_test_alleleEffectSize <- data.frame(gene_symbol = motifbreakR_HAR14_test_sig$geneSymbol,
                                    alleleEffectSize = motifbreakR_HAR14_test_sig$alleleEffectSize)
encode_RPKML_sig_alleleEffectSize_test <- merge(TF_HAR14_test_alleleEffectSize,encode_RPKML_sig_test)
encode_RPKML_sig_alleleEffectSize_test$alleleEffectSize <- as.numeric(encode_RPKML_sig_alleleEffectSize_test$alleleEffectSize)
encode_RPKML_sig_alleleEffectSize_test_human <- encode_RPKML_sig_alleleEffectSize_test %>% filter(Species == "Human")
encode_RPKML_sig_alleleEffectSize_test_human$alleleEffectSize <- as.character(encode_RPKML_sig_alleleEffectSize_test_human$alleleEffectSize)

median_lines <- encode_RPKML_sig_alleleEffectSize_test %>%
  filter(Species == "Human") %>%
  group_by(gene_symbol) %>%
  summarize(median_RPKM_log = median(RPKM_log), 
            median_alleleEffectSize = median(alleleEffectSize)) %>%
  ungroup()

mean_points <- encode_RPKML_sig_alleleEffectSize_test %>%
  filter(Species == "Human") %>%
  group_by(gene_symbol) %>%
  summarize(mean_RPKM_log = mean(RPKM_log), 
            mean_alleleEffectSize = mean(alleleEffectSize)) %>%
  ungroup()

ggplot()+
  geom_rect(aes(xmin = mean(log2(encode_RPKML_gene_symbol$average_expression)) - sd(log2(encode_RPKML_gene_symbol$average_expression)),
                xmax = mean(log2(encode_RPKML_gene_symbol$average_expression)) + sd(log2(encode_RPKML_gene_symbol$average_expression)),
                ymin = -Inf, ymax = Inf),
            fill = "gray", alpha = 0.5) +
  geom_vline(xintercept = mean(log2(encode_RPKML_gene_symbol$average_expression)), linetype = "solid", color = "black", size = 1)+
  geom_boxplot(encode_RPKML_sig_alleleEffectSize_test_human,
             mapping=aes(x=(RPKM_log),y=reorder(gene_symbol,RPKM_log)),outlier.shape = "-")+
  geom_col(median_lines,
           mapping = aes(x = -0.25, y = as.numeric(factor(gene_symbol)), 
               fill = median_alleleEffectSize),
           width = 0.5, alpha = 0.8) +
  labs(y="TF",x="scaled RPKM")+
  scale_fill_viridis(option="B",direction=-1)+
  theme_linedraw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

ggplot()+
  geom_rect(aes(xmin = mean(log2(encode_RPKML_gene_symbol$average_expression)) - sd(log2(encode_RPKML_gene_symbol$average_expression)),
                xmax = mean(log2(encode_RPKML_gene_symbol$average_expression)) + sd(log2(encode_RPKML_gene_symbol$average_expression)),
                ymin = -Inf, ymax = Inf),
            fill = "gray", alpha = 0.5) +
  geom_vline(xintercept = mean(log2(encode_RPKML_gene_symbol$average_expression)), linetype = "solid", color = "black", size = 1)+
  geom_point(mean_points, mapping=aes(x=mean_RPKM_log,y=reorder(gene_symbol,mean_RPKM_log), color=mean_alleleEffectSize),size = 3)+
  labs(y="TF",x="scaled RPKM")+
  scale_color_viridis(option="B",direction=-1)+
  theme_linedraw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

# TFBS analysis - motifbreakR final ---------------------------------------------

motifs_jaspar2022 <- query(MotifDb,
                           andStrings=c("hsapiens"),
                           orStrings=c("jaspar2022"))

HAR_snp <- file.path(sys_dir,"HAR/HAR_maf_vcf/merged_HAR_snp_singles.bed")
HAR_snp <- read.table(HAR_snp, header = FALSE)
remove_close_entries_df <- function(df, bp_threshold = 5) {
  df <- df[order(df$V2), ]
  start_diff <- diff(df$V2)
  close_entries <- which(start_diff <= bp_threshold)
  to_remove <- unique(c(close_entries, close_entries + 1))
  return(df[-to_remove, ])
}
HAR_snp_singles <- remove_close_entries_df(HAR_snp, 9)
write.table(HAR_snp_singles, file = file.path("merged_HAR_snp_singles.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
HAR_snp_2 <- file.path(sys_dir,"/HAR/merged_HAR_snp_singles.bed")
read.table(HAR_snp_2, header = FALSE)
HAR_snp_gr <- snps.from.file(file = HAR_snp_2,
                             search.genome = BSgenome.Hsapiens.UCSC.hg38,
                             format = "bed")

HAR_snps_order <- c(3, 1, 3, 2, 6)
L1 <- 1
L2 <- 3
L3 <- 4
L4 <- 7
L5 <- 9
HAR14_snp_gr <- HAR_snp_gr[L1:(L1+HAR_snps_order[1]-1)]
HAR16_snp_gr <- HAR_snp_gr[L2:(L2+HAR_snps_order[2]-1)]
HAR17_snp_gr <- HAR_snp_gr[L3:(L3+HAR_snps_order[3]-1)]
HAR18_snp_gr <- HAR_snp_gr[L4:(L4+HAR_snps_order[4]-1)]
HAR20_snp_gr <- HAR_snp_gr[L5:(L5+HAR_snps_order[5]-1)]

motifbreakR_reslist_single <- list()
motifbreakR_reslist_singles <- list()
HAR_res_single <- NA

for (i in c(14,16,17,18,20)) {
  HAR_snp_gr <- get(paste0("HAR", i, "_snp_gr"))
  HAR_res_single <- motifbreakR(snpList = HAR_snp_gr,
                         pwmList = motifs_jaspar2022,
                         threshold = 1e-4,
                         method = "ic",
                         bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                         BPPARAM = BiocParallel::bpparam())

  motifbreakR_reslist_singles[[paste0("motifbreakR_HAR", i, "_res")]] <- HAR_res_single
  HAR_res_sig_single <- HAR_res_single[HAR_res_single$alleleEffectSize > quantile(HAR_res_single$alleleEffectSize, 0.95) | HAR_res_single$alleleEffectSize < quantile(HAR_res_single$alleleEffectSize, 0.05),]
  motifbreakR_reslist_single[[paste0("motifbreakR_HAR", i, "_res_sig")]] <- HAR_res_sig_single
}

motifbreakR_allHARs <- do.call(rbind, lapply(motifbreakR_reslist_singles, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))

motifbreakR_allHARs <- do.call(rbind, lapply(motifbreakR_reslist_singles, function(x) {
  as.data.frame(matrix(x, nrow = 1), stringsAsFactors = FALSE)
}))

motifbreakR_allHARs <- as.data.frame(unlist(motifbreakR_reslist_singles))

TF_sig_HAR_list_single <- list()

for (i in c(14,16,17,18,20)) {
  HAR_res_sig_single <- motifbreakR_reslist_single[[paste0("motifbreakR_HAR", i, "_res_sig")]]
  TF_sig_HAR_list_single[[paste0("TF_sig_HAR", i)]] <- data.frame(
    geneSymbol = HAR_res_sig_single$geneSymbol,
    alleleEffectSize = HAR_res_sig_single$alleleEffectSize
  )
}

TF_allHAR_single <- bind_rows(TF_sig_HAR_list_single, .id = "HAR_ID")
TF_allHAR_brainRPKML_single <- merge(TF_allHAR_single, encode_RPKML_avgexp_human)
nrow(TF_allHAR_single)
nrow(TF_allHAR_brainRPKML_single)
TF_allHAR_brainRPKML_single$Effect <- ifelse(TF_allHAR_brainRPKML_single$alleleEffectSize < 0, "improving",
                                      ifelse(TF_allHAR_brainRPKML_single$alleleEffectSize > 0, "disrupting", "none"))
TF_allHAR_brainRPKML_single$avg_RPKM_log2 <- log2(TF_allHAR_brainRPKML_single$avg_RPKM+1)
encode_RPKML_avgexp_human$avg_RPKM_log2 <- log2(encode_RPKML_avgexp_human$avg_RPKM+1)
TF_allHAR_brainRPKML_highexp_single <- TF_allHAR_brainRPKML_single %>% filter(avg_RPKM_log2 > mean(TF_allHAR_brainRPKML_single$avg_RPKM_log2))
nrow(TF_allHAR_brainRPKML_single)
nrow(TF_allHAR_brainRPKML_highexp_single)
TF_allHAR_brainRPKML_highexp_single_unique <- TF_allHAR_brainRPKML_highexp_single[!duplicated(TF_allHAR_brainRPKML_highexp_single$geneSymbol) &
                                                                      !duplicated(TF_allHAR_brainRPKML_highexp_single$geneSymbol, fromLast = TRUE), ]
nrow(TF_allHAR_brainRPKML_highexp_single)
nrow(TF_allHAR_brainRPKML_highexp_single_unique)

ggplot()+
  geom_rect(aes(xmin = mean(encode_RPKML_avgexp_human$avg_RPKM_log2) - sd(encode_RPKML_avgexp_human$avg_RPKM_log2), 
                xmax = mean(encode_RPKML_avgexp_human$avg_RPKM_log2) + sd(encode_RPKML_avgexp_human$avg_RPKM_log2), 
                ymin = -Inf, ymax = Inf), 
            fill = "gray", alpha = 0.8) + 
  geom_vline(xintercept = mean(encode_RPKML_avgexp_human$avg_RPKM_log2), linetype = "solid", color = "black", size = 1)+
  geom_point(data=subset(TF_allHAR_brainRPKML_highexp_single_unique, Effect == "improving"), mapping=aes(x=avg_RPKM_log2,y=reorder(geneSymbol,avg_RPKM_log2),col=Effect))+
  geom_errorbar(data=subset(TF_allHAR_brainRPKML_highexp_single_unique, Effect == "improving"), mapping=aes(xmin=avg_RPKM_log2-sd_RPKM, xmax=avg_RPKM_log2+sd_RPKM,y=geneSymbol,col=Effect), width=0.3)+
  geom_point(data=subset(encode_RPKML_avgexp_human, geneSymbol == "GRIK2"), mapping=aes(x=avg_RPKM_log2,y=reorder(geneSymbol,avg_RPKM_log2),color=geneSymbol))+
  geom_errorbar(data=subset(encode_RPKML_avgexp_human, geneSymbol == "GRIK2"), mapping=aes(xmin=avg_RPKM_log2-sd_RPKM, xmax=avg_RPKM_log2+sd_RPKM,y=geneSymbol,color=geneSymbol), width=0.3)+
  # scale_color_manual(values = c("#FC766AFF","black","#5B84B1FF"))+
  scale_color_manual(values = c("black","#5B84B1FF"))+
  labs(y="TF",x="Expression levels")+
  theme_linedraw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

ggsave(file.path(sys_dir,"thesis_fig/fig2_ba.png"), plot=ggplot2::last_plot())

ggplot()+
  geom_rect(aes(xmin = mean(encode_RPKML_avgexp_human$avg_RPKM_log2) - sd(encode_RPKML_avgexp_human$avg_RPKM_log2), 
                xmax = mean(encode_RPKML_avgexp_human$avg_RPKM_log2) + sd(encode_RPKML_avgexp_human$avg_RPKM_log2), 
                ymin = -Inf, ymax = Inf), 
            fill = "gray", alpha = 0.8) + 
  geom_vline(xintercept = mean(encode_RPKML_avgexp_human$avg_RPKM_log2), linetype = "solid", color = "black", size = 1)+
  geom_point(data=subset(TF_allHAR_brainRPKML_highexp_single_unique, Effect == "disrupting"), mapping=aes(x=avg_RPKM,y=reorder(geneSymbol,avg_RPKM),col=Effect))+
  geom_errorbar(data=subset(TF_allHAR_brainRPKML_highexp_single_unique, Effect == "disrupting"), mapping=aes(xmin=avg_RPKM-sd_RPKM, xmax=avg_RPKM+sd_RPKM,y=geneSymbol,col=Effect), width=0.3)+
  geom_point(data=subset(encode_RPKML_avgexp_human, geneSymbol == "GRIK2"), mapping=aes(x=avg_RPKM,y=reorder(geneSymbol,avg_RPKM),color=geneSymbol))+
  geom_errorbar(data=subset(encode_RPKML_avgexp_human, geneSymbol == "GRIK2"), mapping=aes(xmin=avg_RPKM-sd_RPKM, xmax=avg_RPKM+sd_RPKM,y=geneSymbol,color=geneSymbol), width=0.3)+
  # scale_color_manual(values = c("#FC766AFF","black","#5B84B1FF"))+
  scale_color_manual(values = c("#FC766AFF","black"))+ 
  labs(y="TF",x="Expression levels")+
  theme_linedraw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

ggsave(file.path(sys_dir,"thesis_fig/fig2_bb.png"), plot=ggplot2::last_plot())

write_xlsx(TF_allHAR_brainRPKML_highexp_single_unique, file.path(sys_dir,"Figure_4B.xlsx"))

## combine?

ggplot()+
  geom_rect(aes(xmin=mean(encode_RPKML_avgexp_human$avg_RPKM_log2)-sd(encode_RPKML_avgexp_human$avg_RPKM_log2),xmax=mean(encode_RPKML_avgexp_human$avg_RPKM_log2)+sd(encode_RPKML_avgexp_human$avg_RPKM_log2),ymin=-Inf,ymax=Inf),fill="gray",alpha=0.8)+
  geom_vline(xintercept=mean(encode_RPKML_avgexp_human$avg_RPKM_log2),linetype="solid",color="black",size=1)+
  geom_point(data=subset(TF_allHAR_brainRPKML_highexp_single_unique,Effect=="improving"),mapping=aes(x=avg_RPKM,y=reorder(geneSymbol,avg_RPKM),col=Effect))+
  geom_errorbar(data=subset(TF_allHAR_brainRPKML_highexp_single_unique,Effect=="improving"),mapping=aes(xmin=avg_RPKM-sd_RPKM,xmax=avg_RPKM+sd_RPKM,y=geneSymbol,col=Effect),width=0.3)+
  geom_point(data=subset(TF_allHAR_brainRPKML_highexp_single_unique,Effect=="disrupting"),mapping=aes(x=avg_RPKM,y=reorder(geneSymbol,avg_RPKM),col=Effect))+
  geom_errorbar(data=subset(TF_allHAR_brainRPKML_highexp_single_unique,Effect=="disrupting"),mapping=aes(xmin=avg_RPKM-sd_RPKM,xmax=avg_RPKM+sd_RPKM,y=geneSymbol,col=Effect),width=0.3)+
  geom_point(data=subset(encode_RPKML_avgexp_human,geneSymbol=="GRIK2"),mapping=aes(x=avg_RPKM,y=reorder(geneSymbol,avg_RPKM),color=geneSymbol))+
  geom_errorbar(data=subset(encode_RPKML_avgexp_human,geneSymbol=="GRIK2"),mapping=aes(xmin=avg_RPKM-sd_RPKM,xmax=avg_RPKM+sd_RPKM,y=geneSymbol,color=geneSymbol),width=0.3)+
  scale_color_manual(values=c("#FC766AFF","black","#5B84B1FF"))+
  labs(y="TF",x="Expression levels")+
  theme_linedraw()+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
  theme(plot.title=element_text(size=16,face="bold"),axis.title.x=element_text(size=14,face="bold"),axis.title.y=element_text(size=14,face="bold"))

ggsave(file.path(sys_dir,"thesis_fig/fig2_bcombined.png"), plot=ggplot2::last_plot())

ggplot()+
  geom_rect(aes(xmin=mean(encode_RPKML_avgexp_human$avg_RPKM)-sd(encode_RPKML_avgexp_human$avg_RPKM),xmax=mean(encode_RPKML_avgexp_human$avg_RPKM)+sd(encode_RPKML_avgexp_human$avg_RPKM),ymin=-Inf,ymax=Inf),fill="gray",alpha=0.8)+
  geom_vline(xintercept=mean(encode_RPKML_avgexp_human$avg_RPKM),linetype="solid",color="black",size=1)+
  geom_point(data=TF_allHAR_brainRPKML_highexp_single_unique,mapping=aes(x=avg_RPKM,y=reorder(geneSymbol,avg_RPKM),col=Effect))+
  geom_errorbar(data=TF_allHAR_brainRPKML_highexp_single_unique,mapping=aes(xmin=avg_RPKM-sd_RPKM,xmax=avg_RPKM+sd_RPKM,y=geneSymbol,col=Effect),width=0.3)+
  geom_point(data=subset(encode_RPKML_avgexp_human,geneSymbol=="GRIK2"),mapping=aes(x=avg_RPKM,y=reorder(geneSymbol,avg_RPKM),color=geneSymbol))+
  geom_errorbar(data=subset(encode_RPKML_avgexp_human,geneSymbol=="GRIK2"),mapping=aes(xmin=avg_RPKM-sd_RPKM,xmax=avg_RPKM+sd_RPKM,y=geneSymbol,color=geneSymbol),width=0.3)+
  scale_color_manual(values=c("black","black","black"))+
  labs(y="TF",x="RPKM")+
  theme_linedraw()+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
  theme(plot.title=element_text(size=16,face="bold"),axis.title.x=element_text(size=14,face="bold"),axis.title.y=element_text(size=14,face="bold"))+
  facet_grid(rows=vars(Effect), scales="free_y")

TF_sig_motifbreakrate <- list()

for (i in c(14,16,17,18,20)) {
  HAR_motifbreakrate <- motifbreakR_reslist_single[[paste0("motifbreakR_HAR", i, "_res_sig")]]
  TF_sig_motifbreakrate[[paste0("TF_sig_HAR", i)]] <- data.frame(
    geneSymbol = HAR_motifbreakrate$geneSymbol,
    alleleEffectSize = HAR_motifbreakrate$alleleEffectSize,
    scoreRef = HAR_motifbreakrate$scoreRef,
    scoreAlt = HAR_motifbreakrate$scoreAlt,
    effect = HAR_motifbreakrate$effect
  )
}

HAR_motifbreakrate_df <- bind_rows(TF_sig_motifbreakrate, .id = "HAR_ID")
HAR_motifbreakrate_df_brain <- merge(HAR_motifbreakrate_df, encode_RPKML_avgexp_human)
HAR_motifbreakrate_df_brain$Effect <- ifelse(HAR_motifbreakrate_df_brain$alleleEffectSize < 0, "improving",
                                             ifelse(HAR_motifbreakrate_df_brain$alleleEffectSize > 0, "disrupting", "none"))
HAR_motifbreakrate_df_brain$avg_RPKM_log2 <- log2(HAR_motifbreakrate_df_brain$avg_RPKM+1)
HAR_motifbreakrate_df_brain <- HAR_motifbreakrate_df_brain %>% filter(avg_RPKM_log2 > mean(HAR_motifbreakrate_df_brain$avg_RPKM_log2))
HAR_motifbreakrate_df_brain <- HAR_motifbreakrate_df_brain[!duplicated(HAR_motifbreakrate_df_brain$geneSymbol) &
                                                                                    !duplicated(HAR_motifbreakrate_df_brain$geneSymbol, fromLast = TRUE), ]

HAR_motifbreakrate_df_brain$changerate <- ((abs(HAR_motifbreakrate_df_brain$scoreRef-HAR_motifbreakrate_df_brain$scoreAlt))/((HAR_motifbreakrate_df_brain$scoreAlt+HAR_motifbreakrate_df_brain$scoreRef)/2))*100
# HAR_motifbreakrate_df_brain$changerate_chimp <- (abs(HAR_motifbreakrate_df_brain$scoreAlt-HAR_motifbreakrate_df_brain$scoreRef))/HAR_motifbreakrate_df_brain$scoreRef
# HAR_motifbreakrate_df_brain$changerate <- (abs(HAR_motifbreakrate_df_brain$scoreRef-HAR_motifbreakrate_df_brain$scoreAlt))/HAR_motifbreakrate_df_brain$scoreAlt
HAR_motifbreakrate_df_brain$number <- seq(1:nrow(HAR_motifbreakrate_df_brain))

HAR_motifbreakrate_df_brain <- HAR_motifbreakrate_df_brain %>% as.data.frame()

HAR_motifbreakrate_df_brain_long <- HAR_motifbreakrate_df_brain %>%
  pivot_longer(cols = c(scoreRef, scoreAlt),
               names_to = "score_type",
               values_to = "score")

test <- unique(HAR_motifbreakrate_df_brain_long$HAR_ID)
motifbreakrate_perHAR <- HAR_motifbreakrate_df_brain_long %>% filter(HAR_ID == test[5])
shapiro.test(motifbreakrate_perHAR$score)
t.test(score ~ score_type, motifbreakrate_perHAR)
wilcox.test(score ~ score_type, motifbreakrate_perHAR)

custom_colors <- c("TF_sig_HAR17" = "#AB337CFF",
                   "TF_sig_HAR14" = "#FEB37BFF","TF_sig_HAR16" = "#D6456CFF",
                   "TF_sig_HAR18" = "#952C80FF","TF_sig_HAR20" = "#29115AFF")
custom_colors <- c("HAR_ctrl" = "#E64B35FF",
                   "TF_sig_HAR14" = "#4DBBD5FF","TF_sig_HAR15" = "#00A087FF",
                   "TF_sig_HAR16" = "#7E6148FF","TF_sig_HAR17" = "#91D1C2FF",
                   "TF_sig_HAR18" = "#3C5488FF","TF_sig_HAR19" = "#F39B7FFF",
                   "TF_sig_HAR20" = "#8491B4FF")

dummy_legend <- data.frame(
  HAR = c("Negative Control", "HAR14", "HAR15", "HAR16", "HAR17", "HAR18", "HAR19", "HAR20"),
  Color = c("#FCFDBF", "#FEB37B", "#F4685C", "#D6456C", "#AB337C", "#952C80", "#6B1D81", "#29115A")
)
dummy_legend$Frequency <- c(10, 20, 15, 25, 30, 12, 18, 22)
dummy_legend$HAR <- factor(dummy_legend$HAR, levels = dummy_legend$HAR)
ggplot(dummy_legend, aes(x = HAR, y = Frequency, fill = HAR)) +
  geom_col() +
  scale_fill_manual(values = dummy_legend$Color) +
  theme_minimal() +
  labs(title = "HAR Frequency Bar Plot",
       x = "HAR Category",
       y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(sys_dir,"thesis_fig/fig3_HARlegand.png"), plot=ggplot2::last_plot())

ggplot(HAR_motifbreakrate_df_brain, aes(x=HAR_ID, y = changerate, fill=HAR_ID)) +
  geom_boxplot()+
  labs(x = "HAR id",
       y = "Motif Change Rate") +
  theme_bw()+
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values=custom_colors)+
  # scale_fill_npg()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

ggsave(file.path(sys_dir,"thesis_fig/fig3_motifchangerate_1.png"), plot=ggplot2::last_plot())

write_xlsx(HAR_motifbreakrate_df_brain, file.path(sys_dir,"Figure_4C.xlsx"))

HAR_snps_order <- c(3, 1, 3, 2, 6)
divergence_perbp <- table(HAR_motifbreakrate_df_brain$HAR_ID) %>% as.data.frame()
divergence_perbp$SNPcount <- HAR_snps_order
divergence_perbp$FreqNorm <- divergence_perbp$Freq/divergence_perbp$SNPcount
divergence_perbp$FreqNorm_sum <- (divergence_perbp$Freq/divergence_perbp$SNPcount)/sum(divergence_perbp$FreqNorm)

ggplot(HAR_motifbreakrate_df_brain_long, aes(x = HAR_ID, y = score, fill = score_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "HAR ID",
       y = "Mean Score") +
  theme_classic()
  
HAR_motifbreakrate_perHAR <- HAR_motifbreakrate_df_brain %>%
  group_by(HAR_ID) %>%
  summarise(
    mean_changerate_chimp = mean(changerate_chimp, na.rm = TRUE),
    mean_changerate = mean(changerate, na.rm = TRUE),
    mean_scoreRef = mean(scoreRef, na.rm = TRUE),
    mean_scoreAlt = mean(scoreAlt, na.rm = TRUE)) %>%
  ungroup()
HAR_motifbreakrate_perHAR$number <- "All_HARs"

HAR_motifbreakrate_long <- HAR_motifbreakrate_perHAR %>%
  pivot_longer(cols = c(mean_scoreRef, mean_scoreAlt),
               names_to = "score_type",
               values_to = "mean_score")

ggplot(HAR_motifbreakrate_long, aes(x = HAR_ID, y = mean_score, fill = score_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Comparison of TF Binding Scores for Each HAR",
       x = "HAR ID",
       y = "Score",
       fill = "Score Type") + 
  theme_minimal()

ggplot(HAR_motifbreakrate_long, aes(x = HAR_ID, y = mean_score, fill = score_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Comparison of TF Binding Scores for Each HAR",
       x = "HAR ID",
       y = "Score",
       fill = "Score Type") + 
  theme_classic()

HAR_motifbreakrate_perHAR$HAR_ID

ggplot(HAR_motifbreakrate_perHAR, aes(x=HAR_ID,y=mean_changerate)) +
  geom_col(fill="steelblue")+
  labs(title = "Comparison of TF Binding Scores for Each HAR",
       x = "HAR ID",
       y = "Score",
       fill = "Score Type") + 
  theme_classic()

TF_motifbreakrate_overall <- list()

for (i in c(14,16,17,18,20)) {
  HAR_res_sig_single <- motifbreakR_reslist_singles[[paste0("motifbreakR_HAR", i, "_res")]]
  TF_motifbreakrate_overall[[paste0("TF_sig_HAR", i)]] <- data.frame(
    geneSymbol = HAR_res_sig_single$geneSymbol,
    alleleEffectSize = HAR_res_sig_single$alleleEffectSize
  )
}

HAR_motifbreakrate_overall <- bind_rows(TF_motifbreakrate_overall, .id = "HAR_ID")
HAR_motifbreakrate_overall_brain <- merge(HAR_motifbreakrate_overall, encode_RPKML_avgexp_human)
HAR_motifbreakrate_overall_brain$Effect <- ifelse(HAR_motifbreakrate_overall_brain$alleleEffectSize < 0, "improving",
                                             ifelse(HAR_motifbreakrate_overall_brain$alleleEffectSize > 0, "disrupting", "none"))
HAR_motifbreakrate_overall_brain$avg_RPKM_log2 <- log2(HAR_motifbreakrate_overall_brain$avg_RPKM+1)
HAR_motifbreakrate_overall_brain <- HAR_motifbreakrate_overall_brain %>% filter(avg_RPKM_log2 > mean(HAR_motifbreakrate_overall_brain$avg_RPKM_log2))
# HAR_motifbreakrate_overall_brain <- HAR_motifbreakrate_overall_brain[!duplicated(HAR_motifbreakrate_overall_brain$geneSymbol) &
#                                                              !duplicated(HAR_motifbreakrate_overall_brain$geneSymbol, fromLast = TRUE), ]

HAR_snps_order <- c(3, 1, 3, 2, 6)
divergence_perbp <- table(HAR_motifbreakrate_overall_brain$HAR_ID) %>% as.data.frame()
divergence_perbp$SNPcount <- HAR_snps_order
divergence_perbp$FreqNorm <- divergence_perbp$Freq/divergence_perbp$SNPcount
divergence_perbp$FreqNorm_sum <- (divergence_perbp$Freq/divergence_perbp$SNPcount)/sum(divergence_perbp$FreqNorm)

ggplot(divergence_perbp, aes(x=Freq)) +
  geom_histogram(binwidth=0.01, color="black", alpha=0.7) +
  # geom_vline(xintercept=grik2_count, color="#EE353E", linetype="dashed", size=0.5) +
  labs(x="HAR count", y="Frequency") +
  theme_classic()

TF_motifbreakrate_overall_df <- as.data.frame(unlist(TF_motifbreakrate_overall)) 

# GREAT HAR  --------------------------------------------------------------

great_HAR <- HAR_3100_grik2 %>% dplyr::select(c("chr","start","end","HAR_ID","size","grik2_distance"))
rownames(great_HAR) <- great_HAR$HAR_ID
great_HAR_gr <- makeGRangesFromDataFrame(df=great_HAR,
                                              keep.extra.columns=T,
                                              ignore.strand=T,
                                              seqnames.field=c("chr"),
                                              start.field="start",
                                              end.field="end",
                                              starts.in.df.are.0based=FALSE)

set.seed(321)
great_HAR_job = submitGreatJob(great_HAR_gr,species="hg38")
great_table = getEnrichmentTables(great_HAR_job)
str(great_table)
head(great_table)
plotRegionGeneAssociations(great_HAR_job)
availableOntologies(great_HAR_job)
great_results = great(great_HAR_gr, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene")
plotRegionGeneAssociations(great_results)
plotVolcano(great_results)
# shinyReport(great_HAR_job)

# no significant enrichment; very small number of region set (7). Do the analysis on all HARs

great_all_HAR_job = submitGreatJob(HAR_3100_gr,species="hg38")
great_all_table = getEnrichmentTables(great_all_HAR_job)
head(great_all_table)
great_all_results = great(HAR_3100_gr, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene")
plotVolcano(great_all_results)
# shinyReport(great_all_HAR_job)

# a lot of significant enrichment; impossible to relate to grik2 HARs. Take all HARs 5Mb around grik2

grik2_start = 100962701
grik2_end = 102081622

HAR_100_grik2 <- HAR_3100_grik2 %>% dplyr::select(c("chr","start","end","HAR_ID","size","grik2_distance")) %>%
  filter(start > grik2_start-5000000) %>% 
  filter(end < grik2_end+5000000)
nrow(HAR_100_grik2)
rownames(great_HAR) <- great_HAR$HAR_ID
great_HAR_gr <- makeGRangesFromDataFrame(df=great_HAR,
                                         keep.extra.columns=T,
                                         ignore.strand=T,
                                         seqnames.field=c("chr"),
                                         start.field="start",
                                         end.field="end",
                                         starts.in.df.are.0based=FALSE)

great_100_HAR_gr <- makeGRangesFromDataFrame(df=HAR_100_grik2,
                                         keep.extra.columns=T,
                                         ignore.strand=T,
                                         seqnames.field=c("chr"),
                                         start.field="start",
                                         end.field="end",
                                         starts.in.df.are.0based=FALSE)
great_subset_HAR_job = submitGreatJob(great_100_HAR_gr,species="hg38")
great_subset_table = getEnrichmentTables(great_subset_HAR_job)
head(great_subset_table)
great_subset_results = great(great_100_HAR_gr, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene")
plotVolcano(great_subset_results)
# shinyReport(great_subset_HAR_job)

# significant enrichment, but including other HARs might not be the best approach. Include 5Kb flanking regions instead

## define a function to add flanking_regions of a specific size
## tried 5Kb, 1Kb, and 500bp. Keep 500bp since its more stringent (shorter, no overlaps, a total of 1Kb per HAR is added)

add_flanking_regions <- function(df, flank_size = 500) {
  flanking_regions <- df %>%
    rowwise() %>%
    mutate(
      upstream_start = start - flank_size,
      upstream_end = start - 1,
      downstream_start = end + 1,
      downstream_end = end + flank_size,
      upstream_HAR_ID = paste0(HAR_ID, "_upstream"),
      downstream_HAR_ID = paste0(HAR_ID, "_downstream"),
      size = flank_size,
    ) %>%
    select(chr, upstream_start, upstream_end, upstream_HAR_ID, downstream_start, downstream_end, downstream_HAR_ID, size)
  upstream_df <- flanking_regions %>%
    select(chr, start = upstream_start, end = upstream_end, HAR_ID = upstream_HAR_ID, size) 
  downstream_df <- flanking_regions %>%
    select(chr, start = downstream_start, end = downstream_end, HAR_ID = downstream_HAR_ID, size)
  combined_df <- bind_rows(df, upstream_df, downstream_df)
    return(combined_df)
}

## remove grik2_distance and add flanking region

HAR_flank_grik2 <- great_HAR %>%
  dplyr::select(-grik2_distance) %>% 
  add_flanking_regions() %>%
  arrange(chr, start)

HAR_flank_grik2_gr <- makeGRangesFromDataFrame(df=HAR_flank_grik2,
                                             keep.extra.columns=T,
                                             ignore.strand=T,
                                             seqnames.field=c("chr"),
                                             start.field="start",
                                             end.field="end",
                                             starts.in.df.are.0based=FALSE)
great_HAR_flank_job = submitGreatJob(HAR_flank_grik2_gr,species="hg38")
great_HAR_flank_table = getEnrichmentTables(great_HAR_flank_job)
head(great_HAR_flank_table)
great_HAR_flank_results = great(HAR_flank_grik2_gr, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene")
plotVolcano(great_HAR_flank_results)
# shinyReport(great_HAR_flank_job)

great_HAR_flank_table_MF_df <- as.data.frame((great_HAR_flank_table$`GO Molecular Function`))
great_HAR_flank_table_BP_df <- as.data.frame((great_HAR_flank_table$`GO Biological Process`)) 
great_HAR_flank_table_CC_df <- as.data.frame((great_HAR_flank_table$`GO Cellular Component`)) 

Kp <- plotKaryotype(genome="hg38", plot.type = 2, chromosomes = c("chr6"))
kpPlotMarkers(Kp, data=grik2_gr, labels=grik2_gr$hgnc_symbol, text.orientation = "horizontal",
              r1=0.3, cex=0.8, adjust.label.position = FALSE)
kpPoints(Kp, data=HAR_flank_grik2_gr, y=y1, cex=2, pch=6 ,col="#EE353E")

zoom.region_grik2 <- toGRanges(data.frame("chr6", grik2_start-500000, grik2_end+500000))
HAR_flank_grik2_kp <- plotKaryotype(chromosomes="chr6", zoom=zoom.region_grik2)
kpDataBackground(HAR_flank_grik2_kp)
kpAddBaseNumbers(HAR_flank_grik2_kp)
kpAddCytobandLabels(HAR_flank_grik2_kp)
kpPlotMarkers(HAR_flank_grik2_kp, data=grik2_gr, labels=grik2_gr$hgnc_symbol, text.orientation = "horizontal",r1=0.3, cex=0.8, adjust.label.position = FALSE)
kpPoints(HAR_flank_grik2_kp, data=HAR_flank_grik2_gr, y=y1,  r0=0.1, r1=1, cex=2, pch=6, col="#EE353E")

# do the new genomic ranges of HAR flanking regions overlap with the nearest genes (GRIK2, HACE1, ASCC3)?

overlaps_GRIK2 <- subsetByOverlaps(HAR_flank_grik2_gr,grik2_gr)

HACE1_ENSG00000085382 <- data.frame(chr = "chr6",
                    start = 104728094,
                    end = 104859919,
                    hgnc_symbol = "HACE1")
HACE1_gr <- makeGRangesFromDataFrame(df=HACE1_ENSG00000085382,
                                     keep.extra.columns=T,
                                     ignore.strand=T,
                                     seqnames.field="chr",
                                     start.field="start",
                                     end.field="end",
                                     starts.in.df.are.0based=FALSE)
ASCC3_ENSG00000112249 <- data.frame(chr = "chr6",
                                    start = 100508194,
                                    end = 100881372,
                                    hgnc_symbol = "ASCC3") 
ASCC3_gr <- makeGRangesFromDataFrame(df=ASCC3_ENSG00000112249,
                                     keep.extra.columns=T,
                                     ignore.strand=T,
                                     seqnames.field="chr",
                                     start.field="start",
                                     end.field="end",
                                     starts.in.df.are.0based=FALSE)
overlaps_HACE1 <- subsetByOverlaps(HAR_flank_grik2_gr,HACE1_gr)
# output: no for 5Kb, 1Kb and 500b
overlaps_ASCC3 <- subsetByOverlaps(HAR_flank_grik2_gr,ASCC3_gr)
# output: no for 5Kb, 1Kb and 500b

# do the new genomic ranges of intronic HAR flanking regions overlap with grik2 exons?

test_grik2_exons <- read.table(file.path(sys_dir,"HAR/test_grik2_exons.txt")) %>%
  dplyr::select(c(V1,V2,V3,V4)) %>%
  filter(V1 == "chr6")
colnames(test_grik2_exons) <- c("chr","start","end","id") 
test_grik2_exons_gr <- makeGRangesFromDataFrame(df=test_grik2_exons,
                                     keep.extra.columns=T,
                                     ignore.strand=T,
                                     seqnames.field="chr",
                                     start.field="start",
                                     end.field="end",
                                     starts.in.df.are.0based=FALSE)
overlaps_GRIK2_exons <- subsetByOverlaps(HAR_flank_grik2_gr,test_grik2_exons_gr)
# output: yes for 5Kb, no for 1Kb and 500bp

# primers design for HARs ----------------------------------------------------------

gc_content <- function(dna_sequence) {
  # Convert the sequence to uppercase to ensure consistency
  dna_sequence <- toupper(dna_sequence)
  # Count the number of G and C
  g_count <- sum(strsplit(dna_sequence, NULL)[[1]] == "G")
  c_count <- sum(strsplit(dna_sequence, NULL)[[1]] == "C")
  # Total number of bases
  total_bases <- nchar(dna_sequence)
  # Calculate GC content
  gc_content_percentage <- ((g_count + c_count) / total_bases) * 100
  return(gc_content_percentage)
}
annealing_temp <- function(primer_sequence) {
  # Convert the sequence to uppercase to ensure consistency
  primer_sequence <- toupper(primer_sequence)
  # Count the number of A, T, G, and C
  a_count <- sum(strsplit(primer_sequence, NULL)[[1]] == "A")
  t_count <- sum(strsplit(primer_sequence, NULL)[[1]] == "T")
  g_count <- sum(strsplit(primer_sequence, NULL)[[1]] == "G")
  c_count <- sum(strsplit(primer_sequence, NULL)[[1]] == "C")
  # Calculate the annealing temperature
  tm <- 2 * (a_count + t_count) + 4 * (g_count + c_count)
  return(tm)
}
grep_fasta_pattern <- function(fasta_file, seq_pattern) {
  # Read the sequences from the FASTA file
  sequences <- read.fasta(file = fasta_file, seqtype = "AA")
  # Convert the sequences to a list of strings
  seq_strings <- sapply(sequences, function(seq) paste(seq, collapse = ""))
  # Grep for the pattern in the sequences
  matching_seqs <- seq_strings[grepl(seq_pattern, seq_strings)]
  return(matching_seqs)
}

## cloning by restriction enzymes - HARs 14,15,16,18,19

for_primer <- "GTTGCGTTTGAGACAGGCGACA"
rev_primer <- "AGTTCTGGACCAGCGAGCTGT"

nchar(for_primer)
nchar(rev_primer)

gc_content(for_primer)
gc_content(rev_primer)

annealing_temp(for_primer)
annealing_temp(rev_primer)

## check primers in your HAR sequences

HAR_path <- file.path(sys_dir,"HAR/HAR_flanking_human_seq.fa")
grik2_promoter_path <- file.path(sys_dir,"grik2_promoter_UTRs/grik2_promoter1000_HumRhesMou.fa")

fasta_file <- HAR_path
seq_pattern <- for_primer
seq_pattern <- rev_primer
grep_fasta_pattern(fasta_file, seq_pattern)

## add RE site to your forward primers
## add sitting sequence to your forward primers "ttat"

SpeI <- "ACTAGT"

sit_seq <- "ttat"

grik2HAR_SpeI_for <- paste0(sit_seq, SpeI, for_primer)
grik2HAR_rev <- reverseComplement(DNAStringSet(c(rev_primer))) %>% as.character()

## cloning by Gibson assembly - HARs 17,20

HARsv2_2517_dna <- "TTTAACCTAGGTAGTCTTACCTGCAATATAAAATGATACAGTTGCAAAAACATGAATTTAAAAACACTTGGAACAGACTTAGATAGACATTTTTGTATTCAGAGTGTGGTGCAGCTGACCATATTGAATGAGACCTCATTTACTCTCTTTAATTACAGTTTAACTTTAAGTTAAAAGTCATTTTTAAATGGTTGGTCTTATATATTCATATGCTTTTATAATACATTCTAAAAATTCAAATTCATTCTTGATTGAAGTCACTACATAAAGTTCAAGTGTTATCTCATGGACACTAATTTTCTCCTCCATAACTAGTTGTGTGAGTCTACATTTCAATTCTAAAGCCTTTTCCTTTTCACTGTAATGTATCTTGAAAGACTTAAATCACACATTCAAAAATAATTGCCATTTTAAATTGAACTCAGAGCAGAACAAACCAATTACTGTAGTTGCATGAGTGTAAAGCATTTTACTCTTCTATCAGCTGTAGAAGTGTGCATACTTCCTctatcacaattattgttacattactattgcatttgtattCACCTCCTGAAGCAGACTAAAAAGCTTC"
HARsv2_2520_dna <- "agtaaatccacaattcaaacctcagctatctaatttctaaaatgtgttcttaacttctGAATTACTTGGTAAGATATCCCAATGTTTAAGAACGAAAAAAAAATAGTACTTTTGAAGAAATCTAGATAATATCCAACCTCGTTAAATTTACTCCCTATAGAGATATTCTTCAGATGCTGAAAGAAATGTTTCCTTTGCCTTTATTTTCTTACCTATATCAATTTTGGCATTCAGGGCCATGTTAGGGCTTAGGATATAAATTGTATGTCTGTAAAAATTTCCATTCTGTTAATGTCTTTTATGTGGCAAATTACTCTTGTTCCTTATTATGGTCTATGAGTGTAATGATTGGGTAATGATATCCAATATTCTGAATCATTTTGTGTATATACACCATATGCATGCAAATTCAATTTGCCTCATACTTTGCTCTTTAATTTGCATTTTCATAGAACTAGTTGAGTATTATGAGGTAAAATTATATTGACTGATTTCAATATCTTAAAAAGTTAAATCGTCATTATATGGGAAATGTTTTAAGCATTAATTATAAATGAATCATAAAAAAATTCTCAAGGCACACATCTTGAGTTATAATGTGATAATAGCTACATGTACTAAAGATTTATATTAGAGTTGCAGTACATAACCCATGCTTTGATTAAAAAATTAAACAGTGATATCAAAAGGAAATAATTCATTTTTTAACTGGAAATACCAATGTAAATATATTAATGTAATTCTAAGCCAGTACA"

for_primer_17   <- "TTTAACCTAGGTAGTCTTACCTGCAATATA"
rev_primer_17   <- "CACCTCCTGAAGCAGACTAAAAAGCTTC"

nchar(for_primer_17)
nchar(rev_primer_17)
gc_content(for_primer_17)
gc_content(rev_primer_17)
annealing_temp(for_primer_17)
annealing_temp(rev_primer_17)

IDTplasmid_path <- file.path(sys_dir,"/HAR/pucidt_amp_goldengate_IDT.fa")
fasta_file <- IDTplasmid_path
seq_pattern <- for_primer_17
seq_pattern <- rev_primer_17
grep_fasta_pattern(fasta_file, seq_pattern)

for_primer_20   <- "agtaaatccacaattcaaacctcagctatc"
rev_primer_20   <- "GTAAATATATTAATGTAATTCTAAGCCAGTACA"

nchar(for_primer_20)
nchar(rev_primer_20)
gc_content(for_primer_20)
gc_content(rev_primer_20)
annealing_temp(for_primer_20)
annealing_temp(rev_primer_20)

IDTplasmid_path <- file.path(sys_dir,"/HAR/pucidt_amp_goldengate_IDT.fa")
fasta_file <- IDTplasmid_path
seq_pattern <- for_primer_20
seq_pattern <- rev_primer_20
grep_fasta_pattern(fasta_file, seq_pattern)

## add vector specific regions to 17,20 forward and reverse primers

vectorcc185_5prime <- "GCCACCTGGGTCGACATTGATTATTGACTAGT"
vectorcc185_3prime <- "ATGGTCGAGGTGAGCCCCACGTTCTGCTTC"

GA_grik2HAR17_cc185_for <- paste0(vectorcc185_5prime, for_primer_17)
GA_grik2HAR17_cc185_rev <- paste0(vectorcc185_3prime, reverseComplement(DNAStringSet(c(rev_primer_17))))
nchar(GA_grik2HAR17_cc185_for)
nchar(GA_grik2HAR17_cc185_rev)
GA_grik2HAR20_cc185_for <- paste0(vectorcc185_5prime, for_primer_20)
GA_grik2HAR20_cc185_rev <- paste0(vectorcc185_3prime, reverseComplement(DNAStringSet(c(rev_primer_20))))
GA_grik2HAR20_cc185_rev <- reverseComplement(DNAStringSet(c(vectorcc185_3prime,rev_primer_20))) %>% as.character()
nchar(GA_grik2HAR20_cc185_for)
nchar(GA_grik2HAR20_cc185_rev)

## inverse PCR for cc185 vector

for_primer_cc185 <- "ATGGTCGAGGTGAGCCCCACGTT"
rev_primer_cc185 <- "CTGGGTCGACATTGATTATTGACTAGT"

nchar(for_primer_cc185)
nchar(rev_primer_cc185)
gc_content(for_primer_cc185)
gc_content(rev_primer_cc185)
annealing_temp(for_primer_cc185)
annealing_temp(rev_primer_cc185)

cc185_inverPCR_noCMV_for <- for_primer_cc185
cc185_inverPCR_noCMV_rev <- reverseComplement(DNAStringSet(c(rev_primer_cc185))) %>% as.character()

# primers design - mutagenesis for grik2 ----------------------------------------------------------

gc_content <- function(dna_sequence) {
  # Convert the sequence to uppercase to ensure consistency
  dna_sequence <- toupper(dna_sequence)
  # Count the number of G and C
  g_count <- sum(strsplit(dna_sequence, NULL)[[1]] == "G")
  c_count <- sum(strsplit(dna_sequence, NULL)[[1]] == "C")
  # Total number of bases
  total_bases <- nchar(dna_sequence)
  # Calculate GC content
  gc_content_percentage <- ((g_count + c_count) / total_bases) * 100
  return(gc_content_percentage)
}
annealing_temp <- function(primer_sequence) {
  # Convert the sequence to uppercase to ensure consistency
  primer_sequence <- toupper(primer_sequence)
  # Count the number of A, T, G, and C
  a_count <- sum(strsplit(primer_sequence, NULL)[[1]] == "A")
  t_count <- sum(strsplit(primer_sequence, NULL)[[1]] == "T")
  g_count <- sum(strsplit(primer_sequence, NULL)[[1]] == "G")
  c_count <- sum(strsplit(primer_sequence, NULL)[[1]] == "C")
  # Calculate the annealing temperature
  tm <- 2 * (a_count + t_count) + 4 * (g_count + c_count)
  return(tm)
}
grep_fasta_pattern <- function(fasta_file, seq_pattern) {
  # Read the sequences from the FASTA file
  sequences <- read.fasta(file = fasta_file, seqtype = "AA")
  # Convert the sequences to a list of strings
  seq_strings <- sapply(sequences, function(seq) paste(seq, collapse = ""))
  # Grep for the pattern in the sequences
  matching_seqs <- seq_strings[grepl(seq_pattern, seq_strings)]
  return(matching_seqs)
}

## protocol A - Thermo Fisher

core_primer       <- "ctaaatagtttctggtttg"
core_primer_mut   <- "ctaaatTCGttTtggtttg"

overhang_3prime <- "gagttggagc"
overhang_5prime <- "caattttaccttg"

for_primer_mut <- paste0(core_primer_mut, overhang_3prime) 
rev_primer_mut <- paste0(overhang_5prime, core_primer_mut) 

nchar(core_primer_mut)
nchar(for_primer_mut)
nchar(rev_primer_mut)

gc_content(core_primer_mut)
gc_content(for_primer_mut)
gc_content(rev_primer_mut)

annealing_temp(core_primer_mut)
annealing_temp(for_primer_mut)
annealing_temp(rev_primer_mut)

grik2mut_overhang_for <- for_primer_mut
grik2mut_overhang_rev <- reverseComplement(DNAStringSet(c(rev_primer_mut))) %>% as.character()

nchar(grik2mut_overhang_for)
nchar(grik2mut_overhang_rev)

## protocol B - Thermo Fisher

for_primer_wildtype <-  "ccttgctaaatagtttctggtttggagttgg"
for_primer_mut <-       "ccttgctaaatTCGttTtggtttggagttgg"
rev_primer_wildtype <-  "ccctgactcagacgtggtggaaaacaatttta"

nchar(for_primer_wildtype)
nchar(for_primer_mut)
nchar(rev_primer_wildtype)

gc_content(for_primer_wildtype)
gc_content(for_primer_mut)
gc_content(rev_primer_wildtype)

annealing_temp(for_primer_wildtype)
annealing_temp(for_primer_mut)
annealing_temp(rev_primer_wildtype)

grik2mut_for_protocolB <- for_primer_mut
grik2_rev_protocolB <- reverseComplement(DNAStringSet(c(rev_primer_wildtype))) %>% as.character()

nchar(grik2mut_for_protocolB)
nchar(grik2_rev_protocolB)

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


