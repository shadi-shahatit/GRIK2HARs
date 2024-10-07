#### The Roles of GRIK2 and HARs in the Evolution and Development of Synapses
#### Shadi Shahatit, Master's Thesis, IBENS, Paris 2024
# libraries ---------------------------------------------------------------

library("GenomicDistributions")
library("GenomicDistributionsData")
library("BSgenome.Hsapiens.UCSC.hg38.masked")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
library("GenomicDistributions")
library("topGO")
library("biomaRt")
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

## note: replace your system's directory in sys_dir
sys_dir <- "/home/shadi/Desktop/S4_Project"

# Promoter and 5' 3' UTRs of grik2 ----------------------------------------

## grik2 species coordinates

## human - GRCH38

human_ncbi_NC_000006.12_2898 <- "chr6:101393708-102070083"
human_ensembl_ENSG00000164418 <- "chr6:100962701-102081622"
human_ensembl_canonical_ENST00000369134.9 <- "chr6:101393708-102070083"

## chimp - Pantro3.0

chimp_ncbi_NC_072403.2 <- "chr5:112823618-114006919"
chimp_ensembl_ENSPTRG00000018449 <- "chr6:104757582-105448610"
chimp_ensembl_canonical_ENSPTRT00000034097.3 <- "chr6:104757582-105448610"

## mouse - GRCm39

mouse_ncbi_NC_000076.7 <- "chr10:48969776-49666523"
mouse_ensembl_ENSMUSG00000056073 <- "chr10:48970929-49664862"
mouse_ensembl_canonical_ENSMUST00000218823.2 <- "chr10:48970929-49664862"

## retrieve sequences - biomaRt

human = useMart("ensembl",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
chimp = useMart("ensembl",dataset="ptroglodytes_gene_ensembl", host="www.ensembl.org")
mouse = useMart("ensembl",dataset="mmusculus_gene_ensembl", host="www.ensembl.org")

getLDS(attributes = c("hgnc_symbol","chromosome_name", "start_position"),
       filters = "hgnc_symbol", values = "grik2", mart = human,
       attributesL = c("refseq_mrna","chromosome_name","start_position"), martL = mouse)

seqH_promoter <- getSequence(id = "grik2", 
                        type = "hgnc_symbol", 
                        seqType="gene_flank",
                        upstream=1000,
                        mart = human)

seqH_5UTR <- getSequence(id = "GRIK2", 
                        type = "hgnc_symbol", 
                        seqType = "5utr", 
                        mart = chimp)

## server error revisit again
## https://www.bioconductor.org/packages/3.7/bioc/vignettes/biomaRt/inst/doc/biomaRt.html#introduction 

seq = getSequence(id = "GRIK2", 
                  type = "hgnc_symbol", 
                  seqType = "coding_gene_flank",
                  downstream = 1000,
                  mart = human)
show(seq)

# Transcription factor binding site analysis ------------------------------

## fimo

## Human
Fimo_human_promoter <- importFimo(file.path(sys_dir,"grik2_promoter_UTRs/TFBS/fimo_out_1/fimo.tsv"))
Fimo_human_promoter_df <- Fimo_human_promoter %>% as.data.frame()
Fimo_human_promoter_df_top1 <- Fimo_human_promoter_df %>% filter(qvalue <= quantile(qvalue, 0.01))
nrow(Fimo_human_promoter_df)
nrow(Fimo_human_promoter_df_top1)

## Macaque
Fimo_macaque_promoter <- importFimo(file.path(sys_dir,"grik2_promoter_UTRs/TFBS/fimo_mac/fimo.tsv"))
Fimo_macaque_promoter_df <- Fimo_macaque_promoter %>% as.data.frame()
Fimo_macaque_promoter_df_top1 <- Fimo_macaque_promoter_df %>% filter(qvalue <= quantile(qvalue, 0.01))
nrow(Fimo_macaque_promoter_df)
nrow(Fimo_macaque_promoter_df_top1)

## Mouse
Fimo_mouse_promoter <- importFimo(file.path(sys_dir,"grik2_promoter_UTRs/TFBS/fimo_mouse/fimo.tsv"))
Fimo_mouse_promoter_df <- Fimo_mouse_promoter %>% as.data.frame()
Fimo_mouse_promoter_df_top1 <- Fimo_mouse_promoter_df %>% filter(qvalue <= quantile(qvalue, 0.01))
nrow(Fimo_mouse_promoter_df)
nrow(Fimo_mouse_promoter_df_top1)

unique_to_human_promoter <- anti_join(Fimo_human_promoter_df, Fimo_macaque_promoter_df, by = c("motif_alt_id"))
unique_to_macaque_promoter <- anti_join(Fimo_macaque_promoter_df, Fimo_human_promoter_df, by = c("motif_alt_id"))
unique_to_mouse_promoter <- anti_join(Fimo_mouse_promoter_df, Fimo_human_promoter_df, by = c("motif_alt_id"))

nrow(unique_to_human_promoter)
length(unique(unique_to_human_promoter$motif_alt_id))
nrow(unique_to_macaque_promoter)
length(unique(unique_to_macaque_promoter$motif_alt_id))
nrow(unique_to_mouse_promoter)
length(unique(unique_to_mouse_promoter$motif_alt_id))

cat(unique(unique_to_human_promoter$motif_alt_id), sep = "\n")
cat(unique(unique_to_macaque_promoter$motif_alt_id), sep = "\n")
cat(unique(unique_to_mouse_promoter$motif_alt_id), sep = "\n")

## cluster buster

## Human

cbust_humanpro <- read.table(file.path(sys_dir,"grik2_promoter_UTRs/TFBS/cbustout_humanpro_1_gene_symbols.txt"), header=T)
cbust_humanpro <- cbust_humanpro %>% dplyr::select(-c(Start,End))
cbust_humanpro_t <- as.data.frame(t(cbust_humanpro[-1]))
colnames(cbust_humanpro_t) <- c("A_hum","B","C","D","E","G")
cbust_humanpro_t$gene_sym <- rownames(cbust_humanpro_t)
rownames(cbust_humanpro_t) <- seq(1:nrow(cbust_humanpro_t))
cbust_humanpro_t_sig <- cbust_humanpro_t %>% filter(A_hum > quantile(A_hum,0.95))
cat(cbust_humanpro_t_sig$gene_sym, sep = "\n")
cat(cbust_humanpro_t$gene_sym, sep = "\n")

## Macaque

cbust_macaquepro <- read.table(file.path(sys_dir,"grik2_promoter_UTRs/TFBS/cbustout_macaquepro_1_gene_symbols.txt"), header=T)
cbust_macaquepro <- cbust_macaquepro %>% dplyr::select(-c(Start,End))
cbust_macaquepro_t <- as.data.frame(t(cbust_macaquepro[-1]))
colnames(cbust_macaquepro_t) <- c("A_mac","B","C","D","E","G","H")
cbust_macaquepro_t$gene_sym <- rownames(cbust_macaquepro_t)
rownames(cbust_macaquepro_t) <- seq(1:nrow(cbust_macaquepro_t))
cbust_macaquepro_t_sig <- cbust_macaquepro_t %>% filter(A_mac > quantile(A_mac,0.99))
cat(cbust_macaquepro_t_sig$gene_sym, sep = "\n")
cat(cbust_macaquepro_t$gene_sym, sep = "\n")

# cbust_macaquepro_t_highcluster <- cbust_macaquepro_t %>% dplyr::select(A_mac,gene_sym)
# cbust_macaquepro_t_human_sig <- merge(cbust_macaquepro_t_highcluster,cbust_humanpro_t_sig) %>% 
#   dplyr::select(gene_sym,A_mac,A_hum)
# options(digits = 15)
# cbust_hummac <- cbust_macaquepro_t_human_sig %>%
#   rowwise() %>%
#   mutate(
#     mean = mean(c(A_mac, A_hum)),
#     sd = sd(c(A_mac, A_hum)))
# 
# ggplot()+
#   geom_col(cbust_hummac,mapping=aes(x=gene_sym,y=(A_hum)))
# ggplot()+
#   geom_col(cbust_hummac,mapping=aes(x=(gene_sym),y=(A_mac)))
# 
# heatmap_cbust_hummac <- cbust_macaquepro_t_human_sig %>%
#   melt(id.vars = "gene_sym", variable.name = "species", value.name = "binding_score")
# 
# ggplot(heatmap_cbust_hummac, aes(x=gene_sym, y=species, fill=scale(binding_score)))+
#   geom_tile()+
#   scale_fill_viridis_c(name="Binding Score")+
#   theme_minimal()+
#   labs(x="Genes", y="Species") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))

## GO enrichment

expressed_TF_human <- cbust_humanpro_t$gene_sym
TF_highscore_human <- cbust_humanpro_t_sig$gene_sym
gene_list_human <- factor(as.integer(cbust_humanpro_t$gene_sym %in% TF_highscore_human))
names(gene_list_human) <- expressed_TF_human
geneID2GO <- mapIds(org.Hs.eg.db, keys=expressed_TF_human, column="GO", keytype="SYMBOL", multiVals="first")
geneID2GO <- split(names(geneID2GO), geneID2GO)
sampleGOdata_human <- new("topGOdata", ontology = "BP", allGenes = gene_list_human, annot = annFUN.GO2genes, GO2genes = geneID2GO, nodeSize = 20)
result_topGO_human <- runTest(sampleGOdata_human, algorithm = "classic", statistic = "fisher")
genetable_human <- GenTable(sampleGOdata_human, classicFisher = result_topGO_human, topNodes = 10) %>% as.data.frame()
write.table(genetable_human, file = file.path("sampleGOdata_human.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

top_terms_human <- GenTable(sampleGOdata_human, classicFisher = result_topGO_human, topNodes = 10)
top_terms_human$classicFisher <- as.numeric(as.character(top_terms_human$classicFisher))
top_terms_human$log_pval <- -log10(top_terms_human$classicFisher)
ggplot(top_terms_human, aes(x = reorder(Term, -log_pval), y = log_pval)) +
  geom_point(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top GO Terms", x = "GO Term", y = "-log10(p-value)") +
  theme_minimal()

expressed_TF_macaque <- cbust_macaquepro_t$gene_sym
TF_highscore_macaque <- cbust_macaquepro_t_sig
gene_list_macaque <- factor(as.integer(cbust_macaquepro_t$gene_sym %in% TF_highscore_macaque))
names(gene_list_macaque) <- expressed_TF_macaque
geneID2GO <- mapIds(org.Hs.eg.db, keys=expressed_TF_macaque, column="GO", keytype="SYMBOL", multiVals="first")
geneID2GO <- split(names(geneID2GO), geneID2GO)
sampleGOdata_macaque <- new("topGOdata", ontology = "BP", allGenes = gene_list_macaque, annot = annFUN.GO2genes, GO2genes = geneID2GO, nodeSize = 20)
result_topGO_macaque <- runTest(sampleGOdata_macaque, algorithm = "classic", statistic = "fisher")
GenTable(sampleGOdata_macaque, classicFisher = result_topGO_macaque, topNodes = 10)

top_terms_macaque <- GenTable(sampleGOdata_macaque, classicFisher = result_topGO_macaque, topNodes = 10)
top_terms_macaque$classicFisher <- as.numeric(as.character(top_terms_macaque$classicFisher))
top_terms_macaque$log_pval <- -log10(top_terms_macaque$classicFisher)
ggplot(top_terms_macaque, aes(x = reorder(Term, -log_pval), y = log_pval)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top GO Terms", x = "GO Term", y = "-log10(p-value)") +
  theme_minimal()

unique_to_human <- setdiff(TF_highscore_human, TF_highscore_macaque)
unique_to_macaque <- setdiff(TF_highscore_macaque, TF_highscore_human)

# Promoter genomic features  ----------------------------------------------

promoter_PhyloP <- read.table(file.path(sys_dir,"grik2_promoter_UTRs/grik2promoter_cons30primates_hg38.txt"))
colnames(promoter_PhyloP) <- c("bp","PhyloP")
promoter_PhyloP$chr <- "chr6"
promoter_PhyloP$start <- promoter_PhyloP$bp-1
promoter_PhyloP$end <- promoter_PhyloP$bp

summary(promoter_PhyloP$PhyloP)

ggplot()+
  geom_line(promoter_PhyloP,mapping=aes(x=bp,y=PhyloP),color="steelblue")+
  theme_bw()

bin          <- 10
df2          <- promoter_PhyloP %>% 
  mutate(window = start %/% bin) %>% 
  group_by(window,chr) %>%
  summarise(PhyloP_mean = mean(PhyloP)) %>%
  mutate(start = (window*bin)+1, stop=(window+1)*bin) 
df_out       <- c()
for(i in unique(df2$chr)) {
  df2_tmp      <- df2 %>% filter(chr == i)
  missing_win  <- setdiff(0:max(df2_tmp$window),unique(df2_tmp$window))
  df_tmp       <- data.frame(window = missing_win) %>% 
    mutate(chr    = i,
           PhyloP_mean = NA,
           start  = (window*bin)+1, 
           stop   = (window+1)*bin) %>% 
    bind_rows(df2_tmp) %>% 
    arrange(chr,window)
  df_out       <- bind_rows(df_out,df_tmp)
}
promoter_PhyloP_binned <- df_out %>% arrange(chr, window) %>% drop_na()

promoter_PhyloP_gr <- makeGRangesFromDataFrame(df=promoter_PhyloP_binned,
                                     keep.extra.columns=T,
                                     ignore.strand=T,
                                     seqnames.field="chr",
                                     start.field="start",
                                     end.field="stop",
                                     starts.in.df.are.0based=FALSE)

start <- 239
end <- 344
rangeA <- 101392706+start
rangeB <- 101392706+end

promoter_PhyloP_binned$label <- lapply(promoter_PhyloP_binned$start, function(x) if(x>=rangeA && x<=rangeB){
  promoter_PhyloP_binned$label="InRange"
}else {
  promoter_PhyloP_binned$label="outRange"
})
promoter_PhyloP_binned$label <- as.character(unlist(promoter_PhyloP_binned$label))
mygray <- rgb(0, 0, 0, alpha = 0.7)

ggplot() +
  geom_col(promoter_PhyloP_binned,mapping=aes(x=window,y=PhyloP_mean,fill=label)) +
  geom_hline(yintercept=0, color="black") +
  scale_fill_manual(values=c("#EE353E",mygray))+
  labs(x="Genomic position (bp)", y="Conservation rate (PhyloP)") +
  theme_classic()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

ggsave(file.path(sys_dir,"thesis_fig/fig2_c.png"), plot=ggplot2::last_plot())

write_xlsx(promoter_PhyloP_binned, file.path(sys_dir,"Figure_4G.xlsx"))

ggplot()+
  # geom_col(promoter_PhyloP_binned,mapping=aes(x=window,y=PhyloP_mean),color="steelblue",fill="steelblue")+
  geom_line(promoter_PhyloP_binned,mapping=aes(x=window,y=PhyloP_mean),color="steelblue")+
  geom_hline(yintercept=0, linetype="dashed", color="black")+
  labs(x="Genomic window",y="Evolutionary conservation rate - PhyloP")+
  theme_classic()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

ggplot()+
  geom_col(promoter_PhyloP_binned,mapping=aes(x=window,y=PhyloP_mean,color=label))+
  geom_hline(yintercept=0, linetype="dashed", color="black")+
  labs(x="Genomic window",y="Evolutionary conservation rate - PhyloP")+
  theme_classic()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

cbust_humanpro <- read.table(file.path(sys_dir,"grik2_promoter_UTRs/TFBS/cbustout_humanpro_1_gene_symbols.txt"), header=T)
cbust_humanpro_start_end <- cbust_humanpro %>% dplyr::select(c(Start,End)) %>%
  filter(cbust_humanpro$Score > 1)
cbust_humanpro_start_end$Start <- cbust_humanpro_start_end$Start+101392707
cbust_humanpro_start_end$End <- cbust_humanpro_start_end$End+101392707

check_within_ranges <- function(start, end, cbust_humanpro_start_end) {
  any((start >= cbust_humanpro_start_end$Start & start <= cbust_humanpro_start_end$End) |
        (end >= cbust_humanpro_start_end$Start & end <= cbust_humanpro_start_end$End) |
        (start <= cbust_humanpro_start_end$Start & end >= cbust_humanpro_start_end$End))
}
promoter_PhyloP_motif_labeled <- promoter_PhyloP %>%
  rowwise() %>%
  mutate(label = ifelse(check_within_ranges(start, end, cbust_humanpro_start_end), "InRange", "OutRange")) %>%
  ungroup()

ggplot() +
  geom_col(promoter_PhyloP_motif_labeled,mapping=aes(x=bp,y=PhyloP,color=label)) +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  labs(x="Genomic Position (bp)", y="Evolutionary conservation rate - PhyloP") +
  theme_classic()

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




