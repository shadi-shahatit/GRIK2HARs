#### The Roles of GRIK2 and HARs in the Evolution and Development of Synapses
#### Shadi Shahatit, Master's Thesis, IBENS, Paris 2024
# libraries ---------------------------------------------------------------

library("writexl")
library("effsize")
library("data.table")
library("magrittr")
library("WGCNA")
library("scales")
library("ggplot2")
library("forcats")
library("wesanderson")
library("stringr")
library("tidyr")
library("dplyr")
library("splitstackshape")
library("ggridges")
library("ggthemes")
library("readxl")
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
library("readr")
library("readxl")
library("maftools")
library("scMethrix")
library("BRGenomics")

## note: replace your system's directory in sys_dir
sys_dir <- "/home/shadi/Desktop/S4_Project"

# Developmental rhesus and human data for grik2 -----------------------------

## download data from http://evolution.psychencode.org/#

grik2_RPKM <- read.table(file.path(sys_dir,"PsychENCODE/grik2_nhp_development_RPKM_rmTechRep.txt"),header=T)
grik2_RPKML <- pivot_longer(grik2_RPKM, everything(), names_to = "Variable", values_to = "Value")
colnames(grik2_RPKML) <- c("Sample","RPKM")
grik2_RPKML[,3:4] <- stringr::str_split_fixed(grik2_RPKML$Sample, "\\.", 2) 
colnames(grik2_RPKML) <- c("Sample","RPKM","ID","Tissue")
# grik2_RPKML <- grik2_RPKML %>% dplyr::select(c("ID","Tissue","RPKM"))
sample_metadata <- read_excel(file.path(sys_dir,"PsychENCODE/aat8077_tables-s1-s3.xlsx"),sheet=1) %>%
  dplyr::select("Sample","Species","Age","Predicted age..PC.Days.","Period")
grik2_RPKML_age <- merge(grik2_RPKML,sample_metadata)
colnames(grik2_RPKML_age) <- c( "Sample","RPKM","ID","Tissue","Species","Age","Predicted_Age","Period")
grik2_RPKML_age$RPKM_log <- log2((grik2_RPKML_age$RPKM+1))
grik2_RPKML_age$Predicted_Age_log <- log2((grik2_RPKML_age$Predicted_Age))

# grik2_RPKML_age$Age <- gsub("1Y","1 Y",grik2_RPKML_age$Age)
# grik2_RPKML_age$Age <- gsub("4 Y","4Y",grik2_RPKML_age$Age)
# grik2_RPKML_age$Age <- gsub("11 Y","11Y",grik2_RPKML_age$Age)
# grik2_RPKML_age$Age <- gsub("11 Y","11Y",grik2_RPKML_age$Age)
# 
# Age_order <- fct_relevel(grik2_RPKML_age$Age,
#                          "E60","E80","E81","E82","E110","E111",
#                          "P0","P2",
#                          "8 PCW","9 PCW","12 PCW","13 PCW","16 PCW","17 PCW","19 PCW","21 PCW","22 PCW", "37 PCW",
#                          "4 M","7M","10 M",
#                          "1 Y","2Y","3 Y","3.5Y","4Y","5Y","7Y","8 Y","11Y","13 Y","15 Y","19 Y","21 Y","23 Y","30 Y","36 Y","37 Y","40 Y")  

ggplot()+
  geom_point(grik2_RPKML_age,mapping=aes(x=Predicted_Age_log,y=RPKM_log,fill=Species,col=Species),size=0.1,shape=23)+
  geom_smooth(grik2_RPKML_age,mapping=aes(x=Predicted_Age_log,y=RPKM_log,col=Species),method = "lm")+
  labs(x="Predicted_Age","log2(RPKM+1)",labs="Species")+
  theme_classic()+
  # facet_wrap(~Period)
  facet_wrap(~Tissue)

# grik2_RPKML_age_combined <- grik2_RPKML_age
# grik2_RPKML_age_combined$Tissue <- gsub("VFC_001", "VFC", grik2_RPKML_age_combined$Tissue)
# grik2_RPKML_age_combined$Tissue <- gsub("STC_GAIIx", "STC", grik2_RPKML_age_combined$Tissue)
# grik2_RPKML_age_combined$Tissue <- gsub("A1C_GAIIx", "M1C", grik2_RPKML_age_combined$Tissue)
# grik2_RPKML_age_combined$Tissue <- gsub("M1C_GAIIx", "M1C", grik2_RPKML_age_combined$Tissue)

## filter out weird tissue

grik2_RPKML_age_filtered <- grik2_RPKML_age %>% 
  filter(Tissue != "CGE", Tissue != "DTH", Tissue != "LGE", Tissue != "MGE", Tissue != "OC", Tissue != "PC",
         Tissue != "TC", Tissue != "DTH", Tissue != "URL", Tissue != "MSC",
         Tissue != "VFC_001", Tissue != "STC_GAIIx", Tissue != "A1C_GAIIx", Tissue != "M1C_GAIIx")

ggplot()+
  geom_point(grik2_RPKML_age_filtered,mapping=aes(x=Predicted_Age_log,y=RPKM_log,fill=Species,col=Species),size=0.1,shape=23)+
  geom_smooth(grik2_RPKML_age_filtered,mapping=aes(x=Predicted_Age_log,y=RPKM_log,col=Species),method = "loess")+
  labs(x="Predicted_Age","log2(RPKM+1)",labs="Species")+
  theme_classic()+
  facet_wrap(~Tissue)

ggplot()+
  geom_violin(grik2_RPKML_age_filtered,mapping=aes(x=Species,y=RPKM_log,fill=Species),trim=FALSE)+
  # geom_boxplot(grik2_RPKML_age_filtered,mapping=aes(x=Species,y=RPKM_log),width=0.1,outlier.shape = "|")+
  stat_summary(grik2_RPKML_age_filtered,mapping=aes(x=Species,y=RPKM_log),fun.data="mean_sdl", geom="pointrange", color="black")+
  labs(x="Species","log2(RPKM+1)",labs="Species")+
  theme_classic()

## focus on specific Predicted_Age_log

grik2_RPKML_age_8_12.15 <- grik2_RPKML_age_filtered %>% 
  filter(Predicted_Age_log > 8, Predicted_Age_log < 12.5)

ggplot()+
  geom_point(grik2_RPKML_age_8_12.15,mapping=aes(x=Predicted_Age_log,y=RPKM_log,fill=Species,col=Species),size=0.1,shape=23)+
  geom_smooth(grik2_RPKML_age_8_12.15,mapping=aes(x=Predicted_Age_log,y=RPKM_log,col=Species),method = "loess")+
  labs(x="Predicted_Age","log2(RPKM+1)",labs="Species")+
  theme_classic()

wilcox_reseults <- wilcox.test(RPKM_log ~ Species, data = grik2_RPKML_age_8_12.15)

ggplot()+
  geom_violin(grik2_RPKML_age_8_12.15,mapping=aes(x=Species,y=RPKM_log,fill=Species),trim=FALSE)+
  stat_summary(grik2_RPKML_age_8_12.15,mapping=aes(x=Species,y=RPKM_log),fun.data="mean_sdl", geom="pointrange", color="black")+
  labs(x="Species","log2(RPKM+1)",labs="Species")+
  theme_classic()+
  annotate("text", x = 1.5, y = max(grik2_RPKML_age_8_12.15$RPKM_log)+0.2, 
           label = paste("Wilcox=", round(wilcox_reseults$statistic, 2), "\n", "p-value=", (wilcox_reseults$p.value)),
           size = 5, color = "black")

## focus on specific Tissue

grik2_RPKML_age_filtered_tissue <- grik2_RPKML_age_filtered

grik2_RPKML_age_filtered_tissue$Tissue <- gsub("MFC", "medial prefrontal cortex", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("OFC", "orbital prefrontal cortex", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("DFC", "dorsolateral prefrontal cortex", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("VFC", "ventrolateral prefrontal cortex", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("M1C", "primary motor cortex", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("S1C", "primary somatosensory cortex", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("IPC", "inferior posterior parietal cortex", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("A1C", "primary auditory cortex", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("A1C", "primary auditory cortex", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("A1C", "primary auditory cortex", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("STC", "superior temporal cortex", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("ITC", "inferior temporal cortex", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("V1C", "primary visual cortex", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("HIP", "hippocampus", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("AMY", "amygdala", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("STR", "striatum", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("MD", "mediodorsal nucleus of the thalamus", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("CBC", "cerebellar cortex", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("AMY", "amygdala", grik2_RPKML_age_filtered_tissue$Tissue)
grik2_RPKML_age_filtered_tissue$Tissue <- gsub("AMY", "amygdala", grik2_RPKML_age_filtered_tissue$Tissue)

ggplot()+
  geom_point(grik2_RPKML_age_filtered_tissue,mapping=aes(x=Predicted_Age_log,y=RPKM_log,fill=Species,col=Species),size=0.1,shape=23)+
  geom_smooth(grik2_RPKML_age_filtered_tissue,mapping=aes(x=Predicted_Age_log,y=RPKM_log,col=Species),method = "loess")+
  labs(x="Predicted_Age","log2(RPKM+1)",labs="Species")+
  theme_classic()+
  facet_wrap(~Tissue)

grik2_expcom <- grik2_RPKML_age_filtered_tissue %>%
  filter(Tissue == "primary auditory cortex" |
          Tissue == "orsolateral prefrontal cortex" |
          Tissue == "inferior posterior parietal cortex" |
           Tissue == "inferior temporal cortex" |
           Tissue == "primary motor cortex" |
           Tissue == "medial prefrontal cortex" |
           Tissue == "primary somatosensory cortex" |
           Tissue == "orbital prefrontal cortex" |
           Tissue == "primary visual cortex" |
           Tissue == "ventrolateral prefrontal cortex" |
           Tissue == "superior temporal cortex")

grik2_expcom
grik2_expcom_human <- subset(grik2_expcom, Species == "Human")
grik2_expcom_human$Period <- as.numeric(grik2_expcom_human$Period)

grik2_expcom_human$timepoint <- sapply(grik2_expcom_human$Period, function(x) {
  if (x > 10) {
    return("Adulthood")
  } else if (7 <= x & x <= 10) {
    return("Late fetal development & Postnatal")
  } else {
    return("Early fetal development")
  }
})
grik2_expcom_human$timepoint <- as.character(unlist(grik2_expcom_human$timepoint))

grik2_expcom_adulthood <- grik2_expcom_human %>% filter(timepoint == "Adulthood")
grik2_expcom_postnatal <- grik2_expcom_human %>% filter(timepoint == "Late fetal development & Postnatal")
grik2_expcom_early <- grik2_expcom_human %>% filter(timepoint == "Early fetal development")
range(grik2_expcom_adulthood$Predicted_Age)
range(grik2_expcom_postnatal$Predicted_Age)
range(grik2_expcom_early$Predicted_Age)

grik2_expcom$timepoint <- sapply(grik2_expcom$Predicted_Age_log, function(x) {
  if (x > min(grik2_expcom_adulthood$Predicted_Age_log)) {
    return("Adulthood")
  } else if (min(grik2_expcom_postnatal$Predicted_Age_log) <= x & x <= max(grik2_expcom_postnatal$Predicted_Age_log)) {
    return("Late fetal development & Postnatal")
  } else {
    return("Early fetal development")
  }
})
timepoint_order <- fct_relevel(grik2_expcom$timepoint,"Early fetal development","Late fetal development & Postnatal","Adulthood")

ggplot()+
  geom_violin(grik2_expcom,mapping=aes(x=timepoint_order,y=RPKM,fill=Species),trim=FALSE)+
  stat_summary(grik2_expcom,mapping=aes(x=timepoint_order,y=RPKM,group=Species),fun.data="mean_sdl", geom="pointrange",
               color="black", position = position_dodge(0.9))+
  labs(x="Age",y="scaled_RPKM",labs="Species")+
  scale_fill_manual(values = c("#E41A1C","#377EB8"))+
  theme_classic()

boxplot_plot <- ggplot()+
  geom_boxplot(grik2_expcom,mapping=aes(x=timepoint_order,y=RPKM_log,fill=Species))+
  labs(x="Age",y="scaled_RPKM",labs="Species")+
  scale_fill_manual(values = c("#E41A1C","#377EB8"))+
  theme_classic()

point_plot <- ggplot()+
  geom_point(grik2_expcom,mapping=aes(x=Predicted_Age_log,y=RPKM_log,fill=Species,col=Species),size=0.1,shape=23)+
  geom_smooth(grik2_expcom,mapping=aes(x=Predicted_Age_log,y=RPKM_log,col=Species),method = "loess")+
  labs(x="Age",y="scaled_RPKM",labs="Species")+
  scale_color_manual(values = c("#E41A1C","#377EB8"))+
  theme_classic()

combined_plot <- boxplot_plot/point_plot +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = 'A')
combined_plot

ggplot()+
  geom_point(grik2_expcom,mapping=aes(x=Predicted_Age_log,y=RPKM_log,fill=Species,col=Species),size=0.1,shape=23)+
  geom_smooth(grik2_expcom,mapping=aes(x=Predicted_Age_log,y=RPKM_log,col=Species),method = "loess")+
  labs(x="Developmental time points",y="Expression level",labs="Species")+
  theme_classic()+
  scale_color_manual(values = c("#E41A1C","#377EB8"))+
  theme_classic()+
  # facet_wrap(~Tissue)
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

ggsave(file.path(sys_dir,"thesis_fig/fig1_expressiondataoverall.png"), plot=ggplot2::last_plot())

write_xlsx(grik2_expcom, file.path(sys_dir,"Figure_2A.xlsx"))

ggplot()+
  geom_violin(data=subset(grik2_expcom,timepoint="Late fetal development & Postnatal"),mapping=aes(x=Species,y=RPKM_log,fill=Species),trim=FALSE)+
  stat_summary(data=subset(grik2_expcom,timepoint="Late fetal development & Postnatal"),mapping=aes(x=Species,y=RPKM_log),fun.data="mean_sdl", geom="pointrange", color="black")+
  labs(x="Species",y="Expression level",labs="Species")+
  scale_fill_manual(values = c("#E41A1C","#377EB8"))+
  theme_classic()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))+
  annotate("text", x=1.5, y=max(grik2_expcom$RPKM_log), 
           label=paste("****"), vjust=-0.5, hjust=0.5)

ggsave(file.path(sys_dir,"thesis_fig/fig1_expressiondata.png"), plot=ggplot2::last_plot())

dummy_df <- grik2_expcom %>% filter(Predicted_Age_log > 6.392317)

ggplot()+
  geom_smooth(dummy_df,mapping=aes(x=Predicted_Age_log,y=RPKM_log,col=Species),method = "loess")+
  labs(x="Age",y="Expression",labs="Species")+
  # scale_color_manual(values = c("#00539CFF","#EEA47FFF"))+
  scale_color_manual(values = c("#E41A1C","#377EB8"))+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

pre_data <- subset(grik2_expcom, timepoint == "Early fetal development")
post_data <- subset(grik2_expcom, timepoint == "Late fetal development & Postnatal")
adult_data <- subset(grik2_expcom, timepoint == "Adulthood")
t_pre <- t.test(RPKM_log ~ Species, data = pre_data)
t_post <- t.test(RPKM_log ~ Species, data = post_data)
t_adult <- t.test(RPKM_log ~ Species, data = adult_data)
cohen_d_pre <- cohen.d(RPKM_log ~ Species, data = pre_data)
cohen_d_post <- cohen.d(RPKM_log ~ Species, data = post_data)
cohen_d_adult <- cohen.d(RPKM_log ~ Species, data = adult_data)
print(cohen_d_pre$estimate)
print(cohen_d_post$estimate)
print(cohen_d_adult$estimate)

post_data_fold <- post_data %>%
  group_by(Species) %>%
  summarize(mean_RPKM = mean((RPKM), na.rm = TRUE))
human_mean <- post_data_fold %>% filter(Species == "Human") %>% pull(mean_RPKM)
macaque_mean <- post_data_fold %>% filter(Species == "Macaque") %>% pull(mean_RPKM)
fold_difference <- log2(human_mean/macaque_mean)
fold_difference

# difference <- (human_mean-macaque_mean)
# mean <- (human_mean+macaque_mean)/2
# relativechange <- (difference/macaque_mean)*100

## do the analysis per tissue

cortical_tissues <- c("medial prefrontal cortex", "orbital prefrontal cortex", "dorsolateral prefrontal cortex", 
                      "ventrolateral prefrontal cortex", "primary motor cortex", "primary somatosensory cortex", 
                      "inferior posterior parietal cortex", "primary auditory cortex", 
                      "superior temporal cortex", "inferior temporal cortex", "primary visual cortex")
non_cortical_tissues <- c("hippocampus","amygdala","striatum","mediodorsal nucleus of the thalamus","cerebellar cortex")

grik2_expcom_noncortical <- grik2_RPKML_age_filtered_tissue %>%
  filter(Tissue %in% non_cortical_tissues)

ggplot()+
  geom_point(grik2_expcom_noncortical,mapping=aes(x=Predicted_Age_log,y=RPKM_log,fill=Species,col=Species),size=0.1,shape=23)+
  geom_smooth(grik2_expcom_noncortical,mapping=aes(x=Predicted_Age_log,y=RPKM_log,col=Species),method = "loess")+
  labs(x="Age",y="scaled_RPKM",labs="Species")+
  scale_color_manual(values = c("#E41A1C","#377EB8"))+
  theme_classic()+
  facet_wrap(~Tissue)

unique_tissues <- unique(grik2_RPKML_age_filtered_tissue$Tissue)

grik2_expcom_human_pertis <- subset(grik2_RPKML_age_filtered_tissue, Species == "Human")
grik2_expcom_human_pertis$Period <- as.numeric(grik2_expcom_human_pertis$Period)
grik2_expcom_human_pertis$timepoint <- sapply(grik2_expcom_human_pertis$Period, function(x) {
  if (x > 10) {
    return("Adulthood")
  } else if (7 <= x & x <= 10) {
    return("Late fetal development & Postnatal")
  } else {
    return("Early fetal development")
  }
})
grik2_expcom_human$timepoint <- as.character(unlist(grik2_expcom_human$timepoint))
grik2_expcom_pertis_adulthood <- grik2_expcom_human_pertis %>% filter(timepoint == "Adulthood")
grik2_expcom_pertis_postnatal <- grik2_expcom_human_pertis %>% filter(timepoint == "Late fetal development & Postnatal")
grik2_expcom_pertis_early <- grik2_expcom_human_pertis %>% filter(timepoint == "Early fetal development")
grik2_RPKML_age_filtered_tissue$timepoint <- sapply(grik2_RPKML_age_filtered_tissue$Predicted_Age_log, function(x) {
  if (x > min(grik2_expcom_pertis_adulthood$Predicted_Age_log)) {
    return("Adulthood")
  } else if (min(grik2_expcom_pertis_postnatal$Predicted_Age_log) <= x & x <= max(grik2_expcom_pertis_postnatal$Predicted_Age_log)) {
    return("Late fetal development & Postnatal")
  } else {
    return("Early fetal development")
  }
})

test_results <- data.frame(Tissue = character(), Stage = character(), 
                           T_Test_P_Value = numeric(), T_Test_Statistic = numeric(), 
                           T_Test_DF = numeric(), Wilcoxon_P_Value = numeric(), 
                           Wilcoxon_Statistic = numeric())
for (tissue in unique_tissues) {
  tissue_data <- grik2_RPKML_age_filtered_tissue %>%
    filter(Tissue == tissue)
  pre_data <- tissue_data %>%
    filter(timepoint == "Early fetal development")
  post_data <- tissue_data %>%
    filter(timepoint == "Late fetal development & Postnatal")
  adult_data <- tissue_data %>%
    filter(timepoint == "Adulthood")
  if (nrow(pre_data) > 0 && length(unique(pre_data$Species)) == 2) {
    shapiro_pre <- shapiro.test(pre_data$RPKM_log)
    if (shapiro_pre$p.value > 0.05) {
      t_test_pre <- t.test(RPKM_log ~ Species, data = pre_data)
      test_results <- rbind(test_results, data.frame(Tissue = tissue, Stage = "pre", 
                                                     Test_Type = "t-test", 
                                                     P_Value = t_test_pre$p.value, 
                                                     Test_Statistic = t_test_pre$statistic, 
                                                     Degrees_Freedom = t_test_pre$parameter))
    } else {
      wilcox_test_pre <- wilcox.test(RPKM_log ~ Species, data = pre_data)
      test_results <- rbind(test_results, data.frame(Tissue = tissue, Stage = "pre", 
                                                     Test_Type = "Wilcoxon", 
                                                     P_Value = wilcox_test_pre$p.value, 
                                                     Test_Statistic = wilcox_test_pre$statistic, 
                                                     Degrees_Freedom = NA))
    }
  }
  if (nrow(post_data) > 0 && length(unique(post_data$Species)) == 2) {
    shapiro_post <- shapiro.test(post_data$RPKM_log)
    if (shapiro_post$p.value > 0.05) {
      t_test_post <- t.test(RPKM_log ~ Species, data = post_data)
      test_results <- rbind(test_results, data.frame(Tissue = tissue, Stage = "post", 
                                                     Test_Type = "t-test", 
                                                     P_Value = t_test_post$p.value, 
                                                     Test_Statistic = t_test_post$statistic, 
                                                     Degrees_Freedom = t_test_post$parameter))
    } else {
      wilcox_test_post <- wilcox.test(RPKM_log ~ Species, data = post_data)
      test_results <- rbind(test_results, data.frame(Tissue = tissue, Stage = "post", 
                                                     Test_Type = "Wilcoxon", 
                                                     P_Value = wilcox_test_post$p.value, 
                                                     Test_Statistic = wilcox_test_post$statistic, 
                                                     Degrees_Freedom = NA))
    }
  }
  if (nrow(adult_data) > 0 && length(unique(adult_data$Species)) == 2) {
    shapiro_adult <- shapiro.test(adult_data$RPKM_log)
    if (shapiro_adult$p.value > 0.05) {
      t_test_adult <- t.test(RPKM_log ~ Species, data = adult_data)
      test_results <- rbind(test_results, data.frame(Tissue = tissue, Stage = "adult", 
                                                     Test_Type = "t-test", 
                                                     P_Value = t_test_adult$p.value, 
                                                     Test_Statistic = t_test_adult$statistic, 
                                                     Degrees_Freedom = t_test_adult$parameter))
    } else {
      wilcox_test_adult <- wilcox.test(RPKM_log ~ Species, data = adult_data)
      test_results <- rbind(test_results, data.frame(Tissue = tissue, Stage = "adult", 
                                                     Test_Type = "Wilcoxon", 
                                                     P_Value = wilcox_test_adult$p.value, 
                                                     Test_Statistic = wilcox_test_adult$statistic, 
                                                     Degrees_Freedom = NA))
    }
  }
}
test_results

test_results$Adjusted_P_Value <- p.adjust(test_results$P_Value, method = "fdr")
test_results_sig <- test_results %>% filter(Adjusted_P_Value < 0.05)

# point_plot <- ggplot()+
#   # geom_rect(aes(xmin = 8, 
#   #               xmax = 12.5, 
#   #               ymin = -Inf, ymax = Inf), 
#   #           fill = "gray", alpha = 0.8) + 
#   geom_point(grik2_expcom,mapping=aes(x=Predicted_Age_log,y=RPKM_log,fill=Color_Group,col=Color_Group),size=0.1,shape=23)+
#   geom_smooth(grik2_expcom,mapping=aes(x=Predicted_Age_log,y=RPKM_log,col=Color_Group),method = "loess")+
#   labs(x="Age","log2(RPKM+1)",labs="Species")+
#   theme_classic()+
#   scale_fill_manual(values = c("#EC9290","#E41A1C", "#981112" , "#9FBED5" , "#377EB8" , "#24547A"))+
#   scale_color_manual(values = c("#EC9290","#E41A1C", "#981112" , "#9FBED5" , "#377EB8" , "#24547A"))

# ggplot()+
#   geom_point(grik2_expcom,mapping=aes(x=Predicted_Age_log,y=RPKM_log,fill=Species,col=Species),size=0.1,shape=23)+
#   geom_smooth(grik2_expcom,mapping=aes(x=Predicted_Age_log,y=RPKM_log,col=Species),method = "loess")+
#   geom_point(grik2_expcom_human,mapping=aes(x=grik2_expcom_human$Period,y=RPKM_log,fill=Species,col=Species),size=0.1,shape=23)+
#   # geom_smooth(grik2_expcom_human,mapping=aes(x=Peroid,y=RPKM_log,col=Species),method = "loess")+
#   labs(x="Predicted_Age","log2(RPKM+1)",labs="Species")+
#   theme_classic()
# 
# ggplot()+
#   geom_violin(grik2_expcom,mapping=aes(x=Species,y=RPKM_log,fill=Species),trim=FALSE)+
#   stat_summary(grik2_expcom,mapping=aes(x=Species,y=RPKM_log),fun.data="mean_sdl", geom="pointrange", color="black")+
#   labs(x="Species","log2(RPKM+1)",labs="Species")+
#   theme_classic()

grik2_expcom_human <- subset(grik2_expcom, Species == "Human")
grik2_expcom_macaque <- subset(grik2_expcom, Species == "Macaque")
unique(grik2_expcom_human$Age)
human_synapto_ref_Pollen <- "PCW26 - 3years"
macaque_synapto_ref_Pollen <- "PCD(55+110)/2 - 1year"

grik2_expcom_human <- grik2_expcom %>%
  filter(Age == "1 Y" |
           Age == "4 M" |
           Age == "10 M" |
           Age == "37 PCW")
grik2_expcom_macaque <- grik2_expcom %>%
  filter(Age == "E81" |
           Age == "E110" |
           Age == "P0" |
           Age == "E82" | 
           Age == "P2" |
           Age == "7M" |
           Age == "1Y" |
           Age == "E80" |
           Age == "E111")

grik2_expcom_synapto <- grik2_expcom %>%
  filter(Age == "1 Y" |
           Age == "4 M" |
           Age == "10 M" |
           Age == "37 PCW" |
           Age == "E81" |
           Age == "E110" |
           Age == "P0" |
           Age == "E82" | 
           Age == "P2" |
           Age == "7M" |
           Age == "1Y" |
           Age == "E80" |
           Age == "E111")

grik2_expcom$timepoint <- lapply(grik2_expcom$Predicted_Age_log, function(x) if(x>12.5){
  grik2_expcom$timepoint="adult"
}else if
  (8 < x && x <= 12.5){grik2_expcom$timepoint="postnatel"
}else 
  grik2_expcom$timepoint="prenatal")
grik2_expcom$timepoint <- as.character(unlist(grik2_expcom$timepoint))
timepoint_order <- fct_relevel(grik2_expcom$timepoint,"prenatal","postnatel","adult")

ggplot()+
  geom_point(grik2_expcom_synapto,mapping=aes(x=Predicted_Age_log,y=RPKM_log,fill=Species,col=Species),size=0.1,shape=23)+
  geom_smooth(grik2_expcom_synapto,mapping=aes(x=Predicted_Age_log,y=RPKM_log,col=Species),method = "loess")+
  labs(x="Predicted_Age","log2(RPKM+1)",labs="Species")+
  theme_classic()

species_colors <- c("#E41A1C","#377EB8")
species_levels <- c("Human", "Macacque")
color_gradients <- lapply(species_colors, function(color) {
  colorRampPalette(c(color, "#F4F3EE"))(5)
}) <- lapply(species_colors, function(color) {
  colorRampPalette(c(color, "#F4F3EE"))(5)
})
names(color_gradients) <- species_levels
color_mapping <- unlist(lapply(seq_along(color_gradients), function(i) {
  setNames(color_gradients[[i]], paste0(names(color_gradients)[i], "_", 1:3))
}))

grik2_expcom$Color_Group <- paste0(grik2_expcom$Species, "_", as.numeric(timepoint_order))

boxplot <- ggplot()+
  geom_boxplot(grik2_expcom,mapping=aes(x=timepoint_order,y=RPKM_log,fill=Color_Group))+
  labs(x="Age","log2(RPKM+1)",labs="Species")+
  theme_classic()+
  theme(axis.line.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  scale_fill_manual(values = c("#EC9290","#E41A1C", "#981112" , "#9FBED5" , "#377EB8" , "#24547A"))

point_plot <- ggplot()+
  # geom_rect(aes(xmin = 8, 
  #               xmax = 12.5, 
  #               ymin = -Inf, ymax = Inf), 
  #           fill = "gray", alpha = 0.8) + 
  geom_point(grik2_expcom,mapping=aes(x=Predicted_Age_log,y=RPKM_log,fill=Color_Group,col=Color_Group),size=0.1,shape=23)+
  geom_smooth(grik2_expcom,mapping=aes(x=Predicted_Age_log,y=RPKM_log,col=Color_Group),method = "loess")+
  labs(x="Age","log2(RPKM+1)",labs="Species")+
  theme_classic()+
  scale_fill_manual(values = c("#EC9290","#E41A1C", "#981112" , "#9FBED5" , "#377EB8" , "#24547A"))+
  scale_color_manual(values = c("#EC9290","#E41A1C", "#981112" , "#9FBED5" , "#377EB8" , "#24547A"))

combined_plot <- boxplot/point_plot +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = 'A')
combined_plot

# Developmental rhesus and human data for all genes -----------------------

memory.limit()

encode_RPKM <- read.table(file.path(sys_dir,"PsychENCODE/nhp_development_RPKM_rmTechRep.txt"),header=T)
encode_RPKM_table <- as.data.table(encode_RPKM)
encode_RPKM_table <- melt(encode_RPKM_table, id.vars = "gene_name", variable.name = "Sample", value.name = "expression")
encode_RPKM_table[, c("Sample", "Tissue") := tstrsplit(tissue, "\\.", perl = TRUE)]
colnames(encode_RPKM_table) <- c("gene_name", "Sample", "RPKM", "ID", "Tissue")
sample_metadata <- read_excel(file.path(sys_dir,"PsychENCODE/aat8077_tables-s1-s3.xlsx"),sheet=1) %>%
  dplyr::select("Sample","Species","Age","Predicted age..PC.Days.","Period")
encode_RPKML_age <- merge(encode_RPKM_table,sample_metadata)
colnames(encode_RPKML_age) <- c("Sample","gene_name","RPKM","ID","Tissue","Species","Age","Predicted_Age","Period")
encode_RPKML_age$RPKM_log <- log2((encode_RPKML_age$RPKM+1))
encode_RPKML_age$Predicted_Age_log <- log2((encode_RPKML_age$Predicted_Age))
encode_RPKML_age_filtered <- encode_RPKML_age %>% 
  filter(Tissue != "CGE", Tissue != "DTH", Tissue != "LGE", Tissue != "MGE", Tissue != "OC", Tissue != "PC",
         Tissue != "TC", Tissue != "DTH", Tissue != "URL", Tissue != "MSC",
         Tissue != "VFC_001", Tissue != "STC_GAIIx", Tissue != "A1C_GAIIx", Tissue != "M1C_GAIIx")
rm(encode_RPKM,encode_RPKM_table)

encode_RPKML_avgexp <- encode_RPKML_age_filtered %>%
  group_by(Species, gene_name) %>%
  summarise(avg_RPKM = mean(RPKM, na.rm = TRUE),
            sd_RPKM = sd(RPKM))
split_gene_name <- strsplit(as.character(encode_RPKML_avgexp$gene_name), split = "\\|")
split_df <- do.call(rbind, split_gene_name)
encode_RPKML_avgexp$ensemblID <- split_df[, 1]
encode_RPKML_avgexp$gene_symbol <- split_df[, 2]
encode_RPKML_avgexp$gene_name <- NULL
colnames(encode_RPKML_avgexp) <- c("Species","avg_RPKM","sd_RPKM","ensemblID","geneSymbol")
encode_RPKML_avgexp_human <- encode_RPKML_avgexp %>% filter(Species == "Human")
encode_RPKML_avgexp_human$RPKM_Zscore <- ((encode_RPKML_avgexp_human$avg_RPKM - (mean(encode_RPKML_avgexp_human$avg_RPKM)))/(sd(encode_RPKML_avgexp_human$avg_RPKM)))
colnames(encode_RPKML_avgexp_human) <- c("Species","avg_RPKM","sd_RPKM","ensemblID","geneSymbol","RPKM_Zscore")

ggplot()+
  geom_point(data=subset(encode_RPKML_age_filtered,gene_name == "ENSG00000129535.12|NRL"),
             mapping=aes(x=Predicted_Age_log,y=RPKM_log,fill=Species,col=Species),size=0.1,shape=23)+
  geom_smooth(data=subset(encode_RPKML_age_filtered,gene_name == "ENSG00000129535.12|NRL"),
              mapping=aes(x=Predicted_Age_log,y=RPKM_log,col=Species),method = "loess")+
  labs(x="Predicted_Age","log2(RPKM+1)",labs="Species")+
  theme_classic()+
  facet_wrap(~Tissue)

# WGCNA - for developmental time points (uncompleted) -------------------------------------------------------------------

## define your data

# memory.limit()
# 
# encode_WGCNA <- read.table(file.path(sys_dir,"PsychENCODE/human_encode_data_filtered.txt"),header=T)
# encode_WGCNA_table <- as.data.table(encode_WGCNA)
# gene_names <- encode_WGCNA$gene_name
# encode_WGCNA <- as.data.frame(lapply(encode_WGCNA, as.numeric))
# encode_WGCNA <- encode_WGCNA[,-1]
# rownames(encode_WGCNA) <- gene_names
# str(encode_WGCNA)
# anyNA(encode_WGCNA)
# encode_WGCNA_t = t(encode_WGCNA)
# 
# encode_WGCNA_longer <- encode_WGCNA %>%
#   pivot_longer(cols=-1, names_to="Sample",values_to="Expression")
# sample_metadata <- read_excel(file.path(sys_dir,"PsychENCODE/aat8077_tables-s1-s3.xlsx"),sheet=1) %>%
#   dplyr::select("Sample","Species","Age","Predicted age..PC.Days.","Period")
# encode_WGCNA_age <- merge(encode_WGCNA_longer,sample_metadata)
# encode_WGCNA_age <- encode_WGCNA_age %>% select(c("Sample","Expression","Age","Predicted age..PC.Days."))
# colnames(encode_WGCNA_age) <- c( "ID","Expression","Age","Predicted_Age")
# encode_WGCNA_age$Predicted_Age_log <- log2((encode_WGCNA_age$Predicted_Age))

memory.limit()

encode_WGCNA <- read.table(file.path(sys_dir,"PsychENCODE/human_encode_data_filtered.txt"),header=T)
encode_WGCNA_table <- as.data.table(encode_WGCNA)
encode_WGCNA_table <- melt(encode_WGCNA_table, id.vars = "gene_name", variable.name = "tissue", value.name = "expression")
encode_WGCNA_table[, c("Sample", "Tissue") := tstrsplit(tissue, "\\.", perl = TRUE)]
colnames(encode_WGCNA_table) <- c("gene_name", "Sample", "RPKM", "ID", "Tissue")
sample_metadata <- read_excel(file.path(sys_dir,"PsychENCODE/aat8077_tables-s1-s3.xlsx"),sheet=1) %>%
  dplyr::select("Sample","Species","Age","Predicted age..PC.Days.","Period")
encode_WGCNA_table_age <- merge(encode_WGCNA_table,sample_metadata)
colnames(encode_WGCNA_table_age) <- c("Sample","gene_name","RPKM","ID","Tissue","Species","Age","Predicted_Age","Period")
encode_WGCNA_table_age$RPKM_log <- log2((encode_WGCNA_table_age$RPKM+1))
encode_WGCNA_table_age$Predicted_Age_log <- log2((encode_WGCNA_table_age$Predicted_Age))
WGCNA_encode <- encode_WGCNA_table_age %>% 
  filter(Species == "Human") %>%
  dplyr::select(c(Predicted_Age_log, gene_name, RPKM))
WGCNA_encode_wide <- pivot_wider(WGCNA_encode, names_from=Predicted_Age_log, values_from=RPKM, values_fn=mean) %>%
  as.data.frame()
rownames(WGCNA_encode_wide) <- WGCNA_encode_wide$gene_name
WGCNA_encode_wide <- WGCNA_encode_wide[, -1]
str(WGCNA_encode_wide)
anyNA(WGCNA_encode_wide)

WGCNA_encode_wide_t = t(WGCNA_encode_wide)

## choose good power

powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(WGCNA_encode_wide_t, powerVector = powers, verbose = 5)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     type = "n", xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2", main = "Scale Independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     labels = powers, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     type = "n", xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", main = "Mean Connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], 
     labels = powers, col = "red")
softPower <- 7  # The best power from the plot (start of the curve)

## run analysis

net_mod <- WGCNA::blockwiseModules(WGCNA_encode_wide_t, power = softPower, TOMType = "unsigned",
                        minModuleSize = 30, reassignThreshold = 0,
                        mergeCutHeight = 0.25, numericLabels = TRUE,
                        pamRespectsDendro = FALSE)

net_mod_sign <- WGCNA::blockwiseModules(WGCNA_encode_wide_t, power = softPower,
                                        TOMType = "signed", ,
                            minModuleSize = 30, reassignThreshold = 0,
                            mergeCutHeight = 0.25, numericLabels = TRUE,
                            pamRespectsDendro = FALSE)
BiocManager::install("WGCNA", force=TRUE)

# moduleColors <- labels2colors(net_mod$colors)
plotDendroAndColors(net_mod$dendrograms[[1]], moduleColors, "Module Colors")

## use https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html#ref-Langfelder2018

bwnet <- blockwiseModules(WGCNA_encode_wide_t,
                          networkType = "signed",
                          TOMType = "signed", # topological overlap matrix
                          power = 7, # soft threshold for network construction
                          numericLabels = TRUE, # Let's use numbers instead of colors for module labels
                          randomSeed = 1234) # there's some randomness associated with this calculation

all.equal(metadata$refinebio_accession_code, rownames(module_eigengenes))

des_mat <- model.matrix(~ metadata$time_point)
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)
fit <- limma::eBayes(fit)

stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

## run nonsense

adjacency <- adjacency(encode_WGCNA_t, power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, 
                             deepSplit = 2, pamRespectsDendro = FALSE, 
                             minClusterSize = 30)
dynamicColors <- labels2colors(dynamicMods)

# Re-cluster genes using the dissimilarity TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, 
                             deepSplit = 2, pamRespectsDendro = FALSE, 
                             minClusterSize = 30)
dynamicColors <- labels2colors(dynamicMods)

mergedColors = labels2colors(net$colors)
plotDendroAndColors(
  net$dendrograms[[1]],
  mergedColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

gene_of_interest <- "ENSG00000164418.19|GRIK2"

# # Convert encode_WGCNA_t to a long format
# encode_WGCNA_long <- as.data.frame(as.table(encode_WGCNA_t))
# colnames(encode_WGCNA_long) <- c("ID", "Gene", "Expression")
# merged_data <- merge(encode_WGCNA_long, encode_WGCNA_age)
# str(merged_data)

## module membership and eigengenes (MM & ME) - Correlation between genes andMEs

MEs <- moduleEigengenes(WGCNA_encode_wide_t, colors = moduleColors)$eigengenes

MEs <- net_mod$MEs
MM <- as.data.frame(cor(WGCNA_encode_wide_t, MEs, use = "p"))
pval_MM <- apply(MM, 2, function(x) corPvalueStudent(x, nSamples = nrow(WGCNA_encode_wide_t))) %>% as.data.frame()

MM <- as.data.frame(cor(WGCNA_encode_wide_t, MEs, use = "p"))
Predicted_Age_log <- as.numeric(rownames(WGCNA_encode_wide_t))

GS <- apply(WGCNA_encode_wide_t, 2, function(gene_expr) {
  abs(cor(gene_expr, Predicted_Age_log, use = "complete.obs"))
})
GS <- as.data.frame(GS)
rownames(GS) <- colnames(WGCNA_encode_wide_t)  # Set gene names as row names
GS$gene_names <- row.names(GS)

## hub genes

hub_genes <- MM[apply(MM, 1, function(x) any(x > 0.8)) & GS$GS > 0.2, , drop = FALSE]
hub_genes <- hub_genes %>% drop_na()
nrow(MM)
nrow(hub_genes)

grik2_module <- net_mod$colors[which(colnames(WGCNA_encode_wide_t) == grik2_id)]
grik2_module <- colnames(WGCNA_encode_wide_t)[net_mod$colors == grik2_module]
length(grik2_module)

# module_colors <- net_mod$colors  
plot_MM_GS <- function(module) {
  module_genes <- grik2_module
  verboseScatterplot(abs(MM[module_genes,]), abs(GS[module_genes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for Age",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}
module_of_interest <- "red"
plot_MM_GS(module_of_interest)

# Proper filtering of MM and GS based on the given thresholds (MM > 0.4 and GS > 0.5)
filtered_genes <- MM[apply(MM, 1, function(x) any(x < 0.05)) & GS$GS < 0.05, , drop = FALSE]
filtered_genes <- filtered_genes %>% drop_na()

module_colors <- net_mod$colors  
MS <- sapply(unique(module_colors), function(mod) {
  mean(GS[module_colors == mod, ])
})

MM_grik2_module <- MM[grik2_module, ]
GS_grik2_module <- GS[grik2_module, ]

MM_GS_grik2 <- merge(MM_grik2_module, GS_grik2_module, by = "row.names") %>%
  mutate(MM_average = mean(c_across(starts_with("ME")), na.rm = TRUE)) %>%
  as.data.frame() %>%
  dplyr::select(gene_names,MM_average,MEred,GS) %>% 
  filter(MEred < 0)

ggplot(MM_GS_grik2, aes(x=MEred,y=GS))+
  geom_point()+
  labs(x="MM",y="GS")+
  theme_bw()

cor.test(MM_GS_grik2$MEred,MM_GS_grik2$GS)

# hub_genes <- which(pval_MM < 0.05, arr.ind=TRUE)
# MM_hub_genes <- MM[hub_genes, ]
# grik2_MM_hub_genes <- hub_genes[grik2_id,]
# nrow(MM)
# nrow(pval_MM)
# nrow(hub_genes)
# nrow(MM_hub_genes)

# expr_data_long <- encode_WGCNA_t %>%
#   as.data.frame() %>%
#   mutate(Sample = rownames(.)) %>%
#   pivot_longer(cols = -Sample, names_to = "Gene", values_to = "Expression")

# GS <- as.data.frame(cor(WGCNA_encode_wide_t, WGCNA_encode$Predicted_Age_log, uas.data.frame()))
# GS <- as.data.frame(cor(encode_WGCNA_long, encode_WGCNA_age$Predicted_Age_log, use = "p"))
# pval_GS <- apply(GS, 2, function(x) corPvalueStudent(x, nSamples = nrow(encode_WGCNA_t)))
# pval_GS <- apply(GS, 2, function(x) corPvalueStudent(x, nSamples = nrow(encode_WGCNA_t)))

## grik2 module and co genes

grik2_id <- "ENSG00000164418.19|GRIK2"

gene_module_df <- data.frame(gene_id = names(net_mod$colors),
                             colors = labels2colors(net_mod$colors))
grik2_module_df <- gene_module_df %>% filter(gene_id == grik2_id)

grik2_MM <- MM[grik2_id,]
grik2_pval_MM <- pval_MM[grik2_id,]

grik2_module <- net_mod$colors[which(colnames(WGCNA_encode_wide_t) == grik2_id)]
grik2_module <- colnames(WGCNA_encode_wide_t)[net_mod$colors == grik2_module]
length(grik2_module)

grik2_module_gene_names <- sapply(strsplit(grik2_module, "\\|"), `[`, 2)
writeLines(grik2_module_gene_names, "grik2_module_gene_names.txt")

plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module Colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)

#Defining the variable Peso10dias containing the column Peso10dias of datTrait
Peso10dias = as.data.frame(datTraits$Peso10dias)
names(Peso10dias) = "Peso10d"

#names (colors) of the modules
modNames = substring(names(MEs), 3)
nSamples = nrow(WGCNA_encode_wide_t)
MM <- as.data.frame(cor(WGCNA_encode_wide_t, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(MM), nSamples))
names(MM) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
Predicted_Age_log <- as.numeric(rownames(WGCNA_encode_wide_t))
GS <- apply(WGCNA_encode_wide_t, 2, function(gene_expr) {
  abs(cor(gene_expr, Predicted_Age_log, use = "complete.obs"))
})
GS <- as.data.frame(GS)
rownames(GS) <- colnames(WGCNA_encode_wide_t)  # Set gene names as row names
GS$gene_names <- row.names(GS)
GSPvalue <- apply(WGCNA_encode_wide_t, 2, function(GS) {
  ((corPvalueStudent(as.matrix(GS), nSamples)))
})
GSPvalue <- GSPvalue %>% as.data.frame()
# names(GS) = paste("GS.", names(Peso10dias), sep="")
# names(GSPvalue) = paste("p.GS.", names(Peso10dias), sep="")

module = "red"
column = match(module, modNames)
moduleGenes = moduleColors==module

#sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(MM[moduleGenes, column]),
                   abs(GS[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Peso 10 dias",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

moduleOfInterest <- moduleColors[which(colnames(WGCNA_encode_wide_t) == grik2_id)]

# Extract genes in the module of interest
genesInModule <- colnames(WGCNA_encode_wide_t)[which(moduleColors == moduleOfInterest)]

# Calculate module membership (MM) for each gene
MM <- as.data.frame(cor(WGCNA_encode_wide_t, MEs[, moduleOfInterest], use = "pairwise.complete.obs"))
MEs <- moduleEigengenes(WGCNA_encode_wide_t, colors = moduleColors)$eigengenes

colnames(MM) <- "ModuleMembership"

# get module eigengenes per cluster

MEs0 <- moduleEigengenes(WGCNA_encode_wide_t, moduleColors)$eigengenes
MEs0 <- orderMEs(MEs0)
module_order <- names(MEs0) %>% gsub("ME", "", .)

MEs0$trait = row.names(MEs0)

mME = MEs0 %>%
  pivot_longer(-trait) %>%
  mutate(name = gsub("ME", "", name),
         name = factor(name, levels = module_order))

mME$trait <- unique(mME$trait)
trait_order <- fct_relevel(mME$trait, "5.8073549220576","5.97727992349992",
                           "6.39231742277876","6.5077946401987","6.8073549220576","6.89481776330794",
                           "7.05528243550119","7.19967234483636","7.2667865406949",
                           "8.01680828768655","8.59991284218713",
                           "9.14720492494223","9.30149619498255",
                           "10.1984450414524","10.3465137331656","10.7540523675289",
                           "11.6384359139905","11.6786001390995",
                           "12.0306671362469","12.2917462808339","12.4880911775335","12.8149829364268","12.9541963103869",
                           "13.0813169892858","13.4542992936199","13.7115594341479","13.7502882675917","13.860699032612")

mME %>%
  ggplot(aes(x=trait_order, y=name, fill=value))+
  geom_tile()+
  theme_bw()+
  # scale_fill_viridis_c(option = "C", limits = c(-1, 1), oob = scales::squish)+
  scale_fill_gradient2(low="blue", high="red", mid="white",
                       midpoint=0,
                       limit=c(-1,1))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"), 
        axis.title.y = element_text(color = "black"))+
labs(title="Expression module over developmental time)",
     x="Developmental time", y="Modules", fill="corr")

mME_red <- mME %>% filter(name == "red")

ggplot(data=mME_red,aes(x=trait,y=name,fill=value))+
  geom_tile()+
  theme_bw()+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.5, 0.5)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(color = "black"))+  
  labs(title = "Expression module over developmental time",
       x = "Developmental time", y = "Grik2 Module", fill = "corr")

length(grik2_module)
expr_grik2_WGCNA <- WGCNA_encode_wide[grik2_module,]

expr_grik2_WGCNA_df = data.frame(expr_grik2_WGCNA) %>%
  mutate(gene_id = row.names(.)) %>%
  pivot_longer(-gene_id)
  
expr_grik2_WGCNA_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = "red"),
            alpha = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(x = "Developmental time", y = "Expression")

summary_expr_grik2_WGCNA_df <- expr_grik2_WGCNA_df %>%
  group_by(name) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
    se = sd(value, na.rm = TRUE) / sqrt(n()),
    lower_ci = mean_value - qt(0.975, df = n() - 1) * se,
    upper_ci = mean_value + qt(0.975, df = n() - 1) * se)

summary_expr_grik2_WGCNA_df %>% ggplot(., aes(x=name, y=mean_value)) +
  geom_point(aes(color = "red"),
            alpha = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(x = "Developmental time", y = "Expression")

# WGCNA - for species (uncompleted) -------------------------------------------------------------------

## define your data

colnames(encode_RPKML_age_filtered)
encode_RPKML_age_filtered_mod <- encode_RPKML_age_filtered %>% filter(encode_RPKML_age_filtered$Predicted_Age_log > 8 &
                                                                      encode_RPKML_age_filtered$Predicted_Age_log < 13.5)
WGCNA_species_wide <- pivot_wider(encode_RPKML_age_filtered_mod, names_from=Species, values_from=RPKM, values_fn=mean) %>%
  as.data.frame()
rownames(WGCNA_species_wide) <- WGCNA_species_wide$gene_name
WGCNA_species_wide <- WGCNA_species_wide[, -1]
WGCNA_species_wide <- WGCNA_species_wide %>% dplyr::select(c(RPKM_log))
str(WGCNA_species_wide)
anyNA(WGCNA_species_wide)
WGCNA_species_wide_t = t(WGCNA_species_wide)

## choose good power

powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(WGCNA_encode_wide_t, powerVector = powers, verbose = 5)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     type = "n", xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2", main = "Scale Independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     labels = powers, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     type = "n", xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", main = "Mean Connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], 
     labels = powers, col = "red")
softPower <- 7  # The best power from the plot (start of the curve)

## run analysis

net_mod_species <- blockwiseModules(WGCNA_encode_wide_t, power = softPower, TOMType = "unsigned", ## should be signed
                            minModuleSize = 30, reassignThreshold = 0,
                            mergeCutHeight = 0.25, numericLabels = TRUE,
                            pamRespectsDendro = FALSE)
moduleColors <- labels2colors(net_mod_species$colors)

## module membership and eigengenes (MM & ME) - Correlation between genes andMEs

MEs_species <- moduleEigengenes(WGCNA_encode_wide_t, colors = moduleColors)$eigengenes
MEs <- net_mod_species$MEs
MM <- as.data.frame(cor(WGCNA_encode_wide_t, MEs, use = "p"))
pval_MM <- apply(MM, 2, function(x) corPvalueStudent(x, nSamples = nrow(WGCNA_encode_wide_t))) %>% as.data.frame()

# hub_genes <- which(pval_MM < 0.05, arr.ind=TRUE)
# MM_hub_genes <- MM[hub_genes, ]
# grik2_MM_hub_genes <- hub_genes[grik2_id,]

Predicted_Age_log <- as.numeric(rownames(WGCNA_encode_wide_t))

GS <- apply(WGCNA_encode_wide_t, 2, function(gene_expr) {
  abs(cor(gene_expr, Predicted_Age_log, use = "complete.obs"))
})
GS <- as.data.frame(GS)
rownames(GS) <- colnames(WGCNA_encode_wide_t)  # Set gene names as row names
GS$gene_names <- row.names(GS)

# Proper filtering of MM and GS based on the given thresholds (MM > 0.4 and GS > 0.5)
filtered_genes <- MM[apply(MM, 1, function(x) any(x < 0.05)) & GS$GS < 0.05, , drop = FALSE]
filtered_genes <- filtered_genes %>% drop_na()

module_colors <- net_mod$colors  
MS <- sapply(unique(module_colors), function(mod) {
  mean(GS[module_colors == mod, ])
})

MM_grik2_module <- MM[grik2_module, ]
GS_grik2_module <- GS[grik2_module, ]

MM_GS_grik2 <- merge(MM_grik2_module, GS_grik2_module, by = "row.names") %>%
  mutate(MM_average = mean(c_across(starts_with("ME")), na.rm = TRUE)) %>%
  as.data.frame() %>%
  dplyr::select(gene_names,MM_average,MEred,GS) %>% 
  filter(MEred < 0)

ggplot(MM_GS_grik2, aes(x=MEred,y=GS))+
  geom_point()+
  labs(x="MM",y="GS")+
  theme_bw()

cor.test(MM_GS_grik2$MEred,MM_GS_grik2$GS)

## grik2 module and co genes

grik2_id <- "ENSG00000164418.19|GRIK2"

gene_module_df <- data.frame(gene_id = names(net_mod$colors),
                             colors = labels2colors(net_mod$colors))
grik2_module_df <- gene_module_df %>% filter(gene_id == grik2_id)

grik2_MM <- MM[grik2_id,]
grik2_pval_MM <- pval_MM[grik2_id,]

grik2_module <- net_mod$colors[which(colnames(WGCNA_encode_wide_t) == grik2_id)]
grik2_module <- colnames(WGCNA_encode_wide_t)[net_mod$colors == grik2_module]
length(grik2_module)

grik2_module_gene_names <- sapply(strsplit(grik2_module, "\\|"), `[`, 2)
writeLines(grik2_module_gene_names, "grik2_module_gene_names.txt")

plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module Colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)

#Defining the variable Peso10dias containing the column Peso10dias of datTrait
Peso10dias = as.data.frame(datTraits$Peso10dias)
names(Peso10dias) = "Peso10d"

#names (colors) of the modules
modNames = substring(names(MEs), 3)
nSamples = nrow(WGCNA_encode_wide_t)
MM <- as.data.frame(cor(WGCNA_encode_wide_t, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(MM), nSamples))
names(MM) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
Predicted_Age_log <- as.numeric(rownames(WGCNA_encode_wide_t))
GS <- apply(WGCNA_encode_wide_t, 2, function(gene_expr) {
  abs(cor(gene_expr, Predicted_Age_log, use = "complete.obs"))
})
GS <- as.data.frame(GS)
rownames(GS) <- colnames(WGCNA_encode_wide_t)  # Set gene names as row names
GS$gene_names <- row.names(GS)
GSPvalue <- apply(WGCNA_encode_wide_t, 2, function(GS) {
  ((corPvalueStudent(as.matrix(GS), nSamples)))
})
GSPvalue <- GSPvalue %>% as.data.frame()
# names(GS) = paste("GS.", names(Peso10dias), sep="")
# names(GSPvalue) = paste("p.GS.", names(Peso10dias), sep="")

module = "red"
column = match(module, modNames)
moduleGenes = moduleColors==module

#sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(MM[moduleGenes, column]),
                   abs(GS[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Peso 10 dias",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

moduleOfInterest <- moduleColors[which(colnames(WGCNA_encode_wide_t) == grik2_id)]

# Extract genes in the module of interest
genesInModule <- colnames(WGCNA_encode_wide_t)[which(moduleColors == moduleOfInterest)]

# Calculate module membership (MM) for each gene
MM <- as.data.frame(cor(WGCNA_encode_wide_t, MEs[, moduleOfInterest], use = "pairwise.complete.obs"))
MEs <- moduleEigengenes(WGCNA_encode_wide_t, colors = moduleColors)$eigengenes

colnames(MM) <- "ModuleMembership"

# get module eigengenes per cluster

MEs0 <- moduleEigengenes(WGCNA_encode_wide_t, moduleColors)$eigengenes
MEs0 <- orderMEs(MEs0)
module_order <- names(MEs0) %>% gsub("ME", "", .)

MEs0$trait = row.names(MEs0)

mME = MEs0 %>%
  pivot_longer(-trait) %>%
  mutate(name = gsub("ME", "", name),
         name = factor(name, levels = module_order))

mME$trait <- unique(mME$trait)
trait_order <- fct_relevel(mME$trait, "5.8073549220576","5.97727992349992",
                           "6.39231742277876","6.5077946401987","6.8073549220576","6.89481776330794",
                           "7.05528243550119","7.19967234483636","7.2667865406949",
                           "8.01680828768655","8.59991284218713",
                           "9.14720492494223","9.30149619498255",
                           "10.1984450414524","10.3465137331656","10.7540523675289",
                           "11.6384359139905","11.6786001390995",
                           "12.0306671362469","12.2917462808339","12.4880911775335","12.8149829364268","12.9541963103869",
                           "13.0813169892858","13.4542992936199","13.7115594341479","13.7502882675917","13.860699032612")

mME %>%
  ggplot(aes(x=trait_order, y=name, fill=value))+
  geom_tile()+
  theme_bw()+
  # scale_fill_viridis_c(option = "C", limits = c(-1, 1), oob = scales::squish)+
  scale_fill_gradient2(low="blue", high="red", mid="white",
                       midpoint=0,
                       limit=c(-1,1))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"), 
        axis.title.y = element_text(color = "black"))+
  labs(title="Expression module over developmental time)",
       x="Developmental time", y="Modules", fill="corr")

mME_red <- mME %>% filter(name == "red")

ggplot(data=mME_red,aes(x=trait,y=name,fill=value))+
  geom_tile()+
  theme_bw()+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.5, 0.5)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(color = "black"))+  
  labs(title = "Expression module over developmental time",
       x = "Developmental time", y = "Grik2 Module", fill = "corr")

length(grik2_module)
expr_grik2_WGCNA <- WGCNA_encode_wide[grik2_module,]

expr_grik2_WGCNA_df = data.frame(expr_grik2_WGCNA) %>%
  mutate(gene_id = row.names(.)) %>%
  pivot_longer(-gene_id)

expr_grik2_WGCNA_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = "red"),
            alpha = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(x = "Developmental time", y = "Expression")

summary_expr_grik2_WGCNA_df <- expr_grik2_WGCNA_df %>%
  group_by(name) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            se = sd(value, na.rm = TRUE) / sqrt(n()),
            lower_ci = mean_value - qt(0.975, df = n() - 1) * se,
            upper_ci = mean_value + qt(0.975, df = n() - 1) * se)

summary_expr_grik2_WGCNA_df %>% ggplot(., aes(x=name, y=mean_value)) +
  geom_point(aes(color = "red"),
             alpha = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(x = "Developmental time", y = "Expression")

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


