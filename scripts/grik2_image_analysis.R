#### The Roles of GRIK2 and HARs in the Evolution and Development of Synapses
#### Shadi Shahatit, Master's Thesis, IBENS, Paris 2024
# libraries ---------------------------------------------------------------

library("ggsci")
library("tidyverse")
library("ggplot2")
library("tidyr")
library("dplyr")
library("readxl")
library("dplyr")
library("tidyverse")
library("purrr")
library("fs")
library("patchwork")
library("viridis")
library("writexl")

## note: replace your system's directory in sys_dir

# Axon branching analysis p21 -------------------------------------------------

sys_dir <- "/Users/Shadi Shahatit/OneDrive/Desktop/"

## shCtrl data 

samples_shCtrl <- c("IUEPau_B1_shCtrl_mVenus_p21_10x_001",
                    "IUEPau_B1_shCtrl_mVenus_p21_10x_004",
                    "IUEPau_B1_shCtrl_mVenus_p21_10x_008",
                    "IUEPau_B4_shCtrl_mVenus_p21_10x_001",
                    "IUEPau_B4_shCtrl_mVenus_p21_10x_002",
                    "IUEPau_B4_shCtrl_mVenus_p21_10x_003",
                    "IUEPau_B6_shCtrl_mVenus_p21_10x_001",
                    "IUEPau_B6_shCtrl_mVenus_p21_10x_002",
                    "IUEPau_B6_shCtrl_mVenus_p21_10x_003")

shCtrl_AxonBranch <- read.table(paste0(sys_dir,"AxonBranch_imageJ_output/AxonBranch_shCtrl_p21_B_1_4_6.csv"), header=T, sep=",") 
shCtrl_AxonBranch$layer <- rep(c("layer_5", "layer_2_3", "layer_4"), length.out = nrow(shCtrl_AxonBranch))
shCtrl_AxonBranch$sample <- rep(samples_shCtrl, each= 3)

## define your background --> every three entries

mean_gray_value_seq <- seq(3, nrow(shCtrl_AxonBranch), by = 3)
shCtrl_AxonBranch$CTCF <- (shCtrl_AxonBranch$IntDen)-(shCtrl_AxonBranch$Mean[mean_gray_value_seq]*shCtrl_AxonBranch$Area)

shCtrl_AxonBranch_filtered <- shCtrl_AxonBranch[-mean_gray_value_seq, ] %>%
  select(c("sample","layer","CTCF"))
row.names(shCtrl_AxonBranch_filtered) <- seq(1:nrow(shCtrl_AxonBranch_filtered))

shCtrl_AxonBranch_CTCF <- shCtrl_AxonBranch_filtered %>%
  pivot_wider(names_from = layer, values_from = CTCF)

shCtrl_AxonBranch_CTCF$condition <- "shCtrl"

## shGrik2 data 

samples_shGrik2 <- c("IUEPau_B7_shGrik2_mVenus_p21_10x_001",
                     "IUEPau_B7_shGrik2_mVenus_p21_10x_002",
                     "IUEPau_B7_shGrik2_mVenus_p21_10x_003",
                     "IUEPau_B8_shGrik2_mVenus_p21_10x_001",
                     "IUEPau_B8_shGrik2_mVenus_p21_10x_002",
                     "IUEPau_B8_shGrik2_mVenus_p21_10x_003",
                     "IUEPau_B9_shGrik2_mVenus_p21_10x_001",
                     "IUEPau_B9_shGrik2_mVenus_p21_10x_002",
                     "IUEPau_B9_shGrik2_mVenus_p21_10x_003")

shGrik2_AxonBranch <- read.table(paste0(sys_dir,"AxonBranch_imageJ_output/AxonBranch_shGrik2_p21_B_7_8_9.csv"), header=T, sep=",") 
shGrik2_AxonBranch$layer <- rep(c("layer_5", "layer_2_3", "layer_4"), length.out = nrow(shGrik2_AxonBranch))
shGrik2_AxonBranch$sample <- rep(samples_shGrik2, each= 3)

## define your background --> every three entries

mean_gray_value_seq <- seq(3, nrow(shGrik2_AxonBranch), by = 3)
shGrik2_AxonBranch$CTCF <- (shGrik2_AxonBranch$IntDen)-(shGrik2_AxonBranch$Mean[mean_gray_value_seq]*shGrik2_AxonBranch$Area)

shGrik2_AxonBranch_filtered <- shGrik2_AxonBranch[-mean_gray_value_seq, ] %>%
  select(c("sample","layer","CTCF"))
row.names(shGrik2_AxonBranch_filtered) <- seq(1:nrow(shGrik2_AxonBranch_filtered))

shGrik2_AxonBranch_CTCF <- shGrik2_AxonBranch_filtered %>%
  pivot_wider(names_from = layer, values_from = CTCF)

shGrik2_AxonBranch_CTCF$condition <- "shGrik2"

## concatenate data 

AxonBranch_CTCF_p21 <- rbind(shCtrl_AxonBranch_CTCF,shGrik2_AxonBranch_CTCF) %>% as.data.frame()

AxonBranch_CTCF_p21$L5overL2_3 <- AxonBranch_CTCF_p21$layer_5/AxonBranch_CTCF_p21$layer_2_3

t.test(layer_5 ~ condition, data = AxonBranch_CTCF_p21)
wilcox.test(layer_5 ~ condition, data = AxonBranch_CTCF_p21)

t.test(L5overL2_3 ~ condition, data = AxonBranch_CTCF_p21)
wilcox.test(L5overL2_3 ~ condition, data = AxonBranch_CTCF_p21)

plot(density(log10(AxonBranch_CTCF_p21$L5overL2_3)))
plot(density(log10(shGrik2_AxonBranch_CTCF$layer_5)))

ggplot()+
  geom_boxplot(AxonBranch_CTCF_p21, mapping=aes(x=condition,y=(L5overL2_3),fill=condition))+
  geom_jitter(AxonBranch_CTCF_p21, mapping=aes(x=condition,y=(L5overL2_3),fill="black"))+
  labs(x="Condition",y="Normalized cell fluorescence")+
  scale_fill_manual(values=c("black","#A2A2A1FF","#195190FF"))+
  theme_classic()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

write_xlsx(AxonBranch_CTCF_p21, file.path(sys_dir,"Figure_3D.xlsx"))

AxonBranch_CTCF_p21$brain <-  c("B1","B1","B1","B4","B4","B4","B6","B6","B6",
                                "B7","B7","B7","B8","B8","B8","B9","B9","B9")

ggplot()+
  geom_boxplot(AxonBranch_CTCF_p21, mapping=aes(x=condition,y=(layer_5),fill=condition))+
  geom_jitter(AxonBranch_CTCF_p21, mapping=aes(x=condition,y=(layer_5),shape=brain,size=1))+
  labs(x="Condition",y="The Corrected Total Cell Fluorescence")+
  scale_fill_manual(values=c("#A2A2A1FF","#195190FF"))+
  scale_shape_manual(values = c(16,8,25,16,8,25))+
  theme_classic()

AxonBranch_CTCF_p21_groupsample <- AxonBranch_CTCF_p21 %>% 
  group_by(brain) %>%
  summarise(across(c(layer_5, layer_2_3), mean),
            condition = first(condition)) %>% as.data.frame()

wilcox.test(layer_5 ~ condition, data = AxonBranch_CTCF_p21_groupsample)

ggplot()+
  geom_col(AxonBranch_CTCF_p21_groupsample,mapping=aes(x=brain,y=(layer_5),fill=condition))+
  # geom_text(,mapping=aes(x=,y=size,label=size), vjust = -0.5, color = "black")+
  scale_fill_manual(values=c("#A2A2A1FF","#195190FF"))+
  labs(x="Sample",y="The Corrected Total Cell Fluorescence (log10)")+
  theme_classic()

# Axon branching analysis p10 ---------------------------------------------



# Reporter assay - Promoters ----------------------------------------------

sys_dir <- "/Users/Shadi Shahatit/OneDrive/Desktop/"

# sample_list <- list.files(path = paste0(sys_dir, "Grik2 promoter Hum. vs Mac"), pattern = "\\.nd2$", full.names = TRUE)
sample_list <- list.files(path = paste0(sys_dir, "promoter_imageJ_output"), pattern = "\\.csv$", full.names = TRUE)
sample_names <- tools::file_path_sans_ext(basename(sample_list))

csv_file_list <- list.files(path = paste0(sys_dir, "promoter_imageJ_output"), pattern = "\\.csv$", full.names = TRUE)

all_promoter_CTCF_results <- list()

for (i in seq_along(csv_file_list)) {
  
  promoter_imageJ_out <- read.table(csv_file_list[i], header = TRUE, sep = ",")
  
  promoter_imageJ_out$channel <- rep(c("tdTomato", "EGFP"), length.out = nrow(promoter_imageJ_out))
  promoter_imageJ_out$sample <- rep(sample_names[i])
  promoter_imageJ_out$cell <- rep(1:nrow(promoter_imageJ_out))
  cell_id <- seq(nrow(promoter_imageJ_out) / 2)
  promoter_imageJ_out$cell <- rep(cell_id, each = 2)
  
  ## define your background --> last two entries
  mean_gray_value_seq <- seq(((nrow(promoter_imageJ_out))-1),nrow(promoter_imageJ_out))
  promoter_imageJ_out$CTCF <- (promoter_imageJ_out$IntDen) - (promoter_imageJ_out$Mean[mean_gray_value_seq] * promoter_imageJ_out$Area)
  promoter_imageJ_out_filtered <- promoter_imageJ_out[-mean_gray_value_seq, ] %>%
    select(c("sample", "cell", "channel", "CTCF"))
  row.names(promoter_imageJ_out_filtered) <- seq(1:nrow(promoter_imageJ_out_filtered))
  promoter_CTCF <- promoter_imageJ_out_filtered %>%
    pivot_wider(names_from = channel, values_from = CTCF)
  ## define your condition
  if (grepl("Human", sample_names[i], ignore.case = TRUE)) {
    promoter_CTCF$condition <- "human_promoter"
  } else if (grepl("Macaque", sample_names[i], ignore.case = TRUE)) {
    promoter_CTCF$condition <- "macaque_promoter"
  }
  ## define your brain
  promoter_CTCF$brain <- str_extract(sample_names[i], "_M[0-9]+_") %>% 
    str_extract("[0-9]+")
  
  promoter_CTCF$sig_int <- promoter_CTCF$EGFP / promoter_CTCF$tdTomato
  
  all_promoter_CTCF_results[[i]] <- promoter_CTCF
}

final_promoter_CTCF <- do.call(rbind, all_promoter_CTCF_results) %>% as.data.frame()

t.test(sig_int ~ condition, data = final_promoter_CTCF)
plot(density(log10(final_promoter_CTCF$sig_int)))

brain_count <- final_promoter_CTCF %>%
  group_by(condition) %>%
  summarise(n_brains = n_distinct(brain)) %>% 
  as.data.frame()
cell_count <- table(final_promoter_CTCF$condition)

wilcox_res <- wilcox.test(sig_int~condition, data=final_promoter_CTCF)
p_value <- wilcox_res$p.value
p_value_text <- ifelse(p_value<0.0001, "****", paste("p=", format(p_value, digits=3)))

boxplot_promoter <- ggplot()+
  geom_boxplot(final_promoter_CTCF, mapping=aes(x=condition,y=(sig_int),fill=condition))+
  geom_jitter(final_promoter_CTCF, mapping=aes(x=condition,y=(sig_int)),size=1)+
  labs(x="Condition",y="Normalized cell fluorescence")+
  scale_fill_manual(values=c("#E41A1C","#377EB8"))+
  theme_classic()+
  annotate("text", x=1, y=max(final_promoter_CTCF$sig_int), 
           label=paste0("B=",brain_count[1,2],",C=", cell_count[1]), vjust=-0.5) +
  annotate("text", x=2, y=max(final_promoter_CTCF$sig_int), 
           label=paste0("B=",brain_count[2,2],",C=", cell_count[2]), vjust=-0.5) +
  annotate("text", x=1.5, y=max(final_promoter_CTCF$sig_int)*1.1, 
           label=paste(p_value_text), vjust=-0.5, hjust=0.5)+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

barplot_promoter <- ggplot(final_promoter_CTCF, aes(x=condition, y=sig_int, fill=condition)) +
  stat_summary(fun=mean, geom="bar", position=position_dodge(), width = 0.7, color="black") +
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2, position=position_dodge(0.7)) +
  labs(x="Condition", y="Normalized cell fluorescence") +
  scale_fill_manual(values=c("#E41A1C","#377EB8"))+
  theme_classic()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

(boxplot_promoter + barplot_promoter) + 
  plot_layout(guides = "collect") &
  labs(x = "Condition", y = "The Corrected Total Cell Fluorescence") &
  theme(legend.position = "bottom")

write_xlsx(final_promoter_CTCF, file.path(sys_dir,"Figure_4F.xlsx"))

final_promoter_CTCF_fold <- final_promoter_CTCF %>%
  group_by(condition) %>%
  summarize(mean_sig_int = mean(sig_int), na.rm = TRUE)
human_mean <- final_promoter_CTCF_fold %>% filter(condition == "human_promoter") %>% pull(mean_sig_int)
macaque_mean <- final_promoter_CTCF_fold %>% filter(condition == "macaque_promoter") %>% pull(mean_sig_int)
fold_difference <- log2(human_mean/macaque_mean)
fold_difference

# Reporter assay - HARs ----------------------------------------------

sys_dir <- "/Users/Shadi Shahatit/OneDrive/Desktop/"
# sample_list <- list.files(path = paste0(sys_dir, "Grik2 HAR14+Ctrl"), pattern = "\\.nd2$", full.names = TRUE)
# sample_list <- list.files(path = paste0(sys_dir, "Grik2 HAR18+19+Ctrl"), pattern = "\\.nd2$", full.names = TRUE)
sample_list <- list.files(path = paste0(sys_dir, "HARs_imageJ_output"), pattern = "\\.csv$", full.names = TRUE)
sample_names <- tools::file_path_sans_ext(basename(sample_list))

csv_file_list <- list.files(path = paste0(sys_dir, "HARs_imageJ_output"), pattern = "\\.csv$", full.names = TRUE)

all_HARs_CTCF_results <- list()

for (i in seq_along(csv_file_list)) {
  
  HARs_imageJ_out <- read.table(csv_file_list[i], header = TRUE, sep = ",")
  
  HARs_imageJ_out$channel <- rep(c("tdTomato", "EGFP"), length.out = nrow(HARs_imageJ_out))
  HARs_imageJ_out$sample <- rep(sample_names[i])
  HARs_imageJ_out$cell <- rep(1:nrow(HARs_imageJ_out))
  cell_id <- seq(nrow(HARs_imageJ_out) / 2)
  HARs_imageJ_out$cell <- rep(cell_id, each = 2)
  
  ## define your background --> last two entries
  mean_gray_value_seq <- seq(((nrow(HARs_imageJ_out))-1),nrow(HARs_imageJ_out))
  HARs_imageJ_out$CTCF <- (HARs_imageJ_out$IntDen) - (HARs_imageJ_out$Mean[mean_gray_value_seq] * HARs_imageJ_out$Area)
  HARs_imageJ_out_filtered <- HARs_imageJ_out[-mean_gray_value_seq, ] %>%
    select(c("sample", "cell", "channel", "CTCF"))
  row.names(HARs_imageJ_out_filtered) <- seq(1:nrow(HARs_imageJ_out_filtered))
  HARs_CTCF <- HARs_imageJ_out_filtered %>%
    pivot_wider(names_from = channel, values_from = CTCF)
  ## define your condition
  if (grepl("HARCtrl", sample_names[i], ignore.case = TRUE)) {
    HARs_CTCF$condition <- "HAR_ctrl"
  } else if (grepl("HAR14", sample_names[i], ignore.case = TRUE)) {
    HARs_CTCF$condition <- "HAR_14"
  } else if (grepl("HAR18", sample_names[i], ignore.case = TRUE)) {
    HARs_CTCF$condition <- "HAR_18"
  } else if (grepl("HAR19", sample_names[i], ignore.case = TRUE)) {
    HARs_CTCF$condition <- "HAR_19"
  } else if (grepl("HAR15", sample_names[i], ignore.case = TRUE)) {
    HARs_CTCF$condition <- "HAR_15"
  } 
  ## define your brain
  HARs_CTCF$brain <- str_extract(sample_names[i], "_M[0-9]+_") %>% 
    str_extract("[0-9]+")
  ## define your postnatal day
  HARs_CTCF$p_day <- str_extract(sample_names[i], "_p[0-9]+_") %>% 
    str_extract("[0-9]+")
  
  HARs_CTCF$sig_int <- HARs_CTCF$EGFP / HARs_CTCF$tdTomato
  
  all_HARs_CTCF_results[[i]] <- HARs_CTCF
}

final_HARs_CTCF <- do.call(rbind, all_HARs_CTCF_results) %>% as.data.frame()

final_HARs_CTCF$condition_pday <- paste0(final_HARs_CTCF$condition,"_p",final_HARs_CTCF$p_day)
final_HARs_CTCF$condition_pday <- factor(final_HARs_CTCF$condition_pday,
                                         levels = c("HAR_ctrl_p14","HAR_ctrl_p11",
                                                    "HAR_14_p14","HAR_15_p11","HAR_18_p11","HAR_19_p11"))
final_HARs_CTCF$condition <- factor(final_HARs_CTCF$condition,
                                         levels = c("HAR_ctrl",
                                                    "HAR_14","HAR_15","HAR_18","HAR_19"))

## remove outliers - optional
final_HARs_CTCF <- final_HARs_CTCF %>% group_by(condition) %>%
  filter(abs(sig_int - mean(sig_int, na.rm=TRUE)) <= 3*sd(sig_int, na.rm=TRUE)) %>%
  ungroup()

brain_count <- final_HARs_CTCF %>%
  group_by(condition_pday) %>%
  summarise(n_brains = n_distinct(brain)) %>% 
  as.data.frame()
cell_count <- table(final_HARs_CTCF$condition_pday)

plot(density(log10(final_HARs_CTCF$sig_int)))

HAR_conditions <- unique(final_HARs_CTCF$condition[final_HARs_CTCF$condition != "HAR_ctrl"])
HAR_conditions <- c("HAR_14","HAR_15","HAR_18","HAR_19")

HARstats_list <- list()
HARp_value_df <- data.frame(condition = character(),
                            wilcox_p_value_sign = numeric(),
                            wilcox_rounded_p_value = numeric(),
                            t_test_p_value_sign = numeric(),
                            t_test_rounded_p_value = numeric(),
                            stringsAsFactors = FALSE)

for (cond in HAR_conditions) {
  p_day_value <- unique(final_HARs_CTCF$p_day[final_HARs_CTCF$condition == cond])
  ctrl_condition <- NA
  if (p_day_value == "14") {
    ctrl_condition <- "HAR_ctrl_p14"
  } else if (p_day_value == "11") {
    ctrl_condition <- "HAR_ctrl_p11"
  }
  subset_data <- final_HARs_CTCF %>% filter(condition == cond | condition_pday == ctrl_condition)
  if (length(unique(subset_data$condition)) == 2) {
    t_test_res <- t.test(sig_int ~ condition, data = subset_data)
    wilcox_res <- wilcox.test(sig_int ~ condition, data = subset_data)
    HARstats_list[[cond]] <- list(wilcox_test = wilcox_res, t_test = t_test_res)
    wilcox_p_plottext <- ifelse(wilcox_res$p.value > 0.05, "ns", 
                                ifelse(wilcox_res$p.value <= 0.05 & wilcox_res$p.value > 0.005, "**", 
                                       ifelse(wilcox_res$p.value <= 0.005 & wilcox_res$p.value > 0.0005, "***", 
                                              ifelse(wilcox_res$p.value <= 0.0005, "****", 
                                                     paste("p=", format(wilcox_res$p.value, digits = 3))))))
    t_test_p_plottext <- ifelse(t_test_res$p.value > 0.05, "ns", 
                                ifelse(t_test_res$p.value <= 0.05 & t_test_res$p.value > 0.005, "**", 
                                       ifelse(t_test_res$p.value <= 0.005 & t_test_res$p.value > 0.0005, "***", 
                                              ifelse(t_test_res$p.value <= 0.0005, "****", 
                                                     paste("p=", format(t_test_res$p.value, digits = 3))))))
    rounded_wilcox_p <- round(wilcox_res$p.value, 5)
    rounded_t_test_p <- round(t_test_res$p.value, 5)
    HARp_value_df <- rbind(HARp_value_df, data.frame(condition = cond,
                                                      wilcox_p_value_sign = wilcox_p_plottext,
                                                      wilcox_rounded_p_value = rounded_wilcox_p,
                                                      t_test_p_value_sign = t_test_p_plottext,
                                                      t_test_rounded_p_value = rounded_t_test_p))
  }
}

HAR14 <- final_HARs_CTCF %>% filter(condition=="HAR_14" | condition_pday=="HAR_ctrl_p14")
HAR15 <- final_HARs_CTCF %>% filter(condition=="HAR_15" | condition_pday=="HAR_ctrl_p11")
HAR18 <- final_HARs_CTCF %>% filter(condition=="HAR_18" | condition_pday=="HAR_ctrl_p11")
HAR19 <- final_HARs_CTCF %>% filter(condition=="HAR_19" | condition_pday=="HAR_ctrl_p11")
HARctrls <- final_HARs_CTCF %>% filter(condition_pday=="HAR_ctrl_p11" | condition_pday=="HAR_ctrl_p14")
HARctrl_p11 <- final_HARs_CTCF %>% filter(condition_pday=="HAR_ctrl_p11")
HARctrl_p14 <- final_HARs_CTCF %>% filter(condition_pday=="HAR_ctrl_p14")

wilcox.test(sig_int ~ condition, data = HAR14)
wilcox.test(sig_int ~ condition, data = HAR15)
wilcox.test(sig_int ~ condition, data = HAR18)
wilcox.test(sig_int ~ condition, data = HAR19)
wilcox.test(sig_int ~ condition_pday, data = HARctrls)

shapiro.test(HAR14$sig_int)
shapiro.test(HAR15$sig_int)
shapiro.test(HAR18$sig_int)
shapiro.test(HAR19$sig_int)
shapiro.test(HARctrl_p11$sig_int)
shapiro.test(HARctrl_p14$sig_int)

# custom_colors <- c("HAR_ctrl_p14" = "#FDE725FF","HAR_ctrl_p11" = "#FDE725FF",
#                    "HAR_14_p14" = "#481A6CFF","HAR_15_p11" = "#23888EFF",
#                    "HAR_18_p11" = "#39568CFF","HAR_19_p11" = "#35B779FF")
# custom_colors <- c("HAR_ctrl_p14" = "#000004FF","HAR_ctrl_p11" = "#000004FF",
#                    "HAR_14_p14" = "#5F187FFF","HAR_15_p11" = "#B63679FF",
#                    "HAR_18_p11" = "#EB5760FF","HAR_19_p11" = "#FCFDBFFF")
custom_colors <- c("HAR_ctrl_p14" = "#FCFDBFFF","HAR_ctrl_p11" = "#FCFDBFFF",
                   "HAR_14_p14" = "#F8765CFF","HAR_15_p11" = "#D3436EFF",
                   "HAR_18_p11" = "#982D80FF","HAR_19_p11" = "#5F187FFF")


boxplot_HARs <- ggplot()+
  geom_boxplot(final_HARs_CTCF, mapping=aes(x=condition_pday,y=(sig_int),fill=condition_pday),outlier.shape="|")+
  geom_jitter(final_HARs_CTCF, mapping=aes(x=condition_pday,y=(sig_int)),size=1)+
  labs(x="Conditions",y="Normalized cell fluorescence")+
  scale_fill_manual(values=custom_colors)+
  theme_classic()+
  annotate("text", x=1, y=max(final_HARs_CTCF$sig_int), 
           label=paste0("B=",brain_count$n_brains[1],",C=", cell_count[1]), vjust=-0.7) +
  annotate("text", x=2, y=max(final_HARs_CTCF$sig_int), 
           label=paste0("B=",brain_count$n_brains[2],",C=", cell_count[2]), vjust=-0.7) +
  annotate("text", x=3, y=max(final_HARs_CTCF$sig_int), 
           label=paste0("B=",brain_count$n_brains[3],",C=", cell_count[3]), vjust=-0.7) +
  annotate("text", x=4, y=max(final_HARs_CTCF$sig_int),
           label=paste0("B=",brain_count$n_brains[4],",C=", cell_count[4]), vjust=-0.7) +
  annotate("text", x=5, y=max(final_HARs_CTCF$sig_int),
           label=paste0("B=",brain_count$n_brains[5],",C=", cell_count[5]), vjust=-0.7) +
  annotate("text", x=6, y=max(final_HARs_CTCF$sig_int),
           label=paste0("B=",brain_count$n_brains[6],",C=", cell_count[6]), vjust=-0.7) +
  annotate("text", x=3, y=max(final_HARs_CTCF$sig_int)*1.1,
           label=paste(HARp_value_df$wilcox_p_value_sign[1]), vjust=-0.7, hjust=0.5)+
  annotate("text", x=4, y=max(final_HARs_CTCF$sig_int)*1.1,
           label=paste(HARp_value_df$wilcox_p_value_sign[2]), vjust=-0.7, hjust=0.5)+
  annotate("text", x=5, y=max(final_HARs_CTCF$sig_int)*1.1,
           label=paste(HARp_value_df$wilcox_p_value_sign[3]), vjust=-0.7, hjust=0.5)+
  annotate("text", x=6, y=max(final_HARs_CTCF$sig_int)*1.1,
           label=paste(HARp_value_df$wilcox_p_value_sign[4]), vjust=-0.7, hjust=0.5)
  # annotate("text", x=3, y=max(final_HARs_CTCF$sig_int)*1.1, 
  #          label=paste(HARp_value_df$wilcox_rounded_p_value[1]), vjust=-0.7, hjust=0.5)+
  # annotate("text", x=4, y=max(final_HARs_CTCF$sig_int)*1.1, 
  #          label=paste(HARp_value_df$wilcox_rounded_p_value[2]), vjust=-0.7, hjust=0.5)+
  # annotate("text", x=5, y=max(final_HARs_CTCF$sig_int)*1.1,
  #          label=paste(HARp_value_df$wilcox_rounded_p_value[3]), vjust=-0.7, hjust=0.5)+
  # annotate("text", x=6, y=max(final_HARs_CTCF$sig_int)*1.1,
  #          label=paste(HARp_value_df$wilcox_rounded_p_value[4]), vjust=-0.7, hjust=0.5)

barplot_HARs <- ggplot(final_HARs_CTCF, aes(x=condition_pday, y=sig_int, fill=condition_pday)) +
  stat_summary(fun=mean, geom="bar", position=position_dodge(), width = 0.7, color="black") +
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2, position=position_dodge(0.7)) +
  labs(x="Conditions",y="Normalized cell fluorescence")+
  scale_fill_manual(values=custom_colors)+
  # scale_fill_viridis(option="A",discrete=TRUE)+
  theme_classic()

(boxplot_HARs + barplot_HARs) + 
  plot_layout(guides = "collect") &
  labs(x="Conditions",y="Normalized Cell Fluorescence")+
  theme(legend.position = "bottom") &
  theme(legend.position = "none")

HARs_CTCF_mod <- final_HARs_CTCF %>%
  group_by(condition) %>%
  summarise(sig_int_avg = mean(sig_int)) %>%
  mutate(sig_int_mod = sig_int_avg-sig_int_avg[1]) %>%
  mutate(sig_int_percentage = (sig_int_mod/sig_int_avg[1])*100) %>%
  as.data.frame()

# custom_colors <- c("HAR_ctrl" = "#FDE725FF",
#                    "HAR_14" = "#481A6CFF","HAR_15" = "#23888EFF",
#                    "HAR_18" = "#39568CFF","HAR_19" = "#35B779FF")
custom_colors <- c("HAR_ctrl" = "#FCFDBFFF",
                   "HAR_14" = "#F8765CFF","HAR_15" = "#D3436EFF",
                   "HAR_18" = "#982D80FF","HAR_19" = "#5F187FFF")

ggplot(HARs_CTCF_mod, aes(y=condition, x=(sig_int_percentage), fill=condition)) +
  geom_bar(stat = "identity")+
  labs(y="Conditions",x="Signal Relative Change")+
  # labs(x="Conditions",y="delta signal intensity")+
  scale_fill_manual(values=custom_colors)+
  theme_classic()

## final figure

# custom_colors <- c("HAR_ctrl" = "#E64B35FF",
#                    "TF_sig_HAR14" = "#4DBBD5FF","TF_sig_HAR15" = "#00A087FF",
#                    "TF_sig_HAR16" = "#7E6148FF","TF_sig_HAR17" = "#91D1C2FF",
#                    "TF_sig_HAR18" = "#3C5488FF","TF_sig_HAR19" = "#F39B7FFF",
#                    "TF_sig_HAR20" = "#8491B4FF")

custom_colors <- c("HAR_ctrl" = "#FCFDBFFF",
                   "HAR_14" = "#FEB37BFF","HAR_15" = "#F4685CFF",
                   "HAR_18" = "#952C80FF","HAR_19" = "#6B1D81FF")

# dummy_legand <- as.data.frame("Negative Control" = "#FCFDBFFF",
#                               "HAR14" = "#FEB37BFF","HAR15" = "#F4685CFF",
#                               "HAR16" = "#D6456CFF","HAR17" = "#AB337CFF",
#                               "HAR18" = "#952C80FF","HAR19" = "#6B1D81FF", "HAR20"="#29115AFF")

brain_count <- final_HARs_CTCF %>%
  group_by(condition) %>%
  summarise(n_brains = n_distinct(brain)) %>% 
  as.data.frame()
cell_count <- table(final_HARs_CTCF$condition)

boxplot_HARs <- ggplot()+
  geom_boxplot(final_HARs_CTCF, mapping=aes(x=condition,y=(sig_int),fill=condition),outlier.shape="|")+
  labs(x="Conditions",y="Normalized Cell Fluorescence")+
  theme_classic()+
  annotate("text", x=1, y=max(final_HARs_CTCF$sig_int), 
           label=paste0("B=",brain_count$n_brains[1],",C=", cell_count[1]), vjust=-0.7) +
  annotate("text", x=2, y=max(final_HARs_CTCF$sig_int), 
           label=paste0("B=",brain_count$n_brains[2],",C=", cell_count[2]), vjust=-0.7) +
  annotate("text", x=3, y=max(final_HARs_CTCF$sig_int), 
           label=paste0("B=",brain_count$n_brains[3],",C=", cell_count[3]), vjust=-0.7) +
  annotate("text", x=4, y=max(final_HARs_CTCF$sig_int),
           label=paste0("B=",brain_count$n_brains[4],",C=", cell_count[4]), vjust=-0.7) +
  annotate("text", x=5, y=max(final_HARs_CTCF$sig_int),
           label=paste0("B=",brain_count$n_brains[5],",C=", cell_count[5]), vjust=-0.7) +
  annotate("text", x=2, y=max(final_HARs_CTCF$sig_int)*1.1,
           label=paste(HARp_value_df$wilcox_p_value_sign[1]), vjust=-0.7, hjust=0.5)+
  annotate("text", x=3, y=max(final_HARs_CTCF$sig_int)*1.1,
           label=paste(HARp_value_df$wilcox_p_value_sign[2]), vjust=-0.7, hjust=0.5)+
  annotate("text", x=4, y=max(final_HARs_CTCF$sig_int)*1.1,
           label=paste(HARp_value_df$wilcox_p_value_sign[3]), vjust=-0.7, hjust=0.5)+
  annotate("text", x=5, y=max(final_HARs_CTCF$sig_int)*1.1,
           label=paste(HARp_value_df$wilcox_p_value_sign[4]), vjust=-0.7, hjust=0.5)+
  scale_fill_npg()+
  scale_fill_manual(values=custom_colors)+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

barplot_HARs <- ggplot(final_HARs_CTCF, aes(x=condition, y=sig_int, fill=condition)) +
  stat_summary(fun=mean, geom="bar", position=position_dodge(), width = 0.7, color="black") +
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2, position=position_dodge(0.7)) +
  labs(x="Conditions",y="Normalized cell fluorescence")+
  theme_classic()+
  scale_fill_npg()+
  scale_fill_manual(values=custom_colors)+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

(boxplot_HARs + barplot_HARs) + 
  plot_layout(guides = "collect") &
  labs(x="Conditions",y="Normalized cell fluorescence")+
  theme(legend.position = "bottom") &
  theme(legend.position = "none")

ggsave(file.path(sys_dir,"fig4_E.png"), plot=ggplot2::last_plot())

write_xlsx(final_HARs_CTCF, file.path(sys_dir,"Figure_4E.xlsx"))

HARs_CTCF_mod <- final_HARs_CTCF %>%
  group_by(condition) %>%
  summarise(sig_int_avg = mean(sig_int)) %>%
  mutate(sig_int_mod = sig_int_avg-sig_int_avg[1]) %>%
  mutate(sig_int_percentage = (sig_int_mod/sig_int_avg[1])*100) %>%
  as.data.frame()

ggplot(HARs_CTCF_mod, aes(y=condition, x=(sig_int_percentage), fill=condition)) +
  geom_bar(stat = "identity")+
  labs(y="Conditions",x="Signal Relative Change")+
  # labs(x="Conditions",y="delta signal intensity")+
  theme_classic()+
  scale_fill_npg()+
  scale_fill_manual(values=custom_colors)+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

ggsave(file.path(sys_dir,"fig4_z.png"), plot=ggplot2::last_plot())

HARctrl_sig_int_avg <- HARs_CTCF_mod$sig_int_avg[1]

HARs_CTCF_modified <- final_HARs_CTCF %>%
  group_by(condition) %>%
  mutate(sig_int_modified = sig_int-HARctrl_sig_int_avg) %>%
  mutate(sig_int_percentage = (sig_int_modified/HARctrl_sig_int_avg)*100)

  #   %>%
  # group_by(condition) %>%
  # summarise(sig_int_percentage_avg = mean(sig_int_percentage)) %>%
  # as.data.frame()

# ggplot(HARs_CTCF_modified, aes(y=condition, x=(sig_int_percentage_avg), fill=condition)) +
#   geom_bar(stat = "identity")+
#   labs(y="Conditions",x="Signal Relative Change")+
#   theme_classic()+
#   scale_fill_manual(values=custom_colors)+
#   theme(
#     plot.title = element_text(size = 16, face = "bold"),
#     axis.title.x = element_text(size = 14, face = "bold"),
#     axis.title.y = element_text(size = 14, face = "bold"))+
#   theme_classic()+
#   scale_fill_npg()+
#   scale_fill_manual(values=custom_colors)

ggplot(HARs_CTCF_modified, aes(x=condition, y=sig_int_percentage, fill=condition)) +
  stat_summary(fun=mean, geom="bar", position=position_dodge(), width = 0.7, color="black") +
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2, position=position_dodge(0.7)) +
  labs(x="Conditions",y="Signal Relative Change")+
  theme_classic()+
  annotate("text", x=1, y=max(HARs_CTCF_modified$sig_int_percentage)/3, 
           label=paste0("n= ", cell_count[1]), vjust=-0.7) +
  annotate("text", x=2, y=max(HARs_CTCF_modified$sig_int_percentage)/3, 
           label=paste0("n= ", cell_count[2]), vjust=-0.7) +
  annotate("text", x=3, y=max(HARs_CTCF_modified$sig_int_percentage)/3, 
           label=paste0("n= ", cell_count[3]), vjust=-0.7) +
  annotate("text", x=4, y=max(HARs_CTCF_modified$sig_int_percentage)/3,
           label=paste0("n= ", cell_count[4]), vjust=-0.7) +
  annotate("text", x=5, y=max(HARs_CTCF_modified$sig_int_percentage)/3,
           label=paste0("n= ", cell_count[5]), vjust=-0.7) +
  annotate("text", x=2, y=max(HARs_CTCF_modified$sig_int_percentage)/3*1.1,
           label=paste(HARp_value_df$wilcox_p_value_sign[1]), vjust=-0.7, hjust=0.5)+
  annotate("text", x=3, y=max(HARs_CTCF_modified$sig_int_percentage)/3*1.1,
           label=paste(HARp_value_df$wilcox_p_value_sign[2]), vjust=-0.7, hjust=0.5)+
  annotate("text", x=4, y=max(HARs_CTCF_modified$sig_int_percentage)/3*1.1,
           label=paste(HARp_value_df$wilcox_p_value_sign[3]), vjust=-0.7, hjust=0.5)+
  annotate("text", x=5, y=max(HARs_CTCF_modified$sig_int_percentage)/3*1.1,
           label=paste(HARp_value_df$wilcox_p_value_sign[4]), vjust=-0.7, hjust=0.5)+
  scale_fill_manual(values=custom_colors)+
  theme(
    plot.title = element_text(size = 11),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12))

ggsave(file.path(sys_dir,"fig4_zz.png"), plot=ggplot2::last_plot())

write_xlsx(HARs_CTCF_modified, file.path(sys_dir,"Figure_4D.xlsx"))

cell_count <- as.data.frame(cell_count)

# Supp_tabels -------------------------------------------------------------

Rdataframes_titles_image_analysis <- data.frame(
  
  Sheet = c("Sheet1",
            "Sheet2",
            "Sheet3",
            "Sheet4",
            "Sheet5",
            "Sheet6"
  ),
  
  DataFrameName = c("AxonBranch_CTCF_p21",
                    "final_promoter_CTCF",
                    "final_HARs_CTCF",
                    "HARs_CTCF_modified",
                    "brain_count" ,
                    "cell_count"
  ),
  
  DataFrameNote = c("image analysis - axon branching CTCF at p21 for shCTRL and shGRIK2",
                    "image analysis - grik2 promoter CTCF at p11 for human and macaque",
                    "image analysis - grik2HARs CTCF at p11-p14", 
                    "image analysis - grik2HARs CTCF at p11-p14 modified with combined ctrl scores and relative change",
                    "brain count for modifed grik2HARs data",
                    "cell count for modifed grik2HARs data"
  ))

write_xlsx(list(
  ContentTable = Rdataframes_titles_image_analysis,
  
  Sheet1 = AxonBranch_CTCF_p21,
  Sheet2 = final_promoter_CTCF,
  Sheet3 = final_HARs_CTCF,
  Sheet4 = HARs_CTCF_modified,
  Sheet5 = brain_count,
  Sheet6 = cell_count
  ),
file.path(sys_dir,"supp_data_Rdataframes_image_analysis.xlsx"))


