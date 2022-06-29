#Libraries
library(ggplot2)
library(tidyverse)
library(cowplot)
library(glue)
library(magick)
library(rstatix)
library(ggpubr)
library(COPASutils)
library(easysorter)
library(RColorBrewer)


#FigureS1

HTA_data <- readr::read_delim("HTA_data.tsv", delim = "\t", 
                       escape_double = FALSE, trim_ws = TRUE)

tissue_rescue <- HTA_data %>%
  dplyr::filter(date == 20210728)

biopruned_tissue <- easysorter::bioprune(tissue_rescue)

outlierpruned_tissue <- easysorter::prune_outliers(biopruned_tissue) %>%
  tidyr::separate(strain, c("strain", "GFP"), sep = "_") %>%
  dplyr::filter(strain != "WASH") %>%
  dplyr::filter(GFP == "+") %>%
  dplyr::filter(strain != "ECA2822") #filter out wash wells and non-GFP worms

regressed_tissue <- easysorter::regress(outlierpruned_tissue)

mean_TOF_tissue <- regressed_tissue %>%
  dplyr::filter(trait == "mean.TOF")


plotstats_tissue <- mean_TOF_tissue %>%
  aov(phenotype ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group2 == "ECA2841")

plotstats2 <- mean_TOF_tissue %>%
  aov(phenotype ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1 == "ECA2841" & group2 == "ECA2844") %>%
  dplyr::rename("group1" = "group2", "group2" = "group1")

plotstats_tissue <- rbind(plotstats_tissue, plotstats2)


cols <- c("ECA2841" = "grey",
          "ECA2823" = "#0072B2",
          "ECA2827" = "#009E73",
          "ECA2832" = "#990000",
          "ECA2835" = "#F0E442",
          "ECA2830" = "#663399",
          "ECA2844" = "#CC79A7")


figS1 <- mean_TOF_tissue%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("ECA2841",
                                                       "ECA2823",
                                                       "ECA2827",
                                                       "ECA2832",
                                                       "ECA2835",
                                                       "ECA2830",
                                                       "ECA2844")))%>%
  ggplot2::ggplot()+
    aes(x=fancy_strain, y = phenotype)+
    scale_x_discrete(breaks = c("ECA2841",
                              "ECA2823",
                              "ECA2827",
                              "ECA2832",
                              "ECA2835",
                              "ECA2830",
                              "ECA2844"),
                   labels = c("Deletion",
                              "Endogenous",
                              "Ubiquitous",
                              "Neurons",
                              "Hypodermis",
                              "Muscles",
                              "Intestine"))+
    geom_jitter(width = 0.1,size =1)+
    geom_boxplot(aes(fill=strain, alpha = 0.4),outlier.shape = NA)+
    stat_pvalue_manual(plotstats_tissue, label = "p.adj.signif",x="group1", y.position = c(125),remove.bracket = TRUE,size = 8)+
    scale_fill_manual(name = "fancy_strain",values = cols)+
    ylab("Normalized Animal Length")+
    cowplot::theme_cowplot(12)+
    theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        legend.position = "none")

ggplot2::ggsave("figS1.jpeg", figS1, width = 7.5, height = 4, units = "in" )




#FigureS2

neurons <- HTA_data %>%
  dplyr::filter(date == 20220418)

biopruned_neurons <- easysorter::bioprune(neurons)

outlierpruned_neurons <- easysorter::prune_outliers(biopruned_neurons) %>%
  tidyr::separate(strain, c("strain", "GFP"), sep = "_") %>%
  dplyr::filter(strain != "WASH") %>%
  dplyr::filter(GFP == "+") %>%
  dplyr::filter(strain != "ECA2822")#filter out wash wells and non-GFP worms

regressed_neurons <- easysorter::regress(outlierpruned_neurons)

mean_TOF_neurons <- regressed_neurons %>%
  dplyr::filter(trait == "mean.TOF")


plotstats_neurons <- mean_TOF_neurons %>%
  aov(phenotype ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group2 == "ECA2841")

plotstats3 <- mean_TOF_neurons %>%
  aov(phenotype ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1 == "ECA2841" & group2 == "ECA3332" |
                  group1 == "ECA2841" & group2 == "ECA3334" |
                  group1 == "ECA2841" & group2 == "ECA3336" |
                  group1 == "ECA2841" & group2 == "ECA3355") %>%
  dplyr::rename("group1" = "group2", "group2" = "group1")

plotstats_neurons <- rbind(plotstats_neurons, plotstats3)


cols2 <- c("ECA2841" = "grey",
           "ECA2823" = "#0072B2",
           "ECA2832" = "#990000",
           "ECA3332" = "#666600",
           "ECA3334" = "#FFCCCC",
           "ECA3336" = "#CCCCFF",
           "ECA3355" = "#663300")


figS2 <- mean_TOF_neurons%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("ECA2841",
                                                       "ECA2823",
                                                       "ECA2832",
                                                       "ECA3334",
                                                       "ECA3355",
                                                       "ECA3332",
                                                       "ECA3336")))%>%
  ggplot2::ggplot()+
    aes(x=fancy_strain, y = phenotype)+
    scale_x_discrete(breaks = c("ECA2841",
                              "ECA2823",
                              "ECA2832",
                              "ECA3334",
                              "ECA3355",
                              "ECA3332",
                              "ECA3336"),
                   labels = c("Deletion",
                              "Endogenous",
                              "All neurons",
                              "Cholinergic",
                              "GABAergic",
                              "Glutamatergic",
                              "Dopaminergic"
                   ))+
    geom_jitter(width = 0.1,size =.5)+
    geom_boxplot(aes(fill=strain, alpha = 0.4),outlier.shape = NA)+
    stat_pvalue_manual(plotstats_neurons, label = "p.adj.signif",x="group1", y.position = c(95),remove.bracket = TRUE,size = 4)+
    scale_fill_manual(name = "fancy_strain",values = cols2)+
    ylab("Normalized Animal Length")+
    xlab("GFP (+/-)") +
    cowplot::theme_cowplot(12)+
    theme(axis.text.x = element_text(size=12, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        legend.position = "none")

ggplot2::ggsave("figS2.jpeg", figS2, width = 3.75, height = 4, units = "in")


#FigureS3

CeNGEN_ben_1_expr <- readr::read_csv("CeNGEN_ben-1_expr.csv")

ben1_neurons <- CeNGEN_ben_1_expr %>%
  dplyr::filter(per_cell != 0)

system <- c("Glu" = "#666600",
            "ACh" = "#FFCCCC",
            "DA" = "#CCCCFF",
            "GABA" = "#663300")

figS3 <- ggplot2::ggplot(ben1_neurons, aes(x=Neurotransmitter, y=TPM, fill=Neurotransmitter)) + 
  scale_fill_manual(values = system) +
  geom_jitter(width = 0.1, size=0.8) +
  ylab("Transcripts Per Million (TPM)")+
  xlab("Neurotransmitter") +
  scale_x_discrete(labels = c("Acetylcholine",
             "Dopamine",
             "GABA",
             "Glutamate",
             "NA")) +
  cowplot::theme_cowplot(12) +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        legend.position = "none")

ggplot2::ggsave("figS3.jpeg", figS3, width = 7, height = 4, units = "in")

  




