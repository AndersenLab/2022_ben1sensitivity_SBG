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

#Figure1

Fig1A <- cowplot::ggdraw()+ cowplot::draw_image("20220518_Fig1.png")
Fig1 <- cowplot::plot_grid(Fig1A, nrow = 1)
ggplot2::ggsave("fig1.jpeg", plot=Fig1, width = 3.75, height = 5)


#Figure2

Fig2A <- cowplot::ggdraw()+ cowplot::draw_image("20220518_Figure2A.png")
Fig2B <- cowplot::ggdraw()+ cowplot::draw_image("20220518_Fig2B.png")

norm_green <- read_csv("norm_green.csv")

density_gfp <- density(norm_green$norm.green)
cutoff <- optimize(approxfun(density_gfp$x,density_gfp$y),interval=c(0,2))$minimum

GFP_dist <-norm_green %>%
  dplyr::mutate(GFP = ifelse(norm_green$norm.green >= cutoff,"GFP","nonGFP"))%>%
  ggplot()+
  aes(x=norm.green)+
  geom_density(aes(fill=GFP))+
  scale_fill_manual(values = c("GFP"="green", "nonGFP"="grey"))+
  geom_vline(xintercept = cutoff)+
  xlab("Normalized GFP expression")+
  ylab("")+
  ylim(0,8)+
  geom_text(x=0.15,y=6,label="GFP-", size=4)+
  geom_text(aes(x=1.75,y=2),label="GFP+", size=4)+
  cowplot::theme_cowplot(12)+
  theme(legend.position = "none")

fig2 <- cowplot::plot_grid(Fig2A,ncol = 2, labels = c("A","B"), cowplot::plot_grid(Fig2B,GFP_dist,nrow = 2,labels = c("B","C")))

ggplot2::ggsave("fig2.jpeg", width = 7.5, height = 8, units = "in")


#Figure 3

HTA_data <- read_delim("HTA_data.tsv", delim = "\t", 
                       escape_double = FALSE, trim_ws = TRUE)

tissue_rescue <- HTA_data %>%
  dplyr::filter(date == 20210824)

biopruned_tissue <- bioprune(tissue_rescue)

outlierpruned_tissue <- prune_outliers(biopruned_tissue) %>%
  tidyr::separate(strain, c("strain", "GFP"), sep = "_") %>%
  dplyr::filter(strain != "WASH") %>%
  dplyr::filter(GFP == "+") #filter out wash wells and non-GFP worms

#3A Control strains

control <- outlierpruned_tissue %>%
  dplyr::filter(strain == "ECA2842"| 
                  strain == "ECA2821")

regressed_control <- regress(control)

mean_TOF_control <- regressed_control %>%
  dplyr::filter(trait == "mean.TOF")

control_stats <-  mean_TOF_control %>%
  aov(phenotype ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1 == "ECA2821" & group2 == "ECA2842")


cols <- c("ECA2821" = "orange",
          "ECA2842" = "grey",
          "ECA2824" = "#0072B2",
          "ECA2828" = "#009E73",
          "ECA2833" = "#990000",
          "ECA2836" = "#F0E442",
          "ECA2846" = "#663399",
          "ECA2840" = "#CC79A7")


Fig3A <- mean_TOF_control%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("ECA2821",
                                                       "ECA2842")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  scale_x_discrete(breaks = c("ECA2821",
                              "ECA2842"),
                   labels = c("Wild-type",
                              "Deletion"))+
  geom_jitter(width = 0.1,size =1)+
  geom_boxplot(aes(fill=strain, alpha = 0.4),outlier.shape = NA)+
  stat_pvalue_manual(control_stats, label = "p.adj.signif",x="group1", y.position = c(105),remove.bracket = TRUE,size = 8)+
  scale_fill_manual(name = "fancy_strain",values = cols)+
  ylab("Normalized Animal Length")+
  cowplot::theme_cowplot(12)+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        legend.position = "none")

tissue_results <- mean_TOF_tissue %>%
  dplyr::filter(strain != "ECA2821")

stats_tissue <- plotstats_tissue %>%
  dplyr::filter(group1 != "ECA2821")




#Fig3B

no_n2 <- outlierpruned_tissue %>%
  dplyr::filter(strain != "ECA2821")

regressed <- regress(no_n2)

mean_TOF_tissue <- regressed %>%
  dplyr::filter(trait == "mean.TOF")


plotstats_tissue <- mean_TOF_tissue %>%
  aov(phenotype ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1 == "ECA2824" & group2 == "ECA2842" |
                  group1 == "ECA2828" & group2 == "ECA2842" |
                  group1 == "ECA2833" & group2 == "ECA2842" |
                  group1 == "ECA2836" & group2 == "ECA2842" |
                  group1 == "ECA2840" & group2 == "ECA2842")

plotstats2 <- mean_TOF_tissue %>%
  aov(phenotype ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1 == "ECA2842" & group2 == "ECA2846") %>%
  dplyr::rename("group1" = "group2", "group2" = "group1")

plotstats_tissue <- rbind(plotstats_tissue, plotstats2)


Fig3B <- mean_TOF_tissue%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("ECA2842",
                                                       "ECA2824",
                                                       "ECA2828",
                                                       "ECA2833",
                                                       "ECA2836",
                                                       "ECA2846",
                                                       "ECA2840")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  scale_x_discrete(breaks = c("ECA2842",
                              "ECA2824",
                              "ECA2828",
                              "ECA2833",
                              "ECA2836",
                              "ECA2846",
                              "ECA2840"),
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

fig3 <- cowplot::plot_grid(Fig3A, Fig3B, nrow = 2, rel_heights = c(2,3), labels = c("A","B"))

ggplot2::ggsave("fig3.jpeg", fig3, width = 7.5, height = 7, units = "in" )



#DMSO data

DMSO <- outlierpruned_tissue %>%
  dplyr::filter(condition == "DMSO") %>%
  dplyr::filter(trait == "mean.TOF")


DMSO %>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("ECA2821",
                                                        "ECA2842",
                                                       "ECA2824",
                                                       "ECA2828",
                                                       "ECA2833",
                                                       "ECA2836",
                                                       "ECA2846",
                                                       "ECA2840")))%>%
  ggplot() +
  aes(x=fancy_strain, y=phenotype) +
  geom_boxplot()



#Figure 4

neurons <- HTA_data %>%
  dplyr::filter(date == 20220426)

biopruned_neurons <- bioprune(neurons)

outlierpruned_neurons <- prune_outliers(biopruned_neurons) %>%
  tidyr::separate(strain, c("strain", "GFP"), sep = "_") %>%
  dplyr::filter(strain != "WASH") %>%
  dplyr::filter(GFP == "+")  %>%
  dplyr::filter(strain != "ECA2821") #filter out wash wells and non-GFP worms

regressed_neurons <- regress(outlierpruned_neurons)

mean_TOF_neurons <- regressed_neurons %>%
  dplyr::filter(trait == "mean.TOF")


plotstats_neurons <- mean_TOF_neurons %>%
  aov(phenotype ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1 == "ECA2824" & group2 == "ECA2842" |
                  group1 == "ECA2833" & group2 == "ECA2842")

plotstats3 <- mean_TOF_neurons %>%
  aov(phenotype ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1 == "ECA2842" & group2 == "ECA3333" |
                  group1 == "ECA2842" & group2 == "ECA3335" |
                  group1 == "ECA2842" & group2 == "ECA3354" |
                  group1 == "ECA2842" & group2 == "ECA3356") %>%
  dplyr::rename("group1" = "group2", "group2" = "group1")

plotstats_neurons <- rbind(plotstats_neurons, plotstats3)


cols2 <- c("ECA2842" = "grey",
          "ECA2824" = "#0072B2",
          "ECA2833" = "#990000",
          "ECA3333" = "#666600",
          "ECA3335" = "#FFCCCC",
          "ECA3354" = "#CCCCFF",
          "ECA3356" = "#663300")


fig4 <- mean_TOF_neurons%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("ECA2842",
                                                       "ECA2824",
                                                       "ECA2833",
                                                       "ECA3335",
                                                       "ECA3356",
                                                       "ECA3333",
                                                       "ECA3354")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  scale_x_discrete(breaks = c("ECA2842",
                              "ECA2824",
                              "ECA2833",
                              "ECA3335",
                              "ECA3356",
                              "ECA3333",
                              "ECA3354"),
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

ggplot2::ggsave("fig4.jpeg", fig4, width = 3.75, height = 4, units = "in")


#DMSO

all_strains_pruned <- prune_outliers(biopruned_neurons) %>%
  tidyr::separate(strain, c("strain", "GFP"), sep = "_") %>%
  dplyr::filter(strain != "WASH") %>%
  dplyr::filter(GFP == "+")
  


DMSO_neurons <- all_strains_pruned %>%
  dplyr::filter(condition == "DMSO") %>%
  dplyr::filter(trait == "mean.TOF")


DMSO_neurons %>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("ECA2821",
                                                       "ECA2842",
                                                       "ECA2824",
                                                       "ECA2833",
                                                       "ECA3335",
                                                       "ECA3356",
                                                       "ECA3333",
                                                       "ECA3354")))%>%
  ggplot() +
  aes(x=fancy_strain, y=phenotype) +
  geom_boxplot()












