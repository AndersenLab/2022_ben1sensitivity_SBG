#Libraries
library(COPASutils)
library(easysorter)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(glue)

#Combine datasets from each assay to 1 file for analysis and figures for paper

#Assay1A
dir1A <- c("/projects/b1059/projects/Sophie/ben1_rescue/20210728_ben1rescue1a")

raw1 <- easysorter::read_data(dir1A)

raw_nocontam1 <- easysorter::remove_contamination(raw1)

#Configure +/- GFP
raw_score1 <- raw_nocontam1[[1]]

density_gfp1 <- density(raw_score1$norm.green)
cutoff1 <- optimize(approxfun(density_gfp1$x,density_gfp1$y),interval=c(0,2))$minimum

gfp_score1 <- raw_score1 %>%
  dplyr::mutate("GFP" = ifelse(norm.green >= cutoff1, "+", "null")) %>%
  tidyr::unite("strain", strain, GFP, sep="_")

raw_nocontam_1 <- list(gfp_score1,raw_nocontam1[[2]])

summedraw1 <- easysorter::sumplate(raw_nocontam_1, directories = FALSE, quantiles = TRUE)


#Assay1B

dir2 <- c("/projects/b1059/projects/Sophie/ben1_rescue/20210824_ben1rescue2")


raw2 <- easysorter::read_data(dir2)

raw_nocontam2 <- easysorter::remove_contamination(raw2)

#Configure +/- GFP
raw_score2 <- raw_nocontam2[[1]]

norm_green <- raw_score2 %>%
  dplyr::select(date, experiment, round, assay, plate, row, col, sort, norm.green)

write_csv(norm_green, "norm_green.csv")

density2_gfp <- density(raw_score2$norm.green)
cutoff2 <- optimize(approxfun(density2_gfp$x,density2_gfp$y),interval=c(0,2))$minimum

gfp_score2 <- raw_score2 %>%
  dplyr::mutate("GFP" = ifelse(norm.green >= cutoff2, "+", "null")) %>%
  tidyr::unite("strain", strain, GFP, sep="_")

raw_nocontam_2 <- list(gfp_score2,raw_nocontam2[[2]])

summedraw2 <- easysorter::sumplate(raw_nocontam_2, directories = FALSE, quantiles = TRUE)

#Assay2A

dir3 <- c("/projects/b1059/projects/Sophie/ben1_rescue/20220418_ben1neuronsA")

raw3 <- easysorter::read_data(dir3)

raw_nocontam3 <- easysorter::remove_contamination(raw3)

#Configure +/- GFP
raw_score3 <- raw_nocontam3[[1]]

density_gfp3 <- density(raw_score3$norm.green)
cutoff3 <- optimize(approxfun(density_gfp3$x,density_gfp3$y),interval=c(0,2))$minimum

gfp_score3 <- raw_score3 %>%
  dplyr::mutate("GFP" = ifelse(norm.green >= cutoff3, "+", "null")) %>%
  tidyr::unite("strain", strain, GFP, sep="_")

raw_nocontam_3 <- list(gfp_score3,raw_nocontam3[[2]])

summedraw3 <- easysorter::sumplate(raw_nocontam_3, directories = FALSE, quantiles = TRUE)





#Assay2B

dir4 <- c("/projects/b1059/projects/Sophie/ben1_rescue/20220426_ben1neuronsB")

raw4 <- easysorter::read_data(dir4)

raw_nocontam4 <- easysorter::remove_contamination(raw4)

#Configure +/- GFP
raw_score4 <- raw_nocontam4[[1]]

density_gfp4 <- density(raw_score4$norm.green)
cutoff4 <- optimize(approxfun(density_gfp4$x,density_gfp4$y),interval=c(0,2))$minimum

gfp_score4 <- raw_score4 %>%
  dplyr::mutate("GFP" = ifelse(norm.green >= cutoff4, "+", "null")) %>%
  tidyr::unite("strain", strain, GFP, sep="_")

raw_nocontam_4 <- list(gfp_score4,raw_nocontam4[[2]])

summedraw4 <- easysorter::sumplate(raw_nocontam_4, directories = FALSE, quantiles = TRUE)



#Combine tables


HTA_data <- rbind(summedraw1, summedraw2, summedraw3, summedraw4)

readr::write_tsv(HTA_data, "/projects/b1059/projects/Sophie/ben1_rescue/HTA_data.tsv")

#Download and run on local RStudio to analyze and make figures

