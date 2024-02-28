# Program: covail_plot.R
# Author : Cindy Molitor, Bioinformatics
# Purpose: Create a line plot of the COVAIL Force Of Infection (FOI) score
#          over the date range between the first Stage 1 enrollment date and
#          last D181 + 7 date.

library(patchwork)
library(readr)
library(tidyverse)

setwd("/Users/cindy/projects/covail_plot/code")

# Read Youyi's COVAIL file to generate an enrollment stage column based on
# the studydose1date and the ranges given in the analysis plan. Note that
# for stages 2 and 3, the dates for some participants were outside the
# range given in the analysis plan and those ranges were expanded per
# Peter Gilbert.
enr_df <- read_csv("../input_files/covail_data_processed_20240103.csv")

enr_stage_df <- enr_df %>% mutate(enr_stage = case_when(
                       studydose1date >= as.Date("2022-03-30") & studydose1date <= as.Date("2022-05-06") ~ "Stage1",
                       studydose1date >= as.Date("2022-05-09") & studydose1date <= as.Date("2022-05-27") ~ "Stage2",
                       studydose1date >= as.Date("2022-06-06") & studydose1date <= as.Date("2022-06-17") ~ "Stage3",
                       studydose1date >= as.Date("2022-10-04") & studydose1date <= as.Date("2022-10-28") ~ "Stage4"
                       )) %>%
                       select(Ptid, enr_stage, studydose1date) %>%
                       rename(ptid = Ptid)

# Craig Magaret gave me the following file that contains a subset of
# columns from Youyi's file.
foi_df <- read_csv("../input_files/covail_foi_v2.csv")

foi_enr_stage_df <- left_join(foi_df, enr_stage_df, by = join_by(ptid))

# Expand the FOI dataset such that a new observation with the FOI score
# is generated for each date between the studydose1date and D181 + 7.
foi_expanded_df = foi_enr_stage_df %>%
                  filter(!is.na(foi)) %>%
                  group_by(r=row_number()) %>%
                  mutate(stage_date = list(studydose1date:date.d188)) %>%
                  ungroup %>% select(-r) %>%
                  unnest(cols=c(stage_date))

foi_expanded_df$stage_date <- as.Date(foi_expanded_df$stage_date)

# Get the last (max) D181 + 7 date for each stage.
foi_max_df <- foi_expanded_df %>%
                  group_by(enr_stage) %>%
                  summarise(max_date = max(stage_date))

foi_max_vector <- pull(foi_max_df, max_date)

foi_stats_df <- foi_expanded_df %>% group_by(stage_date) %>%
                                    summarise(foi_median = median(foi),
                                              foi_q25 = quantile(foi, probs=0.25),
                                              foi_q75 = quantile(foi, probs=0.75),
                                              foi_min = min(foi),
                                              foi_max = max(foi),
                                              foi_n = n())

png(filename="../output_files/COVAIL_FOI_combined.png", width=700, height=700)

plot1 <- foi_stats_df %>% ggplot(aes(stage_date, foi_median)) +
                 geom_ribbon(aes(ymin = foi_q25,
                                 ymax = foi_q75),
                             fill = "steelblue2") +
                 geom_line(color = "firebrick",
                           linewidth = 1) +
                 ylim(0, 16) +
                 geom_vline(aes(xintercept=as.Date("2022-03-30"), color="Stage 1 Start/End"), linetype="solid") +
                 geom_vline(xintercept=as.Date(foi_max_vector[1]), linetype="solid", color="royalblue") +
                 geom_vline(aes(xintercept=as.Date("2022-05-09"), color="Stage 2 Start/End"), linetype="solid") +
                 geom_vline(xintercept=as.Date(foi_max_vector[2]), linetype="solid", color="cyan2") +
                 geom_vline(aes(xintercept=as.Date("2022-06-06"), color="Stage 3 Start/End"), linetype="solid") +
                 geom_vline(xintercept=as.Date(foi_max_vector[3]), linetype="solid", color="deeppink") +
                 geom_vline(aes(xintercept=as.Date("2022-10-04"), color="Stage 4 Start/End"), linetype="solid") +
                 geom_vline(xintercept=as.Date(foi_max_vector[4]), linetype="solid", color="green4") +
                 labs(title = "COVAIL FOI By Date",
                      subtitle = "Shaded area represents 25th and 75th quantile FOI",
                      x = "Date",
                      y = "Median FOI") +
                 scale_color_manual(name="Enrollment Stages", values=c("Stage 1 Start/End" = "royalblue",
                                                                       "Stage 2 Start/End" = "cyan2",
                                                                       "Stage 3 Start/End" = "deeppink",
                                                                       "Stage 4 Start/End" = "green4")) +
                  theme(legend.position = "bottom")

plot2 <- foi_stats_df %>% ggplot(aes(stage_date, foi_n)) +
                 geom_line(color = "navy",
                           linewidth = 1) +
                 labs(title = "Number of FOI Observations For Each Date",
                      x = "Date",
                      y = "Number of Observations")

# Stack the plots (patchwork).
plot1 / plot2

dev.off()
