rm(list = ls())

# load packages
library(igraph)
library(tidyverse)

# Load data
degree_est0 <- read.csv("C:\\Users\\raven\\Documents\\degree_est5.csv") %>%
  mutate(label= "(a) individually variable rates of switching roost, cluster, and partner")
degree_est1 <- read.csv("C:\\Users\\raven\\Documents\\degree_est6.csv") %>%
  mutate(label= "(b) individually variable rates of roost switching")
degree_est2 <- read.csv("C:\\Users\\raven\\Documents\\degree_est7.csv") %>%
  mutate(label= "(c) individually variable rates of cluster switching")
degree_est3 <- read.csv("C:\\Users\\raven\\Documents\\degree_est8.csv") %>%
  mutate(label= "(d) individually variable rates of partner switching")
degree_est4 <- read.csv("C:\\Users\\raven\\Documents\\degree_est9.csv") %>%
  mutate(label= "(e) individually variable rates of switching that are correlated across types")

# compile data
d <- rbind(degree_est0,degree_est1,degree_est2,degree_est3,degree_est4)

# get baselines
base_deg_roost1 <-
  median(degree_est2$scale.rs.)
base_deg_roost2 <-
  median(degree_est3$scale.rs.)
base_deg_cluster1 <-
  median(degree_est1$scale.cs.)
base_deg_cluster2 <-
  median(degree_est3$scale.cs.)
base_deg_partner1 <-
  median(degree_est1$scale.ps.)
base_deg_partner2 <-
  median(degree_est2$scale.ps.)

# plot
d %>%
  pivot_longer(scale.rs.:scale.ps.,
               names_to = "t",
               values_to = "coefficient") %>%
  # put them in order from big to small
  mutate(
    type = case_when(
      t == "scale.rs." ~ "1. roost switching",
      t == "scale.cs." ~ "2. cluster switching",
      t == "scale.ps." ~ "3. partner switching")) %>%
  mutate(type= fct_rev(factor(type))) %>% 
  ggplot(aes(x = type, y = coefficient, color = type)) +
  facet_wrap(~label, ncol=1)+
  geom_violin(width = 0.4) +
  geom_boxplot(width = 0.1) +
  geom_hline(yintercept = base_deg_roost1, linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = base_deg_roost2,
             linetype = 'dashed',
             color = 'red') +
  geom_hline(yintercept = base_deg_cluster1, linetype = 'dashed', color = 'dark green') +
  geom_hline(yintercept = base_deg_cluster2,
             linetype = 'dashed',
             color = 'dark green') +
  geom_hline(yintercept = base_deg_partner1, linetype = 'dashed', color = 'blue') +
  geom_hline(yintercept = base_deg_partner2,
             linetype = 'dashed',
             color = 'blue') +
  coord_flip()+
  theme_bw() +
  theme(strip.text.x = element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("effect on grooming degree (standardized coefficient)") +
  scale_color_manual(values= c("blue", "dark green", "red"))

