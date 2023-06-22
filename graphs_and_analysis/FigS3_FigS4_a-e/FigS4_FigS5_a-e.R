rm(list = ls())

# load packages
library(igraph)
library(tidyverse)

# Load data
pagerank_est0 <- read.csv("C:\\Users\\raven\\Documents\\pagerank_est0.csv") %>%
  mutate(label= "(a) individually variable rates of switching roost, cluster, and partner")
pagerank_est1 <- read.csv("C:\\Users\\raven\\Documents\\pagerank_est1.csv") %>%
  mutate(label= "(b) individually variable rates of partner switching")
pagerank_est2 <- read.csv("C:\\Users\\raven\\Documents\\pagerank_est2.csv") %>%
  mutate(label= "(c) individually variable rates of cluster switching")
pagerank_est3 <- read.csv("C:\\Users\\raven\\Documents\\pagerank_est3.csv") %>%
  mutate(label= "(d) individually variable rates of roost switching")
pagerank_est4 <- read.csv("C:\\Users\\raven\\Documents\\pagerank_est4.csv") %>%
  mutate(label= "(e) individually variable rates of switching that are correlated across types")

# compile data
f <- rbind(pagerank_est0,pagerank_est1,pagerank_est2,pagerank_est3,pagerank_est4)

# get baselines
base_pr_roost1 <-
  median(pagerank_est2$scale.rs.)
base_pr_roost2 <-
  median(pagerank_est3$scale.rs.)
base_pr_cluster1 <-
  median(pagerank_est1$scale.cs.)
base_pr_cluster2 <-
  median(pagerank_est3$scale.cs.)
base_pr_partner1 <-
  median(pagerank_est1$scale.ps.)
base_pr_partner2 <-
  median(pagerank_est2$scale.ps.)

# plot
f %>%
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
  ggplot(aes(x = type, y = coefficient, color = type), hjust) +
  facet_wrap(~label, ncol=1)+
  geom_violin(width = 0.4) +
  geom_boxplot(width = 0.1) +
  geom_hline(yintercept = base_pr_roost1, linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = base_pr_roost2,
             linetype = 'dashed',
             color = 'red') +
  geom_hline(yintercept = base_pr_cluster1, linetype = 'dashed', color = 'dark green') +
  geom_hline(yintercept = base_pr_cluster2,
             linetype = 'dashed',
             color = 'dark green') +
  geom_hline(yintercept = base_pr_partner1, linetype = 'dashed', color = 'blue') +
  geom_hline(yintercept = base_pr_partner2,
             linetype = 'dashed',
             color = 'blue') +
  coord_flip()+
  theme_bw() +
  theme(strip.text.x = element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("effect on pagerank (standardized coefficient)") +
  scale_color_manual(values= c("blue", "dark green", "red"))

