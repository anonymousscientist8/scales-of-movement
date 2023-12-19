rm(list = ls())

# load packages
library(igraph)
library(tidyverse)

# Load data
degree_est0a <- read.csv("Filepath\\degree_est0.csv") %>%
  mutate(label= "(a) individually variable rates of switching roost, cluster, and partner across 200 bats")
degree_est0 <- read.csv("Filepath\\degree_est0.csv") %>%
  mutate(label= "individually variable rates of switching roost, cluster, and partner across 200 bats")
degree_est1 <- read.csv("Filepath\\degree_est1.csv") %>%
  mutate(label= "(a) individually variable rates of roost switching with 200 bats")
degree_est2 <- read.csv("Filepath\\degree_est2.csv") %>%
  mutate(label= "(b) individually variable rates of cluster switching with 200 bats")
degree_est3 <- read.csv("Filepath\\degree_est3.csv") %>%
  mutate(label= "(c) individually variable rates of partner switching with 200 bats")
degree_est4 <- read.csv("Filepath\\degree_est4.csv") %>%
  mutate(label= "individually variable switching rates correlated across types for 200 bats")
degree_est5 <- read.csv("Filepath\\degree_est5.csv") %>%
  mutate(label= "(a) individually variable rates of switching roost, cluster, and partner across 100 bats")
degree_est6 <- read.csv("Filepath\\degree_est6.csv") %>%
  mutate(label= "(b) individually variable rates of roost switching with 100 bats")
degree_est7 <- read.csv("Filepath\\degree_est7.csv") %>%
  mutate(label= "(c) individually variable rates of cluster switching with 100 bats")
degree_est8 <- read.csv("Filepath\\degree_est8.csv") %>%
  mutate(label= "(d) individually variable rates of partner switching with 100 bats")
degree_est9 <- read.csv("Filepath\\degree_est9.csv") %>%
  mutate(label= "(b) individually variable switching rates correlated across types for 100 bats")

# compile data
d0 <- degree_est0
d1 <- rbind(degree_est1,degree_est2,degree_est3)
d2 <- rbind(degree_est0a,degree_est1,degree_est2,degree_est3,degree_est4)
d3 <- rbind(degree_est5,degree_est6,degree_est7,degree_est8,degree_est9)

# get reference effects across both population sizes
# 200
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
# 100
base_deg_roost3 <-
  median(degree_est7$scale.rs.)
base_deg_roost4 <-
  median(degree_est8$scale.rs.)
base_deg_cluster3 <-
  median(degree_est6$scale.cs.)
base_deg_cluster4 <-
  median(degree_est8$scale.cs.)
base_deg_partner3 <-
  median(degree_est6$scale.ps.)
base_deg_partner4 <-
  median(degree_est7$scale.ps.)

# plot figure 3
d0 %>% pivot_longer(scale.rs.:scale.ps.,
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
  geom_boxplot(width = 0.1, outlier.size = 0.25) +
  coord_flip()+
  theme_bw() +
  theme(strip.text.x = element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("effect on grooming degree (standardized coefficient)") +
  scale_color_manual(values= c("blue", "dark green", "red"))

# plot figure 4
d1 %>%
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
  geom_boxplot(width = 0.1, outlier.size = 0.25) +
  geom_hline(yintercept = max(c(base_deg_roost1,base_deg_roost2)), linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = max(c(base_deg_cluster1,base_deg_cluster2)), linetype = 'dashed', color = 'dark green') +
  geom_hline(yintercept = max(c(base_deg_partner1,base_deg_partner2)), linetype = 'dashed', color = 'blue') +
  coord_flip()+
  theme_bw() +
  theme(strip.text.x = element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("effect on grooming degree (standardized coefficient)") +
  scale_color_manual(values= c("blue", "dark green", "red"))

# plot Figure S3
d2 %>%
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
  geom_boxplot(width = 0.1, outlier.size = 0.25) +
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

# plot Figure S4
d3 %>%
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
  geom_boxplot(width = 0.1, outlier.size = 0.25) +
  geom_hline(yintercept = base_deg_roost3, linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = base_deg_roost4,
             linetype = 'dashed',
             color = 'red') +
  geom_hline(yintercept = base_deg_cluster3, linetype = 'dashed', color = 'dark green') +
  geom_hline(yintercept = base_deg_cluster4,
             linetype = 'dashed',
             color = 'dark green') +
  geom_hline(yintercept = base_deg_partner3, linetype = 'dashed', color = 'blue') +
  geom_hline(yintercept = base_deg_partner4,
             linetype = 'dashed',
             color = 'blue') +
  coord_flip()+
  theme_bw() +
  theme(strip.text.x = element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("effect on grooming degree (standardized coefficient)") +
  scale_color_manual(values= c("blue", "dark green", "red"))
