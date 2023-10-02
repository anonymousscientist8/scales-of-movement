rm(list = ls())

# load packages
library(igraph)
library(tidyverse)

# Import Data
perms1 <- read.csv("C:\\Users\\raven\\Documents\\exp548.csv")
colnames(perms1) <- c("X","exp1")
perms2 <- read.csv("C:\\Users\\raven\\Documents\\exp2400.csv")
colnames(perms2) <- c("X","exp2")
perms3 <- read.csv("C:\\Users\\raven\\Documents\\exp2401.csv")
colnames(perms3) <- c("X","exp3")
perms4 <- read.csv("C:\\Users\\raven\\Documents\\exp2402.csv")
colnames(perms4) <- c("X","exp4")
perms5 <- read.csv("C:\\Users\\raven\\Documents\\exp2403.csv")
colnames(perms5) <- c("X","exp5")
obs <- read.csv("C:\\Users\\raven\\Documents\\obs.csv")
colnames(obs) <- c("cols","value")
obs$cols <- c("exp1","exp2","exp3","exp4","exp5")

# bind into data frame
d <- cbind(perms1,perms2,perms3,perms4,perms5)
d <- d[,-1]
d <- d[,-2]
d <- d[,-3]
d <- d[,-4]
d <- d[,-5]

# Add label
label <- c(
  "exp1" = "(a) original simulation",
  "exp2" = "(b) removing hierarchically embedded scales of movement (bats stay within one cluster)",
  "exp3" = "(c) removal of individual variation in partner switching (bats also have same switching propensity)",
  "exp4" = "(d) removal of partner fidelity due to location (bats also groom a random available partner every minute)",
  "exp5" = "(e) removal of variation in behavior (bats also groom at same time)"
)

# Theme
theme_new <- theme_set(theme_bw())
theme_new <- theme_update(strip.background = element_blank(),strip.text = element_text(hjust = 0))

# Plot
(plot <- 
    ggplot() + 
    geom_histogram(data = gather(d, cols, value), aes(x = value), binwidth = 0.001, colour = "lightblue", fill = "lightblue") +
    geom_vline(data = obs, aes(xintercept = value), color = "black", linetype = "solid") +
    facet_wrap(.~cols, scales = "fixed", labeller = as_labeller(label), nrow = 5, strip.position = "top") +
    xlab("social differentiation: expected (blue) vs observed (black)") +
    geom_hline(yintercept = 0, color = "light blue") +
    ylab("frequency")+
    labs(subtitle = "Convergence of agent-based model with reference model shows that\nhierarchically embedded scales of movement creates false social differentiation"))
