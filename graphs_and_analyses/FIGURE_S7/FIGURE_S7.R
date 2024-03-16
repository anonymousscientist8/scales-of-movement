rm(list = ls())

# load packages
library(igraph)
library(tidyverse)

# Import Data
perms1 <- read.csv("filepath\\exp2600a.csv")
colnames(perms1) <- c("X","exp1")
perms2 <- read.csv("filepath\\exp2600b.csv")
colnames(perms2) <- c("X","exp2")
perms3 <- read.csv("filepath\\exp2600c.csv")
colnames(perms3) <- c("X","exp3")
obs <- read.csv("filepath\\obs2.csv")
colnames(obs) <- c("cols","value")
obs$cols <- c("exp1","exp2","exp3")

# bind into data frame
d <- cbind(perms1,perms2,perms3)
d <- d[,-1]
d <- d[,-2]
d <- d[,-3]

# Add label
label <- c(
  "exp1" = "(a) unconstrained permutation test",
  "exp2" = "(b) semi-constrained permutation test",
  "exp3" = "(c) fully constrained permutation test"
)

# Theme
theme_new <- theme_set(theme_bw())
theme_new <- theme_update(strip.background = element_blank(),strip.text = element_text(hjust = 0))

# Plot
(plot <-
    ggplot() +
    geom_histogram(data = gather(d, cols, value), aes(x = value), binwidth = 0.001, colour = "light blue", fill = "light blue") +
    geom_hline(yintercept = 0, color = "azure2") +
    geom_vline(data = obs, aes(xintercept = value), color = "black", linetype = "solid") +
    facet_wrap(.~cols, scales = "fixed", labeller = as_labeller(label), nrow = 1, strip.position = "top") +
    xlab("social differentiation: expected (blue) vs observed (black)") +
    ylab("frequency"))
