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
  "exp1" = "A. Original simulation",
  "exp2" = "B. Removing hierarchically embedded scales of movement",
  "exp3" = "C. Removing individual variation in partner switching propensity",
  "exp4" = "D. Bats groom a random available partner every minute",
  "exp5" = "E. All bats groom at same time"
)

# Theme
theme_new <- theme_set(theme_bw())
theme_new <- theme_update(strip.background = element_blank(),strip.text = element_text(hjust = 0))

# Plot
(plot <-
    ggplot() +
    geom_histogram(data = gather(d, cols, value), aes(x = value), binwidth = 0.001, colour = "black", fill = "light blue") +
    geom_vline(data = obs, aes(xintercept = value), color = "red", linetype = "solid") +
    facet_wrap(.~cols, scales = "free", labeller = as_labeller(label), nrow = 5, strip.position = "top") +
    xlab("Social Differentiation: Expected (blue) vs Observed (red)") +
    ylab("Frequency"))
