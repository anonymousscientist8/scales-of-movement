library(asnipe)
library(tidyverse)
library(igraph)
library(stringi)
rm(list = ls())

# Download social differentiation data
roost_soc_diff <- read.csv("C:\\Users\\raven\\Documents\\roost_soc_diff.csv")
cluster_soc_diff <- read.csv("C:\\Users\\raven\\Documents\\cluster_soc_diff.csv")
partner_soc_diff <- read.csv("C:\\Users\\raven\\Documents\\partner_soc_diff.csv")

# Get the number of rows for each dataset and find the max number of rows between
# the three
r_rows <- nrow(roost_soc_diff)
c_rows <- nrow(cluster_soc_diff)
p_rows <- nrow(partner_soc_diff)
max_row <- max(r_rows,c_rows,p_rows)

# Whichever doesn't have the max number of rows, add NA's until the same length
# as the largest dataset
if (r_rows != max_row) {
  roost_soc_diff[r_rows + (max_row-r_rows),] <- NA
}
if (c_rows != max_row) {
  cluster_soc_diff[c_rows + (max_row-c_rows),] <- NA
}
if (p_rows != max_row) {
  partner_soc_diff[p_rows + (max_row-p_rows),] <- NA
}

# Combine and order the data
merged_data <- cbind(roost_soc_diff,cluster_soc_diff,partner_soc_diff)
merged_data <- merged_data[,-c(1,3,4,6,7,9)]
colnames(merged_data) <- c("a","b","c")

# Find the mean of each column and make a dataframe with labels
Mean <- data.frame(colMeans(merged_data, na.rm = T))
colnames(Mean) <- 'value'
cols <- c("a","b","c")
Mean <- cbind(cols,Mean)

# Create graph labels
label <- c(
  "a" = "(a) Roost Association Network Social Diff.",
  "b" = "(b) Cluster Association Network Social Diff.",
  "c" = "(c) Grooming Nework Social Diff."
)

# Setup graph theme
theme_new <- theme_set(theme_bw())
theme_new <- theme_update(strip.background = element_blank(),strip.text = element_text(hjust = 0))

# Plot as histogram
ggplot() + 
  geom_histogram(data = gather(merged_data, cols, value), aes(x = value), bins = 25, colour = "black", fill = "light blue") +
  geom_vline(data = Mean, aes(xintercept = value), color = "red", linetype = "dashed") +
  facet_wrap(.~cols, scales = "free", labeller = as_labeller(label),nrow = 1) +
  xlab("Social Differentiation")

