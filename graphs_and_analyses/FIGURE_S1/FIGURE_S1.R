library(asnipe)
library(tidyverse)
library(igraph)
library(stringi)
rm(list = ls())

# Download average movement rate data
roost_sw <- read.csv("filepath\\roost_sw.csv")
cluster_sw <- read.csv("filepath\\cluster_sw.csv")
partner_sw <- read.csv("filepath\\partner_sw.csv")

# Find the number of rows, and the maximum nnumber of rows between the three
r_rows <- nrow(roost_sw)
c_rows <- nrow(cluster_sw)
p_rows <- nrow(partner_sw)
max_row <- max(r_rows,c_rows,p_rows)

# Add NA rows until all datasets are the same length
if (r_rows != max_row) {
  roost_sw[r_rows + (max_row-r_rows),] <- NA
}
if (c_rows != max_row) {
  cluster_sw[c_rows + (max_row-c_rows),] <- NA
}
if (p_rows != max_row) {
  partner_sw[p_rows + (max_row-p_rows),] <- NA
}

# Merge the data and rearrange
merged_data <- cbind(roost_sw,cluster_sw,partner_sw)
merged_data <- merged_data[,-1]
merged_data <- merged_data[,-2]
merged_data <- merged_data[,-2]
merged_data <- merged_data[,-3]
merged_data <- merged_data[,-3]
merged_data <- merged_data[,-4]
colnames(merged_data) <- c("a","b","c")
merged_data$c <- merged_data$c * 24

# Find the mean of each column and add to data
Mean <- data.frame(colMeans(merged_data, na.rm = T))
colnames(Mean) <- 'value'
cols <- c("a","b","c")
Mean <- cbind(cols,Mean)

# Graph labels
label <- c(
  "a" = "(a) roost switching",
  "b" = "(b) cluster switching",
  "c" = "(c) partner switching"
)

# Graph theme
theme_new <- theme_set(theme_bw())
theme_new <- theme_update(strip.background = element_blank(),strip.text = element_text(hjust = 0))

# Plot
ggplot() + 
  geom_histogram(data = gather(merged_data, cols, value), aes(x = value), bins = 30, colour = "black", fill = "light blue") +
  geom_vline(data = Mean, aes(xintercept = value), color = "red", linetype = "dashed") +
  facet_wrap(.~cols, scales = "free", labeller = as_labeller(label), nrow = 1) +
  xlab("switching rate (switches / day)")
