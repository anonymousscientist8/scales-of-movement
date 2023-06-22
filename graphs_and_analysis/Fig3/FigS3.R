library(asnipe)
library(tidyverse)
library(igraph)
library(stringi)
rm(list = ls())

# Load data
gd <- read.csv("C:\\Users\\raven\\Documents\\gd.csv")
gd <- gd[,c(-1,-2,-3)] # Reorganize

# Pull out cluster switching effects and partner switching effects
cs <- c(gd$cs,gd$cs)
ps <- c(gd$ps,gd$ps2)

# And made a list of what type of switching is being looked at
ps_type <- rep(0,length(cs))

# Then mark as partner switching including or not including switches resulting
# from cluster switching
for (i in 1:length(cs)) {
  if (i > length(gd$cs)) {
    ps_type[i] <- "ps2" 
  } else {
    ps_type[i] <- "ps"
  }
}

# Then make the data into a dataframe, making sure the elements are numeric
gd2 <- data.frame(cbind(ps_type,cs,ps))
gd2$cs <-  as.numeric(gd2$cs)
gd2$ps <-  as.numeric(gd2$ps)

# Create the label
label <- c(
  "ps" = "(a) Total Partner Switching",
  "ps2" = "(b) Within-Cluster Partner Switching"
)

# Set the theme
theme_new <- theme_set(theme_bw())
theme_new <- theme_update(strip.background = element_blank(),strip.text = element_text(hjust = 0))

# Then graph within-cluster parter switching vs cluster switching
ggplot(data = gd, mapping =  aes(x = cs, y = ps2)) + 
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab("cluster switching rate (switches / day)") +
  ylab("partner switching rate (switches / hour)")

