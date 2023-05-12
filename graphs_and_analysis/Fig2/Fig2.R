library(asnipe)
library(tidyverse)
library(igraph)
library(stringi)
rm(list = ls())

# Load centrality and switching rate data
gd <- read.csv("filepath\\gd.csv")
centrality <- read.csv("filepath\\centrality.csv")

# Make two lists of switching rate and grooming degree
switch <- c(centrality$roost_switch,gd$cs,gd$ps2)
groomed <- c(centrality$centrality,gd$n,gd$n)

# Make a third list of the type of movement
type <- rep(0,length(groomed))
for (i in 1:length(type)) {
  if (i <= length(centrality$centrality)) {
    type[i] <- "a"
  } else {
    if (i <= (length(centrality$centrality)+length(gd$cs))) {
      type[i] <- "b"
    } else {
      if (i <= (length(centrality$centrality)+length(gd$cs)*2)) {
        type[i] <- "c"
      }
    }
  }
}

# Combine lists into datafeame
gd2 <- data.frame(cbind(type,switch,groomed))
gd2$switch <- as.numeric(gd2$switch)
gd2$groomed <- as.numeric(gd2$groomed)

# Graph label
label <- c(
  "a" = "(a) Roost Switching (switches / day)",
  "b" = "(b) Cluster Switching (switches / day)",
  "c" = "(c) Within-Cluster Partner Switching (switches / hour)"
)

# Graph theme
theme_new <- theme_set(theme_bw())
theme_new <- theme_update(strip.text = element_text(hjust = 0), strip.placement = "outside", axis.title.x=element_blank(), strip.background = element_blank())

# Plot
ggplot(data = gd2, mapping =  aes(x = switch, y = groomed)) + 
  geom_point() +
  geom_smooth(method = 'lm') +
  ylab("Partners Groomed")  +
  facet_wrap(~type, nrow = 3, scales = "free", strip.position = "bottom", labeller = as_labeller(label)) +
  geom_smooth(data = data.frame(gd2,z="c"), method = '')

# Run permutation tests and get statistic metrics, relationship metrics
model <- lm(centrality~roost_switch, data=centrality)
summary(model)
root <- rep(0,1000)
for (i in 1:1000) {
  data2 <- centrality[sample(nrow(centrality),size = nrow(centrality),replace = T),]
  root[i] <- summary(lm(centrality ~ roost_switch, data = data2))$r.squared
}
quantile(root, probs = c(0.025,0.975))
lm(data.frame(scale(model$model)))
model <- lm(n~cs, data=gd)
summary(model)
coot <- rep(0,1000)
for (i in 1:1000) {
  data2 <- gd[sample(nrow(gd),size = nrow(gd),replace = T),]
  coot[i] <- summary(lm(n ~ cs, data = data2))$r.squared
}
quantile(coot, probs = c(0.025,0.975))
lm(data.frame(scale(model$model)))
model <- lm(n~ps, data=gd)
summary(model)
lm(data.frame(scale(model$model)))
model <- lm(n~ps2, data=gd)
summary(model)
poot <- rep(0,1000)
for (i in 1:1000) {
  data2 <- gd[sample(nrow(gd),size = nrow(gd),replace = T),]
  poot[i] <- summary(lm(n ~ ps2, data = data2))$r.squared
}
quantile(poot, probs = c(0.025,0.975))
lm(data.frame(scale(model$model)))
model <- lm(cs~ps, data=gd)
summary(model)
lm(data.frame(scale(model$model)))
model <- lm(cs~ps2, data=gd)
summary(model)
lm(data.frame(scale(model$model)))

