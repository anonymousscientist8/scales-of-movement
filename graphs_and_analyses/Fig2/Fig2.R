library(asnipe)
library(tidyverse)
library(igraph)
library(stringi)
rm(list = ls())

# Load centrality and switching rate data
gd <- read.csv("C:\\Users\\raven\\Documents\\gd.csv")
centrality <- read.csv("C:\\Users\\raven\\Documents\\centrality.csv")

# Make two lists of switching rate and grooming degree
switch <- c(centrality$roost_switch,gd$cs,gd$ps2)
groomed <- c(centrality$centrality,gd$n,gd$n)

# Make a third list of the type of movement
type <- rep(0,length(groomed))
for (i in 1:length(type)) {
  if (i <= length(centrality$centrality)) {
    type[i] <- "a" # Denote where roost switching data is stored
  } else {
    if (i <= (length(centrality$centrality)+length(gd$cs))) {
      type[i] <- "b" # Denote where cluster switching data is stored
    } else {
      if (i <= (length(centrality$centrality)+length(gd$cs)*2)) {
        type[i] <- "c" # Denote where partner switching data is stored
      }
    }
  }
}

# Combine lists into datafeame and make sure theyh are numeric
gd2 <- data.frame(cbind(type,switch,groomed))
gd2$switch <- as.numeric(gd2$switch)
gd2$groomed <- as.numeric(gd2$groomed)

# Graph label
label <- c(
  "a" = "Roost Switching (switches / day)",
  "b" = "Cluster Switching (switches / day)",
  "c" = "Partner Switching (switches / hour)"
)

# Graph theme
theme_new <- theme_set(theme_bw())
theme_new <- theme_update(legend.position = "none", strip.text = element_text(size = 16), strip.placement = "outside", axis.title.x=element_blank(), strip.background = element_blank(), axis.title.y = element_text(size = 16))

# Plot
ggplot(data = gd2, mapping =  aes(x = switch, y = groomed, color = type)) + 
  geom_point() +
  geom_smooth(method = 'lm') +
  ylab("Count of Partners Groomed")  +
  facet_wrap(~type, nrow = 1, scales = "free", strip.position = "bottom", labeller = as_labeller(label))

# Run permutation tests and get statistic metrics, relationship metrics
model <- lm(centrality~scale(roost_switch), data=centrality) # Get the roost switching linear model
summary(model) # Look at statstics associated with the model
root <- rep(0,1000) # Initialize vector
rootb <- rep(0,1000) # Initialize vector
for (i in 1:1000) { # Run 1000 permutations and find the R^2 value
  data2 <- centrality[sample(nrow(centrality),size = nrow(centrality),replace = T),]
  root[i] <- summary(lm(centrality ~ scale(roost_switch), data = data2))$r.squared
  rootb[i] <- summary(lm(centrality ~ scale(roost_switch), data = data2))$coefficients[2]
}
quantile(root, probs = c(0.025,0.975)) # Find where 95% of the R^2 values lie
quantile(rootb, probs = c(0.025,0.975)) # Find where 95% of the beta values lie
lm(data.frame(scale(model$model))) # And observe more statistics
model <- lm(n~scale(cs), data=gd) # Get the cluster switching linear model
summary(model) # Look at statistics
coot <- rep(0,1000) # Initialize vector
cootb <- rep(0,1000)
for (i in 1:1000) { # Run 1000 permutations and find the R^2 value
  data2 <- gd[sample(nrow(gd),size = nrow(gd),replace = T),]
  coot[i] <- summary(lm(n ~ scale(cs), data = data2))$r.squared
  cootb[i] <- summary(lm(n ~ scale(cs), data = data2))$coefficients[2]
}
quantile(coot, probs = c(0.025,0.975)) # Find where 95% of the R^2 values lie
quantile(cootb, probs = c(0.025,0.975)) # Find where 95% of the beta values lie
lm(data.frame(scale(model$model))) # And observe more statistics
model <- lm(n~ps, data=gd) # Get the partner switching linear model
summary(model) # Look at statistics
lm(data.frame(scale(model$model))) # And view more statistics
model <- lm(n~scale(ps2), data=gd)  # Get the within-cluster partner switching linear model
summary(model) # View statistics
poot <- rep(0,1000) # Initialize vector
pootb <- rep(0,1000)
for (i in 1:1000) { # # Run 1000 permutations and find the R^2 value
  data2 <- gd[sample(nrow(gd),size = nrow(gd),replace = T),]
  poot[i] <- summary(lm(n ~ scale(ps2), data = data2))$r.squared
  pootb[i] <- summary(lm(n ~ scale(ps2), data = data2))$coefficients[2]
}
quantile(poot, probs = c(0.025,0.975)) # Find where 95% of the R^2 values lie
quantile(pootb, probs = c(0.025,0.975)) # Find where 95% of the beta values lie
lm(data.frame(scale(model$model))) # And observe more statistics
# Get statistics for the relation between cluster switching and partner switching
model <- lm(cs~ps, data=gd)
summary(model)
lm(data.frame(scale(model$model)))
model <- lm(scale(cs)~scale(ps2), data=gd)
summary(model)
lm(data.frame(scale(model$model)))

