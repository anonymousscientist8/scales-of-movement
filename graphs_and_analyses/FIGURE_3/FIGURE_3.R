library(asnipe)
library(tidyverse)
library(igraph)
library(stringi)
library(glmmTMB)
library(MASS)
library(performance)
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
  geom_smooth(method = 'glm', method.args = list(family="quasipoisson")) +
  ylab("Count of Partners Groomed")  +
  facet_wrap(~type, nrow = 1, scales = "free", strip.position = "bottom", labeller = as_labeller(label))

# Run permutation tests and get statistic metrics, relationship metrics
#model <- glm(centrality~scale(roost_switch), family = "poisson", data=centrality) # Get the roost switching linear model
overdisp_fun <- function(model) { # Check for dispersion
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
#overdisp_fun(model)
model <- glmmTMB(centrality ~ scale(roost_switch), family = "nbinom1", data = centrality)
summary(model) # Look at statstics associated with the model
#1 - (model$deviance/model$null.deviance)
rootb <- rep(0,1000) # Initialize vector
for (i in 1:1000) { # Run 1000 permutations and find the R^2 value
  data2 <- centrality[sample(nrow(centrality),size = nrow(centrality),replace = T),]
  a <- summary(glmmTMB(centrality ~ scale(roost_switch), family = "nbinom1", data = data2))
  b <- a$coefficients[1]
  c <- data.frame(b)
  rootb[i] <- c$cond.Estimate[2]
}
quantile(rootb, probs = c(0.025,0.975)) # Find where 95% of the beta values lie
#model <- glm(n~scale(cs), family = "poisson", data=gd) # Get the cluster switching linear model
#overdisp_fun(model)
#1 - (model$deviance/model$null.deviance)
gd$scs <- as.numeric(scale(gd$cs))
model <- glmmTMB(n ~ scs, family = "nbinom1", data = gd)
summary(model) # Look at statistics
cootb <- rep(0,1000)
for (i in 1:1000) { # Run 1000 permutations and find the R^2 value
  data2 <- gd[sample(nrow(gd),size = nrow(gd),replace = T),]
  a <- summary(glmmTMB(n ~ scs, family = "nbinom1", data = data2))
  b <- a$coefficients[1]
  c <- data.frame(b)
  cootb[i] <- c$cond.Estimate[2]
}
quantile(cootb, probs = c(0.025,0.975)) # Find where 95% of the beta values lie
#model <- glm(n~scale(ps2), family = "poisson", data=gd)  # Get the within-cluster partner switching linear model
#overdisp_fun(model)
model <- glmmTMB(n ~ scale(ps2), family = "nbinom1", data = gd)
summary(model) # View statistics
pootb <- rep(0,1000)
for (i in 1:1000) { # # Run 1000 permutations and find the R^2 value
  data2 <- gd[sample(nrow(gd),size = nrow(gd),replace = T),]
  a <- summary(glmmTMB(n ~ scale(ps2), family = "nbinom1", data = data2))
  b <- a$coefficients[1]
  c <- data.frame(b)
  pootb[i] <- c$cond.Estimate[2]
}
quantile(pootb, probs = c(0.025,0.975)) # Find where 95% of the beta values lie
# Get statistics for the relation between cluster switching and partner switching
model <- glm(scale(cs)~scale(ps2), data=gd)
summary(model) # View statistics
poot <- rep(0,1000) # Initialize vector
pootb <- rep(0,1000)
for (i in 1:1000) { # # Run 1000 permutations and find the R^2 value
  data2 <- gd[sample(nrow(gd),size = nrow(gd),replace = T),]
  poot[i] <- summary(lm(scale(cs) ~ scale(ps2), data = data2))$r.squared
  pootb[i] <- summary(lm(scale(cs) ~ scale(ps2), data = data2))$coefficients[2]
}
quantile(poot, probs = c(0.025,0.975)) # Find where 95% of the R^2 values lie
quantile(pootb, probs = c(0.025,0.975)) # Find where 95% of the beta values lie
lm(data.frame(scale(model$model)))

