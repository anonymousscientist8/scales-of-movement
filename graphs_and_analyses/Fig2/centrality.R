library(asnipe)
library(tidyverse)
library(igraph)
library(stringi)

###############################################################################
# Predicting grooming centrality from cluster switching and partner switching
rm(list = ls())

# Get cluster switching and partner switching rates
cs <- read.csv("C:\\Users\\raven\\Documents\\cluster_sw.csv")
ps <- read.csv("C:\\Users\\raven\\Documents\\partner_sw.csv")
ps2 <- read.csv("C:\\Users\\raven\\Documents\\partner_sw2.csv")

# get grooming degree
gd <- 
  read.csv("C:\\Users\\raven\\Documents\\partner_switching.csv") %>% 
  filter(duration >0) %>% 
  mutate(bat= actor) %>% 
  group_by(bat, receiver) %>% 
  summarize(n=n()) %>% 
  group_by(bat) %>% 
  summarize(n=n()) 

# add cluster switching data
gd$cs <- cs$cluster_sw[match(gd$bat, cs$bats_i)]

# add partner switching data
gd$ps <- ps$partner_sw[match(gd$bat, ps$bats_i)]

# add within-cluster partner switching
gd$ps2 <- ps2$partner_sw2[match(gd$bat, ps$bats_i)]

# Obtain best fit and relationship metrics
coot <- rep(0,1000) # Initialize vector
for (i in 1:1000) { # Run 1000 permutation tests and extract coefficients
  data2 <- gd[sample(nrow(gd),size = nrow(gd),replace = T),]
  fit <- lm(n ~ cs, data = data2)
  a <- coefficients(fit)
  coot[i] <- a[2]
}
poot <- rep(0,1000) # Initialize vector
for (i in 1:1000) { # Run 1000 permutation tests and extract coefficients
  data2 <- gd[sample(nrow(gd),size = nrow(gd),replace = T),]
  fit <- lm(n ~ ps, data = data2)
  a <- coefficients(fit)
  poot[i] <- a[2]
}
poot2 <- rep(0,1000) # Initialize vector
for (i in 1:1000) { # Run 1000 permutation tests and extract coefficients
  data2 <- gd[sample(nrow(gd),size = nrow(gd),replace = T),]
  fit <- lm(n ~ ps, data = data2)
  a <- coefficients(fit)
  poot2[i] <- a[2]
}

# plot effect of cluster switching
gd %>%
  ggplot(aes(x=cs, y=n))+
  geom_point()+
  geom_smooth(method= "lm")+
  xlab("Proportion of Vampire Bat Observations with Cluster Switch")+
  ylab("Grooming Outdegree Centrality")+
  ggtitle("Cluster Switching Predicts Grooming Outdegree Centrality")+
  theme_bw()

# plot effect of partner switching
gd %>%
  ggplot(aes(x=ps, y=n))+
  geom_point()+
  geom_smooth(method= "lm")+
  xlab("Proportion of Vampire Bat Observations with Partner Switch")+
  ylab("Grooming Outdegree Centrality")+
  ggtitle("Partner Switching Predicts Grooming Outdegree Centrality")+
  theme_bw()

# plot effect of within-cluster partner switching
gd %>%
  ggplot(aes(x=ps2, y=n))+
  geom_point()+
  geom_smooth(method= "lm")+
  xlab("Proportion of Vampire Bat Observations with Partner Switch")+
  ylab("Grooming Outdegree Centrality")+
  ggtitle("Within-Cluster Partner Switching Predicts Degree Centrality")+
  theme_bw()

# Export
write.csv(gd, "C:\\Users\\raven\\Documents\\gd.csv")

################################################################################
# Predicting grooming centrality from roost switching
rm(list = ls())

# load packages
library(tidyverse)
library(stringi)

# get roost switching prob
rs <- read.csv("C:\\Users\\raven\\Documents\\roost_sw.csv")
bats_r <- unique(rs$Bats_w)

  
# Load in transcripted data of grooming events for roost switching data
grooming <- read.csv("C:\\Users\\raven\\Documents\\transcribe.csv")

# Find the grooming degree centrality
bats <- unique(grooming$donor) # Find all unique bats
centrality <- rep(0,length(bats)) # Initialize vector
for (i in 1:length(centrality)) { # For all bats
  temp <- grooming[!(grooming$donor != bats[i]),] # Find all situations where the focal bat is the donor
  centrality[i] <- length(temp$donor) # And find the number of bats groomed
}
centrality <- data.frame(bats,centrality) # Then make a data frame
for (i in 1:length(centrality$bats)) { # Mark these as females
  stri_sub(centrality$bats[i],nchar(centrality$bats[i])-6,nchar(centrality$bats[i])-3) <- "F"
}

# Look at common elements
bats_p <- centrality$bats # Get bats in this dataset
bool <- bats_p %in% bats_r # Find all bats found in both datasets
for (i in 1:length(bool)) { # For all bats
  if (bool[i] == FALSE) { # If bats are not in both datasets
    centrality <- centrality[!(centrality$bats == bats_p[i]),] # Remove them 
  }
}
roost_switch <- rep(0,length(centrality$bats)) # Initialize vector
for (i in 1:length(centrality$bats)) { # For the length of the vector
  for (j in 1:length(rs$Bats_w)) { # For all bats in the roost switching data
    bool <- rs$Bats_w[j] == centrality$bats[i] # Determine whether the bats in both datasets are the same
    if (bool == TRUE) { # And if they are
      roost_switch[i] <- rs$roost_sw[j] # store the roost switching data
    }
  }
}
centrality <- data.frame(cbind(centrality, roost_switch)) # Combine the centrality and roost switching data

root <- rep(0,1000) # Initialize vector
for (i in 1:1000) { # Run 1000 permutations and extract coefficients
  data2 <- centrality[sample(nrow(centrality),size = nrow(centrality),replace = T),]
  fit <- lm(centrality ~ roost_switch, data = data2)
  a <- coefficients(fit)
  root[i] <- a[2]
}

# Plot Degree Centrality against Roost Switching Rate
ggplot(data = centrality) +
  geom_point(mapping = aes(x = roost_switch, y = centrality)) +
  geom_smooth(mapping = aes(x = roost_switch, y = centrality), method = 'lm') +
  ylab("Grooming Outdegree Centrality")+
  xlab("Proportion of Vampire Bat Observations with Roost Switch")+
  ggtitle("Roost Switching Does Not Predict Outdegree Centrality")+
  theme_bw()

# Export
write.csv(centrality, "C:\\Users\\raven\\Documents\\centrality.csv")

