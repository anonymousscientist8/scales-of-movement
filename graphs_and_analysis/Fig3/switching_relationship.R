library(asnipe)
library(tidyverse)
library(igraph)
library(stringi)
rm(list = ls())

###############################################################################
# Predicting grooming centrality from cluster switching
rm(list = ls())

# get cluster switching prob
cs <- 
  read.csv("filepath\\crevice_switching.csv") %>% 
  group_by(bat) %>% 
  summarize(cs.prob= mean(switch, na.rm=T))

# get partner switching prob
ps <- 
  read.csv("filepath\\partner_switching.csv") %>% 
  mutate(bat= actor) %>% 
  group_by(bat) %>% 
  summarize(ps.prob= mean(switch, na.rm=T))

# get partner switching prob WITHIN CLUSTER
ps2 <- 
  read.csv("filepath\\partner_switching.csv") %>% 
  mutate(bat= actor) %>% 
  # remove switches due to bat switching cameras
  arrange(bat, camera, video) %>% 
  # delete first observation for each camera
  filter(camera == lag(camera)) %>% 
  group_by(camera, bat) %>% 
  summarize(switch= mean(switch, na.rm=T)) %>% 
  group_by(bat) %>% 
  summarize(ps2.prob= mean(switch, na.rm=T))

# get grooming degree
gd <- 
  read.csv("filepath\\partner_switching.csv") %>% 
  filter(duration >0) %>% 
  mutate(bat= actor) %>% 
  group_by(bat, receiver) %>% 
  summarize(n=n()) %>% 
  group_by(bat) %>% 
  summarize(n=n()) 

# add cluster switching data
gd$cs <- cs$cs.prob[match(gd$bat, cs$bat)]

# add partner switching data
gd$ps <- ps$ps.prob[match(gd$bat, ps$bat)]

# add within-cluster partner switching
gd$ps2 <- ps2$ps2.prob[match(gd$bat, ps$bat)]

write.csv(gd, "filepath\\gd.csv")
