# clear workspace
rm(list=ls())

# load packages
library(tidyverse)
library(lme4)

# get roost switching data, filtering out all instances where with more than
# 14 days since last switch
d <- 
  read.csv("filepath\\roost_switching.csv") %>% 
  mutate(switch= as.numeric(switch_w>0)) %>% 
  filter(n.days <=14) 

# check number of observations per bat
d %>% 
  group_by(ID) %>% 
  summarize(n= n()) %>% 
  arrange(n) %>% 
  pull(n) %>% 
  range()

# plot logistic regressions with all bats
d %>% 
  ggplot(aes(x=n.days, y=switch))+
  geom_point(size=2)+
  stat_smooth(aes(color= ID),method="glm", se=FALSE, method.args = list(family=binomial))+
  stat_smooth(method="glm", se=FALSE, method.args = list(family=binomial), color= "black")+
  theme(legend.position= "none") +
  xlab("number of days since last observed switch")+
  ylab("probability of switching")

# fit logistic mixed model using all bats as random intercepts
fit <- glmer(switch ~ n.days + (1 | ID), data = d, family = binomial)

# get intercepts
random_effects <- 
  coefficients(fit)$ID %>%
  rownames_to_column() %>% 
  dplyr::select(bat= rowname, intercept = '(Intercept)', coefficient= n.days)

# get max number of days to consider
max.days <-  14

# get dataframe of days for which to estimate a probability
days.to.predict <- data.frame(n.days= 1:max.days)

# get bats
bats <- unique(d$ID)

# make empty matrix of bats (rows) by days (cols)
bat.days <- matrix(data= NA, nrow= length(bats), ncol= max.days)

# fill in matrix with predicted probabilities
for (i in 1:nrow(bat.days)) {
  
  # choose bat
  temp <- bats[i]
  
  # get data from that bat
  d2 <- 
    random_effects %>% 
    filter(bat==temp)
  
  # make empty vector for each hour 0 to 24
  days <- rep(NA, max.days)
  # get predicted probability for each day for that bat
  for (j in 1:length(days)){
    # get response as log odds
    y <- d2$coefficient*j + d2$intercept
    # get odds  
    odds <- exp(y)
    # get probability
    prob <- odds / (1 + odds)
    days[j] <- prob  
  }
  
  # then put those values into matrix of probs by bat and h
  bat.days[i,] <- days
}

# label matrix
colnames(bat.days) <- paste0("day",1:max.days)
rownames(bat.days) <- bats
bat.days <- bat.days[order(bat.days[,1]),]

# save
write.csv(bat.days, file= "roost_switching_probs.csv")

# plot predictions for first ten bats
bat.days %>% 
  as.data.frame() %>% 
  rownames_to_column("bat") %>% 
  pivot_longer(cols= day1:day14, names_to= "day", values_to= "prob") %>% 
  mutate(day= as.numeric(substr(day, start=4, stop=6))) %>% 
  ggplot(aes(x=day, y=prob,color=bat))+
    geom_point(size=2)+
    geom_line()+
  ylab("probability of switch")+
  theme(legend.position= "none")

# clear workspace
rm(list=ls())

# get cluster switching data
d <- 
  read.csv("filepath\\cluster_switching.csv") %>% 
  mutate(latency= 2*latency) %>% 
  # only use latencies within the day
  filter(latency <=24)

# check number of observations per bat
d %>% 
  group_by(bat) %>% 
  summarize(n= n()) %>% 
  arrange(n) %>% 
  pull(n) %>% 
  range()

# plot logistic regressions with all bats
d %>% 
  ggplot(aes(x=latency, y=switch))+
  geom_point(size=2)+
  stat_smooth(aes(color= bat),method="glm", se=FALSE, method.args = list(family=binomial))+
  stat_smooth(method="glm", se=FALSE, method.args = list(family=binomial), color= "black")+
  theme(legend.position= "none") +
  xlab("number of hour samples since last observed switch")+
  ylab("probability of switching")

# fit logistic mixed model using all bats as random intercepts
fit <- glmer(switch ~ latency + (1 | bat), data = d, family = binomial)

# get intercepts
random_effects <- 
  coefficients(fit)$bat %>%
  rownames_to_column() %>% 
  dplyr::select(bat= rowname, intercept = '(Intercept)', coefficient= latency)

# get max number of hours
max.h <-  800

# get bats
bats <- unique(d$bat)

# make empty matrix of bats (rows) by hour (cols)
bat.h <- matrix(data= NA, nrow= length(bats), ncol= max.h)

# fill in matrix with predicted probabilities
for (i in 1:nrow(bat.h)) {
  
  # choose bat
  focal <- bats[i]
  
  # get random effects from that bat
  t <- 
    random_effects %>% 
    filter(bat==focal)
  
  # make empty vector for each hour 0 to 24
  hours <- rep(NA, max.h)
  # get predicted probability for each day for that bat
  for (j in 1:length(hours)){
    # get response as log odds
    y <- t$coefficient*j + t$intercept
    # get odds  
    odds <- exp(y)
    # get probability
    prob <- odds / (1 + odds)
    hours[j] <- prob  
  }
  
  # then put those values into matrix of probs by bat and h
  bat.h[i,] <- hours
}


# label matrix
colnames(bat.h) <- paste0("hour",1:max.h)
rownames(bat.h) <- bats
bat.h <- bat.h[order(bat.h[,1]),]

# plot predictions 
bat.h %>% 
  as.data.frame() %>% 
  rownames_to_column("bat") %>% 
  pivot_longer(cols= hour1:hour800, names_to= "hour", values_to= "prob") %>% 
  mutate(hour= as.numeric(substr(hour, start=5, stop=10))) %>% 
  ggplot(aes(x=hour, y=prob,color=bat))+
  geom_point(size=1)+
  geom_line()+  
  ylab("probability of switch")+
  theme(legend.position= "none")

# The following was used exclusively for correlation of probabilities between
# roost, cluster, and partner switching
# Break first cluster switching probabilities into quantiles
a <- quantile(bat.h[,1], probs = seq(0,1,1/80))
a <- data.frame(a)
# Create a second vector
b <- rep(0,length(a$a))
# For the length of the dataframe
for (i in 1:length(a$a)) {
  # Create a vector for each bat
  c <- rep(0,length(bat.h[,1]))
  # And for each entry in c
  for (j in 1:length(c)) {
  # Find the absolute difference between the quantiles and the actual cluster switching probabilities
    c[j] <- abs(bat.h[j,1] - a$a[i]) 
  }
  # The find where that difference is minimized
  d <- which.min(c)
  # And assign that switching probability to the proper position, starting
  # from least to most likely to switch
  b[i] <- bat.h[d,1]
}
# And combine into a dataframe
a <- cbind(a,b)
a <- data.frame(a)

# save
write.csv(bat.h, file= "cluster_switching_probs.csv")

# clear workspace
rm(list=ls())

# get data 
# start with partner switching data as example
d <- 
  read.csv("filepath\\partner_switching.csv") %>% 
  dplyr:: select(bat=actor, switch, latency) %>% 
  # only use latencies within an hour (latency is minutes)
  filter(latency <= 60)

# check number of observations per bat
d %>% 
  group_by(bat) %>% 
  summarize(n= n()) %>% 
  arrange(n) %>% 
  pull(n) %>% 
  range()

# only include bats with at least 100 obs
bats <- 
  d %>% 
  group_by(bat) %>% 
  summarize(n= n()) %>% 
  filter(n>100) %>% 
  pull(bat) 

d <- d %>% filter(bat %in% bats)

# plot logistic regressions with all bats
d %>% 
  ggplot(aes(x=latency, y=switch))+
  geom_point(size=2)+
  stat_smooth(aes(color= bat),method="glm", se=FALSE, method.args = list(family=binomial))+
  stat_smooth(method="glm", se=FALSE, method.args = list(family=binomial), color= "black")+
  theme(legend.position= "none") +
  xlab("number of 30-min samples since last observed switch")+
  ylab("probability of switching")

# fit logistic mixed model using all bats as random intercepts
fit <- glmer(switch ~ latency + (1 | bat), data = d, family = binomial)

# get intercepts
random_effects <- 
  coefficients(fit)$bat %>%
  rownames_to_column() %>% 
  dplyr::select(bat= rowname, intercept = '(Intercept)', coefficient= latency)

# get probability of partner switching across each bat separately

# get max number of mins
max.min <-  300

# get bats
bats <- unique(d$bat)

# make empty matrix of bats (rows) by min (cols)
bat.min <- matrix(data= NA, nrow= length(bats), ncol= max.min)

# fill in matrix with predicted probabilities
for (i in 1:nrow(bat.min)) {
  
  # choose bat
  focal <- bats[i]
  
  # get random effects from that bat
  t <- 
    random_effects %>% 
    filter(bat==focal)
  
  # make empty vector for each min 0 to 300
  mins <- rep(NA, max.min)
  # get predicted probability for each day for that bat
  for (j in 1:length(mins)){
    # get response as log odds
    y <- t$coefficient*j + t$intercept
    # get odds  
    odds <- exp(y)
    # get probability
    prob <- odds / (1 + odds)
    mins[j] <- prob  
  }
  
  # then put those values into matrix of probs by bat and h
  bat.min[i,] <- mins
}


# label matrix
colnames(bat.min) <- paste0("min",1:max.min)
rownames(bat.min) <- bats
bat.min <- bat.min[order(bat.min[,1]),]

# The following was used exclusively for correlation of probabilities between
# roost, cluster, and partner switching
# Break first cluster switching probabilities into quantiles
a <- quantile(bat.min[,1], probs = seq(0,1,1/80))
a <- data.frame(a)
# Create a second vector
b <- rep(0,length(a$a))
# For the length of the dataframe
for (i in 1:length(a$a)) {
  # Create a vector for each bat
  c <- rep(0,length(bat.min[,1]))
  # And for each entry in c
  for (j in 1:length(c)) {
    # Find the absolute difference between the quantiles and the actual cluster switching probabilities
    c[j] <- abs(bat.min[j,1] - a$a[i]) 
  }
  # The find where that difference is minimized
  d <- which.min(c)
  # And assign that switching probability to the proper position, starting
  # from least to most likely to switch
  b[i] <- bat.min[d,1]
}
# And combine into a dataframe
a <- cbind(a,b)
a <- data.frame(a)


# plot predictions 
bat.min %>% 
  as.data.frame() %>% 
  rownames_to_column("bat") %>% 
  pivot_longer(cols= min1:min300, names_to= "min", values_to= "prob") %>% 
  mutate(min= as.numeric(substr(min, start=4, stop=15))) %>% 
  ggplot(aes(x=min, y=prob,color=bat))+
  geom_point(size=1)+
  geom_line()+  
  ylab("probability of switch")+
  theme(legend.position= "none")

# save
write.csv(bat.min, file= "partner_switching_probs.csv")

