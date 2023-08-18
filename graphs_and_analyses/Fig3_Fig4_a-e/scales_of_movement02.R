# Predicting grooming centrality from scales of movement
rm(list = ls())

# load packages
library(igraph)
library(tidyverse)

# Make 3 x 100 data frames
degree_est <- data.frame(matrix(ncol=100,nrow=3))
pagerank_est <- data.frame(matrix(ncol=100,nrow=3))

for (a in 0000:0099) { # For the specified range
print(a)
# get raw data
t <- read.table(gsub(' ', '', paste('C:\\Users\\raven\\Documents\\interactions',as.character(a),'.csv')), sep = ',', fill = T)

# get bat attributes
bats <- t[,1:4]
sort(bats$V1)
bats[,1] <- paste0("bat", bats$V1)


# label bat attributes
#"So, the first column of both files I sent you should be the ID of the bats. The next column should be that bat's roost switching slope, then intercept, then cluster switching slope, then intercept, then partner switching slope and intercept."
colnames(bats) <- c("id", "rs", "cs", "ps")

# get last column name
lastcol <- colnames(t)[ncol(t)]

# get grooming events
d <- 
  # get grooming partners
  t[,5:ncol(t)] %>% 
  # add groomer
  mutate(id= bats$id) %>% 
  # convert wide to long
  pivot_longer(cols= V5:lastcol, names_to= "order", values_to= "partner") %>% 
  # label order of events
  mutate(order= as.numeric(str_remove(order, "V"))) %>% 
  # start at 1 not 8
  mutate(order= order -4) %>% 
  # sort by bat and order of event
  arrange(id, order) %>% 
  # delete parentheses in partner name
  mutate(partner= str_replace_all(partner, pattern= "[()]", replacement= "")) %>% 
  # relabel turtles as bats
  mutate(partner= str_replace_all(partner, pattern= "turtle ", replacement= "bat")) %>% 
  mutate(duration= 1) %>% 
  # delete NAs
  filter(!is.na(partner))

# Remove blanks
d <- d[!(d$partner == ""),]

# inspect events
unique(d$id)
unique(d$partner)
unique(d$order)

# make grooming network-----
# summarize event data into an edgelist of directed interaction rates
el <- 
  d %>% 
  group_by(id, partner) %>% 
  summarize(duration= sum(duration)) %>% 
  ungroup()

# convert the edgelist into a graph object
g <- graph_from_data_frame(el)

# convert the graph into a sociomatrix
m <- as_adjacency_matrix(g, attr= 'duration', sparse=F)

# inspect sociomatrix
m[1:10, 1:10]

# OPTIONAL
plot_network <- TRUE
# plot the network in a way that looks nice.
if(plot_network){
  plot(g, 
       layout= layout_with_fr(g),
       vertex.shape= 'circle', 
       vertex.size= 15,
       vertex.color= 'orange', 
       vertex.label.color= 'black',
       edge.width=(E(g)$weight),
       edge.arrow.size=0.2,
       vertex.label.cex=0.5,
       vertex.label.dist=0)
}

# update bat attributes with network centrality

# sort bats in order
bats <- bats %>% arrange(id)

# add degree
bats$degree <- 
  degree(g, mode= "out") %>% 
  enframe(name= "id", value= "metric") %>% 
  arrange(id) %>% 
  pull(metric)

# get grooming received matrix
el.r <- 
  d %>% 
  group_by(partner, id) %>% 
  summarize(duration= sum(duration)) %>% 
  ungroup()

# convert the edgelist into a graph object
g.r <- graph_from_data_frame(el.r)

# convert the graph into a sociomatrix
m.r <- as_adjacency_matrix(g.r, attr= 'duration', sparse=F)

# check that it is the transpose
mean(m.r == t(m))

# add reverse page rank
bats$pagerank <- 
  page_rank(g.r, directed = T)$vector %>% 
  enframe(name= "id", value= "metric") %>% 
  arrange(id) %>% 
  pull(metric)

# what predicts degree centrality?
bats %>% 
  dplyr::select(id, rs, cs, ps, degree, pagerank) %>% 
  pivot_longer(rs:ps, names_to= "type", values_to= "switching_rate") %>% 
  mutate(type= case_when(
    type== "rs" ~ "roost switching",
    type== "cs" ~ "cluster switching",
    type== "ps" ~ "partner switching")) %>% 
  ggplot(aes(x=switching_rate, y=degree))+
  facet_wrap(~type, scales= "free")+
  geom_point()+
  geom_smooth(method= "lm")

# Get coefficients from generalized linear model
fit1 <- summary(glm(degree~ scale(rs)+ scale(cs)+ scale(ps), family = "poisson", data= bats))
overdisp_fun <- function(model) { # Check for dispersion
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(fit1)
fit1 <- coefficients(fit1)
fit1 <- data.frame(fit1) # Make into data frame
degree_est[1,a-2100+1] <- fit1[2,1] # Get effects from roost switching 
degree_est[2,a-2100+1] <- fit1[3,1] # Get effects from cluster switching
degree_est[3,a-2100+1] <- fit1[4,1] # Get effects from partner switching
rownames(degree_est) <- c('scale(rs)','scale(cs)','scale(ps)')

# what predicts reverse pagerank centrality?
bats %>% # Make basic test plots
  dplyr::select(id, rs, cs, ps, degree, pagerank) %>% 
  pivot_longer(rs:ps, names_to= "type", values_to= "switching_rate") %>% 
  mutate(type= case_when(
    type== "rs" ~ "roost switching",
    type== "cs" ~ "cluster switching",
    type== "ps" ~ "partner switching")) %>% 
  ggplot(aes(x=switching_rate, y=pagerank))+
  facet_wrap(~type, scales= "free")+
  geom_point()+
  geom_smooth(method= "lm")

# Same as above, now with pagerank
fit2 <- summary(glm(pagerank~ scale(rs)+ scale(cs)+ scale(ps), family = "poisson", data= bats))
fit2 <- coefficients(fit2)
fit2 <- data.frame(fit2)
pagerank_est[1,a-2100+1] <- fit2[2,1]
pagerank_est[2,a-2100+1] <- fit2[3,1]
pagerank_est[3,a-2100+1] <- fit2[4,1]
rownames(pagerank_est) <- c('scale(rs)','scale(cs)','scale(ps)')
}

# Transpose the pagerank and degree estimates
degree_est <- data.frame(t(degree_est))
pagerank_est <- data.frame(t(pagerank_est))

# Export as csv's
write.csv(degree_est, "C:\\Users\\raven\\Documents\\degree_est21.csv")
write.csv(pagerank_est, "C:\\Users\\raven\\Documents\\pagerank_est21.csv")

# 
for (i in 1:3) {
  sdall <- rep(NA,500)
  if (i == 1) {variation_in <- rep(NA,500)}
  for (a in 0500:5099) {
    # get raw data
    t <- read.table(gsub(' ', '', paste('C:\\Users\\raven\\Documents\\interactions',as.character(a),'.csv')), sep = ',', fill = T)
    sdall[a+1-500] <- sd(t[,i+1])
    if (i == 1) {variation_in[a+1-500] <- "a"}
  }
  for (a in 1000:1099) {
    # get raw data
    t <- read.table(gsub(' ', '', paste('C:\\Users\\raven\\Documents\\interactions',as.character(a),'.csv')), sep = ',', fill = T)
    sdall[a+1] <- sd(t[,i+1-500], na.rm = T)
    if (i == 1) {variation_in[a+1-500] <- "b"}
  }
  for (a in 1100:1199) {
    # get raw data
    t <- read.table(gsub(' ', '', paste('C:\\Users\\raven\\Documents\\interactions',as.character(a),'.csv')), sep = ',', fill = T)
    sdall[a+1] <- sd(t[,i+1-500], na.rm = T)
    if (i == 1) {variation_in[a+1-500] <- "c"}
  }
  for (a in 1200:1299) {
    # get raw data
    t <- read.table(gsub(' ', '', paste('C:\\Users\\raven\\Documents\\interactions',as.character(a),'.csv')), sep = ',', fill = T)
    sdall[a+1] <- sd(t[,i+1-500], na.rm = T)
    if (i == 1) {variation_in[a+1-500] <- "d"}
  }
  for (a in 1300:1399) {
    # get raw data
    t <- read.table(gsub(' ', '', paste('C:\\Users\\raven\\Documents\\interactions',as.character(a),'.csv')), sep = ',', fill = T)
    sdall[a+1] <- sd(t[,i+1-500], na.rm = T)
    if (i == 1) {variation_in[a+1-500] <- "e"}
  }
  if (i == 1) {
    sds <- cbind(variation_in,sdall)
  } else {
    sds <- cbind(sds,sdall)
  }
}
sds <- data.frame(sds)
sds$sdall <- as.numeric(as.character(sds$sdall))
sds$sdall.1 <- as.numeric(as.character(sds$sdall.1))
sds$sdall.2 <- as.numeric(as.character(sds$sdall.2))

label <- c(
  "a" = "(a) Variable, Uncorrelated Switching Rates",
  "b" = "(b) Variable Roost Switching Rates",
  "c" = "(c) Variable Cluster Switching Rates",
  "d" = "(d) Variable Partner Switching Rates",
  "e" = "(e) Variable, Correlated Switching Rates"
)

theme_new <- theme_set(theme_bw())
theme_new <- theme_update(strip.background = element_blank(),strip.text = element_text(hjust = 0))


ggplot(data = sds) +
  geom_histogram(mapping = aes(x = sdall),colour = "black", fill = "light blue", position='identity') +
  facet_wrap(~variation_in, ncol = 1, labeller = as_labeller(label)) +
  xlab("standard deviations from mean roost switching count")
