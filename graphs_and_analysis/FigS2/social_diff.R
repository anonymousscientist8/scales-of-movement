library(asnipe)
library(tidyverse)
library(igraph)
library(stringi)
rm(list = ls())

# Import Wilkinson Data, Provides Roost Switching Rate
wilkinson <- read.csv("filepath\\wilkinson2.csv")
wilkinson <- wilkinson[(is.na(wilkinson$Band)==F),]

# All bats in the model
Bats_w <- unique(wilkinson$ID)
Bats_w <- Bats_w[Bats_w != ""]

## Social differentiation (roost association) for common vampire bats

# For each bat, count the number of times the bat has been observed
obs_count <- rep(0,length(Bats_w))
for (i in 1:length(Bats_w)) {
  temp <- wilkinson[!(wilkinson$ID != Bats_w[i]),]
  obs_count[i] <- length(unique(temp$Date))
}

# Count the number of intervals where both were observed in or not in the same group
# For each bat, create a matrix of the bat's position at each interval
for (i in 1:length(Bats_w)) {
  temp <- wilkinson[!(wilkinson$ID != Bats_w[i]),]
  
  # Create a vector that is used to store number of instances where bats were in the same area
  x_vec <- rep(0,length(Bats_w))
  # Same thing, but instances where they aren't the same
  yab_vec <- rep(0,length(Bats_w))
  
  # and at each interval
  for (j in 1:length(temp$Date)) {
    # Use filtered dataset to get roost position at each time
    temp2 <- wilkinson[which(wilkinson$Date == temp$Date[j] & wilkinson$Roost == temp$Roost[j]),]
    temp3 <- wilkinson[which(wilkinson$Date == temp$Date[j] & wilkinson$Roost != temp$Roost[j]),]
    # And for each bat
    for (k in 1:length(Bats_w)) {
      # If looking at a different bat
      if (i != k) {
        # Add the number of observations where they share the same roost
        x_vec[k] <- x_vec[k] + sum(temp2$ID == Bats_w[k])
        # And different roosts
        yab_vec[k] <- yab_vec[k] + sum(temp3$ID == Bats_w[k])
      }
    }
  }
  # Then create a data frame of all the x_vecs and yab_vecs
  if (i == 1) {
    x <- data.frame(x_vec)
    yab <- data.frame(yab_vec)
  } else {
    x <- cbind(x,x_vec)
    yab <- cbind(yab,yab_vec)
  }
}

# Get the number of observatoins of each bat
Observations <- obs_count

# Filter out everything with less than or equal to Y observations
Y <- 25

# Now create SRI association matrix
sri <- data.frame(matrix(0,sum(Observations >  25),sum(Observations >  25)))
m <- 1
n <- 1
for (i in 1:length(Bats_w)) {
  n <- 1
  for (j in 1:length(Bats_w)) {
    if (Observations[i] > 25 & Observations[j] > 25) {
      sri[m,n] <- x[i,j]/(x[i,j]+yab[i,j]+(obs_count[i]-yab[i,j]-x[i,j])+(obs_count[j]-yab[i,j]-x[i,j]))
      n <- n + 1
    }
  }
  if (Observations[i] > 25) {
    m <- m + 1
  }
}

# And find the coefficient of variation for the sri's (social differentiation)
soc_diff <- rep(0,ncol(sri))
for (i in 1:length(sri[1,])) {
  soc_diff[i] <- sd(sri[,i], na.rm = T)/mean(sri[,i], na.rm = T)
}
soc_diff <- data.frame(soc_diff)
soc_diff <- cbind(soc_diff,Observations[Observations > 25])

# Main Graph
soc_diff_avg <- mean(soc_diff$soc_diff, na.rm = T)

roost_soc_diff <- soc_diff
write.csv(roost_soc_diff, "filepath\\roost_soc_diff.csv")

# Now we perform a permutation test to test for significance
# Get group-by-individual matrix
Wilkinson <- wilkinson[(is.na(wilkinson$Band)==F),]
roosts <- unique(Wilkinson$Roost)
times <- unique(Wilkinson$Date)
for (i in 1:length(times)) {
  for (j in 1:length(roosts)) {
    temp2 <- Wilkinson[!(Wilkinson$Date != times[i]),]
    temp2 <- temp2[!(temp2$Roost != roosts[j]),]
    if (length(temp2$ID) > 0) {
      temp <- rep(NA,length(Bats_w)+2)
      temp[1] <- times[i]
      temp[2] <- roosts[j]
      l <- length(temp2$ID)
      for (k in 1:l) {
        temp[k+2] <- temp2$ID[k]
      }
      if (i == 1 & j == 1) {
        mat <- temp
      } else {
        mat <- rbind(mat,temp)
      }
    }
  }
}
mat <- data.frame(mat)
rownames(mat) <- NULL
colnames(mat) <- c("Date","Roost",paste("Bat",1:(ncol(mat)-2),sep = "_"))
mat2 <- apply(mat[,3:ncol(mat)], 1, function(x) match(Bats_w,x))
mat2[is.na(mat2)] <- 0
mat2[mat2>0] <- 1
rownames(mat2) <- Bats_w
colnames(mat2) <- paste('group', 1:ncol(mat2), sep="_")
mat2 <- mat2[which(rowSums(mat2)>25),]
gbi <- t(mat2)
adj <- get_network(gbi, data_format="GBI", association_index = "SRI")
soc_diff2 <- rep(0,ncol(adj))
for (i in 1:length(adj[1,])) {
  soc_diff2[i] <- sd(adj[,i], na.rm = T)/mean(adj[,i], na.rm = T)
}
soc_diff2 <- data.frame(soc_diff2)

# making randomized networks
rand.nets <- network_permutation(gbi, data_format = "GBI", permutations = 10000, association_index = "SRI")

# Vector social differentition values
exp <- c()
get_cv <- function(x){ sd(x,na.rm=T)/mean(x,na.rm=T)}
for (n in 1:10000) {
  exp <- append(exp, get_cv(rand.nets[n,,]))
}

# p-value
obs <- mean(soc_diff2$soc_diff2,na.rm = T)
p <- mean(exp>=obs)
p

# Graph of permutations and social differentiation score, with line showing observed differentiation
data.frame(exp=exp, permutation = 1:length(exp)) %>%
  ggplot(aes(x=permutation, y=exp)) +
  geom_line() +
  geom_hline(aes(yintercept = obs), color="blue", size = 1) +
  ylab("social differentiation") + theme_bw()

# Cluster / crevice information
rm(list = ls())
load("filepath\\imran.RData")
imran <- associations2019

# Only consider the primary cage
imran <- imran[!(imran$cage != 'big_cage'),]

# Obtain all bats and remove any that are NA
bats_i <- unique(c(imran$A,imran$B,imran$C,
                   imran$D,imran$E,imran$F,
                   imran$G,imran$H,imran$I,
                   imran$J,imran$K,imran$L,
                   imran$M,imran$N,imran$O,
                   imran$P,imran$Q,imran$R,
                   imran$S,imran$T,imran$U,
                   imran$V))
bats_i <- bats_i[!is.na(bats_i)]

# List of unique dates
times_i <- unique(imran$period)

# Make my life easier by turning NA's to 0's
imran[is.na(imran)] <- 0

# For each time
for (i in 1:length(times_i)) {
  # And each bat
  for (j in 1:length(bats_i)) {
    # Create a data frame that is just the bats at each time
    if ((i == 1) & (j == 1)) {
      imran2 <- data.frame(time = times_i[i],
                           bat = bats_i[j])
    } else {
      imran2_add <- data.frame(time = times_i[i],
                               bat = bats_i[j])
      imran2 <- rbind(imran2,imran2_add)
    }
  }
}

# Then for each combination of day and time
for (i in 1:length(imran2$time)) {
  # Filter out the original df to only times that matter
  temp <- imran[!(imran$period != imran2$time[i]),]
  # And the bat we're looking at
  r <- which(temp == imran2$bat[i], arr.ind=TRUE)
  temp <- temp[r[1],]
  # And then either create a vector showing the camera that bat was observed at
  if (i == 1) {
    camera <- temp$camera[1]
    # Or add to the vector
  } else {
    if ((length(temp$camera) == 0) || (is.na(temp$camera[1]) == TRUE)) {
      camera_add <- 0
    } else {
      camera_add <- temp$camera[1]
    }
    camera <- rbind(camera,camera_add)
  }
}
imran2 <- cbind(imran2,camera)

# Get rid of row names
rownames(imran2) <- NULL

# Add a spacer between hours and minutes so it can be converted into a number
for (i in 1:length(imran2$time)) {
  stri_sub(imran2$time[i],nchar(imran2$time[i])-1,nchar(imran2$time[i])-2) <- ":"
}

# Convert time into a number divisible by half hours
imran2$time <- as.numeric(as.POSIXlt(imran2$time, format="%Y.%m.%d_%H:%M"))
imran2$time <- imran2$time/(60*30)

## Social differentiation (cluster association) for common vampire bats

# For each bat, count the number of times the bat has been observed
obs_count <- rep(0,length(bats_i))
for (i in 1:length(bats_i)) {
  temp <- imran2[!(imran2$bat != bats_i[i]),]
  temp <- temp[!(temp$camera == 0),]
  obs_count[i] <- length(temp$camera)
}

# Count the number of intervals where both were observed in or not in the same group
# For each bat, create a matrix of the bat's position at each interval
for (i in 1:length(bats_i)) {
  temp <- imran2[!(imran2$bat != bats_i[i]),]
  temp <- temp[!(temp$camera == 0),]
  
  # Create a vector that is used to store number of instances where bats were in the same area
  x_vec <- rep(0,length(bats_i))
  # Same thing, but instances where they aren't the same
  yab_vec <- rep(0,length(bats_i))
  
  # and at each interval
  for (j in 1:length(temp$time)) {
    # Use filtered dataset to get roost position at each time
    temp2 <- imran2[which(imran2$time == temp$time[j] & imran2$camera == temp$camera[j]),]
    temp3 <- imran2[which(imran2$time == temp$time[j] & imran2$camera != temp$camera[j] & imran2$camera != 0),]
    # And for each bat
    for (k in 1:length(bats_i)) {
      # If looking at a different bat
      if (i != k) {
        # Add the number of observations where they share the same roost
        x_vec[k] <- x_vec[k] + sum(temp2$bat == bats_i[k])
        # And different roosts
        yab_vec[k] <- yab_vec[k] + sum(temp3$bat == bats_i[k])
      }
    }
  }
  # Then create a data frame of all the x_vecs and yab_vecs
  if (i == 1) {
    x <- data.frame(x_vec)
    yab <- data.frame(yab_vec)
  } else {
    x <- cbind(x,x_vec)
    yab <- cbind(yab,yab_vec)
  }
}

# Get the number of observatoins of each bat
Observations <- rep(0,length(bats_i))
for (i in 1:length(bats_i)) {
  temp <- imran2[!(imran2$bat != bats_i[i]),]
  temp <- temp[!(temp$camera == 0),]
  Observations[i] <- length(temp$bat)
}

# Filter out everything with less than or equal to Y observations
Y <- 100

# Now create SRI association matrix
sri <- data.frame(matrix(0,sum(Observations >  100),sum(Observations >  100)))
m <- 1
n <- 1
for (i in 1:length(bats_i)) {
  n <- 1
  for (j in 1:length(bats_i)) {
    if (Observations[i] > 100 & Observations[j] > 100) {
      sri[m,n] <- x[i,j]/(x[i,j]+yab[i,j]+(obs_count[i]-yab[i,j]-x[i,j])+(obs_count[j]-yab[i,j]-x[i,j]))
      n <- n + 1
    }
  }
  if (Observations[i] > 100) {
    m <- m + 1
  }
}

# And find the coefficient of variation for the sri's (social differentiation)
soc_diff <- rep(0,nrow(sri))
for (i in 1:length(sri[1,])) {
  soc_diff[i] <- sd(sri[,i], na.rm = T)/mean(sri[,i], na.rm = T)
}
soc_diff <- data.frame(soc_diff)
soc_diff <- cbind(soc_diff,Observations[Observations > 100])

# Plot the social differentiation compared to the number of observations
soc_diff_avg <- mean(soc_diff$soc_diff, na.rm = T)

cluster_soc_diff <- soc_diff
write.csv(cluster_soc_diff, "filepath\\cluster_soc_diff.csv")

# Now we perform a permutation test to test for significance
# Get group-by-individual matrix
imran3 <- imran[!(imran$cluster == 0),]
cameras <- c(1,2,3)
z <- 0
for (i in 1:length(times_i)) {
  for (j in 1:3) {
    temp2 <- imran3[!(imran3$period != times_i[i]),]
    temp2 <- temp2[!(temp2$camera != cameras[j]),]
    if (length(temp2$A) > 0) {
      temp <- rep(NA,length(bats_i[Observations > 100])+2)
      temp[1] <- times_i[i]
      temp[2] <- cameras[j]
      l <- unique(as.vector(as.matrix(temp2[,9:30])),na.rm = T)
      m <- l[l != 0]
      m <- m[m != "alexi"]
      for (n in 1:30) {
        if (length(m) == 30) {
          temp[2+n] <- m[n]
        } else {
          if (n <= length(m)) {
            temp[2+n] <- m[n]
          } else {
            temp[2+n] <- NA
          }
        }
      }
      if (z == 0 & (length(temp2$A) > 0)) {
        mat <- temp
        z <- 1
      } else {
        mat <- rbind(mat,temp)
      }
    }
  }
}
mat <- data.frame(mat)
mat2 <- apply(mat[,3:ncol(mat)], 1, function(x) match(bats_i,x))
mat2[is.na(mat2)] <- 0
mat2[mat2>0] <- 1
rownames(mat2) <- bats_i
colnames(mat2) <- paste('group', 1:ncol(mat2), sep="_")
mat2 <- mat2[which(rowSums(mat2)>100),]
gbi <- t(mat2)
adj <- get_network(gbi, data_format="GBI", association_index = "SRI")
soc_diff2 <- rep(0,ncol(adj))
for (i in 1:length(adj[1,])) {
  soc_diff2[i] <- sd(adj[,i], na.rm = T)/mean(adj[,i], na.rm = T)
}
soc_diff2 <- data.frame(soc_diff2)

# making randomized networks
rand.nets <- network_permutation(gbi, data_format = "GBI", permutations = 50000, association_index = "SRI")

# Vector social differentition values
exp <- c()
get_cv <- function(x){ sd(x,na.rm=T)/mean(x,na.rm=T)}
for (n in 1:50000) {
  exp <- append(exp, get_cv(rand.nets[n,,]))
}

# p-value
obs <- mean(soc_diff2$soc_diff2,na.rm = T)
p <- mean(exp>=obs)
p

# Graph of permutations and social differentiation score, with line showing observed differentiation
data.frame(exp=exp, permutation = 1:length(exp)) %>%
  ggplot(aes(x=permutation, y=exp)) +
  geom_line() +
  geom_hline(aes(yintercept = obs), color="blue", size = 1) +
  ylab("social differentiation") + theme_bw()

rm(list = ls())

# Load data for analysis
load("filepath\\events.RData")
imran <- events2019

# Filter out non-grooming events and events not in the big cage
imran <- imran[!(imran$cage != 'big_cage'),]
imran <- imran[!(imran$behav != 'g'),]

# Obtain all bats and remove any that are NA
bats_i <- unique(c(imran$actor,imran$receiver))
bats_i <- bats_i[!is.na(bats_i)]

# Modify time stamps to include minutes
for (i in 1:length(imran$duration)) {
  # Change the string to include minutes
  if (nchar(imran$min.start[i]) > 1) {
    substring(imran$period[i], nchar(imran$period[i])-1, nchar(imran$period[i])) <- as.character(imran$min.start[i])
  } else {
    substring(imran$period[i], nchar(imran$period[i]), nchar(imran$period[i])) <- as.character(imran$min.start[i])
  }
  # Add ':' between minutes and hours and minutes and seconds
  stri_sub(imran$period[i],nchar(imran$period[i])-1,nchar(imran$period[i])-2) <- ":"
  stri_sub(imran$period[i],nchar(imran$period[i])+1,nchar(imran$period[i])+1) <- ":"
  
  # Add seconds to the end
  if (nchar(imran$sec.start[i]) > 1) {
    stri_sub(imran$period[i],nchar(imran$period[i])+1,nchar(imran$period[i])+2) <-  as.character(imran$sec.start[i])   
  } else {
    stri_sub(imran$period[i],nchar(imran$period[i])+1,nchar(imran$period[i])+1) <-  "0"
    stri_sub(imran$period[i],nchar(imran$period[i])+2,nchar(imran$period[i])+2) <-  as.character(imran$sec.start[i])
  }
}

# Convert to seconds
imran$period <- as.numeric(as.POSIXlt(imran$period, format="%Y.%m.%d_%H:%M:%S"))

# Order by time
imran <- imran[order(imran$period),]

## Social differentiation (grooming interaction) for common vampire bats
# Create a dataframe of actors, order of action, receivers, and grooming duration
# For each bat
for (i in 1:length(bats_i)) {
  # Filter out the non-focal bats
  temp <- imran[!(imran$actor != bats_i[i]),]
  # And for each row
  k <- 0
  if (length(temp$obs) >= 1) {
    for (j in 1:length(temp$obs)) {
      groom_vec <- rep(0,4)
      # Add the actor
      groom_vec[1] <- bats_i[i]
      # Update the observation number
      k <- k + 1
      groom_vec[2] <- k
      # Add the receiver
      groom_vec[3] <- temp$receiver[j]
      # Add the duration
      groom_vec[4] <- temp$duration[j]
      # If it's the first entry to the data frame, make it. Otherwise, add to it
      if (i == 1 & j == 1) {
        groom_df <- data.frame(t(groom_vec))
      } else {
        groom_df <- rbind(groom_df,groom_vec)
      }
    }
  }
}

# Add rownames
colnames(groom_df) <- c('id','order','partner','duration')

# make grooming network-----
# summarize event data into an edgelist of directed interaction rates
el <- 
  groom_df %>% 
  group_by(id, partner) %>% 
  summarize(duration= sum(as.numeric(duration))) %>% 
  ungroup()

# convert the edgelist into a graph object
g <- graph_from_data_frame(el)

# convert the graph into a sociomatrix
m <- as_adjacency_matrix(g, attr= 'duration', sparse=F)


# And find the coefficient of variation for the grooming matrix (social differentiation)
soc_diff <- rep(0,length(bats_i))
n <- t(m)
for (i in 1:length(m[1,])) {
  soc_diff[i] <- sd(n[,i], na.rm = T)/mean(n[,i], na.rm = T)
}
soc_diff_avg <- mean(soc_diff,na.rm = T)
soc_diff <- data.frame(soc_diff)

# Get the number of observatoins of each bat
Observations <- rep(0,length(bats_i))
for (i in 1:length(bats_i)) {
  temp <- imran[!(imran$actor != bats_i[i]),]
  Observations[i] <- length(temp$actor)
}
soc_diff <- cbind(soc_diff,Observations)

# Filter out everything with less than or equal to Y observations
Y <- 0
bool_filter <- soc_diff$Observations <= Y
soc_diff <- soc_diff[!bool_filter,]

# Plot the social differentiation compared to the number of observations
soc_diff_avg <- mean(soc_diff$soc_diff, na.rm = T)

partner_soc_diff <- soc_diff
write.csv(partner_soc_diff, "filepath\\partner_soc_diff.csv")

# Permutation test (randomize within cage, because I don't know who is in what cluster)
groom_df2 <- groom_df
psdm <- rep(0,1000)
for (i in 1:1000) {
  groom_df2$partner <- resample(groom_df2$partner)
  # make grooming network-----
  # summarize event data into an edgelist of directed interaction rates
  el <- 
    groom_df2 %>% 
    group_by(id, partner) %>% 
    summarize(duration= sum(as.numeric(duration))) %>% 
    ungroup()
    
  # convert the edgelist into a graph object
  g <- graph_from_data_frame(el)
  
  # convert the graph into a sociomatrix
  m <- as_adjacency_matrix(g, attr= 'duration', sparse=F)
    
    
  # And find the coefficient of variation for the grooming matrix (social differentiation)
  temp <- rep(0,length(bats_i))
  n <- t(m)
  for (k in 1:length(m[1,])) {
    temp[k] <- sd(n[,k], na.rm = T)/mean(n[,k], na.rm = T)
  }
  temp <- mean(temp,na.rm = T)
  temp <- data.frame(temp)
  psdm[i] <- mean(temp$temp, na.rm = T)
}

# p-value
exp <- psdm
obs <- soc_diff_avg
p <- mean(exp>=obs)
p

# Graph of permutations and social differentiation score, with line showing observed differentiation
data.frame(exp=exp, permutation = 1:length(exp)) %>%
  ggplot(aes(x=permutation, y=exp)) +
  geom_line() +
  geom_hline(aes(yintercept = obs), color="blue", size = 1) +
  ylab("social differentiation") + theme_bw()
