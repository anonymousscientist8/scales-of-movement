library(asnipe)
library(tidyverse)
library(igraph)
library(stringi)
rm(list = ls())

# Import Data, Provides Roost Switching Rate
W <- read.csv("Filepath\\roosting_information.csv")
W <- W[(is.na(W$Band)==F),]

# All bats in the model
Bats_w <- unique(W$ID)
Bats_w <- Bats_w[Bats_w != ""]

## Social differentiation (roost association) for common vampire bats

# For each bat, count the number of times the bat has been observed
obs_count <- rep(0,length(Bats_w))
for (i in 1:length(Bats_w)) {
  temp <- W[!(W$ID != Bats_w[i]),]
  obs_count[i] <- length(unique(temp$Date))
}

# Count the number of intervals where both were observed in or not in the same group
# For each bat, create a matrix of the bat's position at each interval
for (i in 1:length(Bats_w)) {
  temp <- W[!(W$ID != Bats_w[i]),]
  
  # Create a vector that is used to store number of instances where bats were in the same area
  x_vec <- rep(0,length(Bats_w))
  # Same thing, but instances where they aren't the same
  yab_vec <- rep(0,length(Bats_w))
  
  # and at each interval
  for (j in 1:length(temp$Date)) {
    # Use filtered dataset to get roost position at each time
    temp2 <- W[which(W$Date == temp$Date[j] & W$Roost == temp$Roost[j]),]
    temp3 <- W[which(W$Date == temp$Date[j] & W$Roost != temp$Roost[j]),]
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
sri <- data.frame(matrix(0,sum(Observations >  25),sum(Observations >  25))) # Initialize
m <- 1 # Index
n <- 1 # Index
for (i in 1:length(Bats_w)) { # For every bat
  n <- 1 # Reset n to 1
  for (j in 1:length(Bats_w)) { # And every other bat
    if (Observations[i] > 25 & Observations[j] > 25) { # If both bats have been seen more than 25 times
      sri[m,n] <- x[i,j]/(x[i,j]+yab[i,j]+(obs_count[i]-yab[i,j]-x[i,j])+(obs_count[j]-yab[i,j]-x[i,j])) # Find the sri
      n <- n + 1 # And add 1 to n
    }
  }
  if (Observations[i] > 25) { # If only the focal bat was seen more than 25 times
    m <- m + 1 # Add 1 to m
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

# Export the roost social differentiation
roost_soc_diff <- soc_diff
write.csv(roost_soc_diff, "Filepath\\roost_soc_diff.csv")

# Now we perform a permutation test to test for significance
# Get group-by-individual matrix
W <- W[(is.na(W$Band)==F),] # Filter
roosts <- unique(W$Roost) # Get all roosts
times <- unique(W$Date) # Get all dates
for (i in 1:length(times)) { # For all dates
  for (j in 1:length(roosts)) { # And for each roost
    temp2 <- W[!(W$Date != times[i]),] # Filter to only the applicable date
    temp2 <- temp2[!(temp2$Roost != roosts[j]),] # And place
    if (length(temp2$ID) > 0) { # And if there's anything to look at
      temp <- rep(NA,length(Bats_w)+2) # Initialize a vector
      temp[1] <- times[i] # Get the applicable times and roost
      temp[2] <- roosts[j]
      l <- length(temp2$ID) # And find the number of entries applicable
      for (k in 1:l) { # And for each entry
        temp[k+2] <- temp2$ID[k] # Store the ID of the bat
      }
      # Then create a matrix showing where bats are at each time
      if (i == 1 & j == 1) {
        mat <- temp
      } else {
        mat <- rbind(mat,temp)
      }
    }
  }
}
mat <- data.frame(mat) # Turn that matrix into a data frame
rownames(mat) <- NULL # Replace the column names with the date, roost, and bat ID
colnames(mat) <- c("Date","Roost",paste("Bat",1:(ncol(mat)-2),sep = "_"))
mat2 <- apply(mat[,3:ncol(mat)], 1, function(x) match(Bats_w,x)) # Get matrix that shows the presence of each bat
mat2[is.na(mat2)] <- 0 # Turn NA's to 0's
mat2[mat2>0] <- 1 # If there are values greater than 0, mark them as 1
rownames(mat2) <- Bats_w # Name the rows after the bats
colnames(mat2) <- paste('group', 1:ncol(mat2), sep="_") # Name the columns as groups
mat2 <- mat2[which(rowSums(mat2)>25),] # Filter
gbi <- t(mat2) # Transpose to get group-by-individual matrix
adj <- get_network(gbi, data_format="GBI", association_index = "SRI") # Get adjacency matrix
soc_diff2 <- rep(0,ncol(adj)) # Initialize social differentiation vector
for (i in 1:length(adj[1,])) { # Find "Social differentiation" for each bat
  soc_diff2[i] <- sd(adj[,i], na.rm = T)/mean(adj[,i], na.rm = T)
}
soc_diff2 <- data.frame(soc_diff2) # Make a data frame

# making randomized networks
rand.nets <- network_permutation(gbi, data_format = "GBI", permutations = 100, association_index = "SRI")

# Vector social differentition values
exp <- c()
get_cv <- function(x){ sd(x,na.rm=T)/mean(x,na.rm=T)}
for (n in 1:100) {
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
load("Filepath\\clustering_information.RData")
I <- associations2019

# Only consider the primary cage
I <- I[!(I$cage != 'big_cage'),]

# Obtain all bats and remove any that are NA
bats_i <- unique(c(I$A,I$B,I$C,
                   I$D,I$E,I$F,
                   I$G,I$H,I$I,
                   I$J,I$K,I$L,
                   I$M,I$N,I$O,
                   I$P,I$Q,I$R,
                   I$S,I$T,I$U,
                   I$V))
bats_i <- bats_i[!is.na(bats_i)]

# List of unique dates
times_i <- unique(I$period)

# Make my life easier by turning NA's to 0's
I[is.na(I)] <- 0

# For each time
for (i in 1:length(times_i)) {
  # And each bat
  for (j in 1:length(bats_i)) {
    # Create a data frame that is just the bats at each time
    if ((i == 1) & (j == 1)) {
      I2 <- data.frame(time = times_i[i],
                           bat = bats_i[j])
    } else {
      I2_add <- data.frame(time = times_i[i],
                               bat = bats_i[j])
      I2 <- rbind(I2,I2_add)
    }
  }
}

# Then for each combination of day and time
for (i in 1:length(I2$time)) {
  # Filter out the original df to only times that matter
  temp <- I[!(I$period != I2$time[i]),]
  # And the bat we're looking at
  r <- which(temp == I2$bat[i], arr.ind=TRUE)
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
I2 <- cbind(I2,camera)

# Get rid of row names
rownames(I2) <- NULL

# Add a spacer between hours and minutes so it can be converted into a number
for (i in 1:length(I2$time)) {
  stri_sub(I2$time[i],nchar(I2$time[i])-1,nchar(I2$time[i])-2) <- ":"
}

# Convert time into a number divisible by half hours
I2$time <- as.numeric(as.POSIXlt(I2$time, format="%Y.%m.%d_%H:%M"))
I2$time <- I2$time/(60*30)

## Social differentiation (cluster association) for common vampire bats

# For each bat, count the number of times the bat has been observed
obs_count <- rep(0,length(bats_i))
for (i in 1:length(bats_i)) {
  temp <- I2[!(I2$bat != bats_i[i]),]
  temp <- temp[!(temp$camera == 0),]
  obs_count[i] <- length(temp$camera)
}

# Count the number of intervals where both were observed in or not in the same group
# For each bat, create a matrix of the bat's position at each interval
for (i in 1:length(bats_i)) {
  temp <- I2[!(I2$bat != bats_i[i]),]
  temp <- temp[!(temp$camera == 0),]
  
  # Create a vector that is used to store number of instances where bats were in the same area
  x_vec <- rep(0,length(bats_i))
  # Same thing, but instances where they aren't the same
  yab_vec <- rep(0,length(bats_i))
  
  # and at each interval
  for (j in 1:length(temp$time)) {
    # Use filtered dataset to get roost position at each time
    temp2 <- I2[which(I2$time == temp$time[j] & I2$camera == temp$camera[j]),]
    temp3 <- I2[which(I2$time == temp$time[j] & I2$camera != temp$camera[j] & I2$camera != 0),]
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
  temp <- I2[!(I2$bat != bats_i[i]),]
  temp <- temp[!(temp$camera == 0),]
  Observations[i] <- length(temp$bat)
}

# Filter out everything with less than or equal to Y observations
Y <- 100

# Now create SRI association matrix
sri <- data.frame(matrix(0,sum(Observations >  100),sum(Observations >  100))) # Create a matrix with dimensions of bats we care about
m <- 1 # index 1
n <- 1 # index 2
for (i in 1:length(bats_i)) { # For each bat
  n <- 1
  for (j in 1:length(bats_i)) { # Look at each other bat
    if (Observations[i] > 100 & Observations[j] > 100) { # And if both are ones we care about, calculate the sri
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

# Export data frame
cluster_soc_diff <- soc_diff
write.csv(cluster_soc_diff, "Filepath\\cluster_soc_diff.csv")

# Now we perform a permutation test to test for significance
# Get group-by-individual matrix
I3 <- I[!(I$cluster == 0),]
cameras <- c(1,2,3)
z <- 0
for (i in 1:length(times_i)) { # For all times
  for (j in 1:3) { # Check each cluster
    temp2 <- I3[!(I3$period != times_i[i]),] # Filter to each time
    temp2 <- temp2[!(temp2$camera != cameras[j]),] # And camera
    if (length(temp2$A) > 0) { # If entries still exist
      temp <- rep(NA,length(bats_i[Observations > 100])+2) # Make vector for bats presence
      temp[1] <- times_i[i] # With the first entry as the time
      temp[2] <- cameras[j] # and the second as the camera
      l <- unique(as.vector(as.matrix(temp2[,9:30])),na.rm = T) # Find all unique bats in that matrix
      m <- l[l != 0] # Filter out zeros
      m <- m[m != "alexi"] # And alexi
      for (n in 1:30) { # Mark where all bats were found
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
mat <- data.frame(mat) # Make into data frame
mat2 <- apply(mat[,3:ncol(mat)], 1, function(x) match(bats_i,x))  # Get matrix that shows the presence of each bat
mat2[is.na(mat2)] <- 0 # Turn NA's to 0's
mat2[mat2>0] <- 1 # Mark presence as 1
rownames(mat2) <- bats_i # Change rownames to bat names
colnames(mat2) <- paste('group', 1:ncol(mat2), sep="_") # And column names to groups
mat2 <- mat2[which(rowSums(mat2)>100),] # Filter
gbi <- t(mat2) # Get group-by-individual matrix
adj <- get_network(gbi, data_format="GBI", association_index = "SRI") # Get adjacency matrix
soc_diff2 <- rep(0,ncol(adj)) # Initialize social differentiation vector
for (i in 1:length(adj[1,])) { # Find social differentiation
  soc_diff2[i] <- sd(adj[,i], na.rm = T)/mean(adj[,i], na.rm = T)
}
soc_diff2 <- data.frame(soc_diff2) # Make data frame

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
load("Filepath\\events.RData")
I <- events2019

# Filter out non-grooming events and events not in the big cage
I <- I[!(I$cage != 'big_cage'),]
I <- I[!(I$behav != 'g'),]

# Obtain all bats and remove any that are NA
bats_i <- unique(c(I$actor,I$receiver))
bats_i <- bats_i[!is.na(bats_i)]

# Modify time stamps to include minutes
for (i in 1:length(I$duration)) {
  # Change the string to include minutes
  if (nchar(I$min.start[i]) > 1) {
    substring(I$period[i], nchar(I$period[i])-1, nchar(I$period[i])) <- as.character(I$min.start[i])
  } else {
    substring(I$period[i], nchar(I$period[i]), nchar(I$period[i])) <- as.character(I$min.start[i])
  }
  # Add ':' between minutes and hours and minutes and seconds
  stri_sub(I$period[i],nchar(I$period[i])-1,nchar(I$period[i])-2) <- ":"
  stri_sub(I$period[i],nchar(I$period[i])+1,nchar(I$period[i])+1) <- ":"
  
  # Add seconds to the end
  if (nchar(I$sec.start[i]) > 1) {
    stri_sub(I$period[i],nchar(I$period[i])+1,nchar(I$period[i])+2) <-  as.character(I$sec.start[i])   
  } else {
    stri_sub(I$period[i],nchar(I$period[i])+1,nchar(I$period[i])+1) <-  "0"
    stri_sub(I$period[i],nchar(I$period[i])+2,nchar(I$period[i])+2) <-  as.character(I$sec.start[i])
  }
}

# Convert to seconds
I$period <- as.numeric(as.POSIXlt(I$period, format="%Y.%m.%d_%H:%M:%S"))

# Order by time
I <- I[order(I$period),]

## Social differentiation (grooming interaction) for common vampire bats
# Create a dataframe of actors, order of action, receivers, and grooming duration
# For each bat
for (i in 1:length(bats_i)) {
  # Filter out the non-focal bats
  temp <- I[!(I$actor != bats_i[i]),]
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
  temp <- I[!(I$actor != bats_i[i]),]
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
write.csv(partner_soc_diff, "Filepath\\partner_soc_diff.csv")

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



#########################################################
library(asnipe)
library(tidyverse)
library(igraph)
library(stringi)
rm(list = ls())

a_vec <- rep(0,50)
p_vec <- rep(0,50)
obs_vec <- rep(0,50)
meanexp_vec <- rep(0,50)

# Randomly select 5 different simulations of each type to test for social differentiation
for (x in 1:50) {
  if (x <= 5)
  {
    a <- sample(100:199, 1) # so picks a random test between 100 and 199
  } else {
    if (x <= 10) {
      a <- sample(200:299, 1)
    } else {
      if (x <= 15) {
        a <- sample(300:399, 1)
      } else {
        if (x <= 20) {
          a <- sample(400:499, 1)
        } else {
          if (x <= 25) {
            a <- sample(500:599, 1)
          } else {
            if (x <= 30) {
              a <- sample(600:699, 1)
            } else {
              if (x <= 35) {
                a <- sample(700:799, 1)
              } else {
                if (x <= 40) {
                  a <- sample(800:899, 1)
                } else {
                  if (x <= 45) {
                    a <- sample(900:999, 1)
                  } else {
                    a <- sample(0:99, 1)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  # Get the interactions
  t <-
    read.table(gsub(
      ' ',
      '',
      paste(
        'Filepath\\interactions',
        as.character(a),
        '.csv'
      )
    ), sep = ',', fill = T)
  # And associations
  y <-
    read.table(gsub(
      ' ',
      '',
      paste(
        'Filepath\\associations',
        as.character(a),
        '.csv'
      )
    ), sep = ',', fill = T)
  
  # get bat attributes
  bats <- t[, 1:4]
  sort(bats$V1)
  bats[, 1] <- paste0("bat", bats$V1)
  
  
  # label bat attributes
  #"So, the first column of both files I sent you should be the ID of the bats. The next column should be that bat's roost switching slope, then intercept, then cluster switching slope, then intercept, then partner switching slope and intercept."
  colnames(bats) <- c("id", "rs", "cs", "ps")
  
  # get last column name
  lastcol <- colnames(t)[ncol(t)]
  
  # get grooming events
  d <-
    # get grooming partners
    t[, 5:ncol(t)] %>%
    # add groomer
    mutate(id = bats$id) %>%
    # convert wide to long
    pivot_longer(cols = V5:lastcol,
                 names_to = "order",
                 values_to = "partner") %>%
    # label order of events
    mutate(order = as.numeric(str_remove(order, "V"))) %>%
    # start at 1 not 8
    mutate(order = order - 4) %>%
    # sort by bat and order of event
    arrange(id, order) %>%
    # delete parentheses in partner name
    mutate(partner = str_replace_all(partner, pattern = "[()]", replacement = "")) %>%
    # relabel turtles as bats
    mutate(partner = str_replace_all(partner, pattern = "turtle ", replacement = "bat")) %>%
    mutate(duration = 1) %>%
    # delete NAs
    filter(!is.na(partner))
  
  d <- d[!(d$partner == ""), ]
  
  # inspect events
  unique(d$id)
  unique(d$partner)
  unique(d$order)
  
  # make grooming network-----
  # summarize event data into an edgelist of directed interaction rates
  el <-
    d %>%
    group_by(id, partner) %>%
    summarize(duration = sum(duration)) %>%
    ungroup()
  
  # convert the edgelist into a graph object
  g <- graph_from_data_frame(el)
  
  # convert the graph into a sociomatrix
  m <- as_adjacency_matrix(g, attr = 'duration', sparse = F)
  
  
  # And find the coefficient of variation for the grooming matrix (social differentiation)
  soc_diff <- rep(0, length(bats[, 1]))
  n <- t(m)
  for (i in 1:length(m[1, ])) {
    soc_diff[i] <- sd(n[, i], na.rm = T) / mean(n[, i], na.rm = T)
  }
  soc_diff_avg <- mean(soc_diff, na.rm = T)
  soc_diff <- data.frame(soc_diff)
  
  # Plot the social differentiation compared to the number of observations
  soc_diff_avg <- mean(soc_diff$soc_diff, na.rm = T)
  ggplot() +
    geom_histogram(
      data = soc_diff,
      aes(x = soc_diff),
      bins = 15,
      colour = "black",
      fill = "light blue"
    ) +
    geom_vline(aes(xintercept = soc_diff_avg),
               color = "red",
               linetype = "dashed") +
    xlab("Social Differentiation")
  
  # Resample function
  resample <- function(x, ...)
    x[sample.int(length(x), ...)]
  
  psdm <- rep(0, 100)
  for (i in 1:100) { # For each permutation test
    d2 <- d # Same the data file
    
    # And add a column that denotes row number
    d2 <- d2 %>%
      mutate(original_row = row_number())
    
    permutations <- 10000
    for (k in 1:permutations) {
      # Select a random row
      row_num1 <- sample(length(d2$partner), 1) # Select random row number
      chosen_id1 <- d2$id[row_num1] # Pull out donor
      chosen_partner1 <- d2$partner[row_num1] # Pull out recipient
      d3 <- d2
      row_num2 <- sample(d3$original_row, 1) # Select second row number
      while ((d2$partner[row_num2] == chosen_id1) ||
             d2$id[row_num2] == chosen_partner1) {
        # Ensure groomers and donators aren't the same
        row_num2 <- sample(d3$original_row, 1)
      }
      # Swap recipients
      chosen_partner2 <- d2$partner[row_num2]
      chosen_id2 <- d2$id[row_num2]
      d2$partner[row_num1] <- chosen_partner2
      d2$partner[row_num2] <- chosen_partner1
    }
    
    # summarize event data into an edgelist of directed interaction rates
    el <-
      d2 %>%
      group_by(id, partner) %>%
      summarize(duration = sum(duration)) %>%
      ungroup()
    
    # convert the edgelist into a graph object
    g <- graph_from_data_frame(el)
    
    # convert the graph into a sociomatrix
    m <- as_adjacency_matrix(g, attr = 'duration', sparse = F)
    
    # And find the coefficient of variation for the grooming matrix (social differentiation)
    temp <- rep(0, length(bats[, 1]))
    n <- t(m)
    for (k in 1:length(m[1, ])) {
      temp[k] <- sd(n[, k], na.rm = T) / mean(n[, k], na.rm = T)
    }
    temp <- mean(temp, na.rm = T)
    temp <- data.frame(temp)
    psdm[i] <- mean(temp$temp, na.rm = T)
    show(i)
  }
  
  # p-value
  exp <- psdm
  obs <- soc_diff_avg
  p <- mean(exp >= obs)
  a_vec[x] <- a
  p_vec[x] <- p
  obs_vec[x] <- obs
  meanexp_vec[x] <- mean(exp)
}

# Graph of permutations and social differentiation score, with line showing observed differentiation
data.frame(exp=exp, permutation = 1:length(exp)) %>%
ggplot() + 
  geom_histogram(aes(x = exp), bins = 20, colour = "black", fill = "light blue") +
  geom_vline(aes(xintercept = obs), color = "red", linetype = "dashed") +
  xlab("Social Differentiation")

###############################################
# Fig 6 plot
library(asnipe)
library(tidyverse)
library(igraph)
library(stringi)
rm(list = ls())

a <- 2600
t1 <- read.table(gsub(' ', '', paste('Filepath\\interactions',as.character(a),'.csv')), sep = ',', fill = T)
t2 <- read.table(gsub(' ', '', paste('Filepath\\interactions-roost',as.character(a),'.csv')), sep = ',', fill = T)
t3 <- read.table(gsub(' ', '', paste('Filepath\\interactions-time',as.character(a),'.csv')), sep = ',', fill = T)
t1 <- t1[order(t1$V1),]
t2 <- t2[order(t2$V1),]
t3 <- t3[order(t3$V1),]

# get bat attributes
bats <- t1[,1:4]
sort(bats$V1)
bats[,1] <- paste0("bat", bats$V1)


# label bat attributes
#"So, the first column of both files I sent you should be the ID of the bats. The next column should be that bat's roost switching slope, then intercept, then cluster switching slope, then intercept, then partner switching slope and intercept."
colnames(bats) <- c("id", "rs", "cs", "ps")

# get last column name
lastcol <- colnames(t1)[ncol(t1)]

# get grooming events
d <- data.frame()
n <- 0
for (j in 1:length(t1[,1])) {
  for (i in 5:length(t1[1,])) {
    n <- n + 1
    d[n,1] <- t1[j,1]
    d[n,2] <- t1[j,i]
    d[n,3] <- t2[j,i]
    d[n,4] <- t3[j,i]
    d[n,5] <- 1
  }
}
colnames(d) <- c('id','partner','roost','day','duration')
d <- d %>%
  # delete parentheses in partner name
  mutate(partner= str_replace_all(partner, pattern= "[()]", replacement= "")) %>% 
  # relabel turtles as bats
  mutate(partner= str_replace_all(partner, pattern= "turtle ", replacement= "bat"))
d$id <- paste0("bat", d$id)
d <- d[!(d$partner == ""),]

# Create the group column
# if generating c
d$group <- paste(d$day/24, d$roost/4, sep = "_")
# if generating b
#d$group <- paste(floor(d$day), ceiling(d$roost/4), sep = "_")
# if generating a
#d$group <- 1

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


# And find the coefficient of variation for the grooming matrix (social differentiation)
soc_diff <- rep(0,length(bats[,1]))
n <- t(m)
for (i in 1:length(m[1,])) {
  soc_diff[i] <- sd(n[,i], na.rm = T)/mean(n[,i], na.rm = T)
}
soc_diff_avg <- mean(soc_diff,na.rm = T)
soc_diff <- data.frame(soc_diff)

# Plot the social differentiation compared to the number of observations
soc_diff_avg <- mean(soc_diff$soc_diff, na.rm = T)
ggplot() + 
  geom_histogram(data = soc_diff, aes(x = soc_diff), bins = 15, colour = "black", fill = "light blue") +
  geom_vline(aes(xintercept = soc_diff_avg), color = "red", linetype = "dashed") +
  xlab("Social Differentiation")


resample <- function(x, ...)
  x[sample.int(length(x), ...)]

psdm <- rep(0, 100) # For permutation test social differentiation
for (i in 1:100) {
  d2 <- d
  d2 <- d2 %>%
    mutate(original_row = row_number())
  
  all_groups <- unique(d2$group) # Get all unique groups
  
  permutations <- 100000
  for (k in 1:permutations) {
    # Select a random row
    row_num1 <- sample(length(d2$partner),1) # Select row number
    chosen_id1 <- d2$id[row_num1] # Pull out donor
    chosen_partner1 <- d2$partner[row_num1] # pull out recipient
    chosen_group1 <- d2$group[row_num1] # pull out group
    d3 <- d2[d2$group == chosen_group1,] # filter so we only care about the group in question
    row_num2 <- sample(d3$original_row,1) # sample a new row in that subsection
    while ((d2$partner[row_num2] == chosen_id1) || d2$id[row_num2] == chosen_partner1) {
      # Keep resampling until the groomers and recipients are not the same (can't groom self)
      row_num2 <- sample(d3$original_row,1)
    }
    # then swap partners
    chosen_partner2 <- d2$partner[row_num2]
    chosen_id2 <- d2$id[row_num2]
    d2$partner[row_num1] <- chosen_partner2
    d2$partner[row_num2] <- chosen_partner1
  }
  
  # Make edge list
  el <-
    d2 %>%
    group_by(id, partner) %>%
    summarize(duration = sum(duration)) %>%
    ungroup()
  
  # convert the edgelist into a graph object
  g <- graph_from_data_frame(el)
  
  # convert the graph into a sociomatrix
  m <- as_adjacency_matrix(g, attr = 'duration', sparse = F)
  
  # And find the coefficient of variation for the grooming matrix (social differentiation)
  temp <- rep(0, length(bats[, 1]))
  n <- t(m)
  for (k in 1:length(m[1, ])) {
    temp[k] <- sd(n[, k], na.rm = T) / mean(n[, k], na.rm = T)
  }
  temp <- mean(temp, na.rm = T)
  temp <- data.frame(temp)
  psdm[i] <- mean(temp$temp, na.rm = T)
  show(i)
}

# p-value
exp <- psdm
obs <- soc_diff_avg
p <- mean(exp>=obs)
p

exp <- data.frame(exp)
obs2600c <- obs
write.csv(exp, "Filepath\\exp2600c.csv")

# Graph of permutations and social differentiation score, with line showing observed differentiation
data.frame(exp=exp, permutation = 1:length(exp)) %>%
  ggplot() + 
  geom_histogram(aes(x = exp), bins = 100, colour = "black", fill = "light blue") +
  geom_vline(aes(xintercept = obs), color = "red", linetype = "dashed") +
  xlab("Social Differentiation")

# Export observation data
obs <- data.frame(c(obs2600a,obs2600a,obs2600a))
colnames(obs) <- "obs"
write.csv(obs, "Filepath\\obs2.csv")

#########################
# For creating plot 6

## Randomly selected base simulation
library(asnipe)
library(tidyverse)
library(igraph)
library(stringi)
rm(list = ls())

a <- 548
t <- read.table(gsub(' ', '', paste('Filepath\\interactions',as.character(a),'.csv')), sep = ',', fill = T)
y <- read.table(gsub(' ', '', paste('Filepath\\associations',as.character(a),'.csv')), sep = ',', fill = T)

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


# And find the coefficient of variation for the grooming matrix (social differentiation)
soc_diff <- rep(0,length(bats[,1]))
n <- t(m)
for (i in 1:length(m[1,])) {
  soc_diff[i] <- sd(n[,i], na.rm = T)/mean(n[,i], na.rm = T)
}
soc_diff_avg <- mean(soc_diff,na.rm = T)
soc_diff <- data.frame(soc_diff)

# Plot the social differentiation compared to the number of observations
soc_diff_avg <- mean(soc_diff$soc_diff, na.rm = T)
ggplot() + 
  geom_histogram(data = soc_diff, aes(x = soc_diff), bins = 15, colour = "black", fill = "light blue") +
  geom_vline(aes(xintercept = soc_diff_avg), color = "red", linetype = "dashed") +
  xlab("Social Differentiation")

resample <- function(x, ...) x[sample.int(length(x), ...)]

# Permutation test (randomize within cage, because I don't know who is in what cluster)
psdm <- rep(0,100)
for (i in 1:100) {
  d2 <- d
  d2 <- d2 %>%
    mutate(original_row = row_number())
  
  permutations <- 100000
  for (k in 1:permutations) {
    # Select a random row
    row_num1 <- sample(length(d2$partner),1) # Find the row number
    chosen_id1 <- d2$id[row_num1] # Find the ID of the donor
    chosen_partner1 <- d2$partner[row_num1] # Find the ID of the recipient
    d3 <- d2 # Store a backup
    row_num2 <- sample(d3$original_row,1) # Sample a second row
    while ((d2$partner[row_num2] == chosen_id1) || d2$id[row_num2] == chosen_partner1) {
      # Continue to resample until neither partner matches either recipient
      row_num2 <- sample(d3$original_row,1)
    }
    # Then swap partners
    chosen_partner2 <- d2$partner[row_num2]
    chosen_id2 <- d2$id[row_num2]
    d2$partner[row_num1] <- chosen_partner2
    d2$partner[row_num2] <- chosen_partner1
  }
  
  # summarize event data into an edgelist of directed interaction rates
  el <- 
    d2 %>% 
    group_by(id, partner) %>% 
    summarize(duration= sum(duration)) %>% 
    ungroup()
  
  # convert the edgelist into a graph object
  g <- graph_from_data_frame(el)
  
  # convert the graph into a sociomatrix
  m <- as_adjacency_matrix(g, attr= 'duration', sparse=F)
  
  # And find the coefficient of variation for the grooming matrix (social differentiation)
  temp <- rep(0,length(bats[,1]))
  n <- t(m)
  for (k in 1:length(m[1,])) {
    temp[k] <- sd(n[,k], na.rm = T)/mean(n[,k], na.rm = T)
  }
  temp <- mean(temp,na.rm = T)
  temp <- data.frame(temp)
  psdm[i] <- mean(temp$temp, na.rm = T)
  show(i)
}

# p-value
exp <- psdm
obs <- soc_diff_avg
p <- mean(exp>=obs)
p

# Export this as a data frame
exp <- data.frame(exp)
obs548 <- obs
write.csv(exp, "Filepath\\exp548.csv")


## Remove Hierarchically Embedded Scales of Movement
a <- 2400
t <- read.table(gsub(' ', '', paste('Filepath\\interactions',as.character(a),'.csv')), sep = ',', fill = T)
y <- read.table(gsub(' ', '', paste('Filepath\\associations',as.character(a),'.csv')), sep = ',', fill = T)

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


# And find the coefficient of variation for the grooming matrix (social differentiation)
soc_diff <- rep(0,length(bats[,1]))
n <- t(m)
for (i in 1:length(m[1,])) {
  soc_diff[i] <- sd(n[,i], na.rm = T)/mean(n[,i], na.rm = T)
}
soc_diff_avg <- mean(soc_diff,na.rm = T)
soc_diff <- data.frame(soc_diff)

# Plot the social differentiation compared to the number of observations
soc_diff_avg <- mean(soc_diff$soc_diff, na.rm = T)
ggplot() + 
  geom_histogram(data = soc_diff, aes(x = soc_diff), bins = 15, colour = "black", fill = "light blue") +
  geom_vline(aes(xintercept = soc_diff_avg), color = "red", linetype = "dashed") +
  xlab("Social Differentiation")

resample <- function(x, ...) x[sample.int(length(x), ...)]

# Permutation test (randomize within cage, because I don't know who is in what cluster)
psdm <- rep(0,100)
for (i in 1:100) {
  d2 <- d
  d2 <- d2 %>%
    mutate(original_row = row_number())
  
  permutations <- 100000
  for (k in 1:permutations) {
    # Select a random row
    row_num1 <- sample(length(d2$partner),1) # Find the row number
    chosen_id1 <- d2$id[row_num1] # Find the ID of the donor
    chosen_partner1 <- d2$partner[row_num1] # Find the ID of the recipient
    d3 <- d2 # Store a backup
    row_num2 <- sample(d3$original_row,1) # Sample a second row
    while ((d2$partner[row_num2] == chosen_id1) || d2$id[row_num2] == chosen_partner1) {
      # Continue to resample until neither partner matches either recipient
      row_num2 <- sample(d3$original_row,1)
    }
    # Then swap partners
    chosen_partner2 <- d2$partner[row_num2]
    chosen_id2 <- d2$id[row_num2]
    d2$partner[row_num1] <- chosen_partner2
    d2$partner[row_num2] <- chosen_partner1
  }
  
  # summarize event data into an edgelist of directed interaction rates
  el <- 
    d2 %>% 
    group_by(id, partner) %>% 
    summarize(duration= sum(duration)) %>% 
    ungroup()
  
  # convert the edgelist into a graph object
  g <- graph_from_data_frame(el)
  
  # convert the graph into a sociomatrix
  m <- as_adjacency_matrix(g, attr= 'duration', sparse=F)
  
  # And find the coefficient of variation for the grooming matrix (social differentiation)
  temp <- rep(0,length(bats[,1]))
  n <- t(m)
  for (k in 1:length(m[1,])) {
    temp[k] <- sd(n[,k], na.rm = T)/mean(n[,k], na.rm = T)
  }
  temp <- mean(temp,na.rm = T)
  temp <- data.frame(temp)
  psdm[i] <- mean(temp$temp, na.rm = T)
  show(i)
}

# p-value
exp <- psdm
obs <- soc_diff_avg
p <- mean(exp>=obs)
p

# Export this as a data frame
exp <- data.frame(exp)
obs2400 <- obs
write.csv(exp, "Filepath\\exp2400.csv")

## Remove individual variation in partner switching
a <- 2401
t <- read.table(gsub(' ', '', paste('Filepath\\interactions',as.character(a),'.csv')), sep = ',', fill = T)
y <- read.table(gsub(' ', '', paste('Filepath\\associations',as.character(a),'.csv')), sep = ',', fill = T)

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


# And find the coefficient of variation for the grooming matrix (social differentiation)
soc_diff <- rep(0,length(bats[,1]))
n <- t(m)
for (i in 1:length(m[1,])) {
  soc_diff[i] <- sd(n[,i], na.rm = T)/mean(n[,i], na.rm = T)
}
soc_diff_avg <- mean(soc_diff,na.rm = T)
soc_diff <- data.frame(soc_diff)

# Plot the social differentiation compared to the number of observations
soc_diff_avg <- mean(soc_diff$soc_diff, na.rm = T)
ggplot() + 
  geom_histogram(data = soc_diff, aes(x = soc_diff), bins = 15, colour = "black", fill = "light blue") +
  geom_vline(aes(xintercept = soc_diff_avg), color = "red", linetype = "dashed") +
  xlab("Social Differentiation")

resample <- function(x, ...) x[sample.int(length(x), ...)]

# Permutation test
psdm <- rep(0,100)
for (i in 1:100) {
  d2 <- d
  d2 <- d2 %>%
    mutate(original_row = row_number())
  
  permutations <- 100000
  for (k in 1:permutations) {
    # Select a random row
    row_num1 <- sample(length(d2$partner),1) # Find the row number
    chosen_id1 <- d2$id[row_num1] # Find the ID of the donor
    chosen_partner1 <- d2$partner[row_num1] # Find the ID of the recipient
    d3 <- d2 # Store a backup
    row_num2 <- sample(d3$original_row,1) # Sample a second row
    while ((d2$partner[row_num2] == chosen_id1) || d2$id[row_num2] == chosen_partner1) {
      # Continue to resample until neither partner matches either recipient
      row_num2 <- sample(d3$original_row,1)
    }
    # Then swap partners
    chosen_partner2 <- d2$partner[row_num2]
    chosen_id2 <- d2$id[row_num2]
    d2$partner[row_num1] <- chosen_partner2
    d2$partner[row_num2] <- chosen_partner1
  }
  
  # summarize event data into an edgelist of directed interaction rates
  el <- 
    d2 %>% 
    group_by(id, partner) %>% 
    summarize(duration= sum(duration)) %>% 
    ungroup()
  
  # convert the edgelist into a graph object
  g <- graph_from_data_frame(el)
  
  # convert the graph into a sociomatrix
  m <- as_adjacency_matrix(g, attr= 'duration', sparse=F)
  
  # And find the coefficient of variation for the grooming matrix (social differentiation)
  temp <- rep(0,length(bats[,1]))
  n <- t(m)
  for (k in 1:length(m[1,])) {
    temp[k] <- sd(n[,k], na.rm = T)/mean(n[,k], na.rm = T)
  }
  temp <- mean(temp,na.rm = T)
  temp <- data.frame(temp)
  psdm[i] <- mean(temp$temp, na.rm = T)
  show(i)
}

# p-value
exp <- psdm
obs <- soc_diff_avg
p <- mean(exp>=obs)
p

# Export this as a data frame
exp <- data.frame(exp)
obs2401 <- obs
write.csv(exp, "Filepath\\exp2401.csv")

## Remove byproduct partner fidelity
a <- 2402
t <- read.table(gsub(' ', '', paste('Filepath\\interactions',as.character(a),'.csv')), sep = ',', fill = T)
y <- read.table(gsub(' ', '', paste('Filepath\\associations',as.character(a),'.csv')), sep = ',', fill = T)

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


# And find the coefficient of variation for the grooming matrix (social differentiation)
soc_diff <- rep(0,length(bats[,1]))
n <- t(m)
for (i in 1:length(m[1,])) {
  soc_diff[i] <- sd(n[,i], na.rm = T)/mean(n[,i], na.rm = T)
}
soc_diff_avg <- mean(soc_diff,na.rm = T)
soc_diff <- data.frame(soc_diff)

# Plot the social differentiation compared to the number of observations
soc_diff_avg <- mean(soc_diff$soc_diff, na.rm = T)
ggplot() + 
  geom_histogram(data = soc_diff, aes(x = soc_diff), bins = 15, colour = "black", fill = "light blue") +
  geom_vline(aes(xintercept = soc_diff_avg), color = "red", linetype = "dashed") +
  xlab("Social Differentiation")

resample <- function(x, ...) x[sample.int(length(x), ...)]

# Permutation test
psdm <- rep(0,100)
for (i in 1:100) {
  d2 <- d
  d2 <- d2 %>%
    mutate(original_row = row_number())
  
  # Iterate over d2$id
  permutations <- 100000
  for (k in 1:permutations) {
    # Select a random row
    row_num1 <- sample(length(d2$partner),1) # Find the row number
    chosen_id1 <- d2$id[row_num1] # Find the ID of the donor
    chosen_partner1 <- d2$partner[row_num1] # Find the ID of the recipient
    d3 <- d2 # Store a backup
    row_num2 <- sample(d3$original_row,1) # Sample a second row
    while ((d2$partner[row_num2] == chosen_id1) || d2$id[row_num2] == chosen_partner1) {
      # Continue to resample until neither partner matches either recipient
      row_num2 <- sample(d3$original_row,1)
    }
    # Then swap partners
    chosen_partner2 <- d2$partner[row_num2]
    chosen_id2 <- d2$id[row_num2]
    d2$partner[row_num1] <- chosen_partner2
    d2$partner[row_num2] <- chosen_partner1
  }
  
  # summarize event data into an edgelist of directed interaction rates
  el <- 
    d2 %>% 
    group_by(id, partner) %>% 
    summarize(duration= sum(duration)) %>% 
    ungroup()
  
  # convert the edgelist into a graph object
  g <- graph_from_data_frame(el)
  
  # convert the graph into a sociomatrix
  m <- as_adjacency_matrix(g, attr= 'duration', sparse=F)
  
  # And find the coefficient of variation for the grooming matrix (social differentiation)
  temp <- rep(0,length(bats[,1]))
  n <- t(m)
  for (k in 1:length(m[1,])) {
    temp[k] <- sd(n[,k], na.rm = T)/mean(n[,k], na.rm = T)
  }
  temp <- mean(temp,na.rm = T)
  temp <- data.frame(temp)
  psdm[i] <- mean(temp$temp, na.rm = T)
  show(i)
}

# p-value
exp <- psdm
obs <- soc_diff_avg
p <- mean(exp>=obs)
p

# Export this as a data frame
exp <- data.frame(exp)
obs2402 <- obs
write.csv(exp, "Filepath\\exp2402.csv")

## Most simplified simulation
a <- 2403
t <- read.table(gsub(' ', '', paste('Filepath\\interactions',as.character(a),'.csv')), sep = ',', fill = T)
y <- read.table(gsub(' ', '', paste('Filepath\\associations',as.character(a),'.csv')), sep = ',', fill = T)

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


# And find the coefficient of variation for the grooming matrix (social differentiation)
soc_diff <- rep(0,length(bats[,1]))
n <- t(m)
for (i in 1:length(m[1,])) {
  soc_diff[i] <- sd(n[,i], na.rm = T)/mean(n[,i], na.rm = T)
}
soc_diff_avg <- mean(soc_diff,na.rm = T)
soc_diff <- data.frame(soc_diff)

# Plot the social differentiation compared to the number of observations
soc_diff_avg <- mean(soc_diff$soc_diff, na.rm = T)
ggplot() + 
  geom_histogram(data = soc_diff, aes(x = soc_diff), bins = 15, colour = "black", fill = "light blue") +
  geom_vline(aes(xintercept = soc_diff_avg), color = "red", linetype = "dashed") +
  xlab("Social Differentiation")

resample <- function(x, ...) x[sample.int(length(x), ...)]

# Permutation test (randomize within cage, because I don't know who is in what cluster)
psdm <- rep(0,100)
for (i in 1:100) {
  d2 <- d
  d2 <- d2 %>%
    mutate(original_row = row_number())
  
  # Iterate over d2$id
  permutations <- 100000
  for (k in 1:permutations) {
    # Select a random row
    row_num1 <- sample(length(d2$partner),1) # Find the row number
    chosen_id1 <- d2$id[row_num1] # Find the ID of the donor
    chosen_partner1 <- d2$partner[row_num1] # Find the ID of the recipient
    d3 <- d2 # Store a backup
    row_num2 <- sample(d3$original_row,1) # Sample a second row
    while ((d2$partner[row_num2] == chosen_id1) || d2$id[row_num2] == chosen_partner1) {
      # Continue to resample until neither partner matches either recipient
      row_num2 <- sample(d3$original_row,1)
    }
    # Then swap partners
    chosen_partner2 <- d2$partner[row_num2]
    chosen_id2 <- d2$id[row_num2]
    d2$partner[row_num1] <- chosen_partner2
    d2$partner[row_num2] <- chosen_partner1
  }
  
  # summarize event data into an edgelist of directed interaction rates
  el <- 
    d2 %>% 
    group_by(id, partner) %>% 
    summarize(duration= sum(duration)) %>% 
    ungroup()
  
  # convert the edgelist into a graph object
  g <- graph_from_data_frame(el)
  
  # convert the graph into a sociomatrix
  m <- as_adjacency_matrix(g, attr= 'duration', sparse=F)
  
  # And find the coefficient of variation for the grooming matrix (social differentiation)
  temp <- rep(0,length(bats[,1]))
  n <- t(m)
  for (k in 1:length(m[1,])) {
    temp[k] <- sd(n[,k], na.rm = T)/mean(n[,k], na.rm = T)
  }
  temp <- mean(temp,na.rm = T)
  temp <- data.frame(temp)
  psdm[i] <- mean(temp$temp, na.rm = T)
  show(i)
}

# p-value
exp <- psdm
obs <- soc_diff_avg
p <- mean(exp>=obs)
p

# Export this as a data frame
exp <- data.frame(exp)
obs2403 <- obs
write.csv(exp, "Filepath\\exp2403.csv")

# Export observation data
obs <- data.frame(c(obs548,obs2400,obs2401,obs2402,obs2403))
colnames(obs) <- "obs"
write.csv(obs, "Filepath\\obs.csv")
