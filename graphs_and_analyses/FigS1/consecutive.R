library(asnipe)
library(tidyverse)
library(igraph)
library(stringi)
rm(list = ls())

## Distribution of roost switching rates
# Import Wilkinson Data, Provides Roost Switching Rate
wilkinson <- read.csv("filepath\\wilkinson2.csv")
#wilkinson <- wilkinson[(wilkinson$Sex == "F"),]

# All bats in the model
Bats_w <- unique(wilkinson$ID)
Bats_w <- Bats_w[Bats_w != ""]

# For each bat, count the number of times the bat has been observed
obs_count <- rep(0,length(Bats_w))
for (i in 1:length(Bats_w)) {
  temp <- wilkinson[!(wilkinson$ID != Bats_w[i]),]
  obs_count[i] <- length(unique(temp$Date))
}
Bats_w <- Bats_w[obs_count > 25] # Filter out bats with less than 25 observations

# Convert date into something a little easier to work with
wilkinson$Date <- as.numeric(as.POSIXlt(wilkinson$Date, format="%m/%d/%Y")) / 
  (60*60*24)
wilkinson <- wilkinson[!is.na(wilkinson$Date),]

# Count of the number of times a bat was observed switching roosts
# on consecutive days (numerator and denominator)
roost_sw_num <- rep(0,length(obs_count[obs_count > 25]))
roost_sw_den <- rep(0,length(obs_count[obs_count > 25]))

# For each bat
for (i in 1:length(Bats_w)) {
  # Filter out observations of bats
  temp1 <- wilkinson[!(wilkinson$ID != Bats_w[i]),]
  # And if there is more than one occurrence of this bat
  if (length(temp1$Date) > 1) {
    # Use a loop to count the number of roost switches on consecutive days
    for (j in 1:(length(temp1$Date)-1)) {
      # Count number of days between observations
      consecutive <- round(temp1$Date[j+1] - temp1$Date[j])
      # And if the number of days is about 1 (accounting for leap years)
      if (consecutive == 1) {
        # Check whether they switched roosts between days
        # If they are the same, then add 1 to the denominator
        if (temp1$Roost[j+1] == temp1$Roost[j]) {
          roost_sw_den[i] <- roost_sw_den[i] + 1
          # And if they aren't, add 1 to numerator and denominator
        } else {
          roost_sw_den[i] <- roost_sw_den[i] + 1
          roost_sw_num[i] <- roost_sw_num[i] + 1
        }
      }
    }
  }
  
  # Now add instances where the same roost was measured the day after or before
  for (j in 1:length(temp1$Date)) {
    temp2 <- wilkinson[which(abs(wilkinson$Date - (temp1$Date[j]+1)) < 1.5 & wilkinson$Roost == temp1$Roost[j]),]
    roost_sw_den[i] <- roost_sw_den[i] + 1
    roost_sw_num[i] <- roost_sw_num[i] + 1
  }
}

# Find roost switching rates for each bat
roost_sw <- roost_sw_num/roost_sw_den

# Find the mean roost switching rate
roost_sw_avg <- mean(roost_sw)

# Plot roost switching rate
roost_sw <- data.frame(roost_sw)
ggplot(data = roost_sw) +
  geom_histogram(mapping = aes(x = roost_sw), colour = "black", fill = "light blue") +
  xlab("Roost Switches / Day") +
  geom_vline(xintercept = roost_sw_avg, color = 'red', linetype = 'dotted') +
  theme_bw()

# Cluster / crevice information
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
      imran2 <- data.frame(times = times_i[i],
                           bat = bats_i[j])
    } else {
      imran2_add <- data.frame(times = times_i[i],
                               bat = bats_i[j])
      imran2 <- rbind(imran2,imran2_add)
    }
  }
}

# Then for each combination of day and time
for (i in 1:length(imran2$times)) {
  # Filter out the original df to only times that matter
  temp <- imran[!(imran$period != imran2$times[i]),]
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
for (i in 1:length(imran2$times)) {
  stri_sub(imran2$times[i],nchar(imran2$times[i])-1,nchar(imran2$times[i])-2) <- ":"
}

# Convert time into a number divisible by half hours
imran2$times <- as.numeric(as.POSIXlt(imran2$times, format="%Y.%m.%d_%H:%M"))
imran2$times <- imran2$times/(60*30*48)

# Create vectors to determine cluster switching probabilities
cluster_sw_num <- rep(0,length(bats_i))
cluster_sw_den <- rep(0,length(bats_i))

# For each bat
X <- 1/24
for (i in 1:length(bats_i)) {
  # Filter to only known positions of the focal bat
  temp <- imran2[!(imran2$bat != bats_i[i]),]
  temp <- temp[!(temp$camera == 0),]
  # Then, for each row
  if (length(temp$camera) > 1) {
    for (j in 1:(length(temp$camera)-1)) {
      # If the time since last observation is less than X days hours
      if (abs(temp$times[j] - temp$times[j+1]) <= X) {
        # Determine whether the bats switched clusters or not
        if (temp$camera[j] ==  temp$camer[j+1]) {
          cluster_sw_den[i] <- cluster_sw_den[i] + abs(temp$times[j] - temp$times[j+1])
        } else {
          cluster_sw_den[i] <- cluster_sw_den[i] + abs(temp$times[j] - temp$times[j+1])
          cluster_sw_num[i] <- cluster_sw_num[i] + 1
        }
      }
    }
  }
}

# Cluster switching percentage
cluster_sw <- cluster_sw_num/cluster_sw_den

# Find the mean cluster switching rate
cluster_sw_avg <- mean(cluster_sw, na.rm = T)

# Plot cluster switching rate
cluster_sw <- data.frame(cluster_sw)
cluster_sw <- cbind(cluster_sw,bats_i)
ggplot(data = cluster_sw) +
  geom_histogram(mapping = aes(x = cluster_sw), colour = "black", fill = "light blue", binwidth = 0.25) +
  xlab("Cluster Switches / Day") +
  geom_vline(xintercept = cluster_sw_avg, color = 'red', linetype = 'dotted') +
  theme_bw()

## Partner switching between consecutive minutes
# Load data for analysis
load("filepath\\events.RData")
imran <- events2019

# Filter out non-grooming events and events not in the big cage
imran <- imran[!(imran$cage != 'big_cage'),]
imran <- imran[!(imran$behav != 'g'),]
hours <- length(unique(imran$period))

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
imran$period <- as.numeric(as.POSIXlt(imran$period, format="%Y.%m.%d_%H:%M:%S"))/(60*60*24)

# Order by time
imran <- imran[order(imran$period),]

# Create values for determining whether bats switched partners between consecutive bouts (w/in X minutes)
partner_sw_num <- rep(0,length(bats_i))
partner_sw_den <- rep(0,length(bats_i))
partner_sw_num2 <- rep(0,length(bats_i))
partner_sw_den2 <- rep(0,length(bats_i))

# For each bat
X <- 1/24
for (i in 1:length(bats_i)) {
  # Filter out everything except the focal bat
  temp <- imran[!(imran$actor != bats_i[i]),]
  # If there are multiple instances of grooming to look at
  if (length(temp$obs) > 1) {
    # And for each observation
    for (j in 1:(length(temp$actor)-1)) {
      # If the grooming initiation was less than or equal to X minutes from each other
      if (abs(temp$period[j] - temp$period[j+1]) <= X) {
        # Then determine whether the bat switched partners from last time
        if (temp$receiver[j] != temp$receiver[j+1]) {
          partner_sw_den[i] <- partner_sw_den[i] + abs(temp$period[j] - temp$period[j+1])
          partner_sw_num[i] <- partner_sw_num[i] + 1
        } else {
          partner_sw_den[i] <- partner_sw_den[i] + abs(temp$period[j] - temp$period[j+1])
        }
      }
    }
  }
}

# For each bat
X <- 1/24
for (i in 1:length(bats_i)) {
  # Filter out everything except the focal bat
  temp <- imran[!(imran$actor != bats_i[i]),]
  # If there are multiple instances of grooming to look at
  if (length(temp$obs) > 1) {
    # And for each observation
    for (j in 1:(length(temp$actor)-1)) {
      # If the grooming initiation was less than or equal to X minutes from each other
      if (abs(temp$period[j] - temp$period[j+1]) <= X) {
        # And there's not a cluster switch
        if (temp$camera[j] == temp$camera[j+1]) {
          # Then determine whether the bat switched partners from last time
          if (temp$receiver[j] != temp$receiver[j+1]) {
            partner_sw_den2[i] <- partner_sw_den2[i] + abs(temp$period[j] - temp$period[j+1])
            partner_sw_num2[i] <- partner_sw_num2[i] + 1
          } else {
            partner_sw_den2[i] <- partner_sw_den2[i] + abs(temp$period[j] - temp$period[j+1])
          }
        }
      }
    }
  }
}

# Get the partner switching rate, mean, and plot
#partner_sw <- partner_sw_num/partner_sw_den
partner_sw <- partner_sw_num/hours
#partner_sw_avg <- mean(partner_sw, na.rm = T)
partner_sw_avg <- mean(partner_sw)
#partner_sw2 <- partner_sw_num2/partner_sw_den2
partner_sw2 <- partner_sw_num2/hours
#partner_sw_avg2 <- mean(partner_sw2, na.rm = T)
partner_sw_avg2 <- mean(partner_sw2)
partner_sw <- data.frame(partner_sw)
partner_sw2 <- data.frame(partner_sw2)
ggplot(data = partner_sw) +
  geom_histogram(mapping = aes(x = partner_sw), colour = "black", fill = "light blue", bins = 20) +
  xlab("Percentage") +
  geom_vline(xintercept = partner_sw_avg, color = 'red', linetype = 'dotted') +
  theme_bw()
ggplot(data = partner_sw2) +
  geom_histogram(mapping = aes(x = partner_sw2), colour = "black", fill = "light blue", bins = 20) +
  xlab("Percentage") +
  geom_vline(xintercept = partner_sw_avg2, color = 'red', linetype = 'dotted') +
  theme_bw()

# Then export the data
roost_sw <- cbind(roost_sw,Bats_w)
partner_sw <- cbind(partner_sw,bats_i)
write.csv(roost_sw, "filepath\\roost_sw.csv")
write.csv(cluster_sw, "filepath\\cluster_sw.csv")
write.csv(partner_sw, "filepath\partner_sw.csv")
write.csv(partner_sw2, "filepath\\partner_sw2.csv")
