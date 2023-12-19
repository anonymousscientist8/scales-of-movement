
# Load packages and clear workspace
library(tidyverse)
library(stringi)
rm(list = ls())

## I's Dataset, used to find cluster switching rate
load("filepath\\cluster_information.RData")
I <- associations2019

# Remove observations from infrequently used small cage
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

# Change NA's to 0's
I[is.na(I)] <- 0

# Create a new dataset showing the position of every bat at every time
# For every time period
for (i in 1:length(times_i)) {
  # And for every bat
  for (j in 1:length(bats_i)) {
    # If it is the first entry
    if ((i == 1) & (j == 1)) {
      # Create the dataframe
      I2 <- data.frame(time = times_i[i],
                           bat = bats_i[j])
    } else {
      # Otherwise, add new dates and bats to the dataframe
      I2_add <- data.frame(time = times_i[i],
                               bat = bats_i[j])
      I2 <- rbind(I2,I2_add)
    }
  }
}

# And for every entry in the new dataset
for (i in 1:length(I2$time)) {
  # Filter out anything from original dataset not coinciding to the time in question
  temp <- I[!(I$period != I2$time[i]),]
  # Filter out anything not coinciding with the bat in question
  r <- which(temp == I2$bat[i], arr.ind=TRUE)
  temp <- temp[r[1],]
  # If this is the first item
  if (i == 1) {
    # Store the camera that bat was observed at that time
    camera <- temp$camera[1]
  } else { # Otherwise
    # Check to see whether there is an instance where that bat was seen at that time
    if ((length(temp$camera) == 0) || (is.na(temp$camera[1]) == TRUE)) {
      # If it wasn't, say the bat was unobserved (camera = 0)
      camera_add <- 0
    } else { # And if it was
      # Mark what camera made this observation
      camera_add <- temp$camera[1]
    }
    # Then update the list
    camera <- rbind(camera,camera_add)
  }
}
# Add the camera making observations to the list
I2 <- cbind(I2,camera)

# Get rid of row names
rownames(I2) <- NULL

# Add ':' to the period stamp to separate hours from minutes
for (i in 1:length(I2$time)) {
  stri_sub(I2$time[i],nchar(I2$time[i])-1,nchar(I2$time[i])-2) <- ":"
}

# Convert the time into a numeric to find time differences
I2$time <- as.numeric(as.POSIXlt(I2$time, format="%Y.%m.%d_%H:%M"))

# Create columns for whether a switch was found on this observation and time
# since last switch
switch <- rep(0,length(I2$time))
latency <- rep(0,length(I2$time))
I2 <- cbind(I2,switch)
I2 <- cbind(I2,latency)

# For each bat
for (j in 1:length(bats_i)) {
  # Marker used to signify if first entry (0 = first entry)
  k <- 0
  # And for every entry in the new dataframe
  for (i in 1:length(I2$time)) {
    # If this the row coincides with the bat in question
    if (I2$bat[i] == bats_i[j]) {
      # And it's the first entry
      if (k == 0) {
        # Where the bat was actually viewed on a camera
        if (I2$camera[i] != 0) {
          # Note that the bat was seen
          k <- 1
          # And update the first time and camera that bat was observed on
          temp_old <- I2$time[i]
          temp_old_c <- I2$camera[i]
        }
      } else { # But if it's not the first time observed
        # Note the new time that the bat was observed
        temp_new <- I2$time[i]
        # And find the time since last switch
        I2$latency[i] <- temp_new - temp_old
        # And if the camera differs from last time
        if ((temp_old_c != I2$camera[i]) & (I2$camera[i] != 0)) {
          # Mark the switch
          I2$switch[i] <- 1
          # And mark the new camera and time of observation
          temp_old <- temp_new
          temp_old_c <- I2$camera[i]
        }
      }
    }
  }
}

# Convert latencies into number of half hours since last switch
I2$latency <- I2$latency/(60*30)

# Remove the males, if desired
males <- c('pup','noa','xanthe','alexei','ad','ald')
#I2 <- I2[((I2$bat %in% males) == F),]

# Then write the to a csv
write.csv(I2, "filepath\\cluster_switching.csv", row.names = F)
