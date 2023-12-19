# Load necessary packages and clear workspace
library(tidyverse)
library(stringi)
rm(list = ls())

# Load interaction data
load("filepath\\events.RData")
I <- events2019

# Filter out non-grooming events and events not in the big cage
I <- I[!(I$cage != 'big_cage'),]
I <- I[!(I$behav != 'g'),]

# Remove the males
#males <- c('pup','noa','xanthe','alexei','ad','ald')
#I <- I[((I$actor %in% males) == F),]
#I <- I[((I$receiver %in% males) == F),]

# Obtain mean duration of grooming events (55.3 seconds)
duration <- mean(I$duration)

# Obtain all bats and remove any that are NA
bats_i <- unique(c(I$actor,I$receiver))
bats_i <- bats_i[!is.na(bats_i)]

# Obtain all times
times_i <- unique(I$period)

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

# Create columns for whether a partner switch occurred and time since last switch
switch <- rep(0,length(I$actor))
latency <- rep(0,length(I$actor))
I <- cbind(I,switch)
I <- cbind(I,latency)

# For every bat
for (j in 1:length(bats_i)) {
  # Determines whether this is the first entry we've seen this bat groom another bat
  # (0 = new)
  k <- 0
  # And for each row in the dataframe
  for (i in 1:length(I$switch)) {
    # If that row coincides with the bat in question
    if (I$actor[i] == bats_i[j]) {
      # And that bat has not been recorded grooming another yet
      if (k == 0) {
        # Mark that's it's groomed another bat
        k <- 1
        # And record the grooming recipient and time of initiation
        temp_old <- I$period[i]
        temp_old_c <- I$receiver[i]
      } else { # And if it's not the first time observed grooming another bat
        # Record the new time
        temp_new <- I$period[i]
        # And determine time since last partner switch
        I$latency[i] <- temp_new - temp_old
        # And specifically when the grooming recipient is changed
        if (temp_old_c != I$receiver[i]) {
          # Record the new time period
          temp_new <- I$period[i]
          # Mark a switch has occurred
          I$switch[i] <- 1
          # Replace the old time with the new one
          temp_old <- temp_new
          # And note the new recipient
          temp_old_c <- I$receiver[i]
        }
      }
    }
  }
}
# And turn the latency from seconds to minutes
I$latency <- I$latency/60

# Then write the dataframe to a csv
write.csv(I, "filepath\\partner_switching.csv", row.names = F)
