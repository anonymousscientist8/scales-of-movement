## Uses the interactions of captive vampire bats, collected by Imran Razik,
## to determine partner switching duration and create a csv showing partner
## switching events and the amount of time since last partner switch
## Created by C. Raven A. Hartman
## Last Modified: 11/4/2022

# Load necessary packages and clear workspace
library(tidyverse)
library(stringi)
rm(list = ls())

# Load interaction data
load("filepath\\events.RData")
imran <- events2019

# Filter out non-grooming events and events not in the big cage
imran <- imran[!(imran$cage != 'big_cage'),]
imran <- imran[!(imran$behav != 'g'),]

# Remove the males
#males <- c('pup','noa','xanthe','alexei','ad','ald')
#imran <- imran[((imran$actor %in% males) == F),]
#imran <- imran[((imran$receiver %in% males) == F),]

# Obtain mean duration of grooming events (55.3 seconds)
duration <- mean(imran$duration)

# Obtain all bats and remove any that are NA
bats_i <- unique(c(imran$actor,imran$receiver))
bats_i <- bats_i[!is.na(bats_i)]

# Obtain all times
times_i <- unique(imran$period)

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

# Create columns for whether a partner switch occurred and time since last switch
switch <- rep(0,length(imran$actor))
latency <- rep(0,length(imran$actor))
imran <- cbind(imran,switch)
imran <- cbind(imran,latency)

# For every bat
for (j in 1:length(bats_i)) {
  # Determines whether this is the first entry we've seen this bat groom another bat
  # (0 = new)
  k <- 0
  # And for each row in the dataframe
  for (i in 1:length(imran$switch)) {
    # If that row coincides with the bat in question
    if (imran$actor[i] == bats_i[j]) {
      # And that bat has not been recorded grooming another yet
      if (k == 0) {
        # Mark that's it's groomed another bat
        k <- 1
        # And record the grooming recipient and time of initiation
        temp_old <- imran$period[i]
        temp_old_c <- imran$receiver[i]
      } else { # And if it's not the first time observed grooming another bat
        # Record the new time
        temp_new <- imran$period[i]
        # And determine time since last partner switch
        imran$latency[i] <- temp_new - temp_old
        # And specifically when the grooming recipient is changed
        if (temp_old_c != imran$receiver[i]) {
          # Record the new time period
          temp_new <- imran$period[i]
          # Mark a switch has occurred
          imran$switch[i] <- 1
          # Replace the old time with the new one
          temp_old <- temp_new
          # And note the new recipient
          temp_old_c <- imran$receiver[i]
        }
      }
    }
  }
}
# And turn the latency from seconds to minutes
imran$latency <- imran$latency/60

# Then write the dataframe to a csv
write.csv(imran, "filepath\\partner_switching.csv", row.names = F)
