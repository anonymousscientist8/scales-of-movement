# Load packages and clear workspace
library(tidyverse)
rm(list = ls())

# Import dataset
W <- read.csv("filepath\\roosting_information.csv")

# Filter to leave out individuals observed less than Xw times per day
# or Yw times overall
Xw <- 0.000000000001
Yw <- 25

# All bats in the model
Bats_w <- unique(W$ID)
Bats_w <- Bats_w[Bats_w != ""] # Filters out blank entries

# Convert date into something a numeric (used to count days between observations)
W$Date <- as.numeric(as.POSIXlt(W$Date, format="%m/%d/%Y")) / 
  (60*60*24)

# Add a column signify roost-switching events
switch_w <- rep(0,length(W$ID))
W <- cbind(W, switch_w)

# Add a column for number of potential switches and cryptic switch observed
potential <- rep(1,length(W$ID))
W <- cbind(W, potential)

# Now, we mark every switching event for each bat resulting from seeing the same bat in two different roosts
for (n in 1:length(W$Date)) { # For every row
  if (n != 1) { # If we aren't looking at the first entry
    # And if we aren't looking at the same bat (already ordered by bat ID)
    if (W$ID[n] == W$ID[n-1]) {
      # Then check to see if the roosts are the same between time steps
      if (W$Roost[n] != W$Roost[n-1]) {
        # And, if there is a roost change, mark this event
        W$switch_w[n] <- 1
      }
    }
  }
}

# Arrange by date and roost secondarily
W <- arrange(W,Date,Roost)

# Now check for cryptic switches (switches occurring without being seen because of missing observations)
for (n in 1:length(W$Date)) { # For every row
  # Find the current roost the bat is in
  temp_roost_1 <- W$Roost[n]

  m <- 0 # changes to 1 if we see an example of this bat again
  p <- 1 # Counter for checking n through length of columns
  
  # While we haven't checked every row or we haven't encountered the bat again
  while ((m == 0) & (p < (length(W$Date) - n))) {
    if (n != length(W$Date)) {
      # If the bat isn't seen in a particular row
      if (W$ID[n+p] != W$ID[n]) {
        # Add to p (check the next one)
        p <- p + 1
        temp_roost_2 <- 0 # Signify that a second roost has not been found
        # If it is seen
      } else {
        m <- 1 # Break the loop
        temp_roost_2 <- W$Roost[n+p] # And store the second roost ID
      }
    }
    # Check for cryptic switches
    # If the bat is observed later
    if (temp_roost_2 != 0) {
      # Find all dates between observations
      Dates_w <- unique(W$Date[n:(n+p)])
      # And for every date in between, if there is any days in between
      if (length(Dates_w) > 2) {
        r <- 0 # Counts cryptic switches
        s <- 0 # Boolean for whether every row in df is looked at
        t <- 0 # Indexing variable for checking every row
        x <- 0 # Stores date seeing first roost again
        y <- 0 # Stores date seeing second roost again
        for (q in 2:length(Dates_w)-1) { # For every date in between
          # If both roosts were the same, check if there's any day where that
          # roost was observed in between
          if (temp_roost_1 == temp_roost_2)  {
            # Filter out irrelevant rows that aren't between observations
            temp <- W[n:(n+p),]
            # Then filter out anything not on the specified date
            temp <- temp[!(temp$Date != Dates_w[q]),]
            # And if that roost appears in between
            if (temp_roost_1 %in% temp$Roost) {
              # Then there was at least two cryptic switches
              r <- 2
            }
            # And if the observations take place at different roosts
          } else {
            # Filter out irrelevant rows that aren't between observations
            temp <- W[n:(n+p),]
            # Then filter out anything not on the specified date
            temp <- temp[!(temp$Date != Dates_w[q]),]
            # While we haven't looked at every row
            while (s == 0) {
              t <- t + 1 # Updates what row is checked
              # If we've checked every date, break the loop
              if (t == length(temp$Date)) {
                s <- 1
              }
              # If the original roost hasn't been found until this row
              if ((temp$Roost[t] == temp_roost_1) & (x == 0)) {
                # Store the date the roost was found as x
                x <- temp$Date[t]
              }
              # If the second roost hasn't been found until this row
              if ((temp$Roost[t] == temp_roost_2) & (y == 0)) {
                # Store the date the roost was found as y
                y <- temp$Date[t]
              }
            }
            # If measurements were taken of both roosts on the same day
            # and the bat wasn't present, then that's an additional roost
            # switch. If the bat's first roost was measured again and the
            # bat was missing, then the bat should be at the second roost
            # later or there was a roost switch.
            if ((x <= y) & (x != 0)) {
              r <- 1
            }
          }
        }
        # Add cryptic roost switches
        W$switch_w[n+p] <- W$switch_w[n+p] + r
        W$potential[n+p] <- W$potential[n+p] + 1
      }
    }
  }
}

# Add a column for number of days since the bat was last observed
n.days <- rep(0,length(W$Date))
W <- cbind(W, n.days)

# For every bat
for (i in 1:length(Bats_w)) {
  # Used to determine whether we are looking at the first entry (1 == yes, 2 == no)
  k <- 1
  # For every entry in the dataset
  for (j in 1:length(W$Date)) {
    # If the entry is associated with the bat in question
    if (W$ID[j] == Bats_w[i]) {
      # And if it is the first entry
      if (k == 1) {
        # Store the first date of observation
        temp_old <- W$Date[j]
        # And show this is no longer the first entry
        k <- k + 1
        # And if it isn't the first observation of that bat
      } else {
        # Store the new date
        temp_new <- W$Date[j]
        # Find the time since last switch
        W$n.days[j] <- round(temp_new - temp_old)
        # And if a switch occurred
        if (W$switch_w[j] != 0) {
          # Reset the timer by setting the new temp_old to the last obs switch
          temp_old <- W$Date[j]
        }
      }
    }
  }
}

# Convert all switches with 2's to 1's to determine the probability that
# there was at least one roost switch between observations (needs to be binary)
W["switch_w"][W["switch_w"] == 2] <- 1


# In order to filter out data with an observation rate of less than Xw
# times per day, or total observations is less than Yw we need to find
# the rate of observation for each bat
for (i in 1:length(Bats_w)) {
  # Filter out bats that we aren't looking at here
  temp <- W[!(W$ID != Bats_w[i]),]
  # Number of observations over total time observed
  rate_temp <- length(temp$Date) / (max(temp$Date ,na.rm = T) -
                                      min(temp$Date, na.rm = T))
  # If this rate is lower than Xw, get rid of it
  if (rate_temp <= Xw) {
    W <- W[!(W$ID == Bats_w[i]),]
  } else {
    # If observations is lower than Yw, get rid of it
    if (length(unique(temp$Date)) <= Yw) {
      W <- W[!(W$ID == Bats_w[i]),]
    }
  }
}

# Filter out males
# W <- W[(W$Sex == "F"),]

# Filter out entries where there is no bat listed
W <- W[(W$ID != ""),]

# Updates the list of bats to those with more than Yw observations
Bats_w <- unique(W$ID)
Bats_w <- Bats_w[Bats_w != ""]

# Report out final modified dataset
write.csv(W, "filepath\\roost_switching.csv", row.names = F)

# Count the number of observations
dates <- unique(W$Date)
obs <- 0
for (i in 1:length(dates)) {
  temp <- W[(W$Date == dates[i]),]
  obs <- obs + length(unique(temp$Roost))
}

# Find the proportion of observation more or less than 2 weeks since last switch
# with a switch
a <- W[(W$n.days > 14),]
b <- W[(W$n.days <= 14),]
mean(a$switch_w, na.rm = T)
mean(b$switch_w, na.rm = T)
