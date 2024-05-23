

require(mice)
require(lattice)
library(tidyverse)
library(relevent)
library(broom.mixed)
library(rem)
library(ggplot2)
library(survival)
library(remify)
library(remstats)
library(devtools)
library(remstimate)


# set working directory to be the current directory
setwd("~/Applied Data Science (Master)/Thesis 2024/Code/REM")



## 1. Read the data and rename
load("UUsummerschool.Rdata")
Apollo <- PartOfApollo_13  %>% 
  rename(
    actor1 = sender,
    actor2 = receiver
  ) 

str(Apollo)   # 3882 * 3 dataframe
n <- length(unique(Apollo$time)) 
n    # N = 3882 samples

# Remove the rest of the data
rm(Class, PartOfApollo_13, Twitter_data_rem3, WTCPoliceCalls, ClassIntercept, 
   ClassIsFemale, ClassIsTeacher, WTCPoliceIsICR, con, pkg)


### Amputing Apollo dataset

## Test 1 MAR ampute

MAR.1 <-  ampute(Apollo, prop = .3, patterns = c(1,1,0),
                 mech = "MAR")

amputed_1 <- MAR.1$amp
for(col in 1:3){
  print(sum(is.na(amputed_1[col])) / nrow(amputed_1[col]))
}
# 0.3 prop achieved for column 3, 0 for others


## Test multiple patterns
patterns <- matrix(1,nrow = 2,ncol = 3)
patterns[1,3] <- patterns[2,2] <- 0
MAR.2 <-  ampute(Apollo, prop = .3, patterns = patterns, mech = "MAR")

amputed_2 <- MAR.2$amp

for(col in 1:3){
  print(sum(is.na(amputed_2[col])) / nrow(amputed_2[col]))
}
# Both missingness patterns have around 50% frequency, now lets check rows with missingness is actually 30%:

rows_with_na <- amputed_2[apply(
  amputed_2, 
  1, 
  function(x) any(is.na(x))
), ]
nrow(rows_with_na) / nrow(amputed_2)
# 0.3 missingness of rows