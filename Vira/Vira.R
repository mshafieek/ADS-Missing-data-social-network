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
cur_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(cur_dir)



## 1. Read the data
load("UUsummerschool.Rdata")
ls()
str(PartOfApollo_13)
summary(PartOfApollo_13)
length(unique(PartOfApollo_13$time))


## 2. Ampute and then impute the data

apollo.renamed <- PartOfApollo_13
apollo.renamed <- apollo.renamed %>% 
  rename(
    actor1 = sender,
    actor2 = receiver
  )



####### ampuing only one sample ####### 

single.MCAR <- ampute(apollo.renamed, prop = .2, patterns = c(1,0,1),
                      mech = "MCAR")


####### imputing only one set 
single.MCAR.imp <- mice(single.MCAR$amp, m = 5, 
                        maxit = 5,
                        method = "pmm",
                        print = F)

# removing time as a predictor
pred <- single.MCAR.imp$pred
pred[ ,"time"] <- 0

single.MCAR.imp <- mice(single.MCAR$amp, m=5, maxit=5,
                        pred=pred,  method = "pmm", print=F)


#### REM ANALYSIS 



############################## COX MODEL ############################## 
# Define effects
effects <- ~ -1 + reciprocity(scaling = ("std")) +
  indegreeSender() + 
  outdegreeReceiver()

# Prepare event history
reh <- remify(edgelist = apollo.renamed,
                   model = "tie")

apollo.imp1 <- complete(single.MCAR.imp, 1)
apollo.imp2 <- complete(single.MCAR.imp, 2)
apollo.imp3 <- complete(single.MCAR.imp, 3)
apollo.imp4 <- complete(single.MCAR.imp, 4)
apollo.imp5 <- complete(single.MCAR.imp, 5)

reh1 <- remify(edgelist = apollo.imp1, model = "tie")
reh2 <- remify(edgelist = apollo.imp2, model = "tie")
reh3 <- remify(edgelist = apollo.imp3, model = "tie")
reh4 <- remify(edgelist = apollo.imp4, model = "tie")
reh5 <- remify(edgelist = apollo.imp5, model = "tie")

listOfreh = list(reh1, reh2, reh3, reh4, reh5)

# Compute statistics
statsObject <- remstats(tie_effects = effects, reh = reh)

statsObject1 <- remstats(tie_effects = effects, reh = listOfreh[[1]])
statsObject2 <- remstats(tie_effects = effects, reh = listOfreh[[2]])
statsObject3 <- remstats(tie_effects = effects, reh = listOfreh[[3]])
statsObject4 <- remstats(tie_effects = effects, reh = listOfreh[[4]])
statsObject5 <- remstats(tie_effects = effects, reh = listOfreh[[5]])

listOfstats = list(statsObject1, statsObject2, statsObject3, statsObject4, 
                   statsObject5)





# Prepare the risk set for the cox model 
# take the single riskset 
risk_sets.1 <- attr(statsObject1, "riskset")
risk_sets.2 <- attr(statsObject2, "riskset")

# remove the id column
risk_sets.1 <- risk_sets.1 %>% select(-'id')
risk_sets.2 <- risk_sets.2 %>% select(-'id')


# creating one set with all risksets for each time point 
combined.1 <- merge(risk_sets.1, apollo.renamed$time, by=NULL)
combined.2 <- merge(risk_sets.2, apollo.renamed$time, by=NULL)


# rename column time
combined.1 <- combined.1 %>% 
  rename(
    time = y)

combined.2 <- combined.2 %>% 
  rename(
    time = y)

# rename apollo columns back
apollo <- apollo.renamed %>%
  rename(
    sender = actor1,
    receiver = actor2
  )

# ordering combined dataset
col_order <- c("time", "sender", "receiver")
combined.1 <- combined.1[, col_order]
combined.2 <- combined.2[, col_order]

str(apollo)
str(combined)

# recoding classes of the variables to match in both datasets
apollo$sender <- as.character(apollo$sender)
apollo$receiver <- as.character(apollo$receiver)


combined.1$status <- do.call(paste0, combined.1) %in%
  do.call(paste0, apollo.renamed)
combined.2$status <- do.call(paste0, combined.2) %in%
  do.call(paste0, apollo.renamed)

head(combined.1[combined.1$status==TRUE,])
head(combined.2[combined.2$status==TRUE,])
count(combined.1[combined.1$status==TRUE,])
count(combined.2[combined.2$status==TRUE,])


status.true.1 <- combined.1 %>% filter(combined.1$status==TRUE)
status.true.2 <- combined.2 %>% filter(combined.2$status==TRUE)

##### time dupes 
n_occur <- data.frame(table(status.true.2$time))
duplicates <- n_occur[n_occur$Freq > 1,]
count(duplicates)



##### manually fixing entries that should have status FALSE
# 1
status.true[status.true$time == 44258.8,]
apollo[apollo$time == 44258.8,]


combined.1$status[combined.1$time == 44258.8 &
                  combined.1$sender == "1" &
                  combined.1$receiver == "17" ] <- FALSE

combined.2$status[combined.2$time == 44258.8 &
                  combined.2$sender == "1" &
                  combined.2$receiver == "17" ] <- FALSE

# 2
status.true[status.true$time == 44263.8,]
apollo[apollo$time == 44263.8,]


combined.1$status[combined.1$time == 44263.8 &
                  combined.1$sender == "1" &
                  combined.1$receiver == "17" ] <- FALSE

combined.2$status[combined.2$time == 44263.8 &
                  combined.2$sender == "1" &
                  combined.2$receiver == "17" ] <- FALSE

# 3
status.true[status.true$time == 44272.8,]
apollo[apollo$time == 44272.8,]


combined.1$status[combined.1$time == 44272.8 &
                  combined.1$sender == "1" &
                  combined.1$receiver == "17" ] <- FALSE
combined.2$status[combined.2$time == 44272.8 &
                  combined.2$sender == "1" &
                  combined.2$receiver == "17" ] <- FALSE

# 4
status.true[status.true$time == 49753.8,]
apollo[apollo$time == 49753.8,]


combined.1$status[combined.1$time == 49753.8 &
                  combined.1$sender == "1" &
                  combined.1$receiver == "17" ] <- FALSE
combined.2$status[combined.2$time == 49753.8 &
                  combined.2$sender == "1" &
                  combined.2$receiver == "17" ] <- FALSE

# 5

status.true[status.true$time == 49758.8,]
apollo[apollo$time == 49758.8,]


combined.1$status[combined.1$time == 49758.8 &
                  combined.1$sender == "1" &
                  combined.1$receiver == "17" ] <- FALSE

combined.2$status[combined.2$time == 49758.8 &
                  combined.2$sender == "1" &
                  combined.2$receiver == "17" ] <- FALSE


### checking again if everything is correct now 

length(unique(combined.1$time[combined.1$status==TRUE]))
length(unique(combined.2$time[combined.2$status==TRUE]))

tail(combined.1[combined.1$status==TRUE,])
tail(combined.2[combined.2$status==TRUE,])
tail(apollo)



######### combining the dataset with riskset to the statistic ######### 



reciprocity.1 <- statsObject1[,,1]
indegreeSender.1 <- statsObject1[,,2]
outdegreeReceiver.1 <- statsObject1[,,3]

reciprocity.2 <- statsObject2[,,1]
indegreeSender.2 <- statsObject2[,,2]
outdegreeReceiver.2 <- statsObject2[,,3]



recip.vector.1 <- c(reciprocity.1)
combined.1$reciprocity.1 <- recip.vector.1

indegSen.vector.1 <- c(indegreeSender.1)
combined.1$indegreeSender.1 <- indegSen.vector.1

outRec.vector.1 <- c(outdegreeReceiver.1)
combined.1$outdegreeReceiver.1 <- outRec.vector.1

combined.1$status <- as.integer(as.logical(combined.1$status))
count(combined.1[combined.2$status==1,])
tail(combined.1)

#### set 2
recip.vector.2 <- c(reciprocity.2)
combined.2$reciprocity.2 <- recip.vector.2

indegSen.vector.2 <- c(indegreeSender.2)
combined.2$indegreeSender.2 <- indegSen.vector.2

outRec.vector.2 <- c(outdegreeReceiver.2)
combined.2$outdegreeReceiver.2 <- outRec.vector.2

combined.2$status <- as.integer(as.logical(combined.2$status))
count(combined.2[combined.2$status==1,])
tail(combined.2)




cox_t.1 <- coxph(Surv(time, status) ~ reciprocity.1 + 
                 indegreeSender.1 + 
                 outdegreeReceiver.1, 
               data=combined.1)

summary(cox_t.1)



cox_t.2 <- coxph(Surv(time, status) ~ reciprocity.2 + 
                   indegreeSender.2+ 
                   outdegreeReceiver.2, 
                 data=combined.2)
class(cox_t.1)

summary(cox_t.2)


listOffits = list(cox_t.1, cox_t.2)
pool(listOffits)





############################## testing what is wrong ############################## 
filter(combined, combined$sender == apollo.renamed$actor1 & 
         combined$receiver == apollo.renamed$actor2 &
         combined$time == apollo.renamed$time)

test.apollo <- apollo.renamed[apollo.renamed$time == 11849.2,]
test.combined <- combined[combined$time == 11849.2 & combined$sender == 18 & 
                            combined$receiver == 2,]


test.apollo <- test.apollo %>% 
  rename(
    sender = actor1,
    receiver = actor2
  )
test.combined <- test.combined %>% select(-'status')

class(test.apollo)
class(test.combined)
str(test.apollo)
str(test.combined)

do.call(paste0, test.apollo) %in% do.call(paste0, test.combined)


library(arsenal)
summary(comparedf(test.apollo, test.combined))

str(apollo.renamed)
str(combined)


combined %>% mutate(
  status = case_when(
    combined$sender == apollo.renamed$actor1 &
      combined$receiver == apollo.renamed$actor2 &
      combined$time == apollo.renamed$time ~ 1))
combined[combined$status == 1,]


