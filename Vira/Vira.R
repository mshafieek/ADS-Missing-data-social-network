require(mice)
require(lattice)
library(tidyverse)
library(relevent)
library(broom.mixed)
library(rem)
library(ggplot2)


# set working directory to be the current directory
cur_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(cur_dir)



## 1. Read the data
load("UUsummerschool.Rdata")
ls()
str(PartOfApollo_13)
summary(PartOfApollo_13)


## 2. Ampute and then impute the data

apollo_tibble <- as_tibble(PartOfApollo_13)

MCAR_apollo <- furrr::future_map(1:3, ~ { # map 3 sims
    apollo_tibble %>% 
      ampute(prop = .2, 
             patterns = c(1,0,1),
             mech = "MCAR"
            ) %>% .$amp %>% 
      mice(m = 5, 
           maxit = 5,
           method = "pmm",
           print = F,
           pred = c(0,0,1))
  }, .options = furrr_options(seed = 2023))


class(MCAR_apollo) # list
class(MCAR_apollo[[1]]) # mids
class(MCAR_apollo[[1]]$imp) # list

attributes(MCAR_apollo)
####### ampuing only one sample ####### 

single.MCAR <- ampute(PartOfApollo_13, prop = .2, patterns = c(1,0,1),
                      mech = "MCAR")

attributes(single.MCAR)
head(single.MCAR$amp, 20)

single.MCAR.imp <- mice(single.MCAR$amp, m = 5, 
                        maxit = 5,
                        method = "pmm",
                        print = F)

# removing time as a predictor
pred <- single.MCAR.imp$pred
pred[ ,"time"] <- 0

single.MCAR.imp <- mice(single.MCAR$amp, m=5, maxit=5,
                        pred=pred,  method = "pmm", print=F)


## 3. Ran REM analysis 


fit <- with(single.MCAR.imp,
            rem.dyad(complete(single.MCAR.imp, i), n=19, 
                     effects=c("PSAB-BA", "NIDSnd",
                               "NODRec"), hessian= TRUE))
fit.test <-  with(single.MCAR.imp,
                  rem(complete(single.MCAR.imp, i), n=19, 
                           effects=c("PSAB-BA", "NIDSnd",
                                     "NODRec"), hessian= TRUE))
fit
summary(fit$analyses[[2]])

class(fit)
ls(fit)



## 4. Pool results 
pool.fit <- pool(fit, dfcom=3879)
summary(pool.fit)



#################### REM on the fully observed data #################### 
apollo.obs.fit <- rem.dyad(PartOfApollo_13, n=19,
                       effects=c("PSAB-BA", "NIDSnd", "NODRec"),
                       hessian=TRUE) 
summary(apollo.obs.fit)


PartOfApollo_13$inertia <- inertiaStat(data = PartOfApollo_13,
                                       time = PartOfApollo_13$time,
                             sender = PartOfApollo_13$sender, 
                             target = PartOfApollo_13$receiver,
                             halflife = 2)

?inertiaStat()
# plot inertia over time

ggplot(PartOfApollo_13, aes (time, inertia) ) +
  geom_point() + geom_smooth()

# calculate reciprocity statistic
PartOfApollo_13$recip <- reciprocityStat(data = PartOfApollo_13,
                               time = PartOfApollo_13$time,
                               sender = PartOfApollo_13$sender,
                               target = PartOfApollo_13$receiver,
                               halflife = 2)
# plot sender-outdegree over time

ggplot(PartOfApollo_13, aes(time, recip)) +
  geom_point()+ geom_smooth()
