library(mice, warn.conflicts = FALSE)
require(lattice)
library(tidyverse)
library(magrittr)
library(dplyr) 
library(purrr)
library(furrr)  
library(relevent)
library(broom.mixed)
library(rem)
library(ggplot2)
library(survival)
library(remify)
library(remstats)
library(devtools)
library(remstimate)
library(tibble)


# set working directory to be the current directory
cur_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(cur_dir)



## 1. Read the data
load("UUsummerschool.Rdata")
rm(Twitter_data_rem3, WTCPoliceCalls, ClassIntercept, ClassIsFemale,
   ClassIsTeacher, WTCPoliceIsICR, Class)



## 2. Ampute and then impute the data

# renaming dataset for later convenience 
apollo.renamed <- PartOfApollo_13
apollo.renamed <- apollo.renamed %>% 
  rename(
    actor1 = sender,
    actor2 = receiver
  )

# making the dataset a tibble 
apollo.renamed <- as_tibble(apollo.renamed)

# determine which column to condition on
whichcol <- c("", "actor2", "actor1")
names(whichcol) <- colnames(apollo.renamed)


# create predictor matrix for imoutations
pred <- make.predictorMatrix(apollo.renamed)
pred[ ,"time"] <- 0
pred

# use the custom pmm method
method <- make.method(apollo.renamed)
method[c(2,3)] <- "pmm.conditional"

MCAR_finite <- furrr::future_map(1:50, ~ { # map over 1000 sims
  apollo.renamed %>% 
    ampute(prop = .2, 
           patterns = c(1,0,1),
           mech = "MCAR") %>% .$amp %>% 
    mice(m = 5, 
         maxit = 5,
         method = method,
         pred=pred,
         whichcolumn = whichcol,
         print = F)
}, .options = furrr_options(seed = 123))


# check if it is correct
MCAR_finite |> map(~.x %>%
  complete("long") |>
  summarize(all(actor1 != actor2)))



# creating functions 
# Define effects
effects <- ~ -1 + reciprocity(scaling = ("std")) +
  indegreeSender() + outdegreeReceiver()

# function to remify each imputed dataset

rem_function <- function(data) {
  
  # Perform the analysis as before
  reh <- remify::remify(edgelist = data, model = "tie")
  statsObject_imp <- remstats(reh = reh, tie_effects = effects)
  
  # Return the fit
  return(statsObject_imp)
}

# function to create a dataset for the cox model 
cox_sets_function <- function(statsObject_imp, apollo_data) {
  
  # take the single riskset 
  risk_sets <- attr(statsObject_imp, "riskset")
  
  # remove the id column
  risk_sets <- risk_sets %>% select(-'id')
  
  # creating one set with all risksets for each time point 
  combined <- merge(risk_sets, apollo_data$time, by=NULL)
  
  # rename column time
  combined <- combined %>% 
    rename(
      time = y)
  
  
  # ordering combined dataset
  col_order <- c("time", "sender", "receiver")
  combined <- combined[, col_order]
  
  # recoding classes of the variables to match in both datasets
  combined$sender <- as.numeric(combined$sender)
  combined$receiver <- as.numeric(combined$receiver)
  
  
  combined$status <- do.call(paste0, combined) %in%
    do.call(paste0, apollo_data)
  
  ##### manually fixing entries that should have status FALSE
  # 1
  combined$status[combined$time == 44258.8 &
                    combined$sender == "1" &
                    combined$receiver == "17" ] <- FALSE
  
  
  # 2
  combined$status[combined$time == 44263.8 &
                    combined$sender == "1" &
                    combined$receiver == "17" ] <- FALSE
  
  
  
  # 3
  combined$status[combined$time == 44272.8 &
                    combined$sender == "1" &
                    combined$receiver == "17" ] <- FALSE
  
  
  # 4
  combined$status[combined$time == 49753.8 &
                    combined$sender == "1" &
                    combined$receiver == "17" ] <- FALSE
  
  # 5
  combined$status[combined$time == 49758.8 &
                    combined$sender == "1" &
                    combined$receiver == "17" ] <- FALSE
  
  #combining the dataset with riskset to the statistic
  
  
  reciprocity <- statsObject_imp[,,1]
  indegreeSender <- statsObject_imp[,,2]
  outdegreeReceiver <- statsObject_imp[,,3]
  
  
  recip.vector <- c(reciprocity)
  combined$reciprocity <- recip.vector
  
  indegSen.vector <- c(indegreeSender)
  combined$indegreeSender <- indegSen.vector
  
  outRec.vector <- c(outdegreeReceiver)
  combined$outdegreeReceiver <- outRec.vector
  
  combined$status <- as.integer(as.logical(combined$status))
  
  return(combined)
}


# prepare the data for the cox model
MCAR_finite_imp <- MCAR_finite %>% 
        map(~.x %>% # for every simulated multiple imputation....
              complete("all") %>% # create a list of completed data sets
              map(~.x %$% # for every completed data set....
                    rem_function(.x) %>% # remify the imputed set
              cox_sets_function(apollo_data=PartOfApollo_13))) # create the sets for cox model
                
# checking if everything is correct
status.true <- MCAR_finite_imp[[1]]$'1' %>% filter(status==1)
dim(status.true)

# fit cox model 
cox.models.sims <- MCAR_finite_imp %>% 
  map(~.x %>% # for every simulation
  map(~.x %$% # for every dataset
        coxph(Surv(time, status) ~ reciprocity + indegreeSender + 
                outdegreeReceiver, 
              data=.x))) # fit cox model

# pool cox models 
pooled.results <- cox.models.sims %>%
  map(~.x %>%
        pool())



######################### cox model on complete data ######################### 
# Prepare event history
true.reh <- remify(edgelist = apollo.renamed, 
                   model = "tie")
# calculate stats
stats <- remstats(tie_effects = effects, 
                                  reh = true.reh)
# use the function to create the correct format of the dataframe
true.cox.set <- cox_sets_function(stats, PartOfApollo_13)
# fit cox model 
true.cox.fit <- coxph(Surv(time, status) ~ reciprocity + indegreeSender + 
                        outdegreeReceiver, 
                      data=true.cox.set)
true <- coefficients(true.cox.fit)




# average the pools 

pooled.average <- cox.models.sims %>%
  map(~.x %>%
        pool(custom.t = ".data$b + .data$b / .data$m") %>% # pool coefficients
  .$pooled %>% # extract table of pooled coefficients
  mutate(true = true, # add true
         df = m-1,  # correct df
         riv = Inf, # correct riv
         std.error = sqrt(t), # standard error
         statistic = estimate / std.error, # test statistic
         p.value = 2 * (pt(abs(statistic), 
                           pmax(df, 0.001), 
                           lower.tail = FALSE)), # correct p.value
         `2.5 %` = estimate - qt(.975, df) * std.error, # lower bound CI
         `97.5 %` = estimate + qt(.975, df) * std.error, # upper bound CI
         cov = `2.5 %` < true & true < `97.5 %`, # coverage
         bias = estimate - true) %>% # bias
  select(term, m, true, estimate, std.error, statistic, p.value, 
         riv, `2.5 %`, `97.5 %`, cov, bias) %>% 
  column_to_rownames("term") # `term` as rownames
) %>% 
  Reduce("+", .) / length(MCAR_finite)



