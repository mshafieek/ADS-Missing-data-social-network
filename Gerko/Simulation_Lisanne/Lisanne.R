library(mice)     
library(purrr)    
library(furrr)    
library(mvtnorm)  
library(magrittr) 
library(dplyr)    
library(tibble)   
library(relevent)
library(Epi)
library(ggplot2)
library(sna)
library(remify)
library(RSiena)
library(rem)
library(survival)
library(remstats)
library(remify)
library(visNetwork)

load("UUsummerschool.Rdata")
str(PartOfApollo_13)
summary(PartOfApollo_13)
head(PartOfApollo_13)

Apollo_tibble <- as_tibble(PartOfApollo_13) %>% 
  rename(
    actor1 = sender,
    actor2 = receiver
  )


############## missing data 
# definieer predictor relaties
# determine which column to condition on
whichcol <- c("", "actor2", "actor1")
names(whichcol) <- colnames(Apollo_tibble)

# create predictor matrix for imoutations
pred <- make.predictorMatrix(Apollo_tibble)
pred[, "time"] <- 0
pred

# use the custom pmm method
method <- make.method(Apollo_tibble)
method[c(2,3)] <- "pmm.conditional"

#### set with sufficient actors & dyads
set.seed(123) # fix seed to realize a sufficient set
indic <- sample(1:nrow(Apollo_tibble), 1500)
remify(Apollo_tibble[indic, ], model = "tie") %>% dim() 

#### Combine the sufficient set and the incomplete set
make_missing <- function(x, indic){
  sufficient <- x[indic, ]
  miss <- x[-c(indic), ] |>
    ampute(prop = .4, 
           mech = "MCAR",
           patterns = c(1,1,0)) %>% 
    .$amp
  combined <- rbind(sufficient, 
                    miss)
  return(combined[order(combined$time), ]) # sort it all like apollo
}
# runs simulaties
mcar_simulations <- furrr::future_map(1:2, ~ {
  Apollo_tibble %>% 
    make_missing(indic = indic) %>% 
    mice(m = 5, 
         maxit = 5,
         method = method,
         pred=pred,
         whichcolumn = whichcol,
         print = F)
}, .options = furrr_options(seed = 123))

########## REM
# Definieer effects
effects <- ~ -1 + reciprocity(scaling = ("std")) +
  indegreeSender() + outdegreeReceiver()

# Definieer een functie om de analyse uit te voeren op een enkele imputed dataset
analysis_function <- function(data) {
  # Rename columns
  data %>% 
    remify::remify(model = "tie") %>% 
    remstats(tie_effects = effects) # GV hier roep je de effects weer aan
}


############## create cox data 

# Define a function that prepares data for coxph
prepare_coxph_data <- function(statsObject_imp, apollo) {
  risk_sets <- attr(statsObject_imp, "riskset")
  risk_sets <- risk_sets %>% select(-'id')
  
  # Get the times
  time <- apollo$time
  
  # merge riskset with each timepoint
  combined <- merge(risk_sets, time, by = NULL)
  
  combined <- combined %>% rename("time" = "y")
  combined <- lapply(combined, as.numeric)
  combined <- as.data.frame(combined)
  
  # Create matrices for subtraction to make a status column for coxph
  combined_matrix <- data.matrix(combined)
  matrix_rows <- nrow(combined)
  
  repeated_df <- apollo[rep(seq_len(nrow(apollo)), each = 240), ] 
  repeated_df <- repeated_df[, c(2,3,1)]
  apollo_matrix <- data.matrix(repeated_df)
  
  status_matrix <- apollo_matrix - combined_matrix
  
  # create a status column
  status <- as.integer(rowSums(status_matrix == 0) == ncol(status_matrix))
  status <- as.data.frame(status)
  
  # Add status to the combined set
  combined <- cbind(combined, status)
  
  # Extract statistics and add them to the dataframe
  reciprocity <- statsObject_imp[,,1]
  indegreeSender <- statsObject_imp[,,2]
  outdegreeReceiver <- statsObject_imp[,,3]
  
  combined$reciprocity <- c(reciprocity)
  combined$indegreeSender <- c(indegreeSender)
  combined$outdegreeReceiver <- c(outdegreeReceiver)
  
  return(combined)
}

# Pas de analyse toe op elke imputed dataset
evaluations_on_simulation <- 
  mcar_simulations %>% 
  map(~.x %>% # for every simulation
        complete("all") %>% # for every imputation
        map(~.x %>% 
              analysis_function() %>% # do analysis function
              prepare_coxph_data(apollo = PartOfApollo_13) %$% # prepare cox ph
              coxph(Surv(time, status) ~ 
                      reciprocity + 
                      indegreeSender + 
                      outdegreeReceiver))) # run cox ph

# Next step: pooling

# definieer de true analyse
true.reh <- remify(edgelist = Apollo_tibble, 
                   model = "tie")
# calculate stats
stats <- remstats(tie_effects = effects, 
                  reh = true.reh)
# use the function to create the correct format of the dataframe
true.cox.set <- prepare_coxph_data(stats, PartOfApollo_13)
# fit cox model 
true.cox.fit <- coxph(Surv(time, status) ~ reciprocity + indegreeSender + 
                        outdegreeReceiver, 
                      data=true.cox.set)
true <- coefficients(true.cox.fit)

# pool de resultaten
 pool <- 
   evaluations_on_simulation %>% 
   map(~.x %>% pool(custom.t = ".data$b + .data$b / .data$m") %>% 
   .$pooled %>% 
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
   column_to_rownames("term"))
 
 # average the pool
 
      