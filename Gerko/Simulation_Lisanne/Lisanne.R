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

Apollo_tibble <- as_tibble(PartOfApollo_13)

############## missing data 

mcar_simulations <- furrr::future_map(1:3, ~ {
  Apollo_tibble %>% 
    ampute(prop = .2, 
           mech = "MCAR", patterns = c(1, 1, 0)) %>%
    .$amp
}, .options = furrr_options(seed = 123))


############ imputing data

imputed_datasets <- furrr::future_map(mcar_simulations, ~ {
  mids <- mice(.x, m = 5, maxit = 5, method = "pmm", print = FALSE)
  pred <- mids$pred
  pred[ ,"time"] <- 0
  imp <- mice(.x, m = 5, maxit = 5, pred = pred, method = "pmm", print = FALSE)
  imp
}, .options = furrr_options(seed = 123))


########## REM
# Definieer effects
effects <- ~ -1 + reciprocity(scaling = ("std")) +
  indegreeSender() + outdegreeReceiver()

# Definieer een functie om de analyse uit te voeren op een enkele imputed dataset
analysis_function <- function(data) {
  # Rename columns
  data %>% 
    rename(
      actor1 = sender,  # vervang 'sender' door de juiste kolomnaam
      actor2 = receiver  # vervang 'receiver' door de juiste kolomnaam
    ) %>% 
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
statsObject_imp <- 
  imputed_datasets %>% 
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
      
      