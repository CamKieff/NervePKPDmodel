library(RxODE)
library(ggplot2)
library(plyr)
library(dplyr)
source("normalizedDF.R")
source("defineModel.R")
source("runModelFunctions.R")

init_params <- c(KAach = 1, #ach model parameters
                 KEach = 1,
                 DVach = 6,
                 EC50ach = 5.383, #calculated from ACH constriction curves of isolated tracheas
                 
                 m2max = 0.5, #complex model parameters
                 chemax = 0.5,
                 IC50m2 = 7,
                 IC50che = 7,
                 
                 KAunk = 1, #unknown relaxant neurotransmitter parameters
                 KEunk = 1,
                 DVunk = 7,
                 EC50unk = 5,
                 MAXunk = 0)

thismodel <- defineModel(ACH_mod="simple", unk_mod="none", effect_mod = "oneNT") #what model
thismodel[[1]]$model                     #check model diagnostic

#Capsaicin 0.1 Hz
Iteration(con_list=c(1,2,5,7), freq_list = 0.1, n=100, ITmodel = thismodel, bestfit = c("KAach", "KEach", "DVach", "Frequency"),
          init_params=init_params, dataDF = "cap", normDF = "cap", lower = TRUE, hyperparams = c(2, 0.1),
          filename = "FormattedLowerTrachea/capsaicin/Results0.1Hz_")
aggregate_stats(con_list = c(1,2,5,7), input_fileform = "FormattedLowerTrachea/capsaicin/Results0.1Hz_",
                output_filename = "FormattedLowerTrachea/capsaicin/cap_aggregate_0.1Hz.csv")

#Capsaicin 0.3, 1, 3, 10 Hz
Iteration(con_list=c(1,2,5,7), freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("KAach", "KEach", "DVach"),
          init_params=init_params, dataDF = "cap", normDF = "cap", lower = TRUE, hyperparams = c(2, 0.1),
          filename = "FormattedLowerTrachea/capsaicin/ResultsHz_")
aggregate_stats(con_list = c(1,2,5,7), input_fileform = "FormattedLowerTrachea/capsaicin/ResultsHz_",
                output_filename = "FormattedLowerTrachea/capsaicin/cap_aggregate_Hz.csv")

#Control 0.1 Hz
Iteration(con_list=c(1,2,5,7), freq_list = 0.1, n=100, ITmodel = thismodel, bestfit = c("KAach", "KEach", "DVach", "Frequency"),
          init_params=init_params, dataDF = "con", normDF = "cap", lower = TRUE, hyperparams = c(2, 0.1),
          filename = "FormattedLowerTrachea/control/Results0.1Hz_")
aggregate_stats(con_list = c(1,2,5,7), input_fileform = "FormattedLowerTrachea/control/Results0.1Hz_",
                output_filename = "FormattedLowerTrachea/control/con_aggregate_0.1Hz.csv")

#Control 0.3, 1, 3, 10 Hz
Iteration(con_list=c(1,2,5,7), freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("KAach", "KEach", "DVach"),
          init_params=init_params, dataDF = "con", normDF = "cap", lower = TRUE, hyperparams = c(2, 0.1),
          filename = "FormattedLowerTrachea/control/ResultsHz_")
aggregate_stats(con_list = c(1,2,5,7), input_fileform = "FormattedLowerTrachea/control/ResultsHz_",
                output_filename = "FormattedLowerTrachea/control/con_aggregate_Hz.csv")


thismodel <- defineModel(ACH_mod="complex", unk_mod="none", effect_mod = "oneNT") #what model
thismodel[[1]]$model                     #check model diagnostic

#load capsaicin best-fit 0.1 Hz parameters from above
init_df <- read.csv("FormattedLowerTrachea/capsaicin/cap_aggregate_0.1Hz.csv", header=TRUE) %>% 
  filter(var == "median") %>% 
  select(names(init_params))

#Capsaicin - complex
for(y in c(1,2,5,7)){
  Iteration(con_list=y, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
            init_params=unlist(init_df[y,]), dataDF="cap", normDF = "cap", lower = TRUE, 
            hyperparams = c(2, 0.1), filename = "FormattedLowerTrachea/capsaicin/complexResults_")
}

aggregate_stats(con_list = c(1,2,5,7), input_fileform = "FormattedLowerTrachea/capsaicin/complexResults_",
                output_filename = "FormattedLowerTrachea/capsaicin/cap_aggregate_complex.csv")

#load control best-fit 0.1 Hz parameters from above
init_df <- read.csv("FormattedLowerTrachea/control/con_aggregate_0.1Hz.csv", header=TRUE) %>% 
  filter(var == "median") %>% 
  select(names(init_params))

#Control - complex
for(y in c(1,2,5,7)){
  Iteration(con_list=y, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
            init_params=unlist(init_df[y,]), dataDF="con", normDF = "cap", lower = TRUE,
            hyperparams = c(2, 0.1), filename = "FormattedLowerTrachea/control/complexResults_")
}

aggregate_stats(con_list = c(1,2,5,7), input_fileform = "FormattedLowerTrachea/control/complexResults_",
                output_filename = "FormattedLowerTrachea/control/con_aggregate_complex.csv")

#statistics for all frequnecy tissue specific numbers. Creates a single final parameter file.
con_simple <- 
  freq_results <- rbind(read.csv("FormattedLowerTrachea/control/con_aggregate_0.1Hz.csv", header = TRUE), 
                        read.csv("FormattedLowerTrachea/control/con_aggregate_Hz.csv")) %>%
  filter(var =="median") %>%
  aggregate(by=list(.$Tissue), FUN=mean) %>%
  mutate(Capsaicin = 0, Complex = 0) %>%
  select(-Group.1, -X, -Freq, -var)

cap_simple <- 
  rbind(read.csv("FormattedLowerTrachea/capsaicin/cap_aggregate_0.1Hz.csv", header = TRUE), 
        read.csv("FormattedLowerTrachea/capsaicin/cap_aggregate_Hz.csv")) %>%
  filter(var =="median") %>%
  aggregate(by=list(.$Tissue), mean) %>%
  mutate(Capsaicin = 1, Complex = 0) %>%
  select(-Group.1, -X, -Freq, -var)

cap_complex <- 
  read.csv("FormattedLowerTrachea/capsaicin/cap_aggregate_complex.csv", header = TRUE) %>%
  filter(var =="median") %>%
  aggregate(by=list(.$Tissue), mean) %>%
  mutate(Capsaicin = 1, Complex = 1) %>%
  select(-Group.1, -X, -Freq, -var)

con_complex <- 
  read.csv("FormattedLowerTrachea/control/con_aggregate_complex.csv", header = TRUE) %>%
  filter(var =="median") %>%
  aggregate(by=list(.$Tissue), mean) %>%
  mutate(Capsaicin = 0, Complex = 1) %>%
  select(-Group.1, -X, -Freq, -var)

write.csv(rbind(cap_simple, con_simple, cap_complex, con_complex), "FormattedLowerTrachea/finalParameters.csv")

#statistics for frequency specific results for each treatment. Creates a separate file for each treatment
#find frequency specific results - control - simple model
freq_results <- rbind(read.csv("FormattedLowerTrachea/control/con_aggregate_0.1Hz.csv", header = TRUE), 
                      read.csv("FormattedLowerTrachea/control/con_aggregate_Hz.csv")) %>%
  filter(var == "median") %>%
  select(Freq, KAach, KEach, DVach) 

meanRes <- rbind(mutate(aggregate(freq_results, by=list(freq_results$Freq), mean), var ="mean"),
                 mutate(aggregate(freq_results, by=list(freq_results$Freq), function(x) sd(x)/sqrt(4)), var ="SEM"))
write.csv(meanRes, "FormattedLowerTrachea/control/con_frequency_Hz.csv")

#find frequency specific results - capsaicin - simple model
freq_results <- rbind(read.csv("FormattedLowerTrachea/capsaicin/cap_aggregate_0.1Hz.csv", header = TRUE), 
                      read.csv("FormattedLowerTrachea/capsaicin/cap_aggregate_Hz.csv")) %>%
  filter(var == "median") %>%
  select(Freq, KAach, KEach, DVach) 

meanRes <- rbind(mutate(aggregate(freq_results, by=list(freq_results$Freq), mean), var ="mean"),
                 mutate(aggregate(freq_results, by=list(freq_results$Freq), function(x) sd(x)/sqrt(4)), var ="SEM"))
write.csv(meanRes, "FormattedLowerTrachea/capsaicin/cap_frequency_Hz.csv")

#find frequency specific results - capsaicin - complex model
freq_results <- read.csv("FormattedLowerTrachea/capsaicin/cap_aggregate_complex.csv") %>%
  filter(var == "median") %>%
  select(Freq, m2max, chemax, IC50m2, IC50che) 

meanRes <- rbind(mutate(aggregate(freq_results, by=list(freq_results$Freq), mean), var ="mean"),
                 mutate(aggregate(freq_results, by=list(freq_results$Freq), function(x) sd(x)/sqrt(4)), var ="SEM"))
write.csv(meanRes, "FormattedLowerTrachea/capsaicin/cap_frequency_complex.csv")

#find frequency specific results - control - complex model
freq_results <- read.csv("FormattedLowerTrachea/control/con_aggregate_complex.csv") %>%
  filter(var == "median") %>%
  select(Freq, m2max, chemax, IC50m2, IC50che) 

meanRes <- rbind(mutate(aggregate(freq_results, by=list(freq_results$Freq), mean), var ="mean"),
                 mutate(aggregate(freq_results, by=list(freq_results$Freq), function(x) sd(x)/sqrt(4)), var ="SEM"))
write.csv(meanRes, "FormattedLowerTrachea/control/con_frequency_complex.csv")