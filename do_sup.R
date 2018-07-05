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
Iteration(con_list=c(1,2,3,4,5), freq_list = 0.1, n=100, ITmodel = thismodel, bestfit = c("KAach", "KEach", "DVach", "Frequency"),
          init_params=init_params, dataDF = "cap", lower = FALSE, hyperparams = c(2, 0.1),
          filename = "FormattedUpperTrachea/capsaicin/Results0.1Hz_")
aggregate_stats(con_list = c(1,2,3,4,5), input_fileform = "FormattedUpperTrachea/capsaicin/Results0.1Hz_",
                output_filename = "FormattedUpperTrachea/capsaicin/cap_aggregate_0.1Hz.csv")

#Capsaicin 0.3, 1, 3, 10 Hz
Iteration(con_list=c(1,2,3,4,5), freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("KAach", "KEach", "DVach"),
          init_params=init_params, dataDF = "cap", lower = FALSE, hyperparams = c(2, 0.1),
          filename = "FormattedUpperTrachea/capsaicin/ResultsHz_")
aggregate_stats(con_list = c(1,2,3,4,5), input_fileform = "FormattedUpperTrachea/capsaicin/ResultsHz_",
                output_filename = "FormattedUpperTrachea/capsaicin/cap_aggregate_Hz.csv")

#Control 0.1 Hz
Iteration(con_list=c(1,2,3,4,5), freq_list = 0.1, n=100, ITmodel = thismodel, bestfit = c("KAach", "KEach", "DVach", "Frequency"),
          init_params=init_params, dataDF = "con", normDF = "con", lower = FALSE, hyperparams = c(2, 0.1),
          filename = "FormattedUpperTrachea/control/Results0.1Hz_")
aggregate_stats(con_list = c(1,2,3,4,5), input_fileform = "FormattedUpperTrachea/control/Results0.1Hz_",
                output_filename = "FormattedUpperTrachea/control/con_aggregate_0.1Hz.csv")

#Control 0.3, 1, 3, 10 Hz
Iteration(con_list=c(1,2,3,4,5), freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("KAach", "KEach", "DVach"),
          init_params=init_params, dataDF = "con", normDF = "con", lower = FALSE, hyperparams = c(2, 0.1),
          filename = "FormattedUpperTrachea/control/ResultsHz_")
aggregate_stats(con_list = c(1,2,3,4,5), input_fileform = "FormattedUpperTrachea/control/ResultsHz_",
                output_filename = "FormattedUpperTrachea/control/con_aggregate_Hz.csv")


thismodel <- defineModel(ACH_mod="complex", unk_mod="none", effect_mod = "oneNT") #what model
thismodel[[1]]$model                     #check model diagnostic

#load capsaicin best-fit 0.1 Hz parameters from above
init_df <- read.csv("FormattedUpperTrachea/capsaicin/cap_aggregate_0.1Hz.csv", header=TRUE) %>% 
  filter(var == "median") %>% 
  select(names(init_params))

#Capsaicin - complex
for (y in c(1,2,3,4,5)){
  Iteration(con_list=y, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
            init_params=unlist(init_df[y,]), dataDF="cap", lower = FALSE, hyperparams = c(2, 0.1), filename = "FormattedUpperTrachea/capsaicin/complexResults_")
}
aggregate_stats(con_list = c(1,2,3,4,5), input_fileform = "FormattedUpperTrachea/capsaicin/complexResults_",
                output_filename = "FormattedUpperTrachea/capsaicin/cap_aggregate_complex.csv")

#load control best-fit 0.1 Hz parameters from above
init_df <- read.csv("FormattedUpperTrachea/control/con_aggregate_0.1Hz.csv", header=TRUE) %>% 
  filter(var == "median") %>% 
  select(names(init_params))
# Control - complex
for( y in c(1,2,3,4,5)){
  Iteration(con_list=y, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
            init_params=unlist(init_df[y,]), dataDF="con", normDF = "con", lower = FALSE, hyperparams = c(2, 0.1), filename = "FormattedUpperTrachea/control/complexResults_")
}

aggregate_stats(con_list = c(1,2,3,4,5), input_fileform = "FormattedUpperTrachea/control/complexResults_",
                output_filename = "FormattedUpperTrachea/control/con_aggregate_complex.csv")

#statistics for all frequency tissue specific numbers. Creates a single final parameter file.
con_simple <- 
  rbind(read.csv("FormattedUpperTrachea/control/con_aggregate_0.1Hz.csv", header = TRUE), 
                        read.csv("FormattedUpperTrachea/control/con_aggregate_Hz.csv")) %>%
  filter(var =="median") %>%
  plyr::rename(c("Group.1" = "Freq", "Group.2" = "Tissue")) %>%
  aggregate(by=list(.$Tissue), FUN=mean) %>%
  mutate(Capsaicin = 0, Complex = 0) %>%
  select(-Group.1, -X, -Freq, -var)

cap_simple <- 
  rbind(read.csv("FormattedUpperTrachea/capsaicin/cap_aggregate_0.1Hz.csv", header = TRUE), 
        read.csv("FormattedUpperTrachea/capsaicin/cap_aggregate_Hz.csv")) %>%
  filter(var =="median") %>%
  plyr::rename(c("Group.1" = "Freq", "Group.2" = "Tissue")) %>%
  aggregate(by=list(.$Tissue), mean) %>%
  mutate(Capsaicin = 1, Complex = 0) %>%
  select(-Group.1, -X, -Freq, -var)

cap_complex <- 
  read.csv("FormattedUpperTrachea/capsaicin/cap_aggregate_complex.csv", header = TRUE) %>%
  filter(var =="median") %>%
  plyr::rename(c("Group.1" = "Freq", "Group.2" = "Tissue")) %>%
  aggregate(by=list(.$Tissue), mean) %>%
  mutate(Capsaicin = 1, Complex = 1) %>%
  select(-Group.1, -X, -Freq, -var)

con_complex <- 
  read.csv("FormattedUpperTrachea/control/con_aggregate_complex.csv", header = TRUE) %>%
  filter(var =="median") %>%
  plyr::rename(c("Group.1" = "Freq", "Group.2" = "Tissue")) %>%
  aggregate(by=list(.$Tissue), mean) %>%
  mutate(Capsaicin = 0, Complex = 1) %>%
  select(-Group.1, -X, -Freq, -var)

write.csv(rbind(cap_simple, con_simple, cap_complex, con_complex), "FormattedUpperTrachea/finalParameters.csv")

#statistics for frequency specific results for each treatment. Creates a separate file for each treatment
#find frequency specific results - control - simple model
freq_results <- rbind(read.csv("FormattedUpperTrachea/control/con_aggregate_0.1Hz.csv", header = TRUE), 
                      read.csv("FormattedUpperTrachea/control/con_aggregate_Hz.csv")) %>%
  plyr::rename(c("Group.1" = "Freq", "Group.2" = "Tissue")) %>%
  filter(var == "median") %>%
  select(Freq, KAach, KEach, DVach) 

meanRes <- rbind(mutate(aggregate(freq_results, by=list(freq_results$Freq), mean), var ="mean"),
                 mutate(aggregate(freq_results, by=list(freq_results$Freq), function(x) sd(x)/sqrt(4)), var ="SEM"))
write.csv(meanRes, "FormattedUpperTrachea/control/con_frequency_Hz.csv")

#find frequency specific results - capsaicin - simple model
freq_results <- rbind(read.csv("FormattedUpperTrachea/capsaicin/cap_aggregate_0.1Hz.csv", header = TRUE), 
                      read.csv("FormattedUpperTrachea/capsaicin/cap_aggregate_Hz.csv")) %>%
  plyr::rename(c("Group.1" = "Freq", "Group.2" = "Tissue")) %>%
  filter(var == "median") %>%
  select(Freq, KAach, KEach, DVach) 

meanRes <- rbind(mutate(aggregate(freq_results, by=list(freq_results$Freq), mean), var ="mean"),
                 mutate(aggregate(freq_results, by=list(freq_results$Freq), function(x) sd(x)/sqrt(4)), var ="SEM"))
write.csv(meanRes, "FormattedUpperTrachea/capsaicin/cap_frequency_Hz.csv")

#find frequency specific results - capsaicin - complex model
freq_results <- read.csv("FormattedUpperTrachea/capsaicin/cap_aggregate_complex.csv") %>%
  plyr::rename(c("Group.1" = "Freq", "Group.2" = "Tissue")) %>%
  filter(var == "median") %>%
  select(Freq, m2max, chemax, IC50m2, IC50che) 

meanRes <- rbind(mutate(aggregate(freq_results, by=list(freq_results$Freq), mean), var ="mean"),
                 mutate(aggregate(freq_results, by=list(freq_results$Freq), function(x) sd(x)/sqrt(4)), var ="SEM"))
write.csv(meanRes, "FormattedUpperTrachea/capsaicin/cap_frequency_complex.csv")

#find frequency specific results - control - complex model
freq_results <- read.csv("FormattedUpperTrachea/control/con_aggregate_complex.csv") %>%
  plyr::rename(c("Group.1" = "Freq", "Group.2" = "Tissue")) %>%
  filter(var == "median") %>%
  select(Freq, m2max, chemax, IC50m2, IC50che) 

meanRes <- rbind(mutate(aggregate(freq_results, by=list(freq_results$Freq), mean), var ="mean"),
                 mutate(aggregate(freq_results, by=list(freq_results$Freq), function(x) sd(x)/sqrt(4)), var ="SEM"))
write.csv(meanRes, "FormattedUpperTrachea/control/con_frequency_complex.csv")