library(RxODE)
library(ggplot2)
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
          init_params=init_params, dataDF = "cap", hyperparams = c(2, 0.1),
          filename = "FormattedLowerTrachea/capsaicin/Results0.1Hz_")
aggregate_stats(con_list = c(1,2,5,7), input_fileform = "FormattedLowerTrachea/capsaicin/Results0.1Hz_",
                output_filename = "FormattedLowerTrachea/capsaicin/cap_aggregate_0.1Hz.csv")

#Capsaicin 0.3, 1, 3, 10 Hz
Iteration(con_list=c(1,2,5,7), freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("KAach", "KEach", "DVach"),
          init_params=init_params, dataDF = "cap", hyperparams = c(2, 0.1),
          filename = "FormattedLowerTrachea/capsaicin/ResultsHz_")
aggregate_stats(con_list = c(1,2,5,7), input_fileform = "FormattedLowerTrachea/capsaicin/ResultsHz_",
                output_filename = "FormattedLowerTrachea/capsaicin/cap_aggregate_Hz.csv")

#Control 0.1 Hz
Iteration(con_list=c(1,2,5,7), freq_list = 0.1, n=100, ITmodel = thismodel, bestfit = c("KAach", "KEach", "DVach", "Frequency"),
          init_params=init_params, dataDF = "con", hyperparams = c(2, 0.1),
          filename = "FormattedLowerTrachea/control/Results0.1Hz_")
aggregate_stats(con_list = c(1,2,5,7), input_fileform = "FormattedLowerTrachea/control/Results0.1Hz_",
                output_filename = "FormattedLowerTrachea/control/con_aggregate_0.1Hz.csv")

#Control 0.3, 1, 3, 10 Hz
Iteration(con_list=c(1,2,5,7), freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("KAach", "KEach", "DVach"),
          init_params=init_params, dataDF = "con", hyperparams = c(2, 0.1), 
          filename = "FormattedLowerTrachea/control/ResultsHz_")
aggregate_stats(con_list = c(1,2,5,7), input_fileform = "FormattedLowerTrachea/control/ResultsHz_",
                output_filename = "FormattedLowerTrachea/control/con_aggregate_Hz.csv")


thismodel <- defineModel(ACH_mod="complex", unk_mod="none", effect_mod = "oneNT") #what model
thismodel[[1]]$model                     #check model diagnostic

#load capsaicin best-fit 0.1 Hz parameters from above
init_df <- read.csv("FormattedLowerTrachea/capsaicin/cap_aggregate_0.1Hz.csv", header=TRUE) %>% 
  filter(var == "median") %>% 
  select(names(init_params))

#File1 - capsaicin - complex
Iteration(con_list=1, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
          init_params=init_df[1,], dataDF="cap", hyperparams = c(2, 0.1), filename = "FormattedLowerTrachea/capsaicin/complexResults_")
#File2 - capsaicin - complex
Iteration(con_list=2, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
          init_params=init_df[2,], dataDF="cap", hyperparams = c(2, 0.1), filename = "FormattedLowerTrachea/capsaicin/complexResults_")
#File5 - capsaicin - complex
Iteration(con_list=5, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
          init_params=init_df[3,], dataDF="cap", hyperparams = c(2, 0.1), filename = "FormattedLowerTrachea/capsaicin/complexResults_")
#File7 - capsaicin - complex
Iteration(con_list=7, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
          init_params=init_df[4,], dataDF="cap", hyperparams = c(2, 0.1), filename = "FormattedLowerTrachea/capsaicin/complexResults_")

aggregate_stats(con_list = c(1,2,5,7), input_fileform = "FormattedLowerTrachea/capsaicin/complexResults_",
                output_filename = "FormattedLowerTrachea/capsaicin/cap_aggregate_complex.csv")

#load control best-fit 0.1 Hz parameters from above
init_df <- read.csv("FormattedLowerTrachea/control/con_aggregate_0.1Hz.csv", header=TRUE) %>% 
  filter(var == "median") %>% 
  select(names(init_params))

#File1 - control - complex
Iteration(con_list=1, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
          init_params=init_df[1,], dataDF="con", hyperparams = c(2, 0.1), filename = "FormattedLowerTrachea/control/complexResults_")
#File2 - control - complex
Iteration(con_list=2, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
          init_params=init_df[2,], dataDF="con", hyperparams = c(2, 0.1), filename = "FormattedLowerTrachea/control/complexResults_")
#File5 - control - complex
Iteration(con_list=5, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
          init_params=init_df[3,], dataDF="con", hyperparams = c(2, 0.1), filename = "FormattedLowerTrachea/control/complexResults_")
#File7 - control - complex
Iteration(con_list=7, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
          init_params=init_df[4,], dataDF="con", hyperparams = c(2, 0.1), filename = "FormattedLowerTrachea/control/complexResults_")

aggregate_stats(con_list = c(1,2,5,7), input_fileform = "FormattedLowerTrachea/control/complexResults_",
                output_filename = "FormattedLowerTrachea/control/con_aggregate_complex.csv")

#statistics go here.