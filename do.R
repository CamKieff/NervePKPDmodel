library(RxODE)
library(ggplot2)
library(plyr)
library(dplyr)

set.seed(1606)
setwd("~/GitHub/NervePKPDmodel")

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

EFSfileDB <- read.csv("EFSfileDB.csv", header = TRUE) %>% filter(TISSUE == "L")
freq_list <- c(0.1, 0.3, 1, 3, 10)
#EFSfiles$DATE <- as.Date(EFSfiles$DATE, format = "%Y-%m%d")

Lcontrol_EFSdata <- vector("list", length=nrow(EFSfileDB))
for (i in 1:length(Lcontrol_EFSdata)){
  df <- read.csv(as.character(EFSfileDB$CONTROL_FILE[i]), header = TRUE)
  df <- df %>% filter(Time <= 60 & (Time*10) %% 1 == 0) %>%
    select(Time, paste0('X', freq_list, "HZ"))
  df_min <- min(df[, 2:ncol(df)])
  df_max <- max(df[, 2:ncol(df)])
  X0.1HZ_max <- max(df$X0.1HZ)
  #df <- df %>% mutate_at(vars(-Time), function(x){(x-df_min)/(EFSfileDB$KCL[i])}) %>%
  df <- df %>% mutate_at(vars(-Time), function(x){(x-df_min)/(X0.1HZ_max-df_min)}) %>%
    mutate(ID = EFSfileDB$CODE[i])
  Lcontrol_EFSdata[[i]] <- df
}
Lcontrol_EFSdata <- bind_rows(Lcontrol_EFSdata)

g1 <- ggplot(Lcontrol_EFSdata[, c("Time", "X1HZ", "ID")], aes(x=Time, y = X1HZ, color = ID))
g1 <- g1 + geom_line()
g1

g2 <- facetgraph(conDF = test, init_params = init_params, consensus_params = results, bestfit = c("KAach", "KEach", "DVach"), consensus = TRUE, freq_list = freq_list)

# should re-do file using the find_allfreq_params functions
# freq_list = c(0.1, 0.3, 1, 3, 10), m = 50, WconDF, bestfit,
# init_params, model = thismodel, hyper_params = c(2,0.1)
results <- find_allfreq_params(HZdf = test, init_params = init_params, bestfit = c("KAach", "KEach", "DVach"))

# 
# #Capsaicin 0.1 Hz
# Iteration(con_list=c(1,2,5,7), freq_list = 0.1, n=100, ITmodel = thismodel, bestfit = c("KAach", "KEach", "DVach", "Frequency"),
#           init_params=init_params, dataDF = "cap", hyperparams = c(2, 0.1),
#           filename = "FormattedLowerTrachea/capsaicin/Results0.1Hz_")
# aggregate_stats(con_list = c(1,2,5,7), input_fileform = "FormattedLowerTrachea/capsaicin/Results0.1Hz_",
#                 output_filename = "FormattedLowerTrachea/capsaicin/cap_aggregate_0.1Hz.csv")
# 
# #Capsaicin 0.3, 1, 3, 10 Hz
# Iteration(con_list=c(1,2,5,7), freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("KAach", "KEach", "DVach"),
#           init_params=init_params, dataDF = "cap", hyperparams = c(2, 0.1),
#           filename = "FormattedLowerTrachea/capsaicin/ResultsHz_")
# aggregate_stats(con_list = c(1,2,5,7), input_fileform = "FormattedLowerTrachea/capsaicin/ResultsHz_",
#                 output_filename = "FormattedLowerTrachea/capsaicin/cap_aggregate_Hz.csv")
# 
# #Control 0.1 Hz
# Iteration(con_list=c(1,2,5,7), freq_list = 0.1, n=100, ITmodel = thismodel, bestfit = c("KAach", "KEach", "DVach", "Frequency"),
#           init_params=init_params, dataDF = "con", hyperparams = c(2, 0.1),
#           filename = "FormattedLowerTrachea/control/Results0.1Hz_")
# aggregate_stats(con_list = c(1,2,5,7), input_fileform = "FormattedLowerTrachea/control/Results0.1Hz_",
#                 output_filename = "FormattedLowerTrachea/control/con_aggregate_0.1Hz.csv")
# 
# #Control 0.3, 1, 3, 10 Hz
# Iteration(con_list=c(1,2,5,7), freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("KAach", "KEach", "DVach"),
#           init_params=init_params, dataDF = "con", hyperparams = c(2, 0.1),
#           filename = "FormattedLowerTrachea/control/ResultsHz_")
# aggregate_stats(con_list = c(1,2,5,7), input_fileform = "FormattedLowerTrachea/control/ResultsHz_",
#                 output_filename = "FormattedLowerTrachea/control/con_aggregate_Hz.csv")
# 
# 
# thismodel <- defineModel(ACH_mod="complex", unk_mod="none", effect_mod = "oneNT") #what model
# thismodel[[1]]$model                     #check model diagnostic
# 
# #load capsaicin best-fit 0.1 Hz parameters from above
# init_df <- read.csv("FormattedLowerTrachea/capsaicin/cap_aggregate_0.1Hz.csv", header=TRUE) %>% 
#   filter(var == "median") %>% 
#   select(names(init_params))
# 
# #File1 - capsaicin - complex
# Iteration(con_list=1, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
#           init_params=unlist(init_df[1,]), dataDF="cap", hyperparams = c(2, 0.1), filename = "FormattedLowerTrachea/capsaicin/complexResults_")
# #File2 - capsaicin - complex
# Iteration(con_list=2, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
#           init_params=unlist(init_df[2,]), dataDF="cap", hyperparams = c(2, 0.1), filename = "FormattedLowerTrachea/capsaicin/complexResults_")
# #File5 - capsaicin - complex
# Iteration(con_list=5, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
#           init_params=unlist(init_df[3,]), dataDF="cap", hyperparams = c(2, 0.1), filename = "FormattedLowerTrachea/capsaicin/complexResults_")
# #File7 - capsaicin - complex
# Iteration(con_list=7, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
#           init_params=unlist(init_df[4,]), dataDF="cap", hyperparams = c(2, 0.1), filename = "FormattedLowerTrachea/capsaicin/complexResults_")
# 
# aggregate_stats(con_list = c(1,2,5,7), input_fileform = "FormattedLowerTrachea/capsaicin/complexResults_",
#                 output_filename = "FormattedLowerTrachea/capsaicin/cap_aggregate_complex.csv")
# 
# #load control best-fit 0.1 Hz parameters from above
# init_df <- read.csv("FormattedLowerTrachea/control/con_aggregate_0.1Hz.csv", header=TRUE) %>% 
#   filter(var == "median") %>% 
#   select(names(init_params))
# 
# #File1 - control - complex
# Iteration(con_list=1, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
#           init_params=unlist(init_df[1,]), dataDF="con", hyperparams = c(2, 0.1), filename = "FormattedLowerTrachea/control/complexResults_")
# #File2 - control - complex
# Iteration(con_list=2, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
#           init_params=unlist(init_df[2,]), dataDF="con", hyperparams = c(2, 0.1), filename = "FormattedLowerTrachea/control/complexResults_")
# #File5 - control - complex
# Iteration(con_list=5, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
#           init_params=unlist(init_df[3,]), dataDF="con", hyperparams = c(2, 0.1), filename = "FormattedLowerTrachea/control/complexResults_")
# #File7 - control - complex
# Iteration(con_list=7, freq_list = c(0.3,1,3,10), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
#           init_params=unlist(init_df[4,]), dataDF="con", hyperparams = c(2, 0.1), filename = "FormattedLowerTrachea/control/complexResults_")
# 
# aggregate_stats(con_list = c(1,2,5,7), input_fileform = "FormattedLowerTrachea/control/complexResults_",
#                 output_filename = "FormattedLowerTrachea/control/con_aggregate_complex.csv")
# 
# #statistics for all frequnecy tissue specific numbers. Creates a single final parameter file.
# con_simple <- 
#   freq_results <- rbind(read.csv("FormattedLowerTrachea/control/con_aggregate_0.1Hz.csv", header = TRUE), 
#                         read.csv("FormattedLowerTrachea/control/con_aggregate_Hz.csv")) %>%
#   filter(var =="median") %>%
#   #plyr::rename(c("Group.1" = "Freq", "Group.2" = "Tissue")) %>%
#   aggregate(by=list(.$Tissue), FUN=mean) %>%
#   mutate(Capsaicin = 0, Complex = 0) %>%
#   select(-Group.1, -X, -Freq, -var)
# 
# cap_simple <- 
#   rbind(read.csv("FormattedLowerTrachea/capsaicin/cap_aggregate_0.1Hz.csv", header = TRUE), 
#         read.csv("FormattedLowerTrachea/capsaicin/cap_aggregate_Hz.csv")) %>%
#   filter(var =="median") %>%
#   #plyr::rename(c("Group.1" = "Freq", "Group.2" = "Tissue")) %>% #aggregate_stats should do this
#   aggregate(by=list(.$Tissue), mean) %>%
#   mutate(Capsaicin = 1, Complex = 0) %>%
#   select(-Group.1, -X, -Freq, -var)
# 
# cap_complex <- 
#   read.csv("FormattedLowerTrachea/capsaicin/cap_aggregate_complex.csv", header = TRUE) %>%
#   filter(var =="median") %>%
#   #plyr::rename(c("Group.1" = "Freq", "Group.2" = "Tissue")) %>%
#   aggregate(by=list(.$Tissue), mean) %>%
#   mutate(Capsaicin = 1, Complex = 1) %>%
#   select(-Group.1, -X, -Freq, -var)
# 
# con_complex <- 
#   read.csv("FormattedLowerTrachea/control/con_aggregate_complex.csv", header = TRUE) %>%
#   filter(var =="median") %>%
#   #plyr::rename(c("Group.1" = "Freq", "Group.2" = "Tissue")) %>%
#   aggregate(by=list(.$Tissue), mean) %>%
#   mutate(Capsaicin = 0, Complex = 1) %>%
#   select(-Group.1, -X, -Freq, -var)
# 
# write.csv(rbind(cap_simple, con_simple, cap_complex, con_complex), "FormattedLowerTrachea/finalParameters.csv")
# 
# #statistics for frequency specific results for each treatment. Creates a separate file for each treatment
# #find frequency specific results - control - simple model
# freq_results <- rbind(read.csv("FormattedLowerTrachea/control/con_aggregate_0.1Hz.csv", header = TRUE), 
#                       read.csv("FormattedLowerTrachea/control/con_aggregate_Hz.csv")) %>%
#   filter(var == "median") %>%
#   select(Freq, KAach, KEach, DVach) 
# 
# meanRes <- rbind(mutate(aggregate(freq_results, by=list(freq_results$Freq), mean), var ="mean"),
#                  mutate(aggregate(freq_results, by=list(freq_results$Freq), function(x) sd(x)/sqrt(4)), var ="SEM"))
# write.csv(meanRes, "FormattedLowerTrachea/control/con_frequency_Hz.csv")
# 
# #find frequency specific results - capsaicin - simple model
# freq_results <- rbind(read.csv("FormattedLowerTrachea/capsaicin/cap_aggregate_0.1Hz.csv", header = TRUE), 
#                       read.csv("FormattedLowerTrachea/capsaicin/cap_aggregate_Hz.csv")) %>%
#   filter(var == "median") %>%
#   select(Freq, KAach, KEach, DVach) 
# 
# meanRes <- rbind(mutate(aggregate(freq_results, by=list(freq_results$Freq), mean), var ="mean"),
#                  mutate(aggregate(freq_results, by=list(freq_results$Freq), function(x) sd(x)/sqrt(4)), var ="SEM"))
# write.csv(meanRes, "FormattedLowerTrachea/capsaicin/cap_frequency_Hz.csv")
# 
# #find frequency specific results - capsaicin - complex model
# freq_results <- read.csv("FormattedLowerTrachea/capsaicin/cap_aggregate_complex.csv") %>%
#   filter(var == "median") %>%
#   select(Freq, m2max, chemax, IC50m2, IC50che) 
# 
# meanRes <- rbind(mutate(aggregate(freq_results, by=list(freq_results$Freq), mean), var ="mean"),
#                  mutate(aggregate(freq_results, by=list(freq_results$Freq), function(x) sd(x)/sqrt(4)), var ="SEM"))
# write.csv(meanRes, "FormattedLowerTrachea/capsaicin/cap_frequency_complex.csv")
# 
# #find frequency specific results - control - complex model
# freq_results <- read.csv("FormattedLowerTrachea/control/con_aggregate_complex.csv") %>%
#   filter(var == "median") %>%
#   select(Freq, m2max, chemax, IC50m2, IC50che) 
# 
# meanRes <- rbind(mutate(aggregate(freq_results, by=list(freq_results$Freq), mean), var ="mean"),
#                  mutate(aggregate(freq_results, by=list(freq_results$Freq), function(x) sd(x)/sqrt(4)), var ="SEM"))
# write.csv(meanRes, "FormattedLowerTrachea/control/con_frequency_complex.csv")