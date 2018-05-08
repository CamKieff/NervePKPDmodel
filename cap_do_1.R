library(RxODE)
library(ggplot2)
source("normalizedDF.R")
source("defineModel.R")
source("runModelFunctions.R")

thismodel <- defineModel(ACH_mod="complex", unk_mod="none", effect_mod = "oneNT") #what model
thismodel[[1]]$model                     #check model diagnostic
bestfit <- c("m2max", "chemax", "IC50m2","IC50che")  #what unknowns are being solved for?

#File1
init_params <- c(KAach = 3.6532, #ach model parameters
                 KEach = 0.94526,
                 DVach = 5.9202,
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

Iteration(con_list=1, freq_list = c(0.3,0.7,1,3,7,10,15,30), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
          init_params=init_params, dataDF="cap", hyperparams = c(4, 0.5), filename = "FormattedLowerTrachea/capsaicin/complexResults_")
