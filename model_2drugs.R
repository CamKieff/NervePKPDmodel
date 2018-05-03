#This script iteratively solves a system of ODEs and finds parameters
#the model is defined in "defineModel.R" It can work with a single ACH-like neurotransmitter (NT)
#or it can use two neurotransmitters, one constrictor (acetylcholine; ACH) and one relaxant (unknown; unk)

library(RxODE)
library(ggplot2)
source("normalizedDF.R")
source("defineModel.R")
source("runModelFunctions.R")

#Define initial parameters.
#for 2 NT: ACh consensus values here followed by initial unknown parameters
init_params <- c(KAach = 3.0228, #ach model parameters
                    KEach = 1.4223,
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

#Define a single model to run and plot
thismodel <- defineModel(ACH_mod="complex", unk_mod="none", effect_mod = "oneNT") #what model
thismodel[[1]]$model                     #check model diagnostic
#freq0 <- 10                             #what frequency
bestfit <- c("m2max", "chemax", "IC50m2","IC50che")  #what unknowns are being solved for?

Iteration(con_list=5, freq_list = c(0.3,0.7,1,3,7,10,15,30), n=100, ITmodel = thismodel, bestfit = c("m2max", "chemax", "IC50m2","IC50che"),
         init_params=init_params, hyperparams = c(4, 0.5), filename = "FormattedLowerTrachea/control/complexResultsHz_")


# freq_list = c(0.3,0.7,1,3,7,10,15,30)
# WconDF <- loadNormalizedDF(2, lower = TRUE, dataDF = "con", normDF = "cap")
# initialresults <- run_mod1(stim_freq = freq0, init_params, chosenmodel = thismodel)
# finalparams <- final_drug_params(stim_freq = freq0, m = 500, WconDF, bestfit = bestfit, init_params, initialresults, model = thismodel, chosen = FALSE)
# finalresults <- run_mod1(stim_freq = finalparams[nrow(finalparams),][["Frequency"]], finalparams[nrow(finalparams),], chosenmodel = thismodel)
# 
# p30<- (ggplot() #plot
#       + geom_path(aes(x=WconDF$Time, y=WconDF[,paste0("X", freq0, "HZ")]), color="black", alpha = 0.5)
#       + geom_path(aes(x=initialresults[,"time"], y=initialresults[,"eff2"]), color="violet",size = 1)
#       + geom_path(aes(x=finalresults[,"time"], y=finalresults[,"eff2"]), color="green",size = 1)
#       +labs(list(title=paste0(freq0, " HZ Response"), x="Time (s)", y='Normalized Response'))
#       +theme_bw()
#       +xlim(0,60)
# )
# p30
