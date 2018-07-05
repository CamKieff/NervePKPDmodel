#Depricated. This script is primarily used to plot a single curve for figures and diagnostics

#This script iteratively solves a system of ODEs and finds parameters
#the model is defined in "defineModel.R" It can work with a single ACH-like neurotransmitter (NT)
#or it can use two neurotransmitters, one constrictor (acetylcholine; ACH) and one relaxant (unknown; unk)

setwd("~/GitHub/NervePKPDmodel")
library(RxODE)
library(ggplot2)
source("normalizedDF.R")
source("defineModel.R")
source("runModelFunctions.R")

#Define initial parameters.
#for 2 NT: ACh consensus values here followed by initial unknown parameters
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

#Define a single model to run and plot
thismodel <- defineModel(ACH_mod="complex", unk_mod="simple", effect_mod = "twoNT_2") #what model
thismodel[[1]]$model                     #check model diagnostic
freq0 <- 0.1                           #what frequency
tissue_num <- 1
bestfit <- c("KAunk", "KEunk")
# bestfit <- c("KAach", "KEach", "DVach", "Frequency")
# bestfit <- c("m2max", "chemax", "IC50m2","IC50che")  #what unknowns are being solved for?

parameterlist <- read.csv("FormattedLowerTrachea/finalParameters.csv") %>% select(-X)
init_params <- parameterlist[which(parameterlist$Capsaicin == 1 & parameterlist$Complex == 1 & parameterlist$Tissue == tissue_num),]
init_params$EC50unk <- 7.4

#freq_list = c(0.3,0.7,1,3,7,10,15,30)
WconDF <- loadNormalizedDF(tissue_num, lower = TRUE, dataDF = "con", normDF = "cap")
WconDF$X0.1HZ <- WconDF$X0.1HZ - 0.015 #for tissue 1 0.1 Hz
initialresults <- run_mod1(stim_freq = freq0, init_params, chosenmodel = thismodel)
finalparams <- final_drug_params(stim_freq = freq0, m = 2000, WconDF, bestfit = bestfit, init_params, init_model=initialresults, hyper_params = c(2,0.1),model = thismodel, chosen = TRUE)
finalresults <- run_mod1(stim_freq = finalparams[nrow(finalparams),][["Frequency"]], finalparams[nrow(finalparams),], chosenmodel = thismodel)

p30<- (ggplot() #plot
      + geom_path(aes(x=WconDF$Time, y=WconDF[,paste0("X", freq0, "HZ")]), color="black", alpha = 0.85)
      + geom_path(aes(x=initialresults[,"time"], y=initialresults[,"eff2"]), color="violet",alpha = 0.75, size = 1)
      + geom_path(aes(x=finalresults[,"time"], y=finalresults[,"eff2"]), color="blue",alpha = 0.75, size = 1)
      +labs(list(title=paste0(freq0, " HZ Response"), x="Time (sec)", y='Normalized Response'))
      +theme_bw()
      +xlim(0,60)
)
p30
filter(as.data.frame(finalresults), eff2 == max(finalresults[,"eff2"]))