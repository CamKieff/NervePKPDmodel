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
init_params <- list(KAach = 0.2417, #ach model parameters
                    KEach = 0.5268,
                    DVach = 5.902,
                    EC50ach = 5.383, #calculated from ACH constriction curves of isolated tracheas

                    m2max = 0.5, #complex model parameters
                    chemax = 0.5,
                    IC50m2 = 7,
                    IC50che = 7,

                    KAunk = 1, #unknown relaxant neurotransmitter parameters
                    KEunk = 1,
                    DVunk = 7,
                    EC50unk = 5,
                    MAXunk = 1)

#Define a single model to run and plot
thismodel <- defineModel(ACH_mod="simple", unk_mod="none", effect_mod = "oneNT") #what model
thismodel[[1]]$model                     #check model diagnostic
freq0 <- 1                               #what frequency
bestfit <- c("KAach", "KEach", "DVach")  #what unknowns are being solved for?

WconDF <- loadNormalizedDF(1, lower = TRUE, dataDF = "cap", normDF = "cap")
initialresults <- run_mod1(stim_freq = freq0, init_params)
finalparams <- final_drug_params(stim_freq = freq0, m = 500, WconDF, bestfit = bestfit, init_params, initialresults)
finalresults <- run_mod1(stim_freq = freq0, finalparams[nrow(finalparams),])

p30<- (ggplot() #plot
      + geom_path(aes(x=WconDF$Time, y=WconDF[,paste0("X", freq0, "HZ")]), color="black", alpha = 0.5)
      + geom_path(aes(x=initialresults[,"time"], y=initialresults[,"eff2"]), color="green",size = 1.5)
      + geom_path(aes(x=finalresults[,"time"], y=finalresults[,"eff2"]), color="red",size = 1.5)
      +labs(list(title=paste0(freq0, " HZ Upper Responses 2Drugs"), x="Time (s)", y='Constriction'))
      +theme_bw()
      +xlim(0,75)
)
p30

#runs the model for all tissues and frequencies
Iteration <- function(con_list = c(1,2,5,7), dataDF = "con", normDF = "cap", lower = TRUE, filename = "/FormattedLowerTrachea/capsaicin/Results_"){
  freq_list <- c(0.1, 0.3,0.7,1, 3, 7, 10, 15, 30)
  for(r in con_list){                                                 #iterates through each raw file
    WconDF <- loadNormalizedDF(r, lower = lower, dataDF = dataDF, normDF = normDF) #load file
    for (q in freq_list){                                             #iterates through each frequency
      initialresults <- run_mod1(q, init_params1)              #find starting point model results from initial parameters
      finalparams <- final_drug_params(q, m = 500, WconDF, bestfit = bestfit, init_params, initialresults)

      finalparamsDF <- c(q, finalparams[nrow(finalparams),])          #take initial parameters and start vector for final values

      for(k in seq(1:100)){
        finalparams <- final_drug_params(q, m = 500, WconDF, bestfit = bestfit, init_params, initialresults)
        finalparamsDF <- rbind(finalparamsDF, c(q, finalparams[nrow(finalparams),]))
      }
      #append results to a single cvs file per index after each frequency
      finalparamsDF<-as.data.frame(finalparamsDF)
      write.table(finalparamsDF, paste0(filename, r, ".csv"), append = TRUE, sep = ",", dec = ".", qmethod = "double", col.names = FALSE)
    }
  }
}

#run facetgraph function. Plots all the frequencies on one graph.
facetgraph(WconDF, init_params=init_params, bestfit=bestfit, consensus = FALSE)
