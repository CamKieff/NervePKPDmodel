#plot consensus value data against raw data. Consensus data is from model_2drugs.R
require(ggplot2)
source("normalizedDF.R")
source("defineModel.R")

#produces graphable consensus data
freq_list <- c(0.1, 0.3,0.7,1, 3, 7, 10, 15, 30)
consensus_params <- list(KAach = 0.1050,
                         KEach = 0.6771,
                         DVach = 5.8202,
                         EC50ach = 5.383,

                         m2max = 0.9798,
                         chemax = 0.1475,
                         IC50m2 = 6.4064,
                         IC50che = 6.1431,

                         KAunk = 1, KEunk = 1, DVunk = 7, EC50unk = 5,MAXunk = 1
                         )

#Define a model (or more) to run and plot
thismodel <- defineModel(ACH_mod="complex", unk_mod="none", effect_mod = "oneNT") #what model
thismodel[[1]]$model #verifty the model is correct

bestfit <- c("KAach", "KEach", "DVach")  #what unknowns are being solved for?
con_list = c(1,2,5,7) #List of indeces for the files used in this analysis
freq_list <- c(0.1, 0.3,0.7,1, 3, 7, 10, 15, 30)

WconDF <- loadNormalizedDF(7, lower = TRUE, dataDF = "cap", normDF = "cap")
initialresults <- run_mod1(stim_freq = freq_list[1], init_params, chosenmodel = thismodel)

consensus_models <- initialresults[,'time']
for (q in freq_list){
  initialresults <- run_mod1(q, init_params, chosenmodel = thismodel) #run model with consensus values
  consensus_models <- cbind(consensus_models, initialresults[,'eff2']) #append freq results
}

consensus_models<-as.data.frame(consensus_models)
names(consensus_models)<-names(WconDF)

#find the sum of squares differences and the Area under the curve (AUC) percent error
SS1 <- colSums((WconDF[1:3001,] - consensus_models)^2)
AUCerror <- (colSums((WconDF[1:3001,]*0.02))-colSums((consensus_models*0.02)))/colSums((WconDF[1:3001,]*0.02))

facetgraph(conDF = WconDF, init_params = init_params, freq_list = c(0.1, 0.3, 1, 3, 10), consensus = TRUE)
