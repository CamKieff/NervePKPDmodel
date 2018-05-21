#plot consensus value data against raw data. Consensus data is from model_2drugs.R
require(ggplot2)
source("normalizedDF.R")
source("defineModel.R")
source("runModelFunctions.R")

#produces graphable consensus data
bestfit <- c("KAach", "KEach", "DVach")  #what unknowns are being solved for?
con_list = c(1,2,5,7) #List of indeces for the files used in this analysis

#Define a model (or more) to run and plot
thismodel <- defineModel(ACH_mod="complex", unk_mod="none", effect_mod = "oneNT") 
thismodel[[1]]$model 
thismodel2 <- defineModel(ACH_mod="simple", unk_mod="none", effect_mod = "oneNT") 
thismodel2[[1]]$model 

parameterlist <- read.csv("FormattedLowerTrachea/finalParameters.csv") %>% select(-X)
AUCresults <- NULL
SSresults <- NULL

#0.10312	0.103568	0.099246	0.099837
tissue_num <- 7
capsaicin_num <- 0
zeropoint1 <- 0.099837
freq_list <- c(zeropoint1, 0.3,1, 3, 10) #0.1 Hz freq will be tissue specific

#First-order
consensus_params <- parameterlist[which(parameterlist$Capsaicin == capsaicin_num & parameterlist$Complex == 0 & parameterlist$Tissue == tissue_num),]
#Inhibition
init_params <- parameterlist[which(parameterlist$Capsaicin == capsaicin_num & parameterlist$Complex == 1 & parameterlist$Tissue == tissue_num),]

WconDF <- loadNormalizedDF(tissue_num, lower = TRUE, dataDF = "con", normDF = "cap")
facetgraph(conDF = WconDF, init_params= init_params, consensus_params = consensus_params, freq_list = freq_list, consensus = TRUE)

#run the code below this line to find the AUC for the raw data, and both models
freq_list <- c(zeropoint1, 0.3,7,1, 3,7, 10,15,30) #0.1 Hz freq will be tissue specific
initialresults <- run_mod1(stim_freq = freq_list[1], init_params, chosenmodel = thismodel)

inhibition_models <- initialresults[,'time']
consensus_models <- initialresults[,'time']
for (q in freq_list){
  inhibitionresults <- run_mod1(q, init_params, chosenmodel = thismodel) #run model with consensus values
  inhibition_models <- cbind(inhibition_models, inhibitionresults[,'eff2']) #append freq results\

  consensusresults <- run_mod1(q, consensus_params, chosenmodel = thismodel2) #run model with consensus values
  consensus_models <- cbind(consensus_models, consensusresults[,'eff2']) #append freq results
}
inhibition_models<-as.data.frame(inhibition_models)
names(inhibition_models)<-names(WconDF)

consensus_models<-as.data.frame(consensus_models)
names(consensus_models)<-names(WconDF)

# find the AUC results
tissue_sums <- as.data.frame(rbind(colSums(WconDF[1:3001,]*0.02),
colSums(consensus_models*0.02),
colSums(inhibition_models*0.02))) %>%
  select(X0.1HZ, X0.3HZ, X1HZ, X3HZ, X10HZ) %>%
  mutate(Tissue = tissue_num, var = c("Control", "Simple", "Complex"), Capsaicin = capsaicin_num)
AUCresults <- rbind(AUCresults, tissue_sums)

#find the sum of squared differences
SS_sums <- as.data.frame(rbind(colSums((WconDF[1:3001,] - consensus_models)^2),
colSums((WconDF[1:3001,] - inhibition_models)^2))) %>%
  select(X0.1HZ, X0.3HZ, X1HZ, X3HZ, X10HZ) %>%
  mutate(Tissue = tissue_num, var = c("Simple", "Complex"), Capsaicin = capsaicin_num)

SSresults <- rbind(SSresults, SS_sums)

# write.csv(AUCresults, "FormattedLowerTrachea/AUCresults.csv")
# write.csv(SSresults, "FormattedLowerTrachea/SSresults.csv")