#plot consensus value data against raw data. Consensus data is from model_2drugs.R
require(ggplot2)
source("normalizedDF.R")
source("defineModel.R")

#produces graphable consensus data
bestfit <- c("KAach", "KEach", "DVach")  #what unknowns are being solved for?
con_list = c(1,2,5,7) #List of indeces for the files used in this analysis

#Define a model (or more) to run and plot
thismodel <- defineModel(ACH_mod="complex", unk_mod="none", effect_mod = "oneNT") 
thismodel[[1]]$model 
thismodel2 <- defineModel(ACH_mod="simple", unk_mod="none", effect_mod = "oneNT") 
thismodel2[[1]]$model 

parameterlist <- read.csv("FormattedLowerTrachea/finalParameters.csv") %>% select(-X)

tissue_num <- 1
freq_list <- c(0.103, 0.3,1, 3, 10) #0.1 Hz freq will be tissue specific

#First-order
consensus_params <- parameterlist[which(parameterlist$Capsaicin == 0 & parameterlist$Complex == 0 & parameterlist$Tissue == tissue_num),]
#Inhibition
init_params <- parameterlist[which(parameterlist$Capsaicin == 0 & parameterlist$Complex == 1 & parameterlist$Tissue == tissue_num),]

WconDF <- loadNormalizedDF(tissue_num, lower = TRUE, dataDF = "con", normDF = "cap")
facetgraph(conDF = WconDF, init_params= init_params, consensus_params = consensus_params, freq_list = freq_list, consensus = TRUE)

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

colSums(WconDF[1:3001,]*0.02)
colSums(consensus_models*0.02)
colSums(inhibition_models*0.02)

#find the sum of squares differences and the Area under the curve (AUC) percent error
colSums((WconDF[1:3001,] - consensus_models)^2)
colSums((WconDF[1:3001,] - inhibition_models)^2)

AUCerror <- (colSums((WconDF[1:3001,]*0.02))-colSums((consensus_models*0.02)))/colSums((WconDF[1:3001,]*0.02))


facetgraph <- function(conDF, init_params, consensus_params, bestfit, freq_list = c(0.1, 0.3, 0.7, 1, 3, 10, 15, 30), consensus = FALSE, chosenmodel = thismodel){
  facetDF <-NULL
  for (i in freq_list){
    working_freq <- i
    if(working_freq < 0.2){
      initialresults <- run_mod1(stim_freq = working_freq, init_params, chosenmodel = thismodel)
      consensusresults <- run_mod1(stim_freq = working_freq, consensus_params, chosenmodel = thismodel2)
      working_freq <- 0.1
    } else{
      initialresults <- run_mod1(stim_freq = working_freq, init_params, chosenmodel = thismodel)
      consensusresults <- run_mod1(stim_freq = working_freq, consensus_params, chosenmodel = thismodel2)

    }
    if(consensus == TRUE){
      plotDF <- data.frame(rep(working_freq, length(initialresults)), initialresults[,"time"], conDF[1:length(initialresults[,"time"]),paste0("X", working_freq, "HZ")], initialresults[,"eff2"],consensusresults[,"eff2"])
    } else{
      finalparams <- final_drug_params(stim_freq = working_freq, m = 500, conDF, bestfit = bestfit, init_params, initialresults)
      finalresults <- run_mod1(stim_freq = working_freq, finalparams[nrow(finalparams),], chosenmodel = chosenmodel)
      plotDF <- data.frame(rep(working_freq, length(finalresults)), initialresults[,"time"], conDF[1:length(initialresults[,"time"]),paste0("X", working_freq, "HZ")], initialresults[,"eff2"], finalresults[,"eff2"])
    }

    facetDF <- rbind(plotDF, facetDF)
  }
  if(consensus == TRUE){
    colnames(facetDF)<-c("Freq", "Time", "Raw", "Initial", "Consensus")
    p <- (ggplot(facetDF)
          + geom_line(aes(x=Time, y=Raw), color="black", alpha = 0.5)
          + geom_line(aes(x=Time, y=Consensus), color="blue",size = 1)
          + geom_line(aes(x=Time, y=Initial), color="red",size = 1)
          #+ facet_wrap(~ Freq, scales="free", ncol=3)
          + facet_wrap(~ Freq, ncol=3)
          + scale_y_continuous(limits = c(0, 1))
          + theme_bw()
    )
  } else{
    colnames(facetDF)<-c("Freq", "Time", "Raw", "Initial", "Final")
    p <- (ggplot(facetDF)
          + geom_line(aes(x=Time, y=Raw), color="black", alpha = 0.5)
          + geom_line(aes(x=Time, y=Initial), color="green",size = 1)
          + geom_line(aes(x=Time, y=Final), color="red",size = 1)
          + facet_wrap(~ Freq, scales="free", ncol=3)
          + theme_bw()
    )
  }

  p
}
