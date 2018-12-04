#plot consensus value data against raw data. Consensus data is from model_2drugs.R



#                                   ---->>   currently set up to run for 2 drug model
# see the bottom for a modified facetgraph function


#
#
#

#

#


###

#


setwd("~/GitHub/NervePKPDmodel")
require(ggplot2)
source("normalizedDF.R")
source("defineModel.R")
source("runModelFunctions.R")

#produces graphable consensus data
bestfit <- c("KAunk", "KEunk", "DVunk")  #what unknowns are being solved for?
con_list = c(1,2,5,7) #List of indeces for the files used in this analysis

#Define a model (or more) to run and plot
thismodel <- defineModel(ACH_mod="complex", unk_mod="simple", effect_mod = "twoNT_2") 
thismodel[[1]]$model 
thismodel2 <- defineModel(ACH_mod="complex", unk_mod="none", effect_mod = "oneNT") 
thismodel2[[1]]$model 

parameterlist <- read.csv("FormattedLowerTrachea/finalParameters.csv") %>% select(-X)

# parameterlist <- read.csv("FormattedLowerTrachea/finalParameters.csv") %>% select(-X)
AUCresults <- NULL
SSresults <- NULL
maxresults <- NULL
#zeropoint <- c(0.096809,	0.103418,	0.103942,	0.096586) # cap 0.1 hz freq lower
zeropoint <- c(0.10312,	0.103568,	0.099246,	0.099837) # con 0.1 hz freq lower

#zeropoint <- c(0.103613789,	0.103665968,	0.103912957,	0.103762472,	0.104186936) #cap 0.1 Hz freq - Upper
#zeropoint <- c(0.095191713,	0.116669075,	0.116671121,	0.101521438,	0.103764105) # con 0.1 Hz freq - Upper

tissue_num <- 5
capsaicin_num <- 0
zeropoint1 <- zeropoint[3]
freq_list <- c(1, 3, 10) #0.1 Hz freq will be tissue specific

#First-order
consensus_params <- parameterlist[which(parameterlist$Capsaicin == capsaicin_num & parameterlist$Complex == 1 & parameterlist$Tissue == tissue_num),]
#Inhibition
init_params <- parameterlist[which(parameterlist$Capsaicin == 1 & parameterlist$Complex == 1 & parameterlist$Tissue == tissue_num),]
init_params$EC50unk <- 0
consensus_params$EC50unk <- 7.9


#WconDF <- loadNormalizedDF(tissue_num, lower = FALSE, dataDF = "con", normDF = "con")
WconDF <- loadNormalizedDF(tissue_num, lower = TRUE, dataDF = "con", normDF = "cap")
#facetgraph(conDF = WconDF, init_params= init_params, consensus_params = consensus_params, freq_list = freq_list, consensus = TRUE)
facetgraph(conDF = WconDF, init_params= init_params, consensus_params = consensus_params, bestfit = c("KAunk", "KEunk", "DVunk"), freq_list = freq_list, consensus = FALSE)

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

# find the max results
tissue_max <- as.data.frame(rbind(apply(WconDF[1:3001,], 2, max),
                                   apply(consensus_models, 2, max),
                                   apply(inhibition_models, 2, max))) %>%
  select(X0.1HZ, X0.3HZ, X1HZ, X3HZ, X10HZ) %>%
  mutate(Tissue = tissue_num, var = c("Control", "Simple", "Complex"), Capsaicin = capsaicin_num)
maxresults <- rbind(maxresults, tissue_max)

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

# write.csv(maxresults, "FormattedLowerTrachea/maxresults.csv")
# write.csv(AUCresults, "FormattedLowerTrachea/AUCresults.csv")
# write.csv(SSresults, "FormattedLowerTrachea/SSresults.csv")

facetgraph <- function(conDF, init_params, consensus_params, bestfit, 
                       freq_list = c(0.1, 0.3, 0.7, 1, 3, 10, 15, 30), 
                       consensus = FALSE, chosenmodel = thismodel, ncolumns = 3){
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
      finalparams <- final_drug_params(stim_freq = working_freq, m = 1000, conDF, bestfit = bestfit, init_params, initialresults, hyper_params = c(2, 0.1))
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
          + facet_wrap(~ Freq, ncol=ncolumns)
          + scale_y_continuous(limits = c(0, 1))
          + theme_bw()
    )
  } else{
    colnames(facetDF)<-c("Freq", "Time", "Raw", "Initial", "Final")
    p <- (ggplot(facetDF)
          + geom_line(aes(x=Time, y=Raw), color="black", alpha = 0.5)
          #+ geom_line(aes(x=Time, y=Initial), color="blue",size = 1)
          + geom_line(aes(x=Time, y=Final), color="red",size = 1)
          + facet_wrap(~ Freq, ncol=ncolumns)
          + theme_bw()
    )
  }
  
  p
}

