#plot consensus value data against raw data. Consensus data is from model_2drugs.R
require(ggplot2)
source("normalizedDF.R")

#produces graphable consensus data
freq_list <- c(0.1, 0.3,0.7,1, 3, 7, 10, 15, 30)
consensus_models <- initialresults[,'time']
for (q in freq_list){
  init_params <- list(KAach = 0.2417,
                      KEach = 0.5268,
                      DVach = 5.902,
                      EC50ach = 5.383,

                      KAunk = 0.9661,
                      KEunk = 10.8071,
                      DVunk = 7.4434,
                      EC50unk = 4.048,
                      MAXunk = 0.861135)

  initialresults <- run_mod1(NT = 2, q, init_params) #run model with consensus values
  consensus_models <- cbind(consensus_models, initialresults[,'eff2']) #append freq results
}
#write.csv(consensus_models, file = paste0(con_directory_name, "consensus_2drugs.csv"))

#path for consensus file(s), make sure it has a header
consensus_con <- read.csv("FormattedLowerTrachea/control/consensus_2drugs.csv", header = TRUE)
consensus_test <- read.csv("FormattedLowerTrachea/control/consensus_2drugs_test1.csv", header = TRUE)

freq0 <- 30
con_list = c(1,2,5,7) #List of indeces for the files used in this analysis

#loads dataframes, normalizes them, and binds data for that frequency to a single DF for plotting
preparePlots <- function(freq, con_list = c(1,2,5,7), lower = TRUE, dataDF = "con", normDF = "cap"){
  plottime <- loadNormalizedDF(con_list[1], lower = lower, dataDF = dataDF, normDF = normDF)[,"Time"] #this could probably just be a seq
  for (y in con_list){
    plottime <- cbind(plottime, loadNormalizedDF(y, lower = lower, dataDF = dataDF, normDF = normDF)[,paste0("X", freq, "HZ")])
  }
  return(plottime)
}

plotStuff <- preparePlots(freq = freq0)
#
#
#this graph is prettier, but in some ways has been replaced by facetgraph in runModelFunctions.R
#
#
p10<- (ggplot()
       + geom_path(aes(x=plotStuff[,1], y=plotStuff[,2]), color="black", alpha = 0.5)
       + geom_path(aes(x=plotStuff[,1], y=plotStuff[,3]), color="black", alpha = 0.5)
       + geom_path(aes(x=plotStuff[,1], y=plotStuff[,4]), color="black", alpha = 0.5)
       + geom_path(aes(x=plotStuff[,1], y=plotStuff[,5]), color="black", alpha = 0.5)
       #+ geom_path(aes(x=plotStuff[,1], y=plotStuff[,6]), color="black", alpha = 0.5)
       + geom_path(aes(x=consensus_con[,"Time"], y=consensus_con[,paste0("X", freq, "HZ")]), color="green",size = 1.5)
       + geom_path(aes(x=consensus_test[,"Time"], y=consensus_test[,paste0("X", freq, "HZ")]), color="red",size = 1.5)
       +labs(list(title=paste0(freq, " HZ Response"), x="Time (s)", y='Constriction'))
       +theme_bw()
       +theme(title=element_text(size=18, face="bold"),axis.title=element_text(size=14, face="bold"),axis.text=element_text(size=14))
)
p10


#calculates percent error and AUC for one consensus file
testConsensusFit <- function(consensus_file = consensus_con, con_list = c(1,2,5,7), dataDF = "con", normDF = "cap"){
  freq_list <- c(0.1, 0.3, 0.7, 1, 3, 5, 7, 10, 15, 30)
  perror <- NULL
  AUC <- NULL
  for (y in con_list){
    workingDF <- loadNormalizedDF(y, lower = lower, dataDF = dataDF, normDF = normDF)
    for(q in freq_list){
      AUC <- append(AUC, (sum(workingDF[,paste0("X", q, "HZ")]))
      perror <- append(perror, (sum(consensus_file[,paste0("X", q, "HZ")])-sum(workingDF[1:3001]))/sum(workingDF[1:3001])*100)
    }
  }

  for(q in freq_list){
    AUC <- append(AUC, sum(consensus_file[,paste0("X", q, "HZ")])
  }

  perror <- cbind(rep(freq_list, length(con_list)), rep(con_list, each=9), perror)
  AUC <- cbind(rep(freq_list, (length(con_list)+1)), c(rep(con_list, each=9), rep(con, 9)), AUC)

  return(list(perror=perror,AUC=AUC)
}
