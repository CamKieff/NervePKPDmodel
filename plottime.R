#here we plot consensus value data against raw data. Consensus data is from controlModel_2drugs_uppertrach.R
require(ggplot2)
source("NerveModel1/normalizedDF.R")

#path for consensus file, make sure you gave it a header
consensus_con <- read.csv("NerveModel1/Formatted Lower Trachea/control/consensus_2drugs.csv", header = TRUE)
consensus_test <- read.csv("NerveModel1/Formatted Lower Trachea/control/consensus_2drugs_test1.csv", header = TRUE)

freq <- 30
con_list = c(1,2,5,7) #List of indeces for the files used in this analysis

#loads dataframes, normalizes them, and binds data for that frequency to a single DF for plotting
preparePlots <- function(freq, con_list = c(1,2,5,7), lower = TRUE, dataDF = "con", normDF = "cap"){
  plottime <- loadNormalizedDF(con_list[1], lower = lower, dataDF = dataDF, normDF = normDF)[,"Time"] #this could probably just be a seq
  for (y in con_list){
    plottime <- cbind(plottime, loadNormalizedDF(y, lower = lower, dataDF = dataDF, normDF = normDF)[,paste0("X", freq, "HZ")])
  }
  return(plottime)
}

plotStuff <- preparePlots(freq = freq)

p10<- (ggplot()
       + geom_path(aes(x=plotStuff[,1], y=plotStuff[,2]), color="black", alpha = 0.5)
       + geom_path(aes(x=plotStuff[,1], y=plotStuff[,3]), color="black", alpha = 0.5)
       + geom_path(aes(x=plotStuff[,1], y=plotStuff[,4]), color="black", alpha = 0.5)
       + geom_path(aes(x=plotStuff[,1], y=plotStuff[,5]), color="black", alpha = 0.5)
       #+ geom_path(aes(x=plotStuff[,1], y=plotStuff[,6]), color="black", alpha = 0.5)
       + geom_path(aes(x=consensus_con[,"Time"], y=consensus_con[,paste0("X", freq, "HZ")]), color="green",size = 1.5)
       + geom_path(aes(x=consensus_test[,"Time"], y=consensus_test[,paste0("X", freq, "HZ")]), color="red",size = 1.5)
       + geom_path(aes(x=consensus_test[,"Time"], y=(consensus_test[,paste0("X", freq, "HZ")]+consensus_con[,paste0("X", freq, "HZ")])), color="blue",size = 1.5)
       +labs(list(title=paste0(freq, " HZ Responses 2Drugs"), x="Time (s)", y='Constriction'))
       +theme_bw()
       +theme(title=element_text(size=18, face="bold"),axis.title=element_text(size=14, face="bold"),axis.text=element_text(size=14))
)
p10

freq_list <- c(0.1, 0.3, 0.7, 1, 3, 5, 7, 10, 15, 30)

#calculates percent error compared to one model
calculate_perror <- function(consensus_file = consensus_con, con_list = c(1,2,5,7), dataDF = "con", normDF = "cap"){
  perror <- NULL
  for (y in con_list){
    perrorDF <- loadNormalizedDF(y, lower = lower, dataDF = dataDF, normDF = normDF)
    for(q in freq_list){
      perror <- append(perror, (sum(consensus_file[,paste0("X", q, "HZ")])-sum(perrorDF[1:3001]))/sum(perrorDF[1:3001])*100)
    }
  }
  perror <- cbind(rep(freq_list, length(con_list)), rep(con_list, each=9), AUC)

  return(perror)
}

#calculates AUC for raw data and consensus model
calculate_AUC <- function(consensus_file = consensus_con, con_list = c(1,2,5,7), dataDF = "con", normDF = "cap"){
  AUC <- NULL
  for (y in con_list){
    AUCDF <- loadNormalizedDF(y, lower = lower, dataDF = dataDF, normDF = normDF)
    for(q in freq_list){
      AUC <- append(AUC, (sum(AUCDF[,paste0("X", q, "HZ")]))
    }
  }
  for(q in freq_list){
    AUC <- append(AUC, sum(consensus_file[,paste0("X", q, "HZ")])
  }
  AUC <- cbind(rep(freq_list, (length(con_list)+1)), c(rep(con_list, each=9), rep(con, 9)), AUC)
  return(AUC)
}

write.csv(AUCresults, paste0(directory_name, "AUC_results.csv"))
