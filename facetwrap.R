require(ggplot2)
require(reshape2)
source("NerveModel1/normalizedDF.R")

freq_list <- c(0.1, 0.3, 0.7, 1, 3, 7, 10 , 15, 30)

#make sure that Ach and initial parameters are correct for the model being run
facetDF <-NULL
for (i in freq_list){
  working_freq <- i
  WconDF <- loadNormalizedDF(1, lower = TRUE, dataDF = "con", normDF = "cap")
  initialresults <- run_2drugs_mod1(init_params1, stim_freq = working_freq)
  finalparams <- final_2Dparams(stim_freq = working_freq, m = 500, WconDF, init_params1, initialresults)
  finalresults <- run_2drugs_mod1(finalparams[nrow(finalparams),], stim_freq = working_freq)
  plotDF <- data.frame(rep(working_freq, length(finalresults)), initialresults[,"time"], WconDF[1:length(initialresults[,"time"]),paste0("X", working_freq, "HZ")], initialresults[,"eff2"], finalresults[,"eff2"])

  facetDF <- rbind(plotDF, facetDF)
}

colnames(facetDF)<-c("Freq", "Time", "Raw", "Initial", "Final")
#facetDF <- as.data.frame(facetDF)

p <- (ggplot(facetDF)
  + geom_line(aes(x=Time, y=Raw), color="black", alpha = 0.5)
  + geom_line(aes(x=Time, y=Initial), color="green",size = 1)
  + geom_line(aes(x=Time, y=Final), color="red",size = 1)
  + facet_wrap(~ Freq, scales="free", ncol=3)
  + theme_bw()
)
p
