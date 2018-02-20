#This script iteratively solves a system of ODEs and finds parameters
#the model is defined in "defineModel.R" It can work with a single ACH-like neurotransmitter (NT)
#or it can use two neurotransmitters, one constrictor (acetylcholine; ACH) and one relaxant (unknown; unk)

#RxODE
library(RxODE)
library(ggplot2)
source("normalizedDF.R")
source("defineModel.R")

thismodel <- defineModel(ACH_mod="simple", unk_mod="none", effect_mod = "oneNT")

#a function that takes parameters and runs the model
run_mod1 <- function(NT=2, stim_freq, mod_params){
  parameters <- c(KA1 = mod_params$KAunk,
                  KE1 = mod_params$KEunk,
                  KA2 = mod_params$KAach,
                  KE2 = mod_params$KEach,
                  max1 = mod_params$MAXunk,
                  EC501 = 10^-mod_params$EC50unk,
                  EC502 = 10^-mod_params$EC50ach
              )
  time_var <- 60 #length of stimulation in seconds

   if (NT == 2){
    pulse_rate <- 1/stim_freq
    num_doses <- time_var/pulse_rate
    testev <- eventTable(amount.units="mol", time.unit="seconds")
    testev$add.dosing(dose = 10^-mod_params$DVunk, nbr.doses = num_doses, dosing.interval = pulse_rate, dosing.to = 1, start.time = 0)
    testev$add.dosing(dose = 10^-mod_params$DVach, nbr.doses = num_doses, dosing.interval = pulse_rate, dosing.to = 3, start.time = 0)
    testev$add.sampling(seq(from = 0, to = time_var, by = 0.02))
  } else if(NT == 1){

    if(stim_freq == 0.1){
      pulse_rate <- 1/stim_freq
    } else{
      pulse_rate <- 1/stim_freq
    }
    num_doses <- time_var/pulse_rate
    testev <- eventTable(amount.units="mol", time.unit="seconds")
    testev$add.dosing(dose = 10^-mod_params$DVach, nbr.doses = num_doses, dosing.interval = pulse_rate, dosing.to = 1, start.time = 0)
    testev$add.sampling(seq(from = 0, to = time_var, by = 0.02))
  }

  finalres <- thismodel[[1]]$run(parameters, testev, thismodel[[2]])
  return(finalres)
}

#returns a dataframe of chosen/tested parameters from the model
#conDF should be formatted using loadNormalizedDF
final_drug_params <-function(NT=2, stim_freq, m = 50, conDF, bestfit, init_params, init_model, chosen=TRUE){

  testexp <- 4     #sum of "squares" exponenet
  lambda <- 2    #learning rate

  #creates workingdata for the model comparison
  workingdata <- data.frame(init_model[,"time"], init_model[,"eff2"], conDF[1:length(init_model[,"time"]),paste0("X", stim_freq, "HZ")])
  names(workingdata)<- c("Time", "Model", "Data")
  workingdata$SS <- (workingdata$Model - workingdata$Data)^testexp

  newparams <- init_params
  chosen_params <- c(init_params,0)
  tested_params <- init_params

  for (j in 1:m){
    for(i in bestfit){
      testparams <- newparams
      testparams[[i]] <- abs(rnorm(n = 1, mean = newparams[[i]], sd = (sqrt(init_params[[i]]^2)*lambda)))
      tested_params <- rbind(tested_params, testparams)

      #calls run_mod1 to run the model
      iteration <- run_mod1(NT=NT, stim_freq, testparams)

      #calculate Sum of Squares Objective function for the new model
      workingdata$newmodel <- iteration[1:length(workingdata$Model),"eff2"]
      workingdata$newSS <- with(workingdata, (newmodel - Data)^testexp)

      #test if the SS for the new model is better than the old model
      #if it is use that test parameter in future runs; update the SS values to the new ones
      if(mean(workingdata$newSS, na.rm=TRUE) < mean(workingdata$SS, na.rm=TRUE)){
        newparams[i] <- testparams[i]
        workingdata$SS <- workingdata$newSS
        chosen_params <- rbind(chosen_params, c(newparams, mean(workingdata$SS, na.rm=TRUE)))
      } else{
        chosen_params <- rbind(chosen_params, c(newparams, mean(workingdata$SS, na.rm=TRUE)))
      }
    }
  }

  print(chosen_params[nrow(chosen_params),]) #prints the final parameters
  if(chosen){
    return(chosen_params)
  } else {
    return(tested_params)
  }
}

#Define initial parameters.
#for 2 NT: ACh consensus values here followed by initial unknown parameters
init_params <- list(KAach = 0.2417,
                    KEach = 0.5268,
                    DVach = 5.902,
                    EC50ach = 5.383,

                    KAunk = 1,
                    KEunk = 1,
                    DVunk = 7,
                    EC50unk = 5,
                    MAXunk = 1)

#Define what model we are running
freq0 <- 1                               #what frequency
num_NT <- 1                              #how many neurotransmitters (should be informed by the model used)
bestfit <- c("KAach", "KEach", "DVach")  #what unknowns are being solving for?

WconDF <- loadNormalizedDF(1, lower = TRUE, dataDF = "cap", normDF = "cap")
initialresults <- run_mod1(NT = num_NT, stim_freq = freq0, init_params)
finalparams <- final_drug_params(NT = num_NT, stim_freq = freq0, m = 500, WconDF, bestfit = bestfit, init_params, initialresults)
finalresults <- run_mod1(NT = num_NT, stim_freq = freq0, finalparams[nrow(finalparams),])

p30<- (ggplot() #plot
      + geom_path(aes(x=WconDF$Time, y=WconDF[,paste0("X", freq0, "HZ")]), color="black", alpha = 0.5)
      + geom_path(aes(x=initialresults[,"time"], y=initialresults[,"eff2"]), color="green",size = 1.5)
      + geom_path(aes(x=finalresults[,"time"], y=finalresults[,"eff2"]), color="red",size = 1.5)
      +labs(list(title=paste0(freq0, " HZ Upper Responses 2Drugs"), x="Time (s)", y='Constriction'))
      +theme_bw()
      +xlim(0,75)
)
p30

massIteration <- function(con_list = c(1,2,5,7), dataDF = "con", normDF = "cap", lower = TRUE, filename = "Results_"){
  freq_list <- c(0.1, 0.3,0.7,1, 3, 7, 10, 15, 30)
  for(r in con_list){                                                 #iterates through each raw file
    for (q in freq_list){                                             #iterates through each frequency
      WconDF <- loadNormalizedDF(r, lower = lower, dataDF = dataDF, normDF = normDF) #load file
      initialresults <- run_mod1(init_params1, q)              #find starting point model results from initial parameters
      finalparams <- final_drug_params(q, m=500, WconDF, init_params1, initialresults)

      finalparamsDF <- c(q, finalparams[nrow(finalparams),])          #take initial parameters and start vector for final values

      for(k in seq(1:100)){
        finalparams <- final_drug_params(q, m=500, WconDF, init_params1, initialresults)
        finalparamsDF <- rbind(finalparamsDF, c(q, finalparams[nrow(finalparams),]))
      }

      finalparamsDF<-as.data.frame(finalparamsDF)
      write.table(finalparamsDF, paste0(con_directory_name, filename, r, ".csv"), append = TRUE, sep = ",", dec = ".", qmethod = "double", col.names = FALSE)
    }
  }
}

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
write.csv(consensus_models, file = paste0(con_directory_name, "consensus_2drugs.csv"))