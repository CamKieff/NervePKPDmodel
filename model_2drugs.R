#This script iteratively solves a system of ODEs and finds parameters
#the model is for two neurotransmistters, one ACH, another unknown relaxant

#RxODE
library(RxODE)
library(ggplot2)
source("normalizedDF.R")

#choose a set of ODEs that describe the desired model
prepareModel <- function(mod_num = 1){
  if(mod_num == 1){
    #statistical summation 2-drug model
    ODE1 <- "
    d/dt(depot1) = -KA1*depot1;
    d/dt(centr1) = KA1*depot1 - KE1*centr1;

    d/dt(depot2) = -KA2*depot2;
    d/dt(centr2) = KA2*depot2 - KE2*centr2;

    effa = max1*centr1/(EC501 + centr1);
    effb = centr2/(EC502 + centr2);
    eff2 = effb - effa + (effa * effb);
    "
  }  else if(mod_num == 2){
    #decrease maximum 2-drug model
    ODE1 <- "
    d/dt(depot1) = -KA1*depot1;
    d/dt(centr1) = KA1*depot1 - KE1*centr1;

    d/dt(depot2) = -KA2*depot2;
    d/dt(centr2) = KA2*depot2 - KE2*centr2;

    eff1 = max1*centr1/(EC501 + centr1);
    eff2 = (1-eff1)*centr2/(EC502 + centr2);
    "
  }

  mod1 <- RxODE(model = ODE1, modName = "mod1") #Compile Model
  ODEinits <- c(depot1 = 0, centr1 = 0, depot2 = 0, centr2 = 0)  #initial conditions

  return(list(mod1, ODEinits))
}
thismodel <- prepareModel(1)

#a function that takes parameters and runs the model
run_2drugs_mod1 <- function(stim_freq, final_mod_params, NT=2){
  if (NT==2){
    #(KAunk, KEunk, logDVunk, logEC50unk, MAXunk)
    mod_params <- c(KA1 = final_mod_params[1],
                    KE1 = final_mod_params[2],
                    KA2 = KAach,
                    KE2 = KEach,
                    max1 = final_mod_params[5],
                    EC501 = 10^-final_mod_params[4],
                    EC502 = 10^-EC50ach
    )
    time_var <- 60 #length of stimulation in seconds
    pulse_rate <- 1/stim_freq
    num_doses <- time_var/pulse_rate
    testev <- eventTable(amount.units="mol", time.unit="seconds")
    testev$add.dosing(dose = 10^-final_mod_params[3], nbr.doses = num_doses, dosing.interval = pulse_rate, dosing.to = 1, start.time = 0)
    testev$add.dosing(dose = 10^-DVach, nbr.doses = num_doses, dosing.interval = pulse_rate, dosing.to = 3, start.time = 0)
    testev$add.sampling(seq(from = 0, to = time_var, by = 0.02))

    finalres <- thismodel[[1]]$run(mod_params, testev, thismodel[[2]])
  } else if(NT==1){
    mod_params <- c(KA2 = final_mod_params[1],
                    KE2 = KEach1,
                    KE3 = KEach2,
                    EC502 = EC50ach
    )
    time_var <- 60 #length of stimulation in seconds
    if(stim_freq == 0.1){
      pulse_rate <- 1/stim_freq
    } else{
      pulse_rate <- 1/stim_freq
    }
    num_doses <- time_var/pulse_rate
    testev <- eventTable(amount.units="mol", time.unit="seconds")
    testev$add.dosing(dose = 10^-DVach, nbr.doses = num_doses, dosing.interval = pulse_rate, dosing.to = 1, start.time = 0)
    testev$add.sampling(seq(from = 0, to = time_var, by = 0.02))

    finalres <- mod1$run(mod_params, testev, ODEinits)
  }
  return(finalres)
}

#returns a dataframe of chosen/tested parameters from the model
#conDF should be formatted using loadNormalizedDF (see below)
final_2Dparams <-function(stim_freq, m = 50, conDF, init_params, init_model, NT=2, chosen=TRUE){

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
    for(i in 1:length(init_params)){
      testparams <- newparams
      testparams[i] <- abs(rnorm(n = 1, mean = newparams[i], sd = (sqrt(init_params[i]^2)*lambda)))
      tested_params <- rbind(tested_params, testparams)

      #calls run_2drugs_mod1 to run the model
      iteration <- run_2drugs_mod1(stim_freq, testparams, NT=NT)

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

#Define initial parameters - > ACh consensus values here
##Decide on format here
KAach <- 0.2417
KEach <- 0.5268
DVach <- 5.902
EC50ach <- 5.383
init_params_ach <- list(KAach = 0.2417,
                        KEach = 0.5268,
                        DVach = 5.902)

#Relaxation initial parameters
KAunk <- 1
KEunk <- 1
DVunk <- 7
EC50unk <- 5
MAXunk <- 1
init_params_unk <- list(KAunk = 1, KEunk = 1, DVunk = 7, EC50unk = 5, MAXunk = 1)

#Unknown parameters in order (KAunk, KEunk, DVunk, EC50unk, MAXunk)
init_params1 <- c(KAunk, KEunk, DVunk, EC50unk, MAXunk)
freq_list <- c(0.1, 0.3, 0.7, 1, 3, 7, 10 , 15, 30)

working_freq <- 30
WconDF <- loadNormalizedDF(1, lower = FALSE, dataDF = "con", normDF = "con")
initialresults <- run_2drugs_mod1(stim_freq = working_freq, init_params1)
finalparams <- final_2Dparams(stim_freq = working_freq, m = 500, WconDF, init_params1, initialresults)
finalresults <- run_2drugs_mod1(stim_freq = working_freq, finalparams[nrow(finalparams),])

p30<- (ggplot() #plot
      + geom_path(aes(x=WconDF$Time, y=WconDF[,paste0("X", working_freq, "HZ")]), color="black", alpha = 0.5)
      + geom_path(aes(x=initialresults[,"time"], y=initialresults[,"eff2"]), color="green",size = 1.5)
      + geom_path(aes(x=finalresults[,"time"], y=finalresults[,"eff2"]), color="red",size = 1.5)
      +labs(list(title=paste0(working_freq, " HZ Upper Responses 2Drugs"), x="Time (s)", y='Constriction'))
      +theme_bw()
      +xlim(0,75)
)
p30

massIteration <- function(con_list = c(1,2,5,7), dataDF = "con", normDF = "cap", lower = TRUE, filename = "Results_"){
  freq_list <- c(0.1, 0.3,0.7,1, 3, 7, 10, 15, 30)
  for(r in con_list){                                                 #iterates through each raw file
    for (q in freq_list){                                             #iterates through each frequency
      WconDF <- loadNormalizedDF(r, lower = lower, dataDF = dataDF, normDF = normDF) #load file
      initialresults <- run_2drugs_mod1(init_params1, q)              #find starting point model results from initial parameters
      finalparams <- final_2Dparams(q, m=500, WconDF, init_params1, initialresults)

      finalparamsDF <- c(q, finalparams[nrow(finalparams),])          #take initial parameters and start vector for final values

      for(k in seq(1:100)){
        finalparams <- final_2Dparams(q, m=500, WconDF, init_params1, initialresults)
        finalparamsDF <- rbind(finalparamsDF, c(q, finalparams[nrow(finalparams),]))
      }

      finalparamsDF<-as.data.frame(finalparamsDF)
      write.table(finalparamsDF, paste0(con_directory_name, filename, r, ".csv"), append = TRUE, sep = ",", dec = ".", qmethod = "double", col.names = FALSE)
    }
  }
}

#produces graphable consensus data
consensus_models <- initialresults[,'time']
for (q in freq_list){
  #ACh consensus values here
  KAach <- 0.2417
  KEach <- 0.5268
  DVach <- 10^-5.902
  EC50ach <- 10^-5.383

  #Relaxation consensus values
  KAunk <- 0.9661
  KEunk <- 10.8071
  DVunk <- 7.4434
  EC50unk <- 4.048
  MAXunk <- 0.861135
  init_params1 <- c(KAunk, KEunk, DVunk, EC50unk, MAXunk)

  initialresults <- run_2drugs_mod1(init_params1, q) #run model with consensus values
  consensus_models <- cbind(consensus_models, initialresults[,'eff2']) #append freq results

}
write.csv(consensus_models, file = paste0(con_directory_name, "consensus_2drugs.csv"))

