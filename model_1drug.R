#This script iteratively solves a system of ODEs and finds parameters
#the model is for two neurotransmistters, one ACH, another unknown relaxant

#RxODE
library(RxODE)
library(ggplot2)
source("NerveModel1/normalizedDF.R")
#cap_directory_name <- 'NerveModel1/FormattedUpperTrachea/capsaicin/'
#con_directory_name <- 'NerveModel1/FormattedUpperTrachea/control/'

#define model control, two neurotransmitters 1-CM
ODE1 <- "
d/dt(depot2) = -KA2*depot2 + KE3*centr2;
d/dt(centr2) = KA2*depot2 - KE2*centr2 - KE3*centr2;
eff2 = centr2/(EC502 + centr2);
"

#Compile Model
mod1 <- RxODE(model = ODE1, modName = "mod1")
ODEinits <- c(depot2 = 0, centr2 = 0)  #initial conditions

#Define initial parameters - ACh consensus values here
KAach <- 0.5
KEach1 <- 1.023
KEach2 <- 0.177
DVach <- 5.3
EC50ach <- 10^-5.383

#a function that takes parameters and runs the model
run_1drug_mod1 <- function(final_mod_params, stim_freq){
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

  return(finalres)
}

#returns a dataframe of chosen/tested parameters from the model
#conDF should be formatted using loadNormalizedDF (see below)
final_1Dparams <-function(stim_freq, m = 50, conDF, init_params, initModel, chosen=TRUE){

  lambda <- 2

  #creates workingdata for the model comparison
  workingdata <- data.frame(initModel[,"time"], initModel[,"eff2"], conDF[1:length(initModel[,"time"]),paste0("X", stim_freq, "HZ")])
  names(workingdata)<- c("Time", "Model", "Data")
  workingdata$SS <- (workingdata$Model - workingdata$Data)^2

  newparams <- init_params
  chosen_params <- c(init_params,0)
  tested_params <- init_params

  for (j in 1:m){
    for(i in 1:length(init_params)){
      testparams <- newparams
      testparams[i] <- abs(rnorm(n = 1, mean = newparams[i], sd = (sqrt(init_params[i]^2)*lambda)))
      tested_params <- rbind(tested_params, testparams)

      #calls run_2drugs_mod1 to run the model
      iteration <- run_1drug_mod1(testparams, stim_freq)

      #calculate Sum of Squares Objective function for the new model
      workingdata$newmodel <- iteration[1:length(workingdata$Model),"eff2"]
      workingdata$newSS <- with(workingdata, (newmodel - Data)^2)

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

freq0<-10

#Unknown paramters for fitting experimental data (KA, KE, D_V)
init_params1 <- c(KAach, 0.1)

WconDF <- loadNormalizedDF(1, lower = TRUE, dataDF = "cap", normDF = "cap")
initialresults <- run_1drug_mod1(init_params1, freq0)
finalparams <- final_1Dparams(freq0, m=500, WconDF, init_params1, initialresults)
finalresults <- run_1drug_mod1(finalparams[nrow(finalparams),],freq0)

#plot
p10<- (ggplot()
      + geom_point(aes(x=WconDF$Time, y=WconDF[,paste0("X", freq0, "HZ")]), color="blue")
      + geom_point(aes(x=initialresults[,"time"], y=initialresults[,"eff2"]), color="green")
      + geom_point(aes(x=finalresults[,"time"], y=finalresults[,"eff2"]), color="red")
)
p10

freq_list <- c(0.1, 0.3,0.7,1, 3, 7, 10, 15, 30)

#iterates through given freq_list and runs final_2Dparams 100 times for each frequency
for (q in freq_list){

  WconDF <- loadNormalizedDF(7, control = TRUE)
  initialresults <- run_2drugs_mod1(init_params1, q)
  finalparams <- final_2Dparams(q, m=500, WconDF, init_params1, initialresults)

  #take initial parameters and start repository vectors for final values
  finalparamsDF <- c(q, finalparams[nrow(finalparams),])

  for(k in seq(1:100)){
    finalparams <- final_2Dparams(q, m=500, WconDF, init_params1, initialresults)
    finalparamsDF <- rbind(finalparamsDF, c(q, finalparams[nrow(finalparams),]))
  }

  finalparamsDF<-as.data.frame(finalparamsDF)
  names(finalparamsDF) <- c("Freq", "KA","KE","Freq0", "D_V", "meanSS")
  write.table(finalparamsDF, paste0(con_directory_name, "Results_0.1-30_7.csv"), append = TRUE, sep = ",", dec = ".", qmethod = "double", col.names = FALSE)
}
