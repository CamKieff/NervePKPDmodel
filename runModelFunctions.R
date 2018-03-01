#functions to run RxODE and find best fit parameters for model

#incorporating 0.1 Hz variability without creating a new column in the results of the plot.
#In fact, if the model only has one NT, there are a lot of extraneous variables in the output
#need some code to trim output from final_drug_params

#some global parameters need to be defined outside of the functions (init_params)
#variables cannot be controlled from outside some of the function wrappers (m, freq_list?)
#initial parameters, model definition

require(ggplot2)
require(reshape2)
source("normalizedDF.R")
source("defineModel.R")

#a function that takes parameters and runs the model
run_mod1 <- function(stim_freq, mod_params){
  parameters <- c(KA1 = mod_params$KAunk,
                  KE1 = mod_params$KEunk,
                  KA2 = mod_params$KAach,
                  KE2 = mod_params$KEach,
                  KA2f = mod_params$m2max,
                  KE2f = mod_params$chemax,
                  IC501 = 10^-mod_params$IC50m2,
                  IC502 = 10^-mod_params$IC50che,
                  max1 = mod_params$MAXunk,
                  EC501 = 10^-mod_params$EC50unk,
                  EC502 = 10^-mod_params$EC50ach
  )
  time_var <- 60 #length of stimulation in seconds

  if (thismodel[[3]] == 2){
    pulse_rate <- 1/stim_freq
    num_doses <- time_var/pulse_rate
    testev <- eventTable(amount.units="mol", time.unit="seconds")
    testev$add.dosing(dose = 10^-mod_params$DVunk, nbr.doses = num_doses, dosing.interval = pulse_rate, dosing.to = 1, start.time = 0)
    testev$add.dosing(dose = 10^-mod_params$DVach, nbr.doses = num_doses, dosing.interval = pulse_rate, dosing.to = 3, start.time = 0)
    testev$add.sampling(seq(from = 0, to = time_var, by = 0.02))
  } else if(thismodel[[3]] == 1){

    if(stim_freq == 0.1){
      pulse_rate <- 1/stim_freq
    } else{
      pulse_rate <- 1/stim_freq
    }
    num_doses <- time_var/pulse_rate
    testev <- eventTable(amount.units="mol", time.unit="seconds")
    testev$add.dosing(dose = 10^-mod_params$DVach, nbr.doses = num_doses, dosing.interval = pulse_rate, dosing.to = 1, start.time = 0)
    testev$add.sampling(seq(from = 0, to = time_var, by = 0.02))
  } else(
    print("Model does not currently support more than two neurotransmitters (NT). Please set NT equal to either 1 or 2.")
  )
  finalres <- thismodel[[1]]$solve(parameters, testev, thismodel[[2]]) #this actually solves the model
  return(finalres)
}

#returns a dataframe of chosen/tested best-fit parameters from the model
#conDF should be formatted using loadNormalizedDF from "normalizedDF.R"
final_drug_params <-function(stim_freq, m = 50, conDF, bestfit, init_params, init_model, chosen=TRUE){

  testexp <- 2     #sum of "squares" exponent (must be even; 2 or 4 are probably optimal)
  lambda <- 2      #learning rate

  #creates workingdata for the model comparison
  workingdata <- data.frame(init_model[,"time"], init_model[,"eff2"], conDF[1:length(init_model[,"time"]),paste0("X", stim_freq, "HZ")])
  names(workingdata)<- c("Time", "Model", "Data")
  workingdata$SS <- (workingdata$Model - workingdata$Data)^testexp

  newparams <- init_params           #stores current best-fit parameters
  chosen_params <- c(init_params,0)  #stores a list of accepted best-fit parameters
  tested_params <- init_params       #stores a list of all attempted best-fit parameters (primarily for diagnostic purposes)

  for (j in 1:m){      #number of iterations to find best fit
    for(i in bestfit){
      testparams <- newparams
      testparams[[i]] <- abs(rnorm(n = 1, mean = newparams[[i]], sd = (sqrt(init_params[[i]]^2)*lambda)))
      tested_params <- rbind(tested_params, testparams)

      iteration <- run_mod1(stim_freq, testparams)       #calls run_mod1 to run the model

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
    return(chosen_params) #all accepted parameters. Final row are the best-fit parameters
  } else {
    return(tested_params) #all attempted parameters. Final row are the best-fit parameters
  }
}

#runs the model 100x for all tissues and frequencies
#can also run a single model if con_list is only set to a single index
Iteration <- function(con_list = c(1,2,5,7), m = 500, n = 100, dataDF = "con", normDF = "cap", lower = TRUE, filename = "/FormattedLowerTrachea/capsaicin/Results_"){
  freq_list <- c(0.1, 0.3,0.7,1, 3, 7, 10, 15, 30)
  for(r in con_list){                                                 #iterates through each raw file
    WconDF <- loadNormalizedDF(r, lower = lower, dataDF = dataDF, normDF = normDF) #load file
    for (q in freq_list){                                             #iterates through each frequency
      initialresults <- run_mod1(q, init_params1)              #find starting point model results from initial parameters
      finalparams <- final_drug_params(q, m = m, WconDF, bestfit = bestfit, init_params, initialresults)

      finalparamsDF <- c(q, finalparams[nrow(finalparams),])          #take initial parameters and start vector for final values

      for(k in seq(1:n)){
        finalparams <- final_drug_params(q, m = m, WconDF, bestfit = bestfit, init_params, initialresults)
        finalparamsDF <- rbind(finalparamsDF, c(q, finalparams[nrow(finalparams),]))
      }
      #append results to a single cvs file per index after each frequency
      finalparamsDF<-as.data.frame(finalparamsDF)
      write.table(finalparamsDF, paste0(filename, r, ".csv"), append = TRUE, sep = ",", dec = ".", qmethod = "double", col.names = FALSE)
    }
  }
}

#run final_drug_params once for each frequency and create a nice facetwrap graph of the best-fit results
facetgraph <- function(conDF, init_params, bestfit, consensus = FALSE){
  freq_list <- c(0.1, 0.3, 0.7, 1, 3, 7, 10 , 15, 30)
  facetDF <-NULL
  for (i in freq_list){
    working_freq <- i
    initialresults <- run_mod1(stim_freq = working_freq, init_params)
    if(consensus == TRUE){
      plotDF <- data.frame(rep(working_freq, length(initialresults)), initialresults[,"time"], conDF[1:length(initialresults[,"time"]),paste0("X", working_freq, "HZ")], initialresults[,"eff2"])
    } else{
      finalparams <- final_drug_params(stim_freq = working_freq, m = 500, conDF, bestfit = bestfit, init_params, initialresults)
      finalresults <- run_mod1(stim_freq = working_freq, finalparams[nrow(finalparams),])
      plotDF <- data.frame(rep(working_freq, length(finalresults)), initialresults[,"time"], conDF[1:length(initialresults[,"time"]),paste0("X", working_freq, "HZ")], initialresults[,"eff2"], finalresults[,"eff2"])
    }

    facetDF <- rbind(plotDF, facetDF)
  }
  if(consensus == TRUE){
    colnames(facetDF)<-c("Freq", "Time", "Raw", "Consensus")
    p <- (ggplot(facetDF)
          + geom_line(aes(x=Time, y=Raw), color="black", alpha = 0.5)
          + geom_line(aes(x=Time, y=Consensus), color="red",size = 1)
          + facet_wrap(~ Freq, scales="free", ncol=3)
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
