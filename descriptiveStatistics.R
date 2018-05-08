#This script finds descriptive statistics after model runs.

require(dplyr)

#set working directory function
modelWD <- function(lowertrachea = TRUE, capsaicin = FALSE){
  if(lowertrachea){
    trachea_directory_name <- 'FormattedLowerTrachea/'
  }else{
    trachea_directory_name <- 'FormattedUpperTrachea/'
  }

  if(capsaicin){
    directory_name <- paste0(trachea_directory_name, 'control/')
  }else{
    directory_name <- paste0(trachea_directory_name, 'capsaicin/')
  }
  return(directory_name)
}

#calculate descriptive statistics for parameter files
findStatistics <- function(filelist, bestfit=bestfit, directory_name=directory_name){
  myFun <- function(x) {                    #define descriptive statistics function
    c(min = min(x), max = max(x),
      mean = mean(x), median = median(x),
      std = sd(x))
  }
  descrstats <- NULL #set up empty data frame

  for(y in filelist){
    #import CSV list of final parameters
    wdata <- read.csv(paste0(directory_name, y, ".csv"), header=FALSE)
    names(wdata) <- c('','Freq', bestfit, 'meanSS')

    #collects summary statistics for values given in the vector. Rbinds them to our empty data frame
    for (i in c(bestfit, "meanSS")){
      testdata <- do.call(rbind, by(wdata[, i], wdata[,"Freq"], myFun))
      testdata <- cbind(rep(y, nrow(testdata)), rep(i, nrow(testdata)), testdata)
      descrstats <- rbind(descrstats, testdata)
    }
  }
  return(descrstats)
}

#directory_name <- modelWD()
#testlist <- c("Results_1", "Results_2", "Results_5", "Results_7") #example list of file names
#bestfit <- c("KAach", "KEach", "DVach")  #what unknowns were solved for
#write.csv(descrstats, file = paste0(directory_name, "descrstats_2drugs_UT.csv")) #export

df<- data.frame(NULL)
bestfit <- c("KAach", "KEach", "DVach", "test_freq")

for(i in c(1,2,5,7)){
  con1 <- read.csv(paste0("FormattedLowerTrachea/control/Results0.1Hz_", i, ".csv"), header= FALSE)
  names(con1)<-c("","freq", names(init_params), "test_freq", "SS")
  
  df<-
  con1 %>%
    select(bestfit, "freq","SS") %>%
    mutate(tissue = i) %>%
    rbind(df)
}

aggdata <- mutate(aggregate.data.frame(df, by=list(df$freq, df$tissue), mean), var = "mean")
aggdata <- rbind(aggdata, mutate(aggregate.data.frame(df, by=list(df$freq, df$tissue), median), var = "median"))
aggdata <- rbind(aggdata, mutate(aggregate.data.frame(df, by=list(df$freq, df$tissue), sd), var = "sd"))

aggdata <- arrange(aggdata, var, Group.1, Group.2)

write.csv(aggdata, "FormattedLowerTrachea/control/aggregate_results_0.1Hzsimple.csv")