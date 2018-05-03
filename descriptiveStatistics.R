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

con1_1 <- read.csv("FormattedLowerTrachea/Control/Results0.1Hz_1.csv", header= TRUE)
names(con1_1)<-c("freq", names(con1_1)[1:15])
con1_1 <- select(con1_1, bestfit, "freq", "test_freq", "SS")

con1_1 %>%
  select(bestfit, "freq", "test_freq", "SS") %>%
  group_by(freq) %>%
  summarise(myFun(test_freq))

