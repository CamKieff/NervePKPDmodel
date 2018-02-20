#This script will help you find descriptive statistics after model runs.
#At the end is some code to plot diagnostic plots of the mean vs median

lowertrachea <- TRUE
if(lowertrachea){
  cap_directory_name <- 'FormattedLowerTrachea/capsaicin/'
  con_directory_name <- 'FormattedLowerTrachea/control/'
}else{
  cap_directory_name <- 'FormattedUpperTrachea/capsaicin/'
  con_directory_name <- 'FormattedUpperTrachea/control/'
}

#function to find descriptive statstics
myFun <- function(x) {
  c(min = min(x), max = max(x),
  mean = mean(x), median = median(x),
  std = sd(x))
}

descrstats <- NULL #set up empty data frame

filename <- "Results_5.csv" #example

#import final parameters
wdata <- read.csv(paste0(con_directory_name, filename), header=FALSE)
names(wdata) <- c('','Freq','KAunk', 'KEunk', 'DVunk', 'EC50unk', 'MAXunk', 'meanSS')
#names(wdata) <- c('','Freq','KAach', 'KEach', 'DVach', 'meanSS')

#log10 transfrom applicable values
#wdata$logIC501 <- log10(wdata$IC501)
#wdata$logIC502 <- log10(wdata$IC502)
wdata$logDVunk <- log10(wdata$DVunk)
wdata$logEC50unk <- log10(wdata$EC50unk)
#wdata$logDVach <- log10(wdata$DVach)

#collects summary statistics for values given in the vector. Rbinds them to our empty data frame
for (i in c("KAunk", "KEunk", "logDVunk", "logEC50unk", "MAXunk", "meanSS")){
  testdata <- do.call(rbind, by(wdata[, i], wdata[,"Freq"], myFun))
  testdata <- cbind(rep(i, nrow(testdata)), testdata)
  descrstats <- rbind(descrstats, testdata)
}

#repeat the above code for each data set manually before exporting
write.csv(descrstats, file = paste0(con_directory_name, "descrstats_2drugs_UT.csv")) #export

#
#
#
#
#

#Below is the graphing functions and formats. First reimport the descrstats as wdata
#this code is not very polished. Many things may have to be altered manually
wdata <- read.csv(paste0(directory_name, "Results_0.3-30_doubleInhibit_7.csv"), header=TRUE)

#assign number for each individual "each" should be the number of frequencies multiplied by the
#number of parameters you have summary stats for
wdata <-cbind(wdata, rep(1:4, each=32))

colnames(wdata)<- c('Freq', 'Var', 'min', 'max', 'mean', 'median', 'std', 'sample')

#subset off whichever factor you want to graph
KAv <- subset(wdata, Var == 'logIC502')

KAplot <- (ggplot()
+ geom_line(data = KAv, aes(log10(Freq), mean, group =sample), color = 'blue', alpha = 0.75)
+ geom_line(data = KAv, aes(log10(Freq), median, group =sample), color = 'red', alpha = 0.75)
+ geom_line(aes(x = log10(c(0.3, 0.7, 1, 3, 7, 10, 15, 30)), y = as.vector(tapply(KAv$mean, KAv$Freq, mean))), colour="blue", size = 2)
+ geom_line(aes(x = log10(c(0.3, 0.7, 1, 3, 7, 10, 15, 30)), y = as.vector(tapply(KAv$median, KAv$Freq, mean))), colour="red", size = 2)
+labs(list(title="logIC502", x="log(Freq)", y="logIC502"))
+theme(title=element_text(size=18, face="bold"),axis.title=element_text(size=14, face="bold"),axis.text=element_text(size=14))
)
KAplot
