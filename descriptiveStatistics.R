#This script finds descriptive statistics after model runs.
require(dplyr)

df<- data.frame(NULL)
bestfit <- c("KAach", "KEach", "DVach", "test_freq")
#bestfit <- c("m2max", "chemax", "IC50m2","IC50che")

for(i in c(1,2,5,7)){
  con1 <- read.csv(paste0("FormattedLowerTrachea/capsaicin/Results0.1Hz_", i, ".csv"), header= FALSE)
  names(con1)<-c("blank","freq", names(init_params), "test_freq", "SS")
  
  df<-
    con1 %>%
    mutate(tissue = i) %>%
    rbind(df)
}
df$blank <- NULL
aggdata <- mutate(aggregate.data.frame(df, by=list(df$freq, df$tissue), mean, na.action=TRUE), var = "mean")
aggdata <- rbind(aggdata, mutate(aggregate.data.frame(df, by=list(df$freq, df$tissue), median, na.action=TRUE), var = "median"))
aggdata <- rbind(aggdata, mutate(aggregate.data.frame(df, by=list(df$freq, df$tissue), sd), var = "sd"))

aggdata <- arrange(aggdata, var, Group.1, Group.2)

aggdata$freq <- NULL
aggdata$tissue <- NULL
names(aggdata) <- c("Freq", "Tissue", names(init_params),"test_freq", "SS", "var")

write.csv(aggdata, "FormattedLowerTrachea/capsaicin/cap_aggregate_0.1Hz.csv")

aggdata <- read.csv("FormattedLowerTrachea/capsaicin/cap_aggregate_complex.csv", header = TRUE)
#bestfit <- c("KAach", "KEach", "DVach")
bestfit <- c("m2max", "chemax", "IC50m2","IC50che")

agg_tissues <-
aggdata %>%
  select("Group.1", bestfit, "SS", "tissue","var") %>%
  filter(var =="median", Group.1 %in% c(0.3, 1,3,10)) %>%
  plyr::rename(., c("Group.1" = "Frequency")) %>%
  aggregate(by=list(.$tissue), mean) %>%
  mutate(var = "mean") %>%
  select(-tissue, -Frequency)
  
agg_means <- 
  agg_tissues %>%
  select(bestfit) %>%
  apply(2, FUN=mean)

agg_sem <- 
  agg_tissues %>%
  select(bestfit) %>%
  apply(., 2, function(x) sd(x)/sqrt(length(x)))
  
agg_results <- rbind(agg_tissues %>% select(bestfit), agg_means, agg_sem)
agg_results$name <- c(1,2,5,7, "Mean", "SEM")

parameterlist<- read.csv("FormattedLowerTrachea/evenMoreFinalParams.csv")
names(parameterlist) <- c("tissue", "control", "complex", names(init_params))

df <- read.csv("FormattedLowerTrachea/control/ResultsHz_5.csv")
names(df) <- c("blank", "freq1",names(init_params), "freq0","SS")
p0 <-
  df %>%
  #filter(freq1 %in% c(0.3, 1,3,10)) %>%
  mutate(var = seq(1:nrow(.))) %>%
  ggplot(.)
         
p0 + geom_point(aes(x=KAach, y=KEach, color=factor(freq1))) #+scale_x_continuous(name ="Frequency (Hz)",breaks = c(0,100,201,302,403),labels=c("0.3","1","3","10",""))
p0 + geom_point(aes(x=chemax, y=IC50che, color=factor(freq1))) #+scale_x_continuous(name ="Frequency (Hz)",breaks = c(0,100,201,302,403),labels=c("0.3","1","3","10",""))

p0 + geom_point(aes(x=var, y=KEach))
