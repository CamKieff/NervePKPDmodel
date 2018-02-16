#load and normalized Control data.frames. if control is true, max is set from capsaicin; if false, both max and min are capsaicin
#data to be normalized (dataDF) = control ("con") or capsaicin ("cap")
#normalizing data frame (normDF) = control ("con") or capsaicin ("cap")
#lower trachea = TRUE; upper trachea = FALSE
loadNormalizedDF <- function(index = 1, lower = TRUE, dataDF = "con", normDF = "cap"){

  #creates a list of file names for the chosen directory
  if (lower){
    cap_directory_name <- 'FormattedLowerTrachea/capsaicin/'
    con_directory_name <- 'FormattedLowerTrachea/control/'
  } else{
    cap_directory_name <- 'FormattedUpperTrachea/capsaicin/'
    con_directory_name <- 'FormattedUpperTrachea/control/'
  }
  cap_list <-dir(cap_directory_name, "*.csv")
  con_list <-dir(con_directory_name, "*.csv")

  #load DF
  capDF <- read.csv(paste0(cap_directory_name, cap_list[index])) #Working capsaicin data.frame
  conDF <- read.csv(paste0(con_directory_name, con_list[index]))

  if(normDF == "cap"){
    FRmax <- (max(capDF[,2:ncol(capDF)])-min(capDF[,2:ncol(capDF)])) #denominator for capsaicin and control normalization using capsaicin data
  }else if(normDF == "con"){
    FRmax <- (max(conDF[,2:ncol(conDF)])-min(conDF[,2:ncol(conDF)])) #denominator for capsaicin and control normalization using capsaicin data
  }else{
    print("The normalization DF (normDF) must be either control (\"con\") or (\"cap\")")
  }

  if(dataDF == "con"){
    #format control DF
    FRmin <- min(conDF[,2:ncol(conDF)]) #baseline mg tension to subtract to get relative contraction values
    conDF[,2:ncol(conDF)] <- (conDF[,2:ncol(conDF)] - FRmin)/FRmax
    return(conDF)
  } else if (dataDF == "cap"){
    #format capsaicin DF
    FRmin <- min(capDF[,2:ncol(capDF)]) #baseline mg tension to subtract to get relative contraction values
    capDF[,2:ncol(capDF)] <- (capDF[,2:ncol(capDF)] - FRmin)/FRmax
    return(capDF)
  } else{
    print("dataDF must be either control data (\"con\") or capsaicin data (\"cap\")")
  }
}
