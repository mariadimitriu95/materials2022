extractMeanfromMelt <- function(ml, centerRange = 500){
  
  numSamples <- length(table(ml$sample))
  samples <- unique(ml$sample)
  numSplits <- length(table(ml$split))
  splits <-unique(ml$split)
  
  ml.cen <- ml[abs(ml$position) <= centerRange,]
  
  df <- as.data.frame(matrix(NA, nrow = numSamples*numSplits, ncol = 3))
  colnames(df) <- c("Sample", "Split", "Mean")
  
  df$Sample <- sort(rep(samples,numSplits))
  df$Split <- rep(splits,numSamples)
  
  for (i in 1:nrow(df)){
    
    df[i, "Mean"] <- mean(ml[ml$split == df[i, "Split"]  & ml$sample == df[i, "Sample"] ,"mean"])
    
    
  }
  return(df)
}