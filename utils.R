# Author :  Fabian Yii                                                      #
# Email  :  fabian.yii@ed.ac.uk                                             #
# Info   :  Helper functions                                                #
library(quantregGrowth)




## Set-up: read the fitted centile curves (CC) for each of the 10 fundus imaging features
fileNames  <- list.files('fittedCentileCurves', pattern = '.rds')
for(i in 1:length(fileNames)){
  varName  <- unlist(strsplit(fileNames[i], '.rds'))
  fittedCC <- readRDS(file.path('fittedCentileCurves', fileNames[i]))
  assign(varName, fittedCC)
  }




#### Function 1: compute Fundus Stretch Index (FunSI) given an eye's spherical ####
####             equivalent refraction (SER) and observed imaging features     ####
getFunSI <- function(SEReye,        # SER of the eye
                     arteryFD,      # arterial fractal dimension
                     veinFD,        # venous fractal dimension 
                     AVR,           # arteriovenous ratio (AVR)
                     arteryTort,    # arterial tortuosity
                     veinTort,      # venous tortuosity
                     arteryConc,    # arterial concavity
                     veinConc,      # venous concavity
                     DFD_DML_ratio, # disc-fovea distance to disc major axis length ratio
                     discTilt,      # disc tilt
                     absDiscTors,   # absolute disc torsion
                     leaveOneOut = T)   
{  
  
  ## Determine the centile position (CP) of the eye on each of the 10 fundus centile charts:
  ## 1. For a given imaging feature, get the SER-specific fitted centiles (5th–95th).
  ## 2. Calculate the absolute difference between each fitted centile and the observed value of the imaging feature.
  ## 3. The CP is the index of the fitted centile with the smallest difference.
  
  # Note: these imaging features generally decrease as myopia increases (5th centile = worst centile), 
  #       so subtract 1 from CP so a higher value is worse.
  # Arterial fractal dimension
  dif     <- abs(arteryFD - predict(m1, data.frame(SER = SEReye)))
  CP1     <- 1 - as.numeric(names(which(dif == min(dif))))  
  dif_SA  <- abs(arteryFD - predict(m1_sensitivity, data.frame(SER = SEReye)))
  CP1_SA  <- 1 - as.numeric(names(which(dif_SA == min(dif_SA))))  
  # Venous fractal dimension
  dif     <- abs(veinFD - predict(m2, data.frame(SER = SEReye)))
  CP2     <- 1 - as.numeric(names(which(dif == min(dif)))) 
  dif_SA  <- abs(veinFD - predict(m2_sensitivity, data.frame(SER = SEReye)))
  CP2_SA  <- 1 - as.numeric(names(which(dif_SA == min(dif_SA)))) 
  # AVR
  dif     <- abs(AVR - predict(m3, data.frame(SER = SEReye)))
  CP3     <- 1 - as.numeric(names(which(dif == min(dif)))) 
  dif_SA  <- abs(AVR - predict(m3_sensitivity, data.frame(SER = SEReye)))
  CP3_SA  <- 1 - as.numeric(names(which(dif_SA == min(dif_SA)))) 
  # Arterial tortuosity 
  dif     <- abs(arteryTort - predict(m4, data.frame(SER = SEReye)))
  CP4     <- 1 - as.numeric(names(which(dif == min(dif)))) 
  dif_SA  <- abs(arteryTort - predict(m4_sensitivity, data.frame(SER = SEReye)))
  CP4_SA  <- 1 - as.numeric(names(which(dif_SA == min(dif_SA)))) 
  # Venous tortuosity 
  dif     <- abs(veinTort - predict(m5, data.frame(SER = SEReye)))
  CP5     <- 1 - as.numeric(names(which(dif == min(dif)))) 
  dif_SA  <- abs(veinTort - predict(m5_sensitivity, data.frame(SER = SEReye)))
  CP5_SA  <- 1 - as.numeric(names(which(dif_SA == min(dif_SA)))) 
  # Absolute disc torsion
  dif     <- abs(absDiscTors - predict(m10, data.frame(SER = SEReye)))
  CP10    <- 1 - as.numeric(names(which(dif == min(dif)))) 
  dif_SA  <- abs(absDiscTors - predict(m10_sensitivity, data.frame(SER = SEReye)))
  CP10_SA <- 1 - as.numeric(names(which(dif_SA == min(dif_SA)))) 
  
  # Note: these imaging features generally increase as myopia increases (95th centile = worst centile)
  # Arterial concavity
  dif     <- abs(arteryConc - predict(m6, data.frame(SER = SEReye)))
  CP6     <- as.numeric(names(which(dif == min(dif))))
  dif_SA  <- abs(arteryConc - predict(m6_sensitivity, data.frame(SER = SEReye)))
  CP6_SA  <- as.numeric(names(which(dif_SA == min(dif_SA))))
  # Venous concavity
  dif     <- abs(veinConc - predict(m7, data.frame(SER = SEReye)))
  CP7     <- as.numeric(names(which(dif == min(dif))))
  dif_SA  <- abs(veinConc - predict(m7_sensitivity, data.frame(SER = SEReye)))
  CP7_SA  <- as.numeric(names(which(dif_SA == min(dif_SA))))
  # Disc-fovea distance to disc major length ratio
  dif     <- abs(DFD_DML_ratio - predict(m8, data.frame(SER = SEReye)))
  CP8     <- 1 - as.numeric(names(which(dif == min(dif)))) 
  dif_SA  <- abs(DFD_DML_ratio - predict(m8_sensitivity, data.frame(SER = SEReye)))
  CP8_SA  <- 1 - as.numeric(names(which(dif_SA == min(dif_SA))))
  # Disc tilt
  dif     <- abs(discTilt - predict(m9, data.frame(SER = SEReye)))
  CP9     <- as.numeric(names(which(dif == min(dif))))
  dif_SA  <- abs(discTilt - predict(m9_sensitivity, data.frame(SER = SEReye)))
  CP9_SA  <- as.numeric(names(which(dif_SA == min(dif_SA))))
  
  ## Compute FunSI by summing the 10 CPs and min-max normalise the summed CP to 0-1 
  CPsum     <- CP1 + CP2 + CP3 + CP4 + CP5 + CP6 + CP7 + CP8 + CP9 + CP10
  FunSI     <- (CPsum - 0.5) / (9.5 - 0.5) 
  CPsum_SA  <- CP1_SA + CP2_SA + CP3_SA + CP4_SA + CP5_SA + CP6_SA + CP7_SA + CP8_SA + CP9_SA + CP10_SA
  FunSI_SA  <- (CPsum_SA - 0.5) / (9.5 - 0.5) 
  
  ## Leave one feature out and compute FunSI by summing the remaining 9 CPs and min-max normalise the summed CP to 0-1 
  CPsum     <- CP2 + CP3 + CP4 + CP5 + CP6 + CP7 + CP8 + CP9 + CP10
  FunSIwo1  <- (CPsum - 0.45) / (8.55 - 0.45) 
  CPsum     <- CP1 + CP3 + CP4 + CP5 + CP6 + CP7 + CP8 + CP9 + CP10
  FunSIwo2  <- (CPsum - 0.45) / (8.55 - 0.45) 
  CPsum     <- CP1 + CP2 + CP4 + CP5 + CP6 + CP7 + CP8 + CP9 + CP10
  FunSIwo3  <- (CPsum - 0.45) / (8.55 - 0.45) 
  CPsum     <- CP1 + CP2 + CP3 + CP5 + CP6 + CP7 + CP8 + CP9 + CP10
  FunSIwo4  <- (CPsum - 0.45) / (8.55 - 0.45) 
  CPsum     <- CP1 + CP2 + CP3 + CP4 + CP6 + CP7 + CP8 + CP9 + CP10
  FunSIwo5  <- (CPsum - 0.45) / (8.55 - 0.45) 
  CPsum     <- CP1 + CP2 + CP3 + CP4 + CP5 + CP7 + CP8 + CP9 + CP10
  FunSIwo6  <- (CPsum - 0.45) / (8.55 - 0.45) 
  CPsum     <- CP1 + CP2 + CP3 + CP4 + CP5 + CP6 + CP8 + CP9 + CP10
  FunSIwo7  <- (CPsum - 0.45) / (8.55 - 0.45) 
  CPsum     <- CP1 + CP2 + CP3 + CP4 + CP5 + CP6 + CP7 + CP9 + CP10
  FunSIwo8  <- (CPsum - 0.45) / (8.55 - 0.45) 
  CPsum     <- CP1 + CP2 + CP3 + CP4 + CP5 + CP6 + CP7 + CP8 + CP10
  FunSIwo9  <- (CPsum - 0.45) / (8.55 - 0.45) 
  CPsum     <- CP1 + CP2 + CP3 + CP4 + CP5 + CP6 + CP7 + CP8 + CP9
  FunSIwo10 <- (CPsum - 0.45) / (8.55 - 0.45) 
  
  if(leaveOneOut){
    return(list(FunSI, FunSI_SA, FunSIwo1, FunSIwo2, FunSIwo3, FunSIwo4, FunSIwo5, FunSIwo6, FunSIwo7, FunSIwo8, FunSIwo9, FunSIwo10)) 
    } else{
      return(list(FunSI, FunSI_SA))
    }
  }




#### Function 2: display centile curves given a fitted a fitted quantile regression model (m) ####
plotCC <- function(m, name){
  plot.gcrq(m, res = F, col = 'gray', conf.level = 0, shade = T, main = name, xlab = '', ylab = '', axes = F)
  plot.gcrq(m, res = F, col = 'red', add = T, select.tau = 0.5, conf.level = 0, shade = T, xlab = '', ylab = '', axes = F)
  axis(1, lwd = 0, lwd.ticks = 1); axis(2, lwd = 0, lwd.ticks = 1)
  }




#### Function 3: return a Boolean indicating if a value lies within the top and bottom q% of a given imaging feature ####
inlier <- function(feature, q = 0.1){
  return(feature > quantile(feature, q/100) & feature < quantile(feature, 1-(q/100)))
  }


