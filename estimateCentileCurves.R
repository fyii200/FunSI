# Author :  Fabian Yii                                                      #
# Email  :  fabian.yii@ed.ac.uk                                             #
# Info   :  R script for constructing 10 fundus centile charts, each based  #
#           on a different dimensionless (unitless) fundus imaging feature  #

# install.packages('quantregGrowth')
library(quantregGrowth)
library(vctrs)
library(dplyr)
library(tidyr)
rm(list=rm(list=ls()))
# setwd('..')




#### Read relevant datasets ####

# Participant dataset containing baseline data from all 89216 eyes of 51086 IDs 
# (out of 68508) with refractive error data and passed fundus image quality assessment
d                <- read.csv(file.path('data', 'participantData.csv'))
nrow(d); length(unique(d$id))

# Hospital procedural data (long format; each episode has its own row), coded based on OPCS-4
operations       <- read.csv(file.path('data', 'proceduralData.csv'))




#### Preprocess data ####

# Make sure baseline assessment date is coded as datetime
d$visitDate      <- as.Date(d$visitDate, '%d/%m/%Y')

# Rename columns in the procedural dataframe and make sure operation date is coded as datetime
operations       <- operations %>% dplyr::rename(id            = Participant.ID.participant...eid.,  
                                                 OPCS4         = Operative.procedures...OPCS4,
                                                 operationDate = Date.of.operation) %>% mutate(operationDate = as.Date(operationDate, '%d/%m/%Y'))

# Create new columns indicating the date an OPCS-4 code related to cataract surgery appears 
# for each ID (NA if no such history), one per timepoint (e.g. cataractSurgery1 = first ever code)
# Step 1: from the procedural data, subset IDs with ≥1 cataract surgery entries C74 
catSurgery       <- subset(operations, substring(OPCS4, 1, 3) == 'C71' |  # extracapsular extraction of lens
                                       substring(OPCS4, 1, 3) == 'C72' |  # intracapsular extraction of lens
                                       substring(OPCS4, 1, 3) == 'C74')   # other extraction of lens
# Step 2: convert from long to wide format, and then merge with the main dataframe
catSurgery       <- catSurgery %>% group_by(id) %>% arrange(operationDate) %>% mutate(TP = row_number())
catSurgery       <- catSurgery %>% ungroup(id) %>% pivot_wider(id_cols = id, names_from = TP, names_prefix = 'catSurgery', values_from = operationDate)
d                <- merge(d, catSurgery, by = 'id', all.x = T)

# Calculate disc-fovea distance (DFD) to disc major axis length (DML) ratio (DFD_DML_ratio)
d$DFD_DML_ratio  <- d$dist / d$major_length

# Calculate disc tilt 
d$discTilt       <- d$major_length / d$minor_length

# Calculate arteriovenous ratio (AVR)
d$AVR            <- d$CRAE_Knudtson / d$CRVE_Knudtson

# Calculate absolute disc torsion (aka orientation)
d$absDiscTors    <- abs(d$orientation)

# Factorise id, eye and sex
d$id             <- factor(d$id)
d$eye            <- factor(d$eye)
d$sex            <- factor(d$sex)




#### Apply selection criteria ####

# Remove eyes with SER > 0 (43528 eyes of 27039 IDs remained)
inputD           <- subset(d, SER <= 0)
nrow(inputD); length(unique(inputD$id))

# Remove eyes with NA value for any of the features considered (failed computation)
# 43147 eyes of 26845 IDs remained
inputD           <- inputD[!is.na(inputD$FD_artery) & !is.na(inputD$FD_vein) & !is.na(inputD$AVR) & !is.na(inputD$Tortuosity_density_artery) & !is.na(inputD$Tortuosity_density_vein) & !is.na(inputD$conc_rp_artery) & !is.na(inputD$conc_rp_vein) & !is.na(inputD$DFD_DML_ratio) & !is.na(inputD$discTilt) & !is.na(inputD$absDiscTors), ]
nrow(inputD); length(unique(inputD$id))

# Remove eyes with ≥1 feature in the top and bottom 0.1% of their respective distributions
# 42339 eyes of 26556 IDs remained
arteryFDpassed   <- inputD$FD_artery > quantile(inputD$FD_artery, 0.001) & inputD$FD_artery < quantile(inputD$FD_artery, 0.999) 
veinFDpassed     <- inputD$FD_vein > quantile(inputD$FD_vein, 0.001) & inputD$FD_vein < quantile(inputD$FD_vein, 0.999) 
AVRpassed        <- inputD$AVR > quantile(inputD$AVR, 0.001) & inputD$AVR < quantile(inputD$AVR, 0.999)
arteryTortPassed <- inputD$Tortuosity_density_artery > quantile(inputD$Tortuosity_density_artery, 0.001) & inputD$Tortuosity_density_artery < quantile(inputD$Tortuosity_density_artery, 0.999) 
veinTortPassed   <- inputD$Tortuosity_density_vein > quantile(inputD$Tortuosity_density_vein, 0.001) & inputD$Tortuosity_density_vein < quantile(inputD$Tortuosity_density_vein, 0.999) 
arteryConcPassed <- inputD$conc_rp_artery > quantile(inputD$conc_rp_artery, 0.001) & inputD$conc_rp_artery < quantile(inputD$conc_rp_artery, 0.999)
veinConcPassed   <- inputD$conc_rp_vein > quantile(inputD$conc_rp_vein, 0.001) & inputD$conc_rp_vein < quantile(inputD$conc_rp_vein, 0.999) 
DFD_DMLpassed    <- inputD$DFD_DML_ratio > quantile(inputD$DFD_DML_ratio, 0.001) & inputD$DFD_DML_ratio < quantile(inputD$DFD_DML_ratio, 0.999) 
discTiltPassed   <- inputD$discTilt > quantile(inputD$discTilt, 0.001) & inputD$discTilt < quantile(inputD$discTilt, 0.999)
discTorsPassed   <- inputD$absDiscTors > quantile(inputD$absDiscTors, 0.001) & inputD$absDiscTors < quantile(inputD$absDiscTors, 0.999)
includeBoolean   <- arteryFDpassed & veinFDpassed & AVRpassed & arteryTortPassed & veinTortPassed & arteryConcPassed & veinConcPassed & DFD_DMLpassed & discTiltPassed & discTorsPassed
inputD           <- inputD[includeBoolean, ]
nrow(inputD); length(unique(inputD$id))

# Remove eyes with myopia > 12D (imprecise estimates beyond this level due to sparse data, ~0.5% of eyes)
# 42130 eyes of 26456 IDs remained 
quantile(inputD$SER, 0.005)
inputD           <- subset(inputD, SER >= -12)
nrow(inputD); length(unique(inputD$id))

# Remove participants with a history of lens extraction surgery at baseline
# 41513 eyes of 26014 IDs remained
inputD           <- subset(inputD, catSurgery1 > visitDate | is.na(catSurgery1))
nrow(inputD); length(unique(inputD$id))

# Remove eyes with self-reported cataract surgery 
# 40258 eyes of 25222 IDs remained remained (1255 eyes removed)
inputD           <- subset(inputD, cataractSurgery_selfReported == F & refractiveLaser_selfReported == F)
nrow(inputD); length(unique(inputD$id))

# Save dataframe as 'FunSI'
write.csv(inputD, file.path('data', 'FunSI.csv'), row.names = F)




#### Estimate SER-specific centiles for each of the 10 imaging features using quantile regression ####

# Extract relevant columns
inputD           <- inputD %>% select(name, age, sex, catSurgery1, SER, FD_artery, FD_vein, AVR, Tortuosity_density_artery, Tortuosity_density_vein, conc_rp_artery, conc_rp_vein, DFD_DML_ratio, discTilt, absDiscTors)
 
# Arterial fractal dimension (save as 'm1')
m1               <- gcrq(FD_artery ~ ps(SER), data=inputD, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m1, file.path('fittedCentileCurves', 'm1.rds'))

# Venous fractal dimension (save as 'm2')
m2               <- gcrq(FD_vein ~ ps(SER), data=inputD, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m2, file.path('fittedCentileCurves', 'm2.rds'))

# AVR (save as 'm3')
m3               <- gcrq(AVR ~ ps(SER), data=inputD, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m3, file.path('fittedCentileCurves', 'm3.rds'))

# Arterial tortuosity (save as 'm4')
m4               <- gcrq(Tortuosity_density_artery ~ ps(SER), data=inputD, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m4, file.path('fittedCentileCurves', 'm4.rds'))

# Venous tortuosity (save as 'm5')
m5               <- gcrq(Tortuosity_density_vein ~ ps(SER), data=inputD, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m5, file.path('fittedCentileCurves', 'm5.rds'))

# Arterial concavity (save as 'm6')
m6               <- gcrq(conc_rp_artery ~ ps(SER), data=inputD, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m6, file.path('fittedCentileCurves', 'm6.rds'))

# Venous concavity (save as 'm7')
m7               <- gcrq(conc_rp_vein ~ ps(SER), data=inputD, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m7, file.path('fittedCentileCurves', 'm7.rds'))

# DFD_DML_ratio (save as 'm8')
m8               <- gcrq(DFD_DML_ratio ~ ps(SER), data=inputD, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m8, file.path('fittedCentileCurves', 'm8.rds'))

# Disc tilt (save as 'm9')
m9               <- gcrq(tilt ~ ps(SER), data=inputD, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m9, file.path('fittedCentileCurves', 'm9.rds'))

# Absolute disc torsion (save as 'm10')
m10              <- gcrq(absOrient ~ ps(SER), data=inputD, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m10, file.path('fittedCentileCurves', 'm10.rds'))















