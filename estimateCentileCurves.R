# Author :  Fabian Yii                                                      #
# Email  :  fabian.yii@ed.ac.uk                                             #
# Info   :  R script for constructing 10 fundus centile charts, each based  #
#           on a different dimensionless (unitless) fundus imaging feature  #

library(quantregGrowth)
library(vctrs)
library(dplyr)
library(tidyr)
rm(list=rm(list=ls()))
# setwd('..')
source(file.path('code', 'utils.R')) # get helper functions




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

# Remove eyes with SER > 0 (43528 eyes of 27039 IDs eligible)
d                <- subset(d, SER <= 0)
initialD         <- d       
nrow(d); length(unique(d$id))

# Remove eyes with NA value for any of the features considered (failed computation)
# 43147 eyes of 26845 IDs eligible 
d                <- d[!is.na(d$FD_artery) & !is.na(d$FD_vein) & !is.na(d$AVR) & !is.na(d$Tortuosity_density_artery) & !is.na(d$Tortuosity_density_vein) & !is.na(d$conc_rp_artery) & !is.na(d$conc_rp_vein) & !is.na(d$DFD_DML_ratio) & !is.na(d$discTilt) & !is.na(d$absDiscTors), ]
nrow(d); length(unique(d$id))
removed          <- setdiff(initialD[c('id', 'SER')], d[c('id', 'SER')])
nrow(removed)
table(removed$SER <= -6 & removed$SER >= -12)                       # 77 non-eligible eyes between -6 and -12D
table(removed$SER < -12)                                            # 29 non-eligible eyes < -12D

# Remove eyes with ≥1 feature in the top and bottom 0.1% of their respective distributions (create a Boolean variable indicating outliers)
# 42339 eyes of 26556 IDs eligible 
includeBoolean   <- inlier(d$FD_artery) & inlier(d$FD_vein) & inlier(d$AVR) & inlier(d$Tortuosity_density_artery) & inlier(d$Tortuosity_density_vein) & inlier(d$conc_rp_artery) & inlier(d$conc_rp_vein) & inlier(d$DFD_DML_ratio) & inlier(d$discTilt) & inlier(d$absDiscTors) 
d$outlier        <- ifelse(includeBoolean, F, T)
nrow(subset(d, !outlier)); length(unique(subset(d, !outlier)$id))
table(subset(d, outlier)$SER <= -6 & subset(d, outlier)$SER >= -12) # 125 non-eligible eyes between -6 and -12D
table(subset(d, outlier)$SER < -12)                                 # 28 non-eligible eyes < -12D

################################################################################################
# For SENSITIVITY ANALYSIS 5 in 'survivalAnalyses.R': set threshold at 0.01% rather than 0.1% 
# Create a Boolean variable indicating outliers based on this new threshold
includeBooleanSA <- inlier(d$FD_artery, 0.01) & inlier(d$FD_vein, 0.01) & inlier(d$AVR, 0.01) & inlier(d$Tortuosity_density_artery, 0.01) & inlier(d$Tortuosity_density_vein, 0.01) & inlier(d$conc_rp_artery, 0.01) & inlier(d$conc_rp_vein, 0.01) & inlier(d$DFD_DML_ratio, 0.01) & inlier(d$discTilt, 0.01) & inlier(d$absDiscTors, 0.01) 
d$outlierSA      <- ifelse(includeBooleanSA, F, T)
################################################################################################

# Remove eyes with myopia > 12D (imprecise estimates beyond this level due to sparse data, ~0.5% of eyes)
# 42130 eyes of 26456 IDs eligible 
quantile(subset(d, !outlier)$SER, 0.005)
d                <- subset(d, SER >= -12)
nrow(subset(d, !outlier)); length(unique(subset(d, !outlier)$id))

# Remove participants with a history of lens extraction surgery at baseline
# 41513 eyes of 26014 IDs remained
d                <- subset(d, catSurgery1 > visitDate | is.na(catSurgery1))
nrow(subset(d, !outlier)); length(unique(subset(d, !outlier)$id))

# Remove eyes with self-reported cataract surgery 
# 40258 eyes of 25222 IDs remained remained (1255 eyes removed)
d                <- subset(d, cataractSurgery_selfReported == F & refractiveLaser_selfReported == F)
nrow(subset(d, !outlier)); length(unique(subset(d, !outlier)$id))

# Save dataframe as 'FunSI'
write.csv(d, file.path('data', 'FunSI.csv'), row.names = F)




#### Estimate SER-specific centiles for each of the 10 imaging features using quantile regression ####

# Randomly select one eye if both eyes of an individual were available (25,222 eyes)
set.seed(1)
inputD           <- subset(d, !outlier) %>% group_by(id) %>% slice_sample(n = 1)
nrow(inputD)

# Extract relevant columns
inputD           <- inputD %>% select(name, age, sex, catSurgery1, SER, FD_artery, FD_vein, AVR, Tortuosity_density_artery, Tortuosity_density_vein, conc_rp_artery, conc_rp_vein, DFD_DML_ratio, discTilt, absDiscTors)

# Start estimating 
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
m9               <- gcrq(discTilt ~ ps(SER), data=inputD, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m9, file.path('fittedCentileCurves', 'm9.rds'))
# Absolute disc torsion (save as 'm10')
m10              <- gcrq(absDiscTors ~ ps(SER), data=inputD, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m10, file.path('fittedCentileCurves', 'm10.rds'))




#### For SENSITIVITY ANALYSIS 5 in 'survivalAnalyses.R': reconstruct the centile #### 
####      charts using the data derived based on the 0.01% outlier threshold     ####

# Randomly select one eye (25,447 eyes)
set.seed(1)
inputD_SA        <- subset(d, !outlierSA) %>% group_by(id) %>% slice_sample(n = 1)
nrow(inputD_SA)

# Extract relevant columns
inputD           <- inputD_SA %>% select(name, age, sex, catSurgery1, SER, FD_artery, FD_vein, AVR, Tortuosity_density_artery, Tortuosity_density_vein, conc_rp_artery, conc_rp_vein, DFD_DML_ratio, discTilt, absDiscTors)

# Arterial fractal dimension (save as 'm1_sensitivity')
m1_SA            <- gcrq(FD_artery ~ ps(SER), data=inputD_SA, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m1_SA, file.path('fittedCentileCurves', 'm1_sensitivity.rds'))

# Venous fractal dimension (save as 'm2_sensitivity')
m2_SA            <- gcrq(FD_vein ~ ps(SER), data=inputD_SA, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m2_SA, file.path('fittedCentileCurves', 'm2_sensitivity.rds'))

# AVR (save as 'm3_sensitivity')
m3_SA            <- gcrq(AVR ~ ps(SER), data=inputD_SA, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m3_SA, file.path('fittedCentileCurves', 'm3_sensitivity.rds'))

# Arterial tortuosity (save as 'm4_sensitivity')
m4_SA            <- gcrq(Tortuosity_density_artery ~ ps(SER), data=inputD_SA, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m4_SA, file.path('fittedCentileCurves', 'm4_sensitivity.rds'))

# Venous tortuosity (save as 'm5_sensitivity')
m5_SA            <- gcrq(Tortuosity_density_vein ~ ps(SER), data=inputD_SA, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m5_SA, file.path('fittedCentileCurves', 'm5_sensitivity.rds'))

# Arterial concavity (save as 'm6_sensitivity')
m6_SA            <- gcrq(conc_rp_artery ~ ps(SER), data=inputD_SA, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m6_SA, file.path('fittedCentileCurves', 'm6_sensitivity.rds'))

# Venous concavity (save as 'm7_sensitivity')
m7_SA            <- gcrq(conc_rp_vein ~ ps(SER), data=inputD_SA, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m7_SA, file.path('fittedCentileCurves', 'm7_sensitivity.rds'))

# DFD_DML_ratio (save as 'm8_sensitivity')
m8_SA            <- gcrq(DFD_DML_ratio ~ ps(SER), data=inputD_SA, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m8_SA, file.path('fittedCentileCurves', 'm8_sensitivity.rds'))

# Disc tilt (save as 'm9_sensitivity')
m9_SA            <- gcrq(discTilt ~ ps(SER), data=inputD_SA, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m9_SA, file.path('fittedCentileCurves', 'm9_sensitivity.rds'))

# Absolute disc torsion (save as 'm10_sensitivity')
m10_SA           <- gcrq(absDiscTors ~ ps(SER), data=inputD_SA, tau=seq(0.05, 0.95, 0.05), n.boot=50)
saveRDS(m10_SA, file.path('fittedCentileCurves', 'm10_sensitivity.rds'))
















