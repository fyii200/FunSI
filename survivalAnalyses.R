# Author :  Fabian Yii                                                      #
# Email  :  fabian.yii@ed.ac.uk                                             #
# Info   :  Association of baseline Fundus Stretch Index (FunSI) with the   #
#           risk of developing rhegmatogenous retinal detachment (RD) and   #
#           parimary open-angle glaucoma (POAG)                             #

rm(list=ls())
library(ggplot2)
library(gridExtra)
library(sjPlot)
library(car)
library(rms)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(survival)
library(survminer)
library(naniar)
library(mice)
# setwd("..")
options(scipen = 10000)




#### Read, preprocess and merge relevant datasets ####

# Baseline data including FunSI from 
d                    <- read.csv(file.path('data', 'FunSI.csv'), check.names = F)

# First-occurrence data
firstOccur           <- read.csv(file.path('data', 'firstOccurrenceData.csv')) 

# Select relevant columns
firstOccur           <- firstOccur %>% select(Participant.ID, 
                                              Date.of.death...Instance.0,
                                              Date.lost.to.follow.up,
                                              Reason.lost.to.follow.up,
                                              Age.when.loss.of.vision.due.to.injury.or.trauma.diagnosed...Instance.0,
                                              Date.H33.first.reported..retinal.detachments.and.breaks.,
                                              Source.of.report.of.H33..retinal.detachments.and.breaks.,
                                              Date.H40.first.reported..glaucoma.,
                                              Source.of.report.of.H40..glaucoma.) 

# Rename the columns 
names(firstOccur)    <- c('id', 
                          'dateDeath', 
                          'dateLossFU', 
                          'reasonLossFU', 
                          'ageVisionLossDueToInjuryOrTrauma', 
                          'firstDateRetinalDetachmentBreaks', 
                          'sourceRetinalDetachmentBreaks', 
                          'firstDateGlaucoma',
                          'sourceGlaucoma')

# Fundus-predicted spherical equivalent refraction (SER)
FER                  <- read.csv(file.path('data', 'FERall.csv')) %>% filter(!duplicated(name))

# Merge all datasets
d                    <- merge(d, firstOccur, by = 'id', all.x = T)
d                    <- merge(d, FER, by = 'name', all.x = T)

# Get earliest cataract surgery date per ID
d$catSurgery         <- pmin(d$catSurgery1, d$catSurgery2, d$catSurgery3, d$catSurgery4, na.rm = T)

# Ensure all date variables are correctly coded as datetime
cvtDate              <- function(x) as.Date(x, '%d/%m/%Y')
d                    <- d %>% mutate(dateDeath                        = cvtDate(dateDeath),
                                     dateLossFU                       = cvtDate(dateLossFU),
                                     firstDateRetinalDetachmentBreaks = cvtDate(firstDateRetinalDetachmentBreaks),
                                     firstDateGlaucoma                = cvtDate(firstDateGlaucoma))

# Read hospital diagnosis data (ICD-9 or ICD-10) and rename columns
diag                 <- read.csv(file.path('data', 'hospitalDiagnosisData.csv'))
diag                 <- diag %>% rename(id    = Participant.ID.participant...eid.,
                                        ICD10 = Diagnoses...ICD10,
                                        ICD09 = Diagnoses...ICD9)

# Identify IDs with a record of RD as derived by the UK Biobank team (includes various subtypes)
IDsWithRDFirstOcc    <- unique(subset(d, !is.na(firstDateRetinalDetachmentBreaks))$id)

# For each ID ascertain if it's rhegmatogenous RD
sum(unique(substr(diag$ICD09, 1, 3)) == '361') # no RD-related ICD-9 code, so use ICD-10 code below
d$rhegmaRDprelim     <- F
for(i in IDsWithRDFirstOcc){
  # Check if the given ID has a hospital record of rhegmatogenous RD (H33.0)
  d[d$id == i,]$rhegmaRDprelim <- 'H33.0' %in% substring(subset(diag, id == i)$ICD10, 1, 5)
  }

# Hospital procedural data (OPCS-4) in long format where each operation has its own row
operations           <- read.csv(file.path('data', 'proceduralData.csv'))

# Rename columns
operations           <- operations %>% rename(id            = Participant.ID.participant...eid.,  
                                              OPCS4         = Operative.procedures...OPCS4,
                                              operationDate = Date.of.operation) %>% mutate(operationDate = cvtDate(operationDate))

# Identify IDs with POAG using a backward elimination approach: among IDs with a 
# record of glaucoma, as derived by the UK Biobank team, remove those with any
# OPCS-4 procedural or ICD codes related to primary angle-closure glaucoma (PACG) 
# or other/secondary glaucoma (OSG)
sum(unique(substr(diag$ICD09, 1, 3)) == '365') # no ICD-9 code related to glaucoma, so only use ICD-10 below
angleOps             <- c('C62.1', 'C62.2', 'C62.3')                                      
PACG_ICD10           <- c('H40.2')                                                        
OSG_ICD10            <- c('H40.3', 'H40.4', 'H40.5', 'H40.6', 'H40.8', 'H42.0', 'H42.8')
nonPOAG_IDs          <- unique(c(subset(operations, (substr(OPCS4, 1, 5) %in% angleOps))$id,
                                 subset(diag, (substr(ICD10, 1,5) %in% PACG_ICD10))$id,
                                 subset(diag, (substr(ICD10, 1,5) %in% OSG_ICD10))$id))
d                    <- d %>% mutate(nonPOAG = ifelse(id %in% nonPOAG_IDs, T, F)) 
length(unique(subset(d, nonPOAG)$id))          # 43 IDs have non-POAG glaucoma codes
d$POAG               <- !is.na(d$firstDateGlaucoma) & !d$nonPOAG

# Remove unreasonable intraocular pressure (IOP) and corneal hysteresis (CH) values
# Note 1: max IOP in glaucomatous eyes reported by a population-based study was 45.6mmHg (https://www.bmj.com/content/358/bmj.j3889)
# Note 2: max CH based on https://karger.com/ore/article-pdf/46/4/187/3367017/000326896.pdf and https://www.sciencedirect.com/science/article/pii/S0002939411008531?via%3Dihub#sec2
d                    <- d %>% mutate(IOP = ifelse(IOP == 0 | IOP > 45, NA, IOP),
                                     CH  = ifelse(CH == 0 | CH > 15, NA, CH))



#### Compute survival time ####                       

# Right-censoring (RC) date
RCdate                             <- as.Date('2022-12-31') 

# If baseline assessment date is missing, infer by summing year of birth and age at the baseline visit
missingVisitDate                   <- which(is.na(d$visitDate))
if(length(missingVisitDate) > 0){d[missingVisitDate, ]$visitDate <- format(as.Date(paste0(d$YOB[missingVisitDate] + d$age[missingVisitDate], '-01-01')), '%Y-%m-%d') }

# Create new columns in preparation for survival analysis 
# 'rhegmaRD'         : indicates if a participant had a record of rhegmatogenous RD during follow-up
# 'POAG'             : indicates if a participant had a record of POAG during follow-up
# 'xx' + 'cumYear'   : time from baseline to the first recorded occurrence of the 'xx' event (disease) or to the right-censoring date 
# 'previousTrauma'   : whether there's a self-reported trauma/injury 'resulting in a loss of vision'
d$RDobservedPeriod                 <- difftime(pmin(d$dateDeath, d$dateLossFU, RCdate, d$firstDateRetinalDetachmentBreaks, d$catSurgery, na.rm = T), d$visitDate, units = 'days')/365.25
d$rhegmaRD                         <- ifelse(d$rhegmaRDprelim == T & (is.na(d$catSurgery) | d$catSurgery > d$firstDateRetinalDetachmentBreaks), T, F)
d$rhegmaRDcumYear                  <- ifelse(d$rhegmaRD == T, difftime(d$firstDateRetinalDetachmentBreaks, d$visitDate, units = 'days')/365.25, d$RDobservedPeriod)
d$glaucomaObservedPeriod           <- difftime(pmin(d$dateDeath, d$dateLossFU, RCdate, d$firstDateGlaucoma, na.rm = T), d$visitDate, units = 'days')/365.25
d$POAG                             <- ifelse(!is.na(d$firstDateGlaucoma) & !d$nonPOAG, T, F)
d$POAGcumYear                      <- ifelse(d$POAG == T, difftime(d$firstDateGlaucoma, d$visitDate, units = 'days')/365.25, d$glaucomaObservedPeriod)
d$previousTrauma                   <- d$ageVisionLossDueToInjuryOrTrauma != '' & !is.na(d$ageVisionLossDueToInjuryOrTrauma) & d$ageVisionLossDueToInjuryOrTrauma != 'Do not know'

# Average data across eyes if data from both eyes of an individual are available
# For SENSITIVITY ANALYSIS 5 later
d_SA                               <- subset(d, !outlierSA) %>% group_by(id) %>% mutate(indSER       = ifelse(sum(!is.na(SER)) == 2, (SER[1] + SER[2]) / 2, SER),
                                                                                        indFunSI_SA  = ifelse(sum(!is.na(FunSI_SA)) == 2, (FunSI_SA[1] + FunSI_SA[2]) / 2, FunSI_SA),
                                                                                        indIOP       = ifelse(sum(!is.na(IOP)) == 2, (IOP[1] + IOP[2]) / 2, IOP),
                                                                                        indCH        = ifelse(sum(!is.na(CH)) == 2, (CH[1] + CH[2]) / 2, CH))
# For the PRIMARY ANALYSIS
d                                  <- subset(d, !outlier) %>% group_by(id) %>% mutate(indSER                    = ifelse(sum(!is.na(SER)) == 2, (SER[1] + SER[2]) / 2, SER),
                                                                                      indFunSI                  = ifelse(sum(!is.na(FunSI)) == 2, (FunSI[1] + FunSI[2]) / 2, FunSI),
                                                                                      indFunSI_wo_Aconc         = ifelse(sum(!is.na(FunSI_wo_Aconc)) == 2, (FunSI_wo_Aconc[1] + FunSI_wo_Aconc[2]) / 2, FunSI_wo_Aconc),
                                                                                      indFunSI_wo_AFD           = ifelse(sum(!is.na(FunSI_wo_AFD)) == 2, (FunSI_wo_AFD[1] + FunSI_wo_AFD[2]) / 2, FunSI_wo_AFD),
                                                                                      indFunSI_wo_Atort         = ifelse(sum(!is.na(FunSI_wo_Atort)) == 2, (FunSI_wo_Atort[1] + FunSI_wo_Atort[2]) / 2, FunSI_wo_Atort),
                                                                                      indFunSI_wo_AVR           = ifelse(sum(!is.na(FunSI_wo_AVR)) == 2, (FunSI_wo_AVR[1] + FunSI_wo_AVR[2]) / 2, FunSI_wo_AVR),
                                                                                      indFunSI_wo_DFD_DML_ratio = ifelse(sum(!is.na(FunSI_wo_DFD_DML_ratio)) == 2, (FunSI_wo_DFD_DML_ratio[1] + FunSI_wo_DFD_DML_ratio[2]) / 2, FunSI_wo_DFD_DML_ratio),
                                                                                      indFunSI_wo_Dtilt         = ifelse(sum(!is.na(FunSI_wo_Dtilt)) == 2, (FunSI_wo_Dtilt[1] + FunSI_wo_Dtilt[2]) / 2, FunSI_wo_Dtilt),
                                                                                      indFunSI_wo_Dtorsion      = ifelse(sum(!is.na(FunSI_wo_Dtorsion)) == 2, (FunSI_wo_Dtorsion[1] + FunSI_wo_Dtorsion[2]) / 2, FunSI_wo_Dtorsion),
                                                                                      indFunSI_wo_Vconc         = ifelse(sum(!is.na(FunSI_wo_Vconc)) == 2, (FunSI_wo_Vconc[1] + FunSI_wo_Vconc[2]) / 2, FunSI_wo_Vconc),
                                                                                      indFunSI_wo_VFD           = ifelse(sum(!is.na(FunSI_wo_VFD)) == 2, (FunSI_wo_VFD[1] + FunSI_wo_VFD[2]) / 2, FunSI_wo_VFD),
                                                                                      indFunSI_wo_Vtort         = ifelse(sum(!is.na(FunSI_wo_Vtort)) == 2, (FunSI_wo_Vtort[1] + FunSI_wo_Vtort[2]) / 2, FunSI_wo_Vtort),
                                                                                      indPredSER_TTA            = ifelse(sum(!is.na(predSER_TTA)) == 2, (predSER_TTA[1] + predSER_TTA[2]) / 2, predSER_TTA),
                                                                                      indIOP                    = ifelse(sum(!is.na(IOP)) == 2, (IOP[1] + IOP[2]) / 2, IOP),
                                                                                      indCH                     = ifelse(sum(!is.na(CH)) == 2, (CH[1] + CH[2]) / 2, CH))

# Extract data from the worse (more myopic) eye (for a senstivity analysis later)
d                                  <- d %>% mutate(worseEye                  = which(SER == min(SER)),
                                                   worseSER                  = ifelse(sum(!is.na(SER)) == 2, SER[worseEye], SER),
                                                   worseFunSI                = ifelse(sum(!is.na(FunSI)) == 2, FunSI[worseEye], FunSI),
                                                   worseIOP                  = ifelse(sum(!is.na(IOP)) == 2, IOP[worseEye], IOP),
                                                   worseCH                   = ifelse(sum(!is.na(CH)) == 2, CH[worseEye], CH))




#### Survival analysis for rhegmatogenous RD ####

# Population at risk: 25,067 IDs without any prior history of RD/breaks of any subtype at baseline
RDdata                    <- subset(d, !duplicated(id) & RDobservedPeriod > 0) 
nrow(RDdata)

# Data missing completely at random (MCAR) based on Little's test
variables                 <- c('indSER', 'indFunSI', 'age', 'sex', 'townsend', 'ethnicBinary', 'diabetes', 'hypertension', 'previousTrauma', 'rhegmaRD', 'rhegmaRDcumYear')
mcar_test(RDdata[, variables])

# Complete-case analysis because MCAR assumption was met (37 IDs removed, 25,030 remained)
RDdata                    <- subset(RDdata, !is.na(townsend))
nrow(RDdata)

# Create a new column storing categorised FunSI (for visualisation purposes; Kaplan-Meier survival plot)
RDdata$FunSIgroup         <- cut(RDdata$indFunSI, quantile(RDdata$indFunSI, seq(0,1,0.25)), include.lowest = T)
levels(RDdata$FunSIgroup) <- c('0.19-0.44', '0.44-0.50', '0.50-0.56', '0.56-0.88')

# 148 incident rhegmatogenous RD cases, with an event rate of 4.9 per 10000 person-years
length(unique(RDdata[RDdata$rhegmaRD == T, ]$id))
(sum(RDdata$rhegmaRD) / sum(RDdata$rhegmaRDcumYear))*10000

# 124 participants with a self-reported history of ocular trauma/injury
length(unique(subset(RDdata, previousTrauma == T)$id))

# Plot Kaplan-Meier survival curve stratified by baseline FRO 
rhegRD_KMfit              <- survfit(Surv(rhegmaRDcumYear, rhegmaRD) ~ FunSIgroup, RDdata)
rhegRD_KMplot             <- ggsurvplot(rhegRD_KMfit,
                                        xlim         = c(0, 12),
                                        ylim         = c(0, 0.01),
                                        fun          = 'event',
                                        surv.scale   = 'percent',
                                        break.x.by   = 4,                             
                                        size         = 0.6,
                                        add.all      = T,
                                        censor       = F,
                                        xlab         = 'Years elapsed',
                                        ylab         = 'Cumulative incidence\n',
                                        subtitle     = 'Rhegmatogenous retinal detachment\n',
                                        legend       = 'top',
                                        legend.title = 'FunSI',
                                        legend.labs  = c('Overall', levels(RDdata$FunSIgroup)),
                                        ggtheme      = theme_minimal() + theme(panel.grid.major = element_blank(), 
                                                                               panel.grid.minor = element_blank(),
                                                                               plot.subtitle    = element_text(hjust = 0.5), 
                                                                               axis.title.x     = element_text(face = 'bold', hjust = 1),
                                                                               axis.title.y     = element_text(face = 'bold', hjust = 1)),
                                        palette       = c('gray', brewer.pal(11, 'RdYlBu')[c(8, 10, 4, 1)]),
                                        risk.table    = T,
                                        tables.height = 0.2,
                                        fontsize      = 2.5,
                                        tables.y.text = F,
                                        tables.theme  = theme_cleantable() + theme(plot.title = element_text(size = 9, face = 'italic')))                              
rhegRD_KMplot$plot        <- rhegRD_KMplot$plot + guides(col = guide_legend(override.aes = list(linewidth = 2)))
pdf(file.path('figures', 'rhegmaRD_KMplot.pdf'), width = 6.5, height = 7)
print(grid.arrange(rhegRD_KMplot$plot, rhegRD_KMplot$table, ncol = 1, heights = c(2, 0.5)))
dev.off()

# PRIMARY ANALYSIS: univariable Cox regression
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ scale(indFunSI), data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ age, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ sex, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ ethnicBinary, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ townsend, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ diabetes, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ hypertension, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ previousTrauma, data = RDdata))

# PRIMARY ANALYSIS: multivariable Cox regression
# Baseline model (without FunSI)
baseline                  <- coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + age + sex + townsend + ethnicBinary + diabetes + hypertension + previousTrauma, data = RDdata)
tab_model(baseline); max(vif(baseline))   
# Baseline model + FunSI (linear)
lin                       <- update(baseline, .~. + scale(indFunSI), RDdata)
tab_model(lin); max(vif(lin))   
# Baseline model + FunSI (linear) + FunSI (quadratic)
quad                      <- update(baseline, .~. + scale(poly(indFunSI, 2)), RDdata)
tab_model(quad); max(vif(quad))   
# Baseline model + FunSI (linear) + FunSI (quadratic) + FunSI (cubic)
cubic                     <- update(baseline, .~. + scale(poly(indFunSI, 3)), RDdata)
tab_model(cubic); max(vif(cubic))   
# Concordance index (higher means better discriminative power)
concordance(baseline, lin, quad, cubic) 

# LEAVE-ONE-FEATURE-OUT analysis 
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + scale(indFunSI_wo_Aconc) + age + sex + townsend + ethnicBinary + diabetes + hypertension + previousTrauma, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + scale(indFunSI_wo_AFD) + age + sex + townsend + ethnicBinary + diabetes + hypertension + previousTrauma, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + scale(indFunSI_wo_Atort) + age + sex + townsend + ethnicBinary + diabetes + hypertension + previousTrauma, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + scale(indFunSI_wo_AVR) + age + sex + townsend + ethnicBinary + diabetes + hypertension + previousTrauma, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + scale(indFunSI_wo_DFD_DML_ratio) + age + sex + townsend + ethnicBinary + diabetes + hypertension + previousTrauma, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + scale(indFunSI_wo_Dtilt) + age + sex + townsend + ethnicBinary + diabetes + hypertension + previousTrauma, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + scale(indFunSI_wo_Dtorsion) + age + sex + townsend + ethnicBinary + diabetes + hypertension + previousTrauma, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + scale(indFunSI_wo_Vconc) + age + sex + townsend + ethnicBinary + diabetes + hypertension + previousTrauma, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + scale(indFunSI_wo_VFD) + age + sex + townsend + ethnicBinary + diabetes + hypertension + previousTrauma, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + scale(indFunSI_wo_Vtort) + age + sex + townsend + ethnicBinary + diabetes + hypertension + previousTrauma, data = RDdata))

# SUBGROUP ANALYSIS: FRO included as another covariate
# Compute FRO in a subset of IDs not previously used to train the FRO model (15,979 IDs)
RDdataSub                 <- subset(RDdata, !is.na(indPredSER_TTA)); nrow(RDdataSub)
RDdataSub$indFRO          <- residuals(lm(indPredSER_TTA ~ indSER, RDdataSub))
# Fit Cox regression
baseline                  <- coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + age + sex + age + sex + townsend + ethnicBinary + diabetes + hypertension + previousTrauma, data = RDdataSub)
FRO                       <- update(baseline, .~. + scale(indFRO), RDdataSub)  # baseline model + FRO
FunSI                     <- update(FRO, .~. + scale(indFunSI), RDdataSub)     # baseline model + FRO + FunSI
tab_model(FunSI); max(vif(FunSI)); sd(RDdataSub$indFRO); sd(RDdataSub$indFunSI)
# Check concordance (higher means better discriminative power)
concordance(baseline, FRO, FunSI) 

# SENSITIVITY ANALYSIS 2: only include IDs with cyl ≤ 2D (23,042 analysed)
RDdataSens                <- subset(RDdata, cyl <= 2); nrow(RDdataSens)
sens                      <- coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + scale(indFunSI) + age + sex + townsend + ethnicBinary + diabetes + hypertension + previousTrauma, data = RDdataSens)
tab_model(sens); max(vif(sens)); sd(RDdataSens$indFunSI) 

# SENSITIVITY ANALYSIS 3: remove 8 rhegmatogenous RD cases that occurred within 1 year of baseline (25,022 analysed)
excludeIDs                <- subset(RDdata, rhegmaRD & (rhegmaRDcumYear < 1))$id; length(excludeIDs)
RDdataSens                <- subset(RDdata, ! id %in% excludeIDs); nrow(RDdataSens)
sens                      <- coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + scale(indFunSI) + age + sex + townsend + ethnicBinary + diabetes + hypertension + previousTrauma, data = RDdataSens)
tab_model(sens); max(vif(sens)); sd(RDdataSens$indFunSI) 

# SENSITIVITY ANALYSIS 4: use data from the more myopic eye when both eyes available (25,030 analysed)
sens                      <- coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ worseSER + scale(worseFunSI) + age + sex + townsend + ethnicBinary + diabetes + hypertension + previousTrauma, data = RDdata)
tab_model(sens); max(vif(sens)); sd(RDdata$worseFunSI) 

# SENSITIVITY ANALYSIS 5: outlier threshold of 0.01% rather than 0.1% (25,252 analysed)
RDdataSens                <- subset(d_SA, !duplicated(id) & RDobservedPeriod > 0)
sens                      <- coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + scale(indFunSI_SA) + age + sex + townsend + ethnicBinary + diabetes + hypertension + previousTrauma, data = RDdataSens)
tab_model(sens); max(vif(sens)); sd(RDdataSens$indFunSI_SA) 




#### Survival analysis for POAG ####

# Population at risk: 24,835 IDs without any prior history of glaucoma of any subtype at baseline
POAGdata                    <- subset(d, !duplicated(id) & glaucomaObservedPeriod > 0)
nrow(POAGdata)

# Little's MCAR test assumption not met
variables                   <- c('indSER', 'indFunSI', 'age', 'sex', 'townsend', 'ethnicBinary', 'diabetes', 'hypertension', 'indIOP', 'indCH', 'POAG', 'POAGcumYear')
mcar_test(POAGdata[, variables])

# Multiple imputation because MCAR assumption was not supported
set.seed(1) 
imputedPOAGdata             <- mice(POAGdata[, variables], m = 10, method = 'pmm') 

# 460 out of 24835 IDs had incident POAG (15.0 per 10000 person-years)
table(POAGdata$POAG)
(sum(POAGdata$POAG) / sum(POAGdata$POAGcumYear))*10000

# Create a new column storing categorised SI (for visualisation purposes; Kaplan-Meier survival plot)
POAGdata$FunSIgroup         <- cut(POAGdata$indFunSI, quantile(POAGdata$indFunSI, seq(0,1,0.25)), include.lowest = T)
levels(POAGdata$FunSIgroup) <- c('0.19-0.44', '0.44-0.50', '0.50-0.56', '0.56-0.88')

# Linkage sources: majority sourced from hospital or primary care (only 23 exclusively self-reported cases)
POAGdata                    <- POAGdata %>% mutate(sourceGlaucoma = ifelse(is.na(sourceGlaucoma), '', sourceGlaucoma))
table(subset(POAGdata, !duplicated(id) & POAG == T)$sourceGlaucoma)

# Plot Kaplan-Meier survival curve stratified by baseline SI (any glaucoma as event)
allPOAG_KMfit               <- survfit(Surv(POAGcumYear, POAG) ~ FunSIgroup, POAGdata)
allPOAG_KMplot              <- ggsurvplot(allPOAG_KMfit,
                                          xlim         = c(0, 12),
                                          ylim         = c(0, 0.03),
                                          fun          = 'event',
                                          surv.scale   = 'percent',
                                          break.x.by   = 4,
                                          size         = 0.6,
                                          add.all      = T,
                                          censor       = F,
                                          xlab         = 'Years elapsed',
                                          ylab         = 'Cumulative incidence\n',
                                          subtitle     = 'Primary open-angle glaucoma\n',
                                          legend       = 'top',
                                          legend.title = 'FunSI',
                                          legend.labs  = c('Overall', levels(POAGdata$FunSIgroup)),
                                          ggtheme      = theme_minimal() + theme(panel.grid.major   = element_blank(), 
                                                                                 panel.grid.minor   = element_blank(),
                                                                                 plot.subtitle      = element_text(hjust = 0.5), 
                                                                                 axis.title.x       = element_text(face = 'bold', hjust = 1),
                                                                                 axis.title.y       = element_text(face = 'bold', hjust = 1)), 
                                          palette       = c('gray', brewer.pal(11, 'RdYlBu')[c(8, 10, 4, 1)]),
                                          risk.table    = T,
                                          tables.height = 0.2,
                                          fontsize      = 2.5,
                                          tables.y.text = F,
                                          tables.theme  = theme_cleantable() + theme(plot.title = element_text(size = 8, face = 'italic'))) 
allPOAG_KMplot$plot         <- allPOAG_KMplot$plot + guides(col = guide_legend(override.aes = list(linewidth = 2))) 
pdf(file.path('figures', 'POAG_KMplot.pdf'), width = 6.5, height = 7)
print(grid.arrange(allPOAG_KMplot$plot, allPOAG_KMplot$table, ncol = 1, heights = c(2, 0.5)))
dev.off()

# PRIMARY ANALYSIS: univariable Cox regression
set.seed(1)
tab_model(coxph(Surv(POAGcumYear, POAG) ~ indSER, data = POAGdata))
tab_model(coxph(Surv(POAGcumYear, POAG) ~ scale(indFunSI), data = POAGdata))
tab_model(pool(with(data = imputedPOAGdata, exp = coxph(Surv(POAGcumYear, POAG) ~ indIOP))))
tab_model(pool(with(data = imputedPOAGdata, exp = coxph(Surv(POAGcumYear, POAG) ~ indCH))))
tab_model(coxph(Surv(POAGcumYear, POAG) ~ age, data = POAGdata))
tab_model(coxph(Surv(POAGcumYear, POAG) ~ sex, data = POAGdata))
tab_model(coxph(Surv(POAGcumYear, POAG) ~ ethnicBinary, data = POAGdata))
tab_model(pool(with(data = imputedPOAGdata, exp = coxph(Surv(POAGcumYear, POAG) ~ townsend))))
tab_model(coxph(Surv(POAGcumYear, POAG) ~ diabetes, data = POAGdata))
tab_model(coxph(Surv(POAGcumYear, POAG) ~ hypertension, data = POAGdata))

# PRIMARY ANALYSIS: multivariable Cox regression
# Baseline model (without FunSI)
baseline                    <- with(data = imputedPOAGdata, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension))
tab_model(pool(baseline)); for(i in 1:10){print(max(vif(baseline$analyses[[i]])))}
# Baseline model + FunSI (linear)
lin                         <- with(data = imputedPOAGdata, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(indFunSI) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension))
tab_model(pool(lin)); for(i in 1:10){print(max(vif(lin$analyses[[i]])))}
# Baseline model + FunSI (linear) + FunSI (quadratic)
quad                        <- with(data = imputedPOAGdata, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(poly(indFunSI,2)) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension))
tab_model(pool(quad)); for(i in 1:10){print(max(vif(quad$analyses[[i]])))}
# Baseline model + FunSI (linear) + FunSI (quadratic) + FunSI (cubic)
cubic                       <- with(data = imputedPOAGdata, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(poly(indFunSI,3)) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension))
tab_model(pool(cubic)); for(i in 1:10){print(max(vif(cubic$analyses[[i]])))}
# Concordance index
mean(unlist(sapply(baseline$analyses, concordance)[1,])) 
mean(unlist(sapply(lin$analyses, concordance)[1,])) 
mean(unlist(sapply(quad$analyses, concordance)[1,])) 
mean(unlist(sapply(cubic$analyses, concordance)[1,]))

# LEAVE-ONE-FEATURE-OUT analysis 
imputedPOAGdata$data        <- replace(imputedPOAGdata$data, 2, POAGdata$indFunSI_wo_Aconc) 
FunSI_wo_Aconc              <- with(data = imputedPOAGdata, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(indFunSI) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension)); tab_model(pool(FunSI_wo_Aconc))
imputedPOAGdata$data        <- replace(imputedPOAGdata$data, 2, POAGdata$indFunSI_wo_AFD) 
FunSI_wo_AFD                <- with(data = imputedPOAGdata, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(indFunSI) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension)); tab_model(pool(FunSI_wo_AFD))
imputedPOAGdata$data        <- replace(imputedPOAGdata$data, 2, POAGdata$indFunSI_wo_Atort) 
FunSI_wo_Atort              <- with(data = imputedPOAGdata, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(indFunSI) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension)); tab_model(pool(FunSI_wo_Atort))
imputedPOAGdata$data        <- replace(imputedPOAGdata$data, 2, POAGdata$indFunSI_wo_AVR) 
FunSI_wo_AVR                <- with(data = imputedPOAGdata, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(indFunSI) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension)); tab_model(pool(FunSI_wo_AVR))
imputedPOAGdata$data        <- replace(imputedPOAGdata$data, 2, POAGdata$indFunSI_wo_DFD_DML_ratio) 
FunSI_wo_DFD_DML_ratio      <- with(data = imputedPOAGdata, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(indFunSI) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension)); tab_model(pool(FunSI_wo_DFD_DML_ratio))
imputedPOAGdata$data        <- replace(imputedPOAGdata$data, 2, POAGdata$indFunSI_wo_Dtilt) 
FunSI_wo_Dtilt              <- with(data = imputedPOAGdata, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(indFunSI) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension)); tab_model(pool(FunSI_wo_Dtilt))
imputedPOAGdata$data        <- replace(imputedPOAGdata$data, 2, POAGdata$indFunSI_wo_Dtorsion) 
FunSI_wo_Dtorsion           <- with(data = imputedPOAGdata, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(indFunSI) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension)); tab_model(pool(FunSI_wo_Dtorsion))
imputedPOAGdata$data        <- replace(imputedPOAGdata$data, 2, POAGdata$indFunSI_wo_Vconc) 
FunSI_wo_Vconc              <- with(data = imputedPOAGdata, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(indFunSI) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension)); tab_model(pool(FunSI_wo_Vconc))
imputedPOAGdata$data        <- replace(imputedPOAGdata$data, 2, POAGdata$indFunSI_wo_VFD) 
FunSI_wo_VFD                <- with(data = imputedPOAGdata, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(indFunSI) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension)); tab_model(pool(FunSI_wo_VFD))
imputedPOAGdata$data        <- replace(imputedPOAGdata$data, 2, POAGdata$indFunSI_wo_Vtort) 
FunSI_wo_Vtort              <- with(data = imputedPOAGdata, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(indFunSI) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension)); tab_model(pool(FunSI_wo_Vtort))

# SUBGROUP ANALYSIS: FRO included as another covariate
# Compute FRO in a subset of IDs not previously used to train the FRO model (15,754 IDs)
POAGdataSub                 <- subset(POAGdata, !is.na(indPredSER_TTA)); nrow(POAGdataSub)
indFRO                      <- residuals(lm(indPredSER_TTA ~ indSER, POAGdataSub))
set.seed(1)
POAGdataSub                 <- mice(POAGdataSub[, variables], m = 10, method = 'pmm') 
POAGdataSub$data$indFRO     <- indFRO
sd(POAGdataSub$data$indFRO); sd(POAGdataSub$data$indFunSI)
# Fit Cox regression
baseline                    <- with(data = POAGdataSub, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension))
tab_model(pool(baseline)); for(i in 1:10){print(max(vif(baseline$analyses[[i]])))}
FRO                         <- with(data = POAGdataSub, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(indFRO) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension))
tab_model(pool(FRO)); for(i in 1:10){print(max(vif(FRO$analyses[[i]])))}
FunSI_lin                   <- with(data = POAGdataSub, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(indFRO) + scale(indFunSI) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension))
tab_model(pool(FunSI_lin)); for(i in 1:10){print(max(vif(FunSI_lin$analyses[[i]])))}
# Concordance index 
mean(unlist(sapply(baseline$analyses, concordance)[1,])) 
mean(unlist(sapply(FRO$analyses, concordance)[1,])) 
mean(unlist(sapply(FunSI_lin$analyses, concordance)[1,])) 

# SENSITIVITY ANALYSIS remove 23 self-reported cases and 118 non-hospital cases (24,694 analysed)
POAGdataSens                <- POAGdata %>% mutate(sourceGlaucoma = ifelse(POAG, sourceGlaucoma, ''))
POAGdataSens                <- subset(POAGdataSens, sourceGlaucoma == '' | sourceGlaucoma == 'Hospital admissions data only' | sourceGlaucoma == 'Hospital admissions data and other source(s)'); nrow(POAGdataSens)
set.seed(1)
POAGdataSens                <- mice(POAGdataSens[, variables], m = 10, method = 'pmm') 
sens                        <- with(data = POAGdataSens, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(indFunSI) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension))
tab_model(pool(sens)); for(i in 1:10){print(max(vif(sens$analyses[[i]])))}; sd(complete(POAGdataSens, action = 1)$indFunSI)

# SENSITIVITY ANALYSIS 2: only include IDs with cyl ≤ 2D (22,882 analysed)
POAGdataSens                <- subset(POAGdata, cyl <= 2); nrow(POAGdataSens)
set.seed(1)
POAGdataSens                <- mice(POAGdataSens[, variables], m = 10, method = 'pmm') 
sens                        <- with(data = POAGdataSens, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(indFunSI) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension))
tab_model(pool(sens)); for(i in 1:10){print(max(vif(sens$analyses[[i]])))}; sd(complete(POAGdataSens, action = 1)$indFunSI)

# SENSITIVITY ANALYSIS 3: remove 3 POAG cases that occurred within 1 year of baseline (24,832 analysed)
excludeIDs                  <- subset(POAGdata, glaucoma & (POAGcumYear < 1))$id; length(excludeIDs)
POAGdataSens                <- subset(POAGdata, ! id %in% excludeIDs)
set.seed(1)
POAGdataSens                <- mice(POAGdataSens[, variables], m = 10, method = 'pmm') 
sens                        <- with(data = POAGdataSens, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(indFunSI) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension))
tab_model(pool(sens)); for(i in 1:10){print(max(vif(sens$analyses[[i]])))}; sd(complete(POAGdataSens, action = 1)$indFunSI)

# SENSITIVITY ANALYSIS 4: use data from the more myopic eye when both eyes available (24,835 analysed)
set.seed(1)
POAGdataSens                <- mice(POAGdata[, c('worseSER', 'worseFunSI', 'age', 'sex', 'townsend', 'ethnicBinary', 'diabetes', 'hypertension', 'worseIOP', 'worseCH', 'POAG', 'POAGcumYear')], m = 10, method = 'pmm') 
sens                        <- with(data = POAGdataSens, exp = coxph(Surv(POAGcumYear, POAG) ~ worseSER + scale(worseFunSI) + age + sex + townsend + ethnicBinary + worseIOP + worseCH + diabetes + hypertension))
tab_model(pool(sens)); for(i in 1:10){print(max(vif(sens$analyses[[i]])))}; sd(complete(POAGdataSens, action = 1)$worseFunSI)

# SENSITIVITY ANALYSIS 5: outlier threshold of 0.01% rather than 0.1% (25,052 analysed)
POAGdataSens                <- subset(d_SA, !duplicated(id) & glaucomaObservedPeriod > 0)
set.seed(1)
POAGdataSens                <- mice(POAGdataSens[, replace(variables, 2, 'indFunSI_SA')], m = 10, method = 'pmm') 
sens                        <- with(data = POAGdataSens, exp = coxph(Surv(POAGcumYear, POAG) ~ indSER + scale(indFunSI_SA) + age + sex + townsend + ethnicBinary + indIOP + indCH + diabetes + hypertension))
tab_model(pool(sens)); for(i in 1:10){print(max(vif(sens$analyses[[i]])))}; sd(complete(POAGdataSens, action = 1)$indFunSI_SA)







