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

# Ensure all date variables are correctly coded as datetime
cvtDate              <- function(x) format(as.Date(x, '%d/%m/%Y'))
d                    <- d %>% mutate(visitDate                        = cvtDate(visitDate),
                                     catSurgery1                      = cvtDate(catSurgery1),
                                     dateDeath                        = cvtDate(dateDeath),
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
d$rhegmaRD           <- F
for(i in IDsWithRDFirstOcc){
  # Check if the given ID has a hospital record of rhegmatogenous RD (H33.0)
  d[d$id == i,]$rhegmaRD <- 'H33.0' %in% substring(subset(diag, id == i)$ICD10, 1, 5)
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

# Remove unreasonable intraocular pressure (IOP) and corneal hysteresis (CH) values
# Note 1: max IOP in glaucomatous eyes reported by a population-based study was 45.6mmHg (https://www.bmj.com/content/358/bmj.j3889)
# Note 2: max CH based on https://karger.com/ore/article-pdf/46/4/187/3367017/000326896.pdf and https://www.sciencedirect.com/science/article/pii/S0002939411008531?via%3Dihub#sec2
d                    <- d %>% mutate(IOP = ifelse(IOP == 0 | IOP > 45, NA, IOP),
                                     CH  = ifelse(CH == 0 | CH > 15, NA, CH))

# Average SER, fundus-predicted SER, FunSI, IOP and CH across eyes if data from both eyes are available
d                    <- d %>% group_by(id) %>% mutate(indSER         = ifelse(sum(!is.na(SER)) == 2, (SER[1] + SER[2]) / 2, SER),
                                                      indFunSI       = ifelse(sum(!is.na(FunSI)) == 2, (FunSI[1] + FunSI[2]) / 2, FunSI),
                                                      indPredSER_TTA = ifelse(sum(!is.na(predSER_TTA)) == 2, (predSER_TTA[1] + predSER_TTA[2]) / 2, predSER_TTA),
                                                      indIOP         = ifelse(sum(!is.na(IOP)) == 2, (IOP[1] + IOP[2]) / 2, IOP),
                                                      indCH          = ifelse(sum(!is.na(CH)) == 2, (CH[1] + CH[2]) / 2, CH))




#### Compute survival time ####                       

# Right-censoring (RC) date
RCdate                             <- as.Date('2022-05-31') 
d$firstDateGlaucoma                <- pmin(as.Date(d$firstDateGlaucoma), as.Date(RCdate))
d$firstDateRetinalDetachmentBreaks <- pmin(as.Date(d$firstDateRetinalDetachmentBreaks), as.Date(RCdate))

# If baseline assessment date is missing, infer by summing year of birth and age at the baseline visit
missingVisitDate                   <- which(is.na(d$visitDate))
if(length(missingVisitDate) > 0){d[missingVisitDate, ]$visitDate <- format(as.Date(paste0(d$YOB[missingVisitDate] + d$age[missingVisitDate], '-01-01')), '%Y-%m-%d') }

# Create new columns in preparation for survival analysis 
# 'RD'               : indicates if a participant had a record of RD of any subtype during follow-up
# 'rhegmaRD'         : indicates if a participant had a record of rhegmatogenous RD during follow-up
# 'glaucoma'         : indicates if a participant had a record of glaucoma of any subtype during follow-up
# 'POAG'             : indicates if a participant had a record of POAG during follow-up
# 'xx' + 'cumYear'   : time from baseline to the first recorded occurrence of the 'xx' event (disease) or to the right-censoring date 
# 'timeToCataractOp' : time from baseline to the first recorded occurrence of cataract surgery
# 'cataractOpDuring' : indicates if a participant had undergone cataract surgery during follow-up and prior to the first-recorded RD or the right-censoring date
# 'ageTrauma'        : age when a participant had a self-reported trauma/injury 'resulting in a loss of vision'
# 'previousTrauma'   : True if 'ageTrauma' is present and smaller than baseline age
d                                 <- d %>% mutate(observedPeriod    = ifelse(!is.na(dateDeath), difftime(dateDeath, visitDate)/365.25,
                                                                             ifelse(!is.na(dateLossFU), difftime(dateLossFU, visitDate)/365.25, difftime(RCdate, visitDate)/365.25)),
                                                  RD                = ifelse(!is.na(firstDateRetinalDetachmentBreaks), T, F),
                                                  RDcumYear         = ifelse(RD == T, difftime(firstDateRetinalDetachmentBreaks, visitDate)/365.25, observedPeriod),
                                                  rhegmaRDcumYear   = ifelse(rhegmaRD == T, difftime(firstDateRetinalDetachmentBreaks, visitDate)/365.25, observedPeriod),
                                                  glaucoma          = ifelse(!is.na(firstDateGlaucoma), T, F),
                                                  glaucomaCumYear   = ifelse(glaucoma == T, difftime(firstDateGlaucoma, visitDate)/365.25, observedPeriod),
                                                  POAG              = ifelse(!is.na(firstDateGlaucoma) & !nonPOAG, T, F),
                                                  POAGcumYear       = ifelse(POAG == T, difftime(firstDateGlaucoma, visitDate)/365.25, observedPeriod),
                                                  timeToCataractOp  = difftime(catSurgery1, visitDate, units = 'days')/365.25,
                                                  cataractOpDuring  = ifelse(timeToCataractOp < RDcumYear & !is.na(catSurgery1), T, F),
                                                  ageTrauma         = ifelse(ageVisionLossDueToInjuryOrTrauma != '' & !is.na(ageVisionLossDueToInjuryOrTrauma) & ageVisionLossDueToInjuryOrTrauma != 'Do not know', as.numeric(ageVisionLossDueToInjuryOrTrauma), NA),
                                                  previousTrauma    = ifelse(!is.na(ageTrauma) & ageTrauma <= age, T, F))




#### Survival analysis for rhegmatogenous RD ####

# Population at risk (25067 IDs): IDs without any history of RD/breaks at baseline
RDdata                    <- subset(d, !duplicated(id) & RDcumYear > 0)
nrow(RDdata); length(unique(RDdata$id))

# Create a new column storing categorised FunSI (for visualisation purposes; Kaplan-Meier survival plot)
RDdata$FunSIgroup         <- cut(RDdata$indFunSI, quantile(RDdata$indFunSI, seq(0,1,0.25)), include.lowest = T)
levels(RDdata$FunSIgroup) <- c('0.19-0.44', '0.44-0.50', '0.50-0.57', '0.57-0.88')

# Summary stats for time to rhegmatogenous RD
summary(subset(RDdata, rhegmaRD == T)$RDcumYear)

# Number of participants with newly onset rhegmatogenous RD (n=164) (n=176)
length(unique(RDdata[RDdata$rhegmaRD == T, ]$id))

# 1696 out of 25067 participants with a history of cataract surgery during follow-up
length(unique(subset(RDdata, cataractOpDuring == T)$id))

# 124 participants with a self-reported history of ocular trauma/injury
length(unique(subset(RDdata, previousTrauma == T)$id))

# Plot Kaplan-Meier survival curve stratified by baseline FRO (rhegmatogenous RD/breaks as event)
rhegRD_KMfit       <- survfit(Surv(rhegmaRDcumYear, rhegmaRD) ~ FunSIgroup, RDdata)
rhegRD_KMplot      <- ggsurvplot(rhegRD_KMfit,
                                 ylim         = c(0.985, 1),
                                 xlim         = c(0, 12),
                                 break.x.by   = 4,                             
                                 size         = 0.6,
                                 add.all      = T,
                                 censor       = F,
                                 xlab         = 'Years elapsed',
                                 ylab         = 'Survival probability\n',
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
rhegRD_KMplot$plot  <- rhegRD_KMplot$plot + guides(col = guide_legend(override.aes = list(linewidth = 2)))
pdf(file.path('figures', 'rhegmaRD_KMplot.pdf'), width = 6.5, height = 7)
print(grid.arrange(rhegRD_KMplot$plot, rhegRD_KMplot$table, ncol = 1, heights = c(2, 0.5)))
dev.off()


## Main analysis 

# Univariable Cox regression
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indFunSI, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ age, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ sex, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ townsend, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ ethnicBinary, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ cataractOpDuring, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ diabetes, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ hypertension, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ BMI, data = RDdata))
tab_model(coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ previousTrauma, data = RDdata))

# Multivariable Cox regression
coxRhegmaRD           <- coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + indFunSI + age + sex + cataractOpDuring, data = RDdata)
vif(coxRhegmaRD); tab_model(coxRhegmaRD)       # unstandardised FunSI
RDdata$indFunSIscaled <- scale(RDdata$indFunSI)
coxRhegmaRDscaled     <- coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + indFunSIscaled + age + sex + cataractOpDuring, data = RDdata)
tab_model(coxRhegmaRDscaled)                   # standardised FunSI (zero mean and unit variance)

# Check concordance (higher means better discriminative power)
fit1                  <- coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ age + sex + cataractOpDuring, data = RDdata)
fit2                  <- update(fit1, .~. + indSER, RDdata)
fit3                  <- update(fit2, .~. + indFunSI, RDdata)
concordance(fit1, fit2, fit3)
AIC(fit2, fit3)

# Diagnostic
cox.zph(coxRhegmaRD)


## Subgroup analysis: Cox regression with FRO as included another independent variable 
## in a subset of IDs not previously used to train the FRO model (15979 IDs)

# Compute FRO
RDdataSub                <- subset(RDdata, !is.na(indPredSER_TTA)); nrow(RDdataSub)
FERmodel                 <- lm(indPredSER_TTA ~ indSER, RDdataSub)
RDdataSub$indFRO         <- residuals(FERmodel)

# Fit Cox regression
coxRhegmaRDsub           <- coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + indFunSI + age + sex + cataractOpDuring + indFRO, data = RDdataSub)
vif(coxRhegmaRDsub); tab_model(coxRhegmaRDsub)        # unstandardised FunSI and FRO
RDdataSub$indFunSIscaled <- scale(RDdataSub$indFunSI); scale(RDdataSub$indFunSI)
RDdataSub$indFROscaled   <- scale(RDdataSub$indFRO); scale(RDdataSub$indFRO)
coxRhegmaRDsubScaled     <- coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + indFunSIscaled + age + sex + cataractOpDuring + indFROscaled, data = RDdataSub)
tab_model(coxRhegmaRDsubScaled)                       # standardised FunSI and FRO (zero mean and unit variance)


## Sensitivity analysis 1: only include eyes with cyl ≤ 2D (23,077 IDs)
RDdataSens1              <- subset(RDdata, cyl <= 2); nrow(RDdataSens1)
coxRhegmaRDSens1         <- coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + indFunSI + age + sex + cataractOpDuring, data = RDdataSens1)
vif(coxRhegmaRDSens1); tab_model(coxRhegmaRDSens1) 


## Sensitivity analysis 2: remove rhegmatogenous RD cases that developed within 1 year of baseline (rule out data processing lag)
# 8 rhegmatogenous RD cases removed
excludeIDs                  <- subset(RDdata, rhegmaRD & (rhegmaRDcumYear < 1))$id; length(excludeIDs)
RDdataSens2                 <- subset(RDdata, ! id %in% excludeIDs)
coxRhegmaRDSens2            <- coxph(Surv(rhegmaRDcumYear, rhegmaRD) ~ indSER + indFunSI + age + sex + cataractOpDuring, data = RDdataSens2)
vif(coxRhegmaRDSens2); tab_model(coxRhegmaRDSens2) 




#### Survival analysis for POAG ####

# Population at risk (24835 IDs): IDs without any history of glaucoma at baseline
POAGdata                    <- subset(d, !duplicated(id) & glaucomaCumYear > 0)
nrow(POAGdata); length(unique(POAGdata$id))

# 460 out of 24835 IDs had incidence POAG
table(POAGdata$POAG)

# Create a new column storing categorised SI (for visualisation purposes; Kaplan-Meier survival plot)
POAGdata$FunSIgroup         <- cut(POAGdata$indFunSI, quantile(POAGdata$indFunSI, seq(0,1,0.25)), include.lowest = T)
levels(POAGdata$FunSIgroup) <- c('0.19-0.44', '0.44-0.50', '0.50-0.57', '0.57-0.88')

# Summary stats for time to RD/breaks
quantile(subset(POAGdata, POAG == T)$POAGcumYear)

# Linkage sources 
# Of the 460 newly onset cases, majority (429) sourced from hospital or primary care (only 23 exclusively self-reported cases)
POAGdata                    <- POAGdata %>% mutate(sourceGlaucoma = ifelse(is.na(sourceGlaucoma), '', sourceGlaucoma))
table(subset(POAGdata, !duplicated(id) & POAG == T)$sourceGlaucoma)

# Plot Kaplan-Meier survival curve stratified by baseline SI (any RD/breaks as event)
allPOAG_KMfit               <- survfit(Surv(POAGcumYear, POAG) ~ FunSIgroup, POAGdata)
allPOAG_KMplot              <- ggsurvplot(allPOAG_KMfit,
                                          ylim         = c(0.965, 1),
                                          xlim         = c(0, 12),
                                          break.x.by   = 4,
                                          size         = 0.6,
                                          add.all      = T,
                                          censor       = F,
                                          xlab         = 'Years elapsed',
                                          ylab         = 'Survival probability\n',
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


## Main analysis

# Univariable Cox regression
tab_model(coxph(Surv(POAGcumYear, POAG) ~ indSER, data = POAGdata))
tab_model(coxph(Surv(POAGcumYear, POAG) ~ indFunSI, data = POAGdata))
tab_model(coxph(Surv(POAGcumYear, POAG) ~ IOP, data = POAGdata))
tab_model(coxph(Surv(POAGcumYear, POAG) ~ CH, data = POAGdata))
tab_model(coxph(Surv(POAGcumYear, POAG) ~ age, data = POAGdata))
tab_model(coxph(Surv(POAGcumYear, POAG) ~ sex, data = POAGdata))
tab_model(coxph(Surv(POAGcumYear, POAG) ~ townsend, data = POAGdata))
tab_model(coxph(Surv(POAGcumYear, POAG) ~ ethnicBinary, data = POAGdata))
tab_model(coxph(Surv(POAGcumYear, POAG) ~ diabetes, data = POAGdata))
tab_model(coxph(Surv(POAGcumYear, POAG) ~ hypertension, data = POAGdata))
tab_model(coxph(Surv(POAGcumYear, POAG) ~ BMI, data = POAGdata))


# Multivariable Cox regression
coxPOAG                     <- coxph(Surv(POAGcumYear, POAG) ~ indSER + indFunSI + age + sex + indIOP + indCH + diabetes + hypertension, data = POAGdata)
vif(coxPOAG); tab_model(coxPOAG)   # Unstandardised FunSI
POAGdata$indFunSIscaled     <- scale(POAGdata$indFunSI) 
coxPOAGscaled               <- coxph(Surv(POAGcumYear, POAG) ~ indSER + indFunSIscaled + age + sex + indIOP + indCH + diabetes + hypertension, data = POAGdata)
tab_model(coxPOAG)   # Standardised FunSI (zero mean and unit variance)

# Check concordance (higher means better discriminative power)
fit1                        <- coxph(Surv(POAGcumYear, POAG) ~ age + sex + IOP + diabetes + hypertension, data = POAGdata)
fit2                        <- update(fit1, .~. + indSER, POAGdata)
fit3                        <- update(fit2, .~. + indFunSI, POAGdata)
concordance(fit1, fit2, fit3)
AIC(fit2, fit3)

# Diagnostic plots
cox.zph(coxPOAG)


## Subgroup analysis: Cox regression with FRO as included another independent variable 
## in a subset of IDs not previously used to train the FRO model (15754 IDs)

# Compute FRO 
POAGdataSub                 <- subset(POAGdata, !is.na(indPredSER_TTA)); nrow(POAGdataSub)
FERmodel                    <- lm(indPredSER_TTA ~ indSER, POAGdataSub)
POAGdataSub$indFRO          <- residuals(FERmodel)

# Fit Cox regression
coxPOAGsub                  <- coxph(Surv(POAGcumYear, POAG) ~ indSER + indFunSI + age + sex + indIOP + indCH + diabetes + hypertension + indFRO, data = POAGdataSub)
vif(coxPOAGsub); tab_model(coxPOAGsub)        # Unstandardised FunSI and FRO
POAGdataSub$indFunSIscaled  <- scale(POAGdataSub$indFunSI)
POAGdataSub$indFROscaled    <- scale(POAGdataSub$indFRO)
coxPOAGsubScaled            <- coxph(Surv(POAGcumYear, POAG) ~ indSER + indFunSIscaled + age + sex + indIOP + indCH + diabetes + hypertension + indFROscaled, data = POAGdataSub)
tab_model(coxPOAGsubScaled)                   # Standardised FunSI and FRO (zero mean and unit variance)


## Sensitivity analysis 1: remove the 23 self-reported cases
POAGdataSens1               <- subset(POAGdata, sourceGlaucoma != 'Self-report only'); nrow(POAGdataSens1)
coxPOAGsens1                <- coxph(Surv(POAGcumYear, POAG) ~ indSER + indFunSI + age + sex + indIOP + indCH + diabetes + hypertension, data = POAGdataSens1)
vif(coxPOAGsens1); tab_model(coxPOAGsens1) 


## Sensitivity analysis 2: only include eyes with cyl ≤ 2D
POAGdataSens2               <- subset(POAGdata, cyl <= 2); nrow(POAGdataSens2)
coxPOAGsens2                <- coxph(Surv(POAGcumYear, POAG) ~ indSER + indFunSI + age + sex + indIOP + indCH + diabetes + hypertension, data = POAGdataSens2)
vif(coxPOAGsens2); tab_model(coxPOAGsens2) 


## Sensitivity analysis 3: remove POAG cases that developed within 1 year of baseline (rule out data processing lag)
# 23 POAG cases removed
excludeIDs                  <- subset(POAGdata, glaucoma & (POAGcumYear < 1))$id; length(excludeIDs)
POAGdataSens3               <- subset(POAGdata, ! id %in% excludeIDs)
coxPOAGsens3                <- coxph(Surv(POAGcumYear, POAG) ~ indSER + indFunSI + age + sex + indIOP + indCH + diabetes + hypertension, data = POAGdataSens3)
vif(coxPOAGsens3); tab_model(coxPOAGsens3) 







