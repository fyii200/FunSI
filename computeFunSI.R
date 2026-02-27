# Author :  Fabian Yii                                                      #
# Email  :  fabian.yii@ed.ac.uk                                             #
# Info   :  R script for computing Fundus Stretch Index (FunSI) for each    #
#           eye based on the fundus centile charts (n=10) created using     #
#           the 'estimateCentileCurves' R script                            #

library(dplyr)
library(ggplot2)
library(scales)
library(sjPlot)
rm(list=rm(list=ls()))
# setwd('..')




#### Set-up ####

# Read the dataset containing data from all eligible eyes
d <- read.csv(file.path('data', 'FunSI.csv'), check.names = F)

# Get the fitted centile curves and helper functions
source(file.path('code', 'utils.R'))




#### Display the fitted centile curves (CC) ####

# Plot and save
pdf(file.path('figures', 'centileCurves_randomOneEye.pdf'), width = 13.5, height = 9) 
layout(matrix(1:10, nrow = 2, byrow = T))
plotCC(m1, 'Arterial fractal dimension')
plotCC(m2, 'Venous fractal dimension')
plotCC(m3, 'Arteriovenous ratio')
plotCC(m4, 'Arterial tortuosity')
plotCC(m5, 'Venous tortuosity')
plotCC(m6, 'Arterial concavity')
plotCC(m7, 'Venous concavity')
plotCC(m8, 'DFD to disc major axis length ratio')
plotCC(m9, 'Disc tilt')
plotCC(m10, 'Absolute disc orientation')
mtext('Spherical equivalent refraction (D)', side = 1, outer = F, line = 3.5, cex = 1, adj = 9)
dev.off()




#### Compute FunSI ####

# Loop through all eyes in the dataset 
d                 <- d %>% rowwise() %>% mutate(FunSI = list(getFunSI(SER, FD_artery, FD_vein, AVR, Tortuosity_density_artery, Tortuosity_density_vein, conc_rp_artery, conc_rp_vein, DFD_DML_ratio, discTilt, absDiscTors))) %>% unnest_wider(FunSI, names_sep = '_')
names(d)[274:285] <- c('FunSI', 'FunSI_SA', 'FunSI_wo_AFD', 'FunSI_wo_VFD', 'FunSI_wo_AVR', 'FunSI_wo_Atort', 'FunSI_wo_Vtort', 'FunSI_wo_Aconc', 'FunSI_wo_Vconc', 'FunSI_wo_DFD_DML_ratio', 'FunSI_wo_Dtilt', 'FunSI_wo_Dtorsion')

# Visualise distribution of FunSI
ggplot(d, aes(x = FunSI))                                                 + 
  geom_histogram(bins = 20, alpha = 0.5) + labs(x = 'FunSI', y = 'Count') + 
  scale_x_continuous(breaks = pretty_breaks(15))                          + 
  scale_y_continuous(breaks = pretty_breaks(10))                          + 
  theme_blank()

# Save FunSI for each eye
write.csv(d, file.path('data', 'FunSI.csv'), row.names = F)


