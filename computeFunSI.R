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




#### Display fitted centile curves (CC) ####

# Plot and save
pdf(file.path('figures', 'centileCurves.pdf'), width = 13.5, height = 9) 
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
d <- d %>% rowwise() %>% mutate(FunSI = getFunSI(SER, FD_artery, FD_vein, AVR, Tortuosity_density_artery, Tortuosity_density_vein, conc_rp_artery, conc_rp_vein, DFD_DML_ratio, discTilt, absDiscTors))

# Visualise distribution of FunSI
ggplot(d, aes(x = FunSI))                                            + 
  geom_histogram(bins = 20, alpha = 0.5) + labs(x = 'FunSI', y = 'Count') + 
  scale_x_continuous(breaks = pretty_breaks(15))                          + 
  scale_y_continuous(breaks = pretty_breaks(10))                          + 
  theme_blank()

# Save FunSI for each eye
write.csv(d, file.path('data', 'FunSI.csv'), row.names = F)




d$FunSI     <- NA
for(i in 1:nrow(d)){
  d$FunSI[i]   <- getFunSI(d$SER[i],
                                d$FD_artery[i],
                                d$FD_vein[i],
                                d$AVR[i],
                                d$Tortuosity_density_artery[i],
                                d$Tortuosity_density_vein[i],
                                d$conc_rp_artery[i],
                                d$conc_rp_vein[i],
                                d$DFD_DML_ratio[i],
                                d$discTilt[i],
                                d$absDiscTors[i]) }






