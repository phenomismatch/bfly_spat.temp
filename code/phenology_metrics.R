#Occurrence data phenology Metrics
## Resident butterfly species groups defined by overwinter stage (OWS)
## Phenometrics Onset (5%), Median (50%), Termination (95%), Duration (95%-5%)
## quantile phenometrics with CIs calculated by M Belitz, U Florida
## Species OWS traits compiled by GU Ries Lab
# Modifications by E Larsen, Georgetown U
# Updated 2021-08

library(tidyverse)

load("data/spatial.domain.RData")


pheno.datafile<-"data/derived/adult_bfly_phenometrics_noCountCircles_withFull2020Data.csv"
pheno.quant<-read_csv(pheno.datafile) %>% rename(cell=HEXcell) %>%
  filter(cell %in% STUDYCELLS, between(q50,152,243), !is.na(q5))



###Calculate durations
pheno.quant<-pheno.quant %>%
  mutate(qdur=q95-q5, qdur_low=q95_low-q5_high, qdur_high=q95_high-q5_low,
         q5_ci=q5_high-q5_low,q50_ci=q50_high-q50_low,
         q95_ci=q95_high-q95_low,qdur_ci=qdur_high-qdur_low)



##Calculate deviations for cells with data 2016-2020
dev.years<-c(2016:2020)
pheno.cell.baseline<-pheno.quant %>%
  filter(year %in% dev.years) %>%
  group_by(cell,code) %>%
  add_tally(name="nyear")  %>%
  filter(nyear==5) %>%
  group_by(cell,code) %>%
  summarize(meanq5=mean(q5, na.rm=T),meanq50=mean(q50, na.rm=T),meanq95=mean(q95, na.rm=T),meanqdur=mean(qdur, na.rm=T)) %>%
  dplyr::select(cell, code, meanq5, meanq50, meanq95, meanqdur)  
  
pheno.dev<-inner_join(pheno.quant, pheno.cell.baseline)  %>%
  mutate(q5_dev=q5-meanq5,q50_dev=q50-meanq50,
         q95_dev=q95-meanq95,qdur_dev=qdur-meanqdur) %>%
  dplyr::select(cell, year, code, onset.dev=q5_dev, onset.ci=q5_ci, median.dev=q50_dev,median.ci=q50_ci, 
         dur.dev=qdur_dev, dur.ci=qdur_ci, cumObsDays,uniqObsDays,totalAbundance)


save(pheno.quant, pheno.dev, file="data/derived/pheno.RData")

#x<-read_csv("data/envir/cpcGDD_hex_2000-2020.csv")


##### END

