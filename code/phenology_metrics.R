#Occurrence data phenology Metrics
## Resident butterfly species groups defined by overwinter stage (OWS)
## Phenometrics Onset (5%), Median (50%), Termination (95%), Duration (95%-5%)
## quantile phenometrics with CIs calculated by M Belitz, U Florida
## Species OWS traits compiled by GU Ries Lab
# Modifications by E Larsen, Georgetown U
# Updated 2023-10

library(tidyverse)
save.output<-F #set to True to write output files

load("data/spatial.domain.RData")

dev.years<-c(2016:2019)

pheno.datafile<-"data/derived/adult_bfly_phenometrics_noCountCircles_withFull2020Data.csv"
pheno.quant<-read_csv(pheno.datafile) %>% rename(cell=HEXcell) %>%
  filter(cell %in% STUDYCELLS, between(q50,152,243), !is.na(q5)) #, between(q50,152,243), !is.na(q5))


###Calculate durations
pheno.quant<-pheno.quant %>%
filter(year<2020) %>%
  mutate(qdur=q95-q5, qdur_low=q95_low-q5_high, qdur_high=q95_high-q5_low,
         q5_ci=q5_high-q5_low,q50_ci=q50_high-q50_low,
         q95_ci=q95_high-q95_low,qdur_ci=qdur_high-qdur_low)


##Calculate deviations for cells with data 2016-2019
pheno.cell.baseline<-pheno.quant %>%
  filter(year %in% dev.years) %>%
  group_by(cell,code) %>%
  add_tally(name="nyear")  %>%
  filter(nyear==4) %>%
  group_by(cell,code) %>%
  summarize(meanq5=mean(q5, na.rm=T),meanq50=mean(q50, na.rm=T),meanq95=mean(q95, na.rm=T),meanqdur=mean(qdur, na.rm=T)) %>%
  dplyr::select(cell, code, meanq5, meanq50, meanq95, meanqdur)  
  
pheno.dev<-left_join(pheno.quant, pheno.cell.baseline, by=c("cell", "code"))  %>%
  mutate(q5_dev=q5-meanq5,q50_dev=q50-meanq50,
         q95_dev=q95-meanq95,qdur_dev=qdur-meanqdur) %>%
  dplyr::select(cell, year, code, q50, onset.dev=q5_dev, onset.ci=q5_ci, median.dev=q50_dev,median.ci=q50_ci, 
         dur.dev=qdur_dev, dur.ci=qdur_ci, cumObsDays,uniqObsDays,totalAbundance)

if(save.output) {
  save(pheno.dev, file="data/derived/phenoDev.RData")
  save(pheno.quant, file="data/derived/phenoQ.RData")
}

##### END