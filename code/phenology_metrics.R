#Occurrence data phenology Metrics
## Resident butterfly species groups defined by overwinter stage (OWS)
## Phenometrics Onset (5%), Median (50%), Termination (95%), Duration (95%-5%)
## quantile phenometrics with CIs calculated by M Belitz, U Florida
## Species OWS traits compiled by GU Ries Lab
# Modifications by E Larsen, Georgetown U
# Updated 2021-08

library(tidyverse)

pheno.datafile<-"data/derived/adult_bfly_phenometrics_quantileCIs.csv"
pheno.quant<-read_csv(pheno.datafile) %>% rename(cell=HEXcell)


###Calculate durations
pheno.quant<-pheno.quant %>%
  mutate(qdur=q95-q5, qdur_low=q95_low-q5_high, qdur_high=q95_high-q5_low,
         q5_ci=q5_high-q5_low,q50_ci=q50_high-q50_low,
         q95_ci=q95_high-q95_low,qdur_ci=qdur_high-qdur_low)

##Calculate deviations for cells with data 2017-2020
dev.codes<-c("RL","RP")
dev.years<-c(2016:2020)
pheno.dev.cells<-pheno.quant %>%
  filter(code %in% dev.codes & year %in% dev.years) %>%
  group_by(cell, year) %>%
  add_tally(name="ncode") %>%
  filter(ncode==2) %>%
  group_by(cell) %>%
  add_tally(name="nyear")  %>%
  filter(nyear==10) %>%
  summarize(meanq5=mean(q5, na.rm=T),meanq50=mean(q50, na.rm=T),meanq95=mean(q95, na.rm=T),meanqdur=mean(qdur, na.rm=T)) %>%
  select(cell, meanq5, meanq50, meanq95, meanqdur)  
  
pheno.dev<-inner_join(pheno.quant, pheno.dev.cells)  %>%
  mutate(q5_dev=q5-meanq5,q50_dev=q50-meanq50,
         q95_dev=q95-meanq95,qdur_dev=qdur-meanqdur) %>%
  select(cell, year, code, q5_dev, q5_ci, q50_dev,q50_ci, qdur_dev, qdur_ci, 
         cumObsDays,uniqObsDays,totalAbundance,gcode)

save(pheno.quant, pheno.dev, file="data/derived/pheno.RData")



##### END

