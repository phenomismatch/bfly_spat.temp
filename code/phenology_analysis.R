#Analysis of Presence-only Phenology Metrics
#eButterfly, NABA Butterflies I've Seen, iNaturalist (research grade)
## Resident butterfly species groups defined by overwinter stage (OWS)
## Phenometrics estimated from Weibull (M Belitz)
## Species OWS traits compiled by GU Ries Lab
#E Larsen, Georgetown U, Updated 2021-06


#libraries
library(tidyverse)
library(ggplot2)
library(ggcorrplot)
library(lme4)
library(lmerTest)
library(MuMIn)
library(sp)
library(sjPlot)

#Parameters
pheno.datafile<-"data/derived/adult_bfly_metrics.csv"
pheno.dev.file<-"data/derived/adult_bfly_metrics_deviations.csv"
dev.analysis<-TRUE
pheno.dev<-TRUE

#Phenology data
pheno.orig<-read_csv(pheno.datafile) %>% rename(cell=HEXcell)
hexyrs<-pheno.orig %>% dplyr::select(year, cell)  %>% group_by(year, cell) %>% tally()
nrow(pheno.orig)

pheno.dev<-read_csv(pheno.dev.file) %>% rename(cell=HEXcell)
hexyrs<-pheno.dev %>% dplyr::select(year, cell)  %>% group_by(year, cell) %>% tally()
nrow(pheno.dev)

#Environmental data
env.var<-read_csv("data/derived/envir.csv")


#Merge data
if(dev.analysis) {
  pheno.input <- pheno.dev 
  } else {
    if(dev.cells) {
      pheno.input <- filter(pheno.orig, HEXcell %in% unique(pheno.dev$HEXcell))
    } else {
      pheno.input <- pheno.orig
    }
  }

pheno.input.dev<-merge(pheno.dev, env.var, by=c("year","cell")) %>%
  filter(code %in% c("RE","RL","RP")) %>%
  mutate(year=year-2000)

pheno.input<-merge(pheno.orig, env.var, by=c("year","cell")) %>%
  filter(code %in% c("RE","RL","RP")) %>%
  mutate(year=year-2000)


### Analysis
library(MASS)
#Mixed effect model
#ONSET DEVIATION
onset.dev.full<-lmer(dev_10~-1+code*(FFD.dev + warmearly + warmlateopen + year) + (1|cell), data=pheno.input.dev)
r.squaredGLMM(onset.dev.full)     
summary(onset.dev.full)
(onset.dev.step <- step(onset.dev.full))
onset.dev.best <- get_model(onset.dev.step) #stepAIC(ab.yr.full)
summary(onset.dev.best)
r.squaredGLMM(onset.dev.best)     
plot_model(onset.dev.best)

#ONSET
onset.full<-lmer(onset~-1+code*(FFD + warmearly + warmlateopen  + year) + (1|cell), data=pheno.input)
onset.best<-step(onset.full)
r.squaredGLMM(onset.full)     
summary(onset.full)
(onset.step <- step(onset.full))
onset.best <- get_model(onset.step) #stepAIC(ab.yr.full)
vif(onset.best)
summary(onset.best)
r.squaredGLMM(onset.best)
plot_model(onset.best)

##MEDIAN DEVIATION

pheno.input.dev<-na.omit(pheno.input.dev)
med.dev.full<-lmer(dev_50~-1+code*(FFD.dev + spring.dev + summer.dev + gr_mn_for + gr_mn_lag + year)+ (1|cell), data=pheno.input.dev)
r.squaredGLMM(med.dev.full)     
summary(med.dev.full)
med.dev.step <- step(med.dev.full)
med.dev.best <- get_model(med.dev.step) #stepAIC(ab.yr.full)
summary(med.dev.best)
r.squaredGLMM(med.dev.best)     
plot_model(med.dev.best)

#med
med.full<-lmer(fiftieth~code*(FFD + warmearly + warmlateopen + summer.dev  + year) + (1|cell), data=pheno.input)
med.best<-step(med.full)
r.squaredGLMM(med.full)     
summary(med.full)
(med.step <- step(med.full))
med.best <- get_model(med.step) #stepAIC(ab.yr.full)
summary(med.best)
r.squaredGLMM(med.best)     


##DURATION DEVIATION
pheno.input.dev<-na.omit(pheno.input.dev) %>%
  mutate(duration<-term-onset, dur.dev=dev_90-dev_10)
dur.dev.full<-lmer(dur.dev~-1+code*(FFD.dev + warmearly + warmlateopen + summer.dev + year) + (1|cell), data=pheno.input.dev)
r.squaredGLMM(dur.dev.full)     
summary(dur.dev.full)
dur.dev.step <- step(dur.dev.full)
dur.dev.best <- get_model(dur.dev.step) #stepAIC(ab.yr.full)
summary(dur.dev.best)
r.squaredGLMM(dur.dev.best)     


#DURATION
pheno.input<-pheno.input %>%
  mutate(duration=term-onset)
dur.full<-lmer(duration~code*(FFD + warmearly + warmlateopen + summer.dev  + year) + (1|cell), data=pheno.input)
dur.best<-step(dur.full)
r.squaredGLMM(dur.full)     
summary(dur.full)
(dur.step <- step(dur.full))
dur.best <- get_model(dur.step) #stepAIC(ab.yr.full)
summary(dur.best)
r.squaredGLMM(dur.best)     
plot_model(dur.best)
vif(dur.best)
