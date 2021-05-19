#Modeling patterns in butterfly abundance
#Elise Larsen, Georgetown U, Updated 2012-05

#libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(r2glmm)
library(sjPlot)
#abundance metrics

abund<-read_csv("data/derived/naba_OWS_abundances.csv") %>%
  select(cell, year=ObsYear, doy, CountID, SurveyID, ObsMonth, ows.grp=group, abund.bph, log.abund, SR)

#explanatory variables
env<-read_csv("data/derived/envir.csv") %>%
  select(cell, year, warmearly, warmlateopen, gr_mn_for, gr_mn_open, gr_mn_lag, spring.dev, summer.dev, FFD.dev, FRD.dev)
pheno<-read_csv("data/derived/simpleton_pheno_pdfs-OutlierDetection.csv") %>%
  select(year, cell=HEXcell,ows.grp=code, doy=x, pdf=y) %>%
  group_by(year, cell, ows.grp) %>% 
  mutate(cumpr=cumsum(pdf)*.5, on=ifelse(cumpr>0.01,1,0),med=ifelse(cumpr>0.5,1,0),term=ifelse(cumpr>0.99,1,0)) %>%
  mutate(onc=cumsum(on),medc=cumsum(med),termc=cumsum(term)) %>%
  summarize(onset=doy[onc==1], median=doy[medc==1], duration=doy[termc==1]-doy[onc==1]) %>%
  group_by(cell, ows.grp) %>% 
  mutate(onset.hm=mean(onset, na.rm=T), duration.hm=mean(duration, na.rm=T), onset.d=onset-onset.hm, dur.d=duration-duration.hm)
  

#previous year abundance
abund.py<-abund %>%
  mutate(year=year+1) %>%
  select(cell, year, CountID, ows.grp, abund.py=abund.bph, logab.py=log.abund)


#combine tables
ab<-merge(x = abund, y = pheno, by = intersect(names(abund), names(pheno)), all.x = TRUE)
ab1<-merge(x = ab, y = env, by = intersect(names(ab), names(env)), all.x = TRUE)
ab.final<-merge(x = ab1, y = abund.py, by.x=c("cell", "year", "CountID", "ows.grp"), by.y=c("cell", "year", "CountID", "ows.grp"), all.x = TRUE)

#abundance Model
naba.1<-na.omit(ab.final) %>% mutate(ows.grp=as.factor(ows.grp),summer.dev1=summer.dev/100,year=year-2002, doy=doy-150, on.dev=onset.d/7, dur.dev=dur.d/7, gr_mn_lag=gr_mn_lag/7, FR.dev=FRD.dev/7)
#naba.1$ows.grp<-factor(naba.1$ows.grp,levels(factor(naba.1$ows.grp))[c(1,2,3)])
ab.yr.full<-lmer(log.abund~0+ows.grp+abund.py+warmearly:ows.grp+warmlateopen+on.dev+dur.dev+year:ows.grp+gr_mn_lag+as.factor(ObsMonth)+doy+FR.dev:ows.grp+(1|cell)+(1|CountID), data=naba.1)
r.squaredGLMM(ab.yr.full)     
summary(ab.yr.full)
(step_resyr <- step(ab.yr.full))
finalyr <- get_model(step_resyr)
anova(finalyr)
summary(finalyr)
r.squaredGLMM(finalyr)
plot_model(finalyr, values=T)


