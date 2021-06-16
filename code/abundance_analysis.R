#Modeling patterns in butterfly abundance
#Elise Larsen, Georgetown U, Updated 2012-05

#libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(r2glmm)
library(sjPlot)
library(MuMIn)
library(MASS)
library(car) #(for VIF function)

#abundance metrics

abund<-read_csv("data/derived/naba_OWS_abundances.csv") %>%
  dplyr::select(cell, year=ObsYear, doy, CountID, SurveyID, ObsMonth, ows.grp=group, abund.bph, log.abund, SR)

#explanatory variables
env<-read_csv("data/derived/envir.csv") %>%
  dplyr::select(cell, year, warmearly, warmlateopen, gr_mn_for, gr_mn_open, gr_mn_lag, spring.dev, summer.dev, FFD.dev, FRD.dev, cell_lat)
pheno<-read_csv("data/derived/simpleton_pheno_pdfs-OutlierDetection.csv") %>%
  dplyr::select(year, cell=HEXcell,ows.grp=code, doy=x, pdf=y) %>%
  group_by(year, cell, ows.grp) %>% 
  mutate(maxpdf=max(pdf, na.rm=T)) %>%
  mutate(cumpr=cumsum(pdf)*.5, on=ifelse(cumpr>0.01,1,0),med=ifelse(cumpr>0.5,1,0),term=ifelse(cumpr>0.99,1,0)) %>%
  mutate(onc=cumsum(on),medc=cumsum(med),termc=cumsum(term),) %>%
  summarize(onset=doy[onc==1], median=doy[medc==1], duration=doy[termc==1]-doy[onc==1], maxdoy=doy[pdf==maxpdf]) %>%
  group_by(cell, ows.grp) %>% 
  mutate(onset.hm=mean(onset, na.rm=T), duration.hm=mean(duration, na.rm=T), onset.d=onset-onset.hm, dur.d=duration-duration.hm)


#previous year abundance
abund.py<-abund %>%
  mutate(year=year+1) %>%
  dplyr::select(cell, year, CountID, ows.grp, abund.py=abund.bph, logab.py=log.abund)


#combine tables
ab<-merge(x = abund, y = pheno, by = intersect(names(abund), names(pheno)), all.x = TRUE)
ab1<-merge(x = ab, y = env, by = intersect(names(ab), names(env)), all.x = TRUE)
ab.final<-merge(x = ab1, y = abund.py, by.x=c("cell", "year", "CountID", "ows.grp"), by.y=c("cell", "year", "CountID", "ows.grp"), all.x = TRUE)

#abundance Model
naba.1<-(ab.final) %>% mutate(ows.grp=as.factor(ows.grp),summer.dev1=summer.dev/100,year=year-2002, on.dev=onset.d/7, dur.dev=dur.d/7, gr_mn_lag=gr_mn_lag/7, FR.dev=FRD.dev/7, daylag=doy-maxdoy, abslag=abs(doy-maxdoy))


#naba.1$ows.grp<-factor(naba.1$ows.grp,levels(factor(naba.1$ows.grp))[c(1,2,3)])

ab.yr.full<-lmer(log.abund~ows.grp+logab.py+warmearly+warmearly:ows.grp+warmlateopen+warmlateopen:ows.grp+on.dev+abslag+abslag:ows.grp+dur.dev+year+year:ows.grp+gr_mn_lag+as.factor(ObsMonth)+doy+FR.dev+FR.dev:ows.grp+(1|cell) + (1|CountID:cell), data=naba.1)
ab.yr.full<-lmer(log.abund~-1+ows.grp*(logab.py+warmearly+warmlateopen+on.dev+abslag+dur.dev+year+gr_mn_lag+as.factor(ObsMonth)+doy+FR.dev)+(1|cell) + (1|CountID:cell), data=naba.1)
r.squaredGLMM(ab.yr.full)     
summary(ab.yr.full)
(step_resyr <- step(ab.yr.full))
finalyr <- get_model(step_resyr) #stepAIC(ab.yr.full)
anova(finalyr)
summary(finalyr)
r.squaredGLMM(finalyr)
AIC(finalyr)
vif(finalyr)

ab.yr.full<-lmer(log.abund~ows.grp*(logab.py+warmearly+warmlateopen+on.dev+abslag+dur.dev+year+gr_mn_lag+as.factor(ObsMonth)+doy+FR.dev)+(1|cell) + (1|CountID:cell), data=naba.1)

## how to interpret VIFs for interactions
plot_model(finalyr)
finalyr<-lmer(log.abund~-1+ows.grp+logab.py+warmearly:ows.grp+warmlateopen:ows.grp+on.dev:ows.grp+abslag+year:ows.grp+as.factor(ObsMonth):ows.grp+doy:ows.grp+(1|cell/CountID), data=naba.1)
(step_resyr <- step(finalyr))
finalyr <- get_model(step_resyr) #stepAIC(ab.yr.full)
anova(finalyr)
ranova(finalyr)
summary(finalyr)
r.squaredGLMM(finalyr)
AIC(finalyr)
vif(finalyr)
plot_model(finalyr)

fixedyr<-lm(log.abund~-1+ows.grp+logab.py+abslag+warmearly:ows.grp+warmlateopen:ows.grp+on.dev:ows.grp+as.factor(ObsMonth):ows.grp+doy:ows.grp, data=naba.1)
phtest_glmer(finalyr, fixedyr)
plotresult<-plot_model(finalyr, values=T)

fixedlabels<-c(rep("",20))
fixedlabels[which(summary(finalyr)$coefficients[,5]<0.05)]<-"*"
plotresult + annotate(geom="text", x=rev(c(1:20)), y=rep(-3,20), label=fixedlabels)

ggplot(data=naba.1, aes(x=ows.grp, y=log.abund)) + geom_boxplot()
#ggplot(data=naba.1, aes(color=ows.grp, fill=ows.grp, x=logab.py, y=log.abund)) + geom_point()
ggplot(data=naba.1, aes(color=as.factor(CountID), fill=CountID, shape=ows.grp, x=logab.py, y=log.abund)) + geom_point() +theme(legend.position="none")
#ggplot(data=naba.1, aes( shape=ows.grp, x=logab.py, y=log.abund)) + geom_point()
ggplot(data=naba.1, aes(color=ows.grp, fill=ows.grp, x=warmearly, y=log.abund)) + geom_point() + geom_smooth(method="lm")
ggplot(data=naba.1, aes(color=ows.grp, fill=ows.grp, x=warmlateopen, y=log.abund)) + geom_point() + geom_smooth(method="lm")
ggplot(data=naba.1, aes(color=ows.grp, fill=ows.grp, x=warmlateopen, y=log.abund)) + geom_point() + geom_smooth(method="lm")
ggplot(data=naba.1, aes(color=ows.grp, fill=ows.grp, x=on.dev, y=log.abund)) + geom_point() + geom_smooth(method="lm")
ggplot(data=naba.1, aes(color=ows.grp, fill=ows.grp, x=doy, y=log.abund)) + geom_point() + geom_smooth(method="lm") + facet_wrap(.~ObsMonth)

ggplot(data=naba.1, aes(x=doy, y=abslag)) + geom_point() + geom_smooth()
ggplot(data=naba.1, aes(x=abslag, y=log.abund)) + geom_point() + geom_smooth()

fixed_only<-lm(log.abund~ows.grp+logab.py+on.dev+as.factor(ObsMonth)+abslag+daylag+warmearly+ows.grp:FR.dev, data=naba.1)
AIC(fixed_only)

#Is there a significant temporal trend in warmearly?
env1<-naba.1 %>% dplyr::select(cell, year, warmearly:abslag) %>% group_by(cell, year, warmearly) %>% tally()
summary(lm(warmearly~year+cell, data=env1))
#spatial trend, not temporal trend


fixed.full<-lm(log.abund~-1+ows.grp+logab.py+warmearly:ows.grp+warmlateopen+on.dev+dur.dev+year:ows.grp+gr_mn_lag+as.factor(ObsMonth)+doy+FR.dev:ows.grp, data=naba.1)
step_fix <- stepAIC(fixed.full, direction = "both", 
                    trace = FALSE)
summary(step_fix)


summary(fixed.full)     
step_fixyr<-step(fixed.full)
fixedyr <- stepAIC(step_fixyr)
anova(fixedyr)
summary(fixedyr)
r.squaredGLMM(fixedyr)
AIC(fixedyr)



phtest_glmer <- function (glmerMod, glmMod, ...)  {  ## changed function call
  coef.wi <- coef(glmMod)
  coef.re <- fixef(glmerMod)  ## changed coef() to fixef() for glmer
  vcov.wi <- vcov(glmMod)
  vcov.re <- vcov(glmerMod)
  names.wi <- names(coef.wi)
  names.re <- names(coef.re)
  coef.h <- names.re[names.re %in% names.wi]
  dbeta <- coef.wi[coef.h] - coef.re[coef.h]
  df <- length(dbeta)
  dvcov <- vcov.re[coef.h, coef.h] - vcov.wi[coef.h, coef.h]
  stat <- abs(t(dbeta) %*% as.matrix(solve(dvcov)) %*% dbeta)  ## added as.matrix()
  pval <- pchisq(stat, df = df, lower.tail = FALSE)
  names(stat) <- "chisq"
  parameter <- df
  names(parameter) <- "df"
  alternative <- ifelse(pval<0.05,"choose fixed effects model","choose mixed effects model")
  res <- list(statistic = stat, p.value = pval, parameter = parameter, 
              method = "Hausman Test",  alternative = alternative,
              data.name=deparse(getCall(glmerMod)$data))  ## changed
  class(res) <- "htest"
  return(res)
}

#plots
plot_model(finalyr, values=T)

phtest_glmer(finalyr, fixedyr)
plotresult<-plot_model(fixedyr, values=T)

fixedlabels<-c(rep("",15))
fixedlabels[which(summary(fixedyr)$coefficients[,4]<0.05)]<-"*"
plotresult + annotate(geom="text", x=rev(c(1:15)), y=rep(-1.5,15), label=fixedlabels)

summary(naba.1)
#####
## Fixed effects model is better

hist(naba.1$year)
hist(naba.1$CountID)



#Plots
ggplot(data=naba.1) + geom_boxplot(aes(factor(ObsMonth),log.abund, fill=ows.grp)) + 
  scale_x_discrete(labels = c('June','July')) + theme_classic()


ggplot(data=naba.1,aes(x=logab.py, y=log.abund, color=ows.grp)) + geom_point() +
  geom_smooth(method="lm") + theme_classic()

ggplot(data=naba.1,aes(x=on.dev, y=log.abund, color=ows.grp)) + geom_point() +
  geom_smooth(method="lm") + theme_classic()

ggplot(data=naba.1,aes(x=FR.dev, y=log.abund, color=ows.grp)) + geom_point() +
  geom_smooth(method="lm") + theme_classic()

ggplot(data=naba.1,aes(x=year, y=log.abund, color=ows.grp)) + geom_point() +
  geom_smooth(method="lm") + theme_classic() + theme(legend.position="none")

ggplot(data=naba.1,aes(x=warmearly, y=log.abund, color=ows.grp)) + geom_point() +
  geom_smooth(method="lm") + theme_classic() 

ggplot(data=env,aes(x=year, y=spring.dev, color=(cell_lat))) + geom_point() +
  geom_smooth(method="lm") + theme_classic() + theme(legend.position="none")
