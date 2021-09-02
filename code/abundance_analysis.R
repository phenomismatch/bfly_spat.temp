#Modeling patterns in butterfly abundance
#Elise Larsen, Georgetown U, Updated 2012-08

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
theme_set(theme_sjplot())
###PHENO DATA
load("data/derived/pheno.RData")
pheno.dev<-left_join(pheno.dev, pheno.quant)

abund<-read_csv("data/derived/naba_OWS_abundances.csv") %>%
  dplyr::select(cell, year=ObsYear, doy, CountID, SurveyID, ObsMonth, code=group, abund.bph, log.abund, SR) %>%
  filter(ObsMonth %in% c(6:8), cell %in% STUDYCELLS)


#explanatory variables
env<-read_csv("data/derived/envir.csv") %>%
  dplyr::select(cell, year, warmearly, warmlateopen, gr_mn_for, gr_mn_open, gr_mn_lag, spring.dev, summer.dev, FFD.dev, FRD.dev, cell_lat)

#previous year abundance
abund.py<-abund %>%
  mutate(year=year+1) %>%
  dplyr::select(cell, year, CountID, code, abund.py=abund.bph, logab.py=log.abund)


#combine tables
ab<-merge(x = abund, y = pheno.dev, by = intersect(names(abund), names(pheno.dev)), all.x = TRUE)
ab1<-merge(x = ab, y = env, by = intersect(names(ab), names(env)), all.x = TRUE)
ab.final<-merge(x = ab1, y = abund.py, by.x=c("cell", "year", "CountID", "code"), by.y=c("cell", "year", "CountID", "code"), all.x = TRUE)
#abundance Model
naba.1<-(ab.final) %>% mutate(code=as.factor(code),summer.dev1=summer.dev/100,year=year-2000, on.dev=onset.dev/7, dur.dev=dur.dev/7, gr_mn_lag=gr_mn_lag/7, FR.dev=FRD.dev/7, daylag=doy-q50, abslag=abs(doy-q50))
naba.1<-naba.1 %>% mutate(MonthF=as.factor(ObsMonth), ObsDay=as.numeric(format(as.Date(doy, origin = paste((year+1999),"12-31",sep="-")),format="%d")))

save(naba.1, file="data/abund.input.RData")

summary(naba.1)
includeNA<-T
if(includeNA==T) {
  ab.yr.full<-lmer(log.abund~-1+code*(warmearly+logab.py+warmlateopen+on.dev+abslag+dur.dev+MonthF*ObsDay+year+FR.dev)+(1|cell) + (1|CountID:cell), data=naba.1)

  extractAIC(ab.yr.full)
  ab.yr.1<-lmer(log.abund~-1+code+code:warmearly+code:logab.py+code:warmlateopen+code:on.dev+code:abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:FR.dev+code:year+(1|cell) + (1|CountID:cell), data=naba.1)
  extractAIC(ab.yr.1)
  t2<-as_tibble(bind_cols(Parameters=row.names(summary(ab.yr.1)$coefficients),summary(ab.yr.1)$coefficients)) %>% arrange(abs(`t value`))
  
  ab.yr.2<-lmer(log.abund~-1+code+code:warmearly+code:logab.py+code:warmlateopen+code:on.dev+code:abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:FR.dev+code:year + (1|cell) + (1|CountID:cell), data=naba.1)
  extractAIC(ab.yr.2)
  r.squaredGLMM(ab.yr.2)
  ab.yr.3<-lmer(log.abund~-1+code+code:warmearly+code:logab.py+code:warmlateopen+code:on.dev+code:abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:FR.dev+code:year +(1|cell) + (1|CountID:cell), data=naba.1)
  extractAIC(ab.yr.3)
  #ab.yr.4<-lmer(log.abund~-1+code+code:warmearly+logab.py+code:warmlateopen+code:on.dev+abslag+dur.dev+code:MonthF + code:MonthF:ObsDay+code:FR.dev+code:year +(1|cell) + (1|CountID:cell), data=naba.1)
  #extractAIC(ab.yr.4)
  #ab.yr.5<-lmer(log.abund~-1+code+code:warmearly+logab.py+code:warmlateopen+code:on.dev+abslag+code:MonthF+code:MonthF:ObsDay+code:FR.dev+code:year + (1|cell) + (1|CountID:cell), data=naba.1)
  #extractAIC(ab.yr.5)
  ab.final<-ab.yr.1
  ab.vif<-lmer(log.abund~code+warmearly+logab.py+warmlateopen+on.dev+abslag+dur.dev+ObsDay+MonthF+FR.dev+year+(1|cell) + (1|CountID:cell), data=naba.1)
  
  vif(ab.vif)
}else {
  naba.1<-na.omit(naba.1)

ab.yr.full<-lmer(log.abund~-1+code*(warmearly+logab.py+warmlateopen+on.dev+abslag+dur.dev+MonthF*ObsDay+year+FR.dev)+(1|cell) + (1|CountID:cell), data=naba.1)
r.squaredGLMM(ab.yr.full)     
summary(ab.yr.full)
(step_resyr <- step(ab.yr.full))
finalyr <- get_model(step_resyr) #stepAIC(ab.yr.full)
anova(finalyr)
summary(finalyr)
ab.final<-lmer(log.abund~-1+code+code:warmearly+code:on.dev+code:poly(on.dev,2)++code:poly(warmearly,2)+logab.py+code:warmlateopen+code:MonthF+code:MonthF:ObsDay+code:year+code:FR.dev+(1|cell) + (1|CountID:cell), data=naba.1)
extractAIC(ab.final)
ab.final1<-lmer(log.abund~-1+code+code:warmearly+code:on.dev+logab.py+code:warmlateopen+code:MonthF+code:MonthF:ObsDay+code:year+code:FR.dev+(1|cell) + (1|CountID:cell), data=naba.1)
extractAIC(ab.final1)

r.squaredGLMM(ab.final)
extractAIC(ab.final)
(step_ab.1 <- step(ab.final))
ab.final <- get_model(step_ab.1) #stepAIC(ab.yr.full)
ab.vif<-lmer(log.abund~code+warmearly+logab.py+warmlateopen+on.dev+ObsDay+MonthF+year+FR.dev+(1|cell) + (1|CountID:cell), data=naba.1)

vif(ab.vif)
#abund.best<-lmer(log.abund~-1+ows.grp+logab.py+abslag+ows.grp:warmearly+ows.grp:warmlateopen+ows.grp:on.dev+ows.grp:as.factor(ObsMonth)+ows.grp:doy+ows.grp:year+FR.dev+(1|cell) + (1|CountID:cell), data=naba.1)
summary(ab.final)
r.squaredGLMM(ab.final)
extractAIC(ab.final)
write.csv(summary(ab.final)$coefficients, file="output/abundance.coefs0831.csv")
plot_model(ab.final, type="eff", terms=c("ObsDay","code","MonthF"))

plot_model(ab.final1, type="eff", terms=c("on.dev","code"))

#plot_model(ab.final, type="eff", terms=c("dur.dev","code"))
plot_model(ab.final, type="eff", terms=c("FR.dev","code"))
plot_model(ab.final1, type="eff", terms=c("warmearly","code"))
plot_model(ab.final, type="eff", terms=c("warmlateopen","code"))
plot_model(ab.final, type="eff", terms=c("logab.py","code"))
plot_model(ab.final, type="eff", terms=c("year","code"))
plot_model(ab.final, type="eff", terms=c("on.dev","code","warmearly"))
plot_model(ab.final, type="eff", terms=c("warmearly","code","on.dev"))


}

names(naba.1)

#### PREDICT FOR VISUALIZATIONS
lims.pc1.ondev<-naba.1 %>%
#  select(code,logab.py,cell, CountID,warmearly,on.dev,warmlateopen,MonthF,year,ObsDay,FR.dev) %>%
  mutate(codelabel=ifelse(code=="RE","egg",ifelse(code=="RL","caterpillar","pupa"))) %>% 
#         logab.py=median(logab.py, na.rm=T),cell=median(cell, na.rm=T),
#         CountID=median(CountID, na.rm=T),) %>%
  group_by(code, codelabel) %>%
  summarize(minPC1=floor(min(warmearly, na.rm=T)*10)/10,maxPC1=ceiling(max(warmearly, na.rm=T)*10)/10,
            minond=floor(min(on.dev, na.rm=T)*10)/10,maxond=ceiling(max(on.dev, na.rm=T)*10)/10)
  
  
  

newcode<-c("RE","RL","RP")
codelabel<-c("egg","caterpillar","pupa")
newwarmearly<-round(c( ((min(naba.1$warmearly, na.rm=T))*2):((max(naba.1$warmearly, na.rm=T))*2)/2),2)
newab.py<-mean(naba.1$logab.py, na.rm=T)
newwarmlateopen<-mean(naba.1$warmlateopen, na.rm=T)
newondev<-round(c( ((round(min(naba.1$on.dev, na.rm=T)))):((round(max(naba.1$on.dev, na.rm=T))))),2)

nrep<-length(newcode)*length(newwarmearly)*length(newondev)
newDat <- data.frame(cell = rep(1,nrep), CountID=rep(1,nrep), 
                     abslag=rep(0,nrep), dur.dev=rep(0,nrep),
                     MonthF=as.factor(rep(7,nrep)),ObsDay=rep(1,nrep),
                     FR.dev=rep(0,nrep),year=rep(10,nrep),
                     logab.py=rep(newab.py,nrep), warmlateopen=rep(newwarmlateopen,nrep),
                     on.dev=rep(newondev, each=length(newwarmearly),3),
                     warmearly=rep(newwarmearly, 3*length(newondev)),
                     code=rep(newcode, each=length(newwarmearly)*length(newondev)),
                     codelabel=rep(codelabel, each=length(newwarmearly)*length(newondev)))

newDat$pred <- predict(ab.final, newDat,allow.new.levels =T)
library(viridis)

newDat1<-inner_join(newDat,lims.pc1.ondev) %>%
  mutate(f1=ifelse(warmearly>=minPC1,ifelse(warmearly<=maxPC1,1,0),0)+ifelse(on.dev>=minond,ifelse(on.dev<=maxond,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel=factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel,as.numeric(code)))


ggplot(data=newDat1, aes(x=on.dev, y=warmearly, fill=pred)) + 
  geom_tile() +  labs(x="Onset delay from baseline", y="GDD + greenup advance") + 
  scale_fill_viridis(name="log abundance") + 
  facet_wrap(~codelabel)


newDat <- data.frame(cell = rep(1,nrep), CountID=rep(1,nrep), 
                     abslag=rep(0,nrep), dur.dev=rep(0,nrep),
                     MonthF=as.factor(rep(7,nrep)),ObsDay=rep(15,nrep),
                     FR.dev=rep(0,nrep),year=rep(10,nrep),
                     logab.py=rep(newab.py,nrep), warmlateopen=rep(newwarmlateopen,nrep),
                     on.dev=rep(newondev, each=length(newwarmearly),3),
                     warmearly=rep(newwarmearly, 3*length(newondev)),
                     code=rep(newcode, each=length(newwarmearly)*length(newondev)),
                     codelabel=rep(codelabel, each=length(newwarmearly)*length(newondev)))

newDat$pred <- predict(ab.final, newDat,allow.new.levels =T)
newDat1<-inner_join(newDat,lims.pc1.ondev) %>%
  mutate(f1=ifelse(warmearly>=minPC1,ifelse(warmearly<=maxPC1,1,0),0)+ifelse(on.dev>=minond,ifelse(on.dev<=maxond,1,0),0)) %>%
  filter(f1==2)

ggplot(data=newDat1, aes(x=on.dev, y=warmearly, fill=pred)) + 
  geom_tile() +  labs(x="Onset delay from baseline", y="GDD + greenup advance") + 
  scale_fill_viridis(name="log abundance7") + 
  facet_wrap(~codelabel)

newDat<-naba.1[,-9]  %>%
  mutate(codelabel=ifelse(code=="RE","egg",ifelse(code=="RL","caterpillar","pupa")),
         logab.py=median(logab.py, na.rm=T), on.dev=median(on.dev, na.rm=T),
         MonthF=as.factor(7),year=10,FR.dev=0,ObsDay=15,cell=622,countID=median(CountID))
newDat$pred <- predict(ab.final, newDat,allow.new.levels =T)

ggplot(data=newDat, aes(x=warmlateopen, y=warmearly, color=pred)) + 
  geom_point(size=3) +  labs(x="PC2", y="PC1") + 
  scale_color_viridis(name="log abundance") + 
  facet_wrap(~codelabel)


newDat<-naba.1[,-9]  %>%
  mutate(codelabel=ifelse(code=="RE","egg",ifelse(code=="RL","caterpillar","pupa")),
         logab.py=median(logab.py, na.rm=T), warmearly=median(warmearly, na.rm=T),
         warmlateopen=median(warmlateopen, na.rm=T), abslag=abs(q50-196),
         MonthF=as.factor(6),year=10,ObsDay=15,cell=622,CountID=median(CountID))
newDat$pred <- predict(ab.final, newDat,allow.new.levels =T)

ggplot(data=newDat, aes(x=FR.dev, y=on.dev, color=pred)) + 
  geom_point(size=3) +  labs(x="# cold days above baseline", y="Onset delay") + 
  scale_color_viridis(name="log abundance") + 
  facet_wrap(~codelabel)


newDat<-naba.1[,-9]  %>%
  mutate(codelabel=ifelse(code=="RE","egg",ifelse(code=="RL","caterpillar","pupa")),
         logab.py=median(logab.py, na.rm=T), warmearly=median(warmearly, na.rm=T),
         warmlateopen=median(warmlateopen, na.rm=T), abslag=abs(q50-196),FR.dev=0,
         MonthF=as.factor(6),ObsDay=15,cell=622,CountID=median(CountID))
newDat$pred <- predict(ab.final, newDat,allow.new.levels =T)

ggplot(data=newDat, aes(x=year, y=on.dev, color=pred)) + 
  geom_point(size=3) +  labs(x="Year", y="Onset delay") + 
  scale_color_viridis(name="log abundance") + 
  facet_wrap(~codelabel)



naba.1 %>% filter(year>15) %>% group_by(cell, code) %>%
  tally()


#finalyr<-lmer(log.abund~-1+ows.grp+logab.py+warmearly:ows.grp+warmlateopen:ows.grp+on.dev:ows.grp+abslag+year:ows.grp+as.factor(ObsMonth):ows.grp+doy:ows.grp+(1|cell/CountID), data=naba.1)
#(step_resyr <- step(finalyr))
#finalyr <- get_model(step_resyr) #stepAIC(ab.yr.full)
#mixedyr<-lmer(log.abund~-1+ows.grp+logab.py+warmearly:ows.grp+warmlateopen:ows.grp+on.dev:ows.grp+abslag+year:ows.grp+year+as.factor(ObsMonth):ows.grp+doy:ows.grp+(1|cell/CountID), data=naba.1)
#mixedyr2 <- get_model(step(mixedyr)) #stepAIC(ab.yr.full)

mixedyr<-abund.best
anova(mixedyr)
ranova(mixedyr)
summary(mixedyr)
r.squaredGLMM(mixedyr)
AIC(mixedyr)
vif(mixedyr)
plot.lmer<-plot_model(mixedyr)
mixedlabels<-c(rep("",nrow(plot.lmer[[1]])))
mixedlabels[which(summary(mixedyr)$coefficients[,5]<0.05)]<-"*"
plot.lmer + annotate(geom="text", x=rev(c(1:nrow(plot.lmer[[1]]))), y=rep(-3,nrow(plot.lmer[[1]])), label=mixedlabels) + ggtitle("Abundance mixed effects model")

test.vif.re<-lmer(log.abund ~ logab.py + abslag + warmearly +  
                 warmlateopen + on.dev + as.factor(ObsMonth) +  
                 doy + year + FR.dev + (1 | cell) + (1 | CountID:cell), data= filter(naba.1, ows.grp=="RE"))
test.vif.rl<-lmer(log.abund ~ logab.py + abslag + warmearly +  
                 warmlateopen + on.dev + as.factor(ObsMonth) +  
                 doy + year + FR.dev + (1 | cell) + (1 | CountID:cell), data= filter(naba.1, ows.grp=="RL"))
test.vif.rp<-lmer(log.abund ~ logab.py + abslag + warmearly +  
                 warmlateopen + on.dev + as.factor(ObsMonth) +  
                 doy + year + FR.dev + (1 | cell) + (1 | CountID:cell), data= filter(naba.1, ows.grp=="RP"))
vif(test.vif.re)
vif(test.vif.rl)
vif(test.vif.rp)


#fixedyr<-lm(log.abund~-1+ows.grp+logab.py+abslag+warmearly:ows.grp+warmlateopen:ows.grp+on.dev:ows.grp+as.factor(ObsMonth):ows.grp+doy:ows.grp, data=naba.1)
fixedyr<-lm(log.abund~-1+ows.grp+logab.py+warmearly:ows.grp+warmlateopen:ows.grp+on.dev:ows.grp+abslag+year:ows.grp+as.factor(ObsMonth):ows.grp+doy:ows.grp, data=naba.1)
phtest_glmer(finalyr, fixedyr)
plotresult<-plot_model(finalyr, values=T)

summary(finalyr)
r.squaredGLMM(finalyr)

plot.lm<-plot_model(fixedyr, values=T)
fixedlabels<-c(rep("",nrow(plot.lm[[1]])))
fixedlabels[which(summary(fixedyr)$coefficients[,4]<0.05)]<-"*"
plot.lm + annotate(geom="text", x=rev(c(1:nrow(plot.lm[[1]]))), y=rep(-3,nrow(plot.lm[[1]])), label=fixedlabels) + ggtitle("Abundance fixed effects model")

ggplot(data=naba.1, aes(x=ows.grp, y=log.abund)) + geom_boxplot()
#ggplot(data=naba.1, aes(color=ows.grp, fill=ows.grp, x=logab.py, y=log.abund)) + geom_point()
ggplot(data=naba.1, aes(color=as.factor(CountID), fill=CountID, shape=ows.grp, x=logab.py, y=log.abund)) + geom_point() +theme(legend.position="none")
#ggplot(data=naba.1, aes( shape=ows.grp, x=logab.py, y=log.abund)) + geom_point()
ggplot(data=naba.1, aes(color=ows.grp, fill=ows.grp, x=warmearly, y=log.abund)) + geom_point() + geom_smooth(method="lm")
ggplot(data=naba.1, aes(color=ows.grp, fill=ows.grp, x=warmlateopen, y=log.abund)) + geom_point() + geom_smooth(method="lm")
ggplot(data=naba.1, aes(color=ows.grp, fill=ows.grp, x=warmlateopen, y=log.abund)) + geom_point() + geom_smooth(method="lm")
ggplot(data=naba.1, aes(color=ows.grp, fill=ows.grp, x=on.dev, y=log.abund)) + geom_point() + geom_smooth(method="lm")
ggplot(data=naba.1, aes(color=ows.grp, fill=ows.grp, x=doy, y=log.abund)) + geom_point() + geom_smooth(method="lm") + facet_wrap(.~ObsMonth)

ggplot(data=naba.1, aes(x=doy, y=abslag)) + geom_point() + geom_smooth(method="lm")
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
plotresult<-plot_model(fixedyr, show.values=T, show.p=T)

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



#all plot_models
models<-list(onset.best, onset.dev.best, med.best, med.dev.best, dur.best, dur.dev.best)

plot.panel<-function(model.i) {
  labels.i<-label.sig(model.i)
  coefs<-summary(model.i)$coefficients
  min.row<-which(coefs[,1]==min(coefs[,1]))
  ylab<-coefs[min.row,1]-2*(coefs[min.row,2])
  plot_model(model.i) + 
    annotate(geom="text", x=rev(c(1:length(labels.i))+.5), y=coefs[,1], label=labels.i) + 
    theme_bw() + 
    labs(subtitle = paste("marginal r-sq.=",round(r.squaredGLMM(onset.best),2)))
}

label.sig<-function(model.i) {
  fixedlabels<-c(rep("",nrow(summary(model.i)$coefficients)))
  fixedlabels[which(summary(model.i)$coefficients[,ncol(summary(model.i)$coefficients)]<0.05)]<-"*"
  fixedlabels
}
plot.panel(onset.best)
fixedlabels<-lapply(models, label.sig)

plots<-plot.panel(fixedyr)
pdf("output/abund.model.pdf")
plots
dev.off()


tempabund<-predict(mixedyr, data=naba.1)
plotdata<-cbind(naba.1, tempabund)
ggplot(data=)

ggplot(data=pheno.dev, aes(x=year, y=q5_dev, color=code)) + geom_point() + geom_smooth(method="lm") + 
  labs(y="Onset deviation")

