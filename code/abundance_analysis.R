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

load("data/spatial.domain.RData")
#abundance metrics
theme_set(theme_sjplot())
###PHENO DATA
load("data/derived/pheno.RData")
#### ENV DATA
load("data/derived/env.input.new.RData")
pheno1<-pheno.quant %>%
  dplyr::select(q50,year,cell,code)

(p1<-pheno.quant %>% group_by(code) %>% summarize(meanonset=mean(q5, na.rm=T), meandur=mean(qdur, na.rm=T), meanmed=mean(q50, na.rm=T)))

env10<-env.var %>%
  group_by(year) %>%
  summarize(pc1=mean(pc1,na.rm=T), pc2=mean(pc2,na.rm=T),wd=mean(warmdays,na.rm=T)/10,cd=mean(colddays,na.rm=T)/10, summer=mean(summer.gdd,na.rm=T),openlag=mean(gr_mn_lag,na.rm=T)/2, lat=mean(cell_lat,na.rm=T))

ggplot(data=filter(env10, year<2018), aes(x=year, y=pc1*2+1)) + geom_line() + 
  geom_line(data=env10, aes(x=year, y=pc2*2), color="purple") + 
  geom_line(data=env10, aes(x=year, y=wd-5), color="blue") + 
  geom_line(data=env10, aes(x=year, y=cd), color="slateblue") + 
  geom_line(data=env10, aes(x=year, y=summer), color="forestgreen") + 
  geom_line(data=env10, aes(x=year, y=openlag), color="springgreen") + 
  theme_minimal()
  


abund<-read_csv("data/derived/naba_OWS_abundances.csv") %>%
  dplyr::select(cell, year=ObsYear, doy, CountID, SurveyID, ObsMonth, code=group, abund.bph, log.abund, SR) %>%
  filter(ObsMonth %in% c(6:8), cell %in% STUDYCELLS)

(a1<-abund %>% filter(year<2018) %>% group_by(code) %>% summarize(meanbph=mean(abund.bph, na.rm=T), sdbph=sd(abund.bph, na.rm=T)))
#previous year abundance
abund.py<-abund %>%
  mutate(year=year+1) %>%
  dplyr::select(cell, year, CountID, code, abund.py=abund.bph, logab.py=log.abund)

a2<-abund %>% filter(year<2018)  %>% group_by(year, code) %>% summarize(ct=mean(log.abund, na.rm=T))
ggplot(data=a2, aes(x=year, y=ct, color=code)) + geom_line() + geom_smooth(method="lm")

a3<-a2 %>% pivot_wider(id_cols=year, names_from=code, values_from=ct)
a3b<-abund %>% pivot_wider(id_cols=c(year,CountID,doy), names_from=code, values_from=log.abund)
cor(a3)

cor(na.omit(a3b[,c(1,4:6)]), method="pearson")

#combine tables
ab0<-merge(abund, pheno1, by=c("year","cell","code"), all.x=TRUE)
ab<-merge(x = ab0, y = pheno.dev, by = intersect(names(ab0), names(pheno.dev)), all.x = TRUE)
ab1<-merge(x = ab, y = env.dev, by = intersect(names(ab), names(env.dev)), all.x = TRUE)
ab.final<-merge(x = ab1, y = abund.py, by.x=c("cell", "year", "CountID", "code"), by.y=c("cell", "year", "CountID", "code"), all.x = TRUE)
#abundance Model
naba.1<-(ab.final) %>% mutate(code=as.factor(code),summer.dev=summer.dev,year=year-2000, on.dev=onset.dev/7, dur.dev=dur.dev/7, lag.dev=lag.dev/7, cold.dev=cold.dev/7, daylag=doy-q50, abslag=abs(doy-q50))
naba.1<-naba.1 %>% mutate(MonthF=as.factor(ObsMonth), ObsDay=as.numeric(format(as.Date(doy, origin = paste((year+1999),"12-31",sep="-")),format="%d")), eggows=ifelse(code=="RE",1,0))

save(naba.1, file="data/abund.input1005.RData")

(n.cellyr.summary<-naba.1 %>%
  group_by(cell, code) %>%
  summarize(x1=1) %>%
  group_by(code) %>%
  summarize(cs1=sum(x1))
)

(pdev.sum<-pheno.dev %>%
  group_by(cell, code) %>%
  summarize(x1=1) %>%
  group_by(code) %>%
  summarize(cs1=sum(x1))
) 

ggplot(data=naba.1, aes(x=onset.dev, y=log.abund, color=code)) + 
  geom_smooth(method="lm")


ab.yr<-lmer(log.abund~-1+code+code:year+(1|cell) + (1|CountID:cell), data=naba.1)
summary(ab.yr)
#ab2.yr<-lmer(abund.bph~-1+code+code:year+(1|cell) + (1|CountID:cell), data=naba.1)
#summary(ab2.yr)


summary(naba.1)
includeNA<-T
if(includeNA==T) {
  ab.yr.full<-lmer(log.abund~-1+code*(pc1.dev+logab.py+pc2.dev+on.dev+abslag+dur.dev+MonthF*ObsDay+year+cold.dev)+(1|cell) + (1|CountID:cell), data=naba.1)

  extractAIC(ab.yr.full)
  ab.yr.1<-lmer(log.abund~-1+code+code:pc1.dev+code:logab.py+code:pc2.dev+code:on.dev+code:abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+code:year+(1|cell) + (1|CountID:cell), data=naba.1)
  extractAIC(ab.yr.1)
  (t2<-as_tibble(bind_cols(Parameters=row.names(summary(ab.yr.1)$coefficients),summary(ab.yr.1)$coefficients)) %>% arrange(abs(`t value`)))
  
  best.model<-ab.yr.1
  ab.yr.noint1<-lmer(log.abund~-1+code+pc1.dev+code:logab.py+code:pc2.dev+code:on.dev+code:abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+code:year+(1|cell) + (1|CountID:cell), data=naba.1)
  extractAIC(ab.yr.noint1)
  if(extractAIC(ab.yr.noint1)[2]<extractAIC(best.model)[2]) {best.model<-ab.yr.noint1}
  
  ab.yr.nointlag<-lmer(log.abund~-1+code+code:pc1.dev+code:logab.py+code:pc2.dev+code:on.dev+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+code:year+(1|cell) + (1|CountID:cell), data=naba.1)
  extractAIC(ab.yr.nointlag)
  if(extractAIC(ab.yr.nointlag)[2]<extractAIC(best.model)[2]) {best.model<-ab.yr.nointlag}
  extractAIC(best.model)
  
  ab.yr.nointyr<-lmer(log.abund~-1+code+code:pc1.dev+code:logab.py+code:pc2.dev+code:on.dev+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+year+(1|cell) + (1|CountID:cell), data=naba.1)
  extractAIC(ab.yr.nointyr)
  if(extractAIC(ab.yr.nointyr)[2]<extractAIC(best.model)[2]) {best.model<-ab.yr.nointyr}
  extractAIC(best.model)
 
  ab.yr.noyr<-lmer(log.abund~-1+code+code:pc1.dev+code:logab.py+code:pc2.dev+code:on.dev+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+(1|cell) + (1|CountID:cell), data=naba.1)
  extractAIC(ab.yr.noyr)
  if(extractAIC(ab.yr.noyr)[2]<extractAIC(best.model)[2]) {best.model<-ab.yr.noyr}
  extractAIC(best.model)

  ab.yr.nopyint<-lmer(log.abund~-1+code+code:pc1.dev+logab.py+code:pc2.dev+code:on.dev+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+(1|cell) + (1|CountID:cell), data=naba.1)
  extractAIC(ab.yr.nopyint)
  if(extractAIC(ab.yr.nopyint)[2]<extractAIC(best.model)[2]) {best.model<-ab.yr.nopyint}
  extractAIC(best.model)
  
  
  
  
  
  (yrpr2<-r.squaredGLMM(best.model)[1]-r.squaredGLMM(ab.yr.noyr)[1])
  ab.yr.nopc1<-lmer(log.abund~-1+code+code:logab.py+code:pc2.dev+code:on.dev+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+year+(1|cell) + (1|CountID:cell), data=naba.1)
  (pc1pr2<-r.squaredGLMM(best.model)[1]-r.squaredGLMM(ab.yr.nopc1)[1])
  ab.yr.nopc2<-lmer(log.abund~-1+code+code:logab.py+code:pc1.dev+code:on.dev+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+year+(1|cell) + (1|CountID:cell), data=naba.1)
  (pc2pr2<-r.squaredGLMM(best.model)[1]-r.squaredGLMM(ab.yr.nopc2)[1])
  ab.yr.noon<-lmer(log.abund~-1+code+code:logab.py+code:pc1.dev+code:pc2.dev+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+year+(1|cell) + (1|CountID:cell), data=naba.1)
  (onpr2<-r.squaredGLMM(best.model)[1]-r.squaredGLMM(ab.yr.noon)[1])
  ab.yr.py<-lmer(log.abund~-1+code+code:on.dev+code:pc1.dev+code:pc2.dev+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+year+(1|cell) + (1|CountID:cell), data=naba.1)
  (pypr2<-r.squaredGLMM(best.model)[1]-r.squaredGLMM(ab.yr.py)[1])
  ab.yr.dd<-lmer(log.abund~-1+code+code:on.dev+code:pc1.dev+code:pc2.dev+abslag+code:logab.py+code:MonthF + code:MonthF:ObsDay+code:cold.dev+year+(1|cell) + (1|CountID:cell), data=naba.1)
  (ddpr2<-r.squaredGLMM(best.model)[1]-r.squaredGLMM(ab.yr.dd)[1])
  ab.yr.cd<-lmer(log.abund~-1+code+code:on.dev+code:pc1.dev+code:pc2.dev+abslag+code:logab.py+code:MonthF + code:MonthF:ObsDay+code:dur.dev+year+(1|cell) + (1|CountID:cell), data=naba.1)
  (cdpr2<-r.squaredGLMM(best.model)[1]-r.squaredGLMM(ab.yr.cd)[1])
  
    ab.vif<-lmer(log.abund~logab.py+code+pc1.dev+pc2.dev+on.dev+abslag+dur.dev+ObsDay:MonthF+cold.dev+year+(1|cell) + (1|CountID:cell), data=naba.1)
    vif(ab.vif)
    
    

    
abund.output<-as_tibble(summary(best.model)$coefficients) %>%
  dplyr::mutate(param=row.names(summary(best.model)$coefficients), sig=ifelse(`Pr(>|t|)`<0.05,1,0), r2m=round(r.squaredGLMM(best.model)[1],2),r2c=round(r.squaredGLMM(best.model)[2],2))
write.csv(abund.output, file="output/abund.model1005.csv")

yr.only<-lmer(log.abund ~ -1+code+code:year + code:logab.py+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay + (1|cell)+(1|CountID:cell), data=naba.1)
summary(yr.only)              


plot_model(best.model, type="eff", terms=c("ObsDay","code","MonthF"))

plot_model(best.model, type="eff", terms=c("on.dev","code"))

plot_model(best.model, type="eff", terms=c("dur.dev","code"))
plot_model(best.model, type="eff", terms=c("cold.dev","code"))
plot_model(best.model, type="eff", terms=c("pc1.dev","code"))
plot_model(best.model, type="eff", terms=c("pc2.dev","code"))
plot_model(best.model, type="eff", terms=c("logab.py","code"))
plot_model(best.model, type="eff", terms=c("year","code"))
plot_model(best.model, type="eff", terms=c("on.dev","code","warmearly"))
plot_model(best.model, type="eff", terms=c("pc1.dev","on.dev","code"))

}

names(naba.1)

#### PREDICT FOR VISUALIZATIONS



lims.pc1.ondev<-naba.1 %>%
#  select(code,logab.py,cell, CountID,warmearly,on.dev,warmlateopen,MonthF,year,ObsDay,FR.dev) %>%
  mutate(codelabel=ifelse(code=="RE","BOE",ifelse(code=="RL","BOL","BOP"))) %>% 
#         logab.py=median(logab.py, na.rm=T),cell=median(cell, na.rm=T),
#         CountID=median(CountID, na.rm=T),) %>%
  group_by(code, codelabel) %>%
  summarize(minPC1=floor(min(pc1.dev, na.rm=T)*10)/10,maxPC1=ceiling(max(pc1.dev, na.rm=T)*10)/10,
            minPC2=floor(min(pc2.dev, na.rm=T)*10)/10,maxPC2=ceiling(max(pc2.dev, na.rm=T)*10)/10,
            minond=floor(min(on.dev, na.rm=T)*10)/10,maxond=ceiling(max(on.dev, na.rm=T)*10)/10,
            mindurd=floor(min(dur.dev, na.rm=T)*10)/10,maxdurd=ceiling(max(dur.dev, na.rm=T)*10)/10,
            minyr=floor(min(year, na.rm=T)*10)/10,maxyr=ceiling(max(year, na.rm=T)*10)/10,
            mincold=floor(min(cold.dev, na.rm=T)*10)/10,maxcold=ceiling(max(cold.dev, na.rm=T)*10)/10  )

  
  

newcode<-c("RE","RL","RP")
codelabel<-c("BOE","BOL","BOP")
nrep<-1200

abund.newDat <- data.frame(cell = rep(507,nrep), CountID=rep(497,nrep),
                         logab.py=rep(median(naba.1$logab.py, na.rm=T),nrep),
                         cold.dev=rep(0,nrep),
                         abslag=rep(10,nrep),
                         pc1.dev=rep(0,nrep),
                         pc2.dev=rep(0,nrep),
                         on.dev=rep(0,nrep), dur.dev=rep(0,nrep),
                         MonthF=as.factor(rep(6,nrep)),ObsDay=rep(30,nrep),
                         year=rep(10, nrep),
                         code=rep(newcode, each=20,20),
                         codelabel=rep(codelabel, each=20,20)) %>%
  mutate(codelabel = factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel, as.numeric(code)))

dev.maxvals<-(as.numeric(apply(na.omit(naba.1[,c(2,12,16,19,24:37)]),2,FUN=max)))
dev.minvals<-(as.numeric(apply(na.omit(naba.1[,c(2,12,16,19,24:37)]),2,FUN=min)))
dev.intervals<-round((dev.maxvals-dev.minvals)/20,3)
dev.names<-names(naba.1)[c(2,12,16,19,24:37)]

#pc1 x year
newpc1<-(-1)+0:19*(0.2)
newyr<-0:19
nrep<-1200
minpc1<-naba.1 %>% filter(year %in% c(5,15)) %>%group_by(year) %>% summarize(minpc1=min(pc1.dev, na.rm=T),maxpc1=max(pc1.dev, na.rm=T))
remove.cols<-which(names(abund.newDat) %in% c("pc1.dev","year"))
abund.newDat1<-cbind(abund.newDat[,-remove.cols], pc1.dev=rep(newpc1,60), year=rep(newyr,each=60))
abund.newDat1$pred <- predict(best.model, abund.newDat1, allow.new.levels =T)

lims.ab<-naba.1 %>%
  group_by(code) %>%
  summarize(minyr=min(year, na.rm=T),maxyr=max(year, na.rm=T),
            minpc1=min(pc1.dev, na.rm=T),maxpc1=max(pc1.dev, na.rm=T),
            minpc2=min(pc2.dev, na.rm=T), maxpc2=max(pc2.dev, na.rm=T))

abund.newDat.F1<-inner_join(abund.newDat1,lims.ab) %>%
  mutate(f1=ifelse(year>=minyr,ifelse(year<=maxyr,1,0),0)+ifelse(pc1.dev>=minpc1,ifelse(pc1.dev<=maxpc1,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel = factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel, as.numeric(code)))


##for fig 
library(viridis)
pc1.yr<-naba.1 %>% filter(year>0) %>%group_by(year) %>% summarize(minpc1=min(pc1.dev, na.rm=T),maxpc1=max(pc1.dev, na.rm=T))

F.abund.1<-ggplot(data=abund.newDat.F1, aes(x=year, y=pc1.dev, fill=pred)) + 
  geom_tile() +  labs(x="Year", y="PC1 Deviation") + 
  scale_fill_viridis(name="Log abundance") +
  geom_line(data=pc1.yr,  inherit.aes = FALSE,aes(x=year, y=maxpc1), color="black") + 
  geom_line(data=pc1.yr,  inherit.aes = FALSE,aes(x=year, y=minpc1), color="white") + 
#  geom_point(data=minpc1, aes(x=year, y=minpc1), shape=24, fill=NA, color="white") + 
#  geom_point(data=minpc1, aes(x=year, y=maxpc1), shape=25,fill=NA, color="white") + 
  facet_wrap(~codelabel)
F.abund.1

save(F.abund.1,file="output/abund.fig.1.png")

#### FIG: pc1 and cold winter days
#pc1 x year
newpc1<-dev.minvals[which(dev.names%in%"pc1.dev")]+0:19*(dev.intervals[which(dev.names%in%"pc1.dev")])
newcold<-dev.minvals[which(dev.names%in%"cold.dev")]+0:19*(dev.intervals[which(dev.names%in%"cold.dev")])
nrep<-1200

remove.cols<-which(names(abund.newDat) %in% c("pc1.dev","cold.dev"))
abund.newDat2<-cbind(abund.newDat[,-remove.cols], pc1.dev=rep(newpc1,60), cold.dev=rep(newcold,each=60))
abund.newDat2$pred <- predict(best.model, abund.newDat2, allow.new.levels =T)

lims.ab<-naba.1 %>%
  group_by(code) %>%
  summarize(minyr=min(year, na.rm=T),maxyr=max(year, na.rm=T),
            minpc1=min(pc1.dev, na.rm=T),maxpc1=max(pc1.dev, na.rm=T),
            minpc2=min(pc2.dev, na.rm=T), maxpc2=max(pc2.dev, na.rm=T),
            mincold=min(cold.dev, na.rm=T), maxcold=max(cold.dev, na.rm=T),
            minondev=min(on.dev, na.rm=T), maxondev=max(on.dev, na.rm=T))

abund.newDat.F2<-inner_join(abund.newDat2,lims.ab) %>%
  mutate(f1=ifelse(cold.dev>=mincold,ifelse(cold.dev<=maxcold,1,0),0)+ifelse(pc1.dev>=minpc1,ifelse(pc1.dev<=maxpc1,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel = factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel, as.numeric(code)))


##for fig 
library(viridis)
F.abund.2<-ggplot(data=abund.newDat.F2, aes(x=cold.dev, y=pc1.dev, fill=pred)) + 
  geom_tile() +  labs(x="Cold winter days", y="PC1 deviation") + 
  geom_point(x=2,y=0.5,shape=24, color="white", size=1.5) + 
  scale_fill_viridis(name="Log abundance") + 
  facet_wrap(~codelabel)
F.abund.2

save(F.abund.2,file="output/abund.fig.2.png")


#### FIG: winter and onset deviation
#pc1 x year
newcold<-dev.minvals[which(dev.names%in%"cold.dev")]+0:19*(dev.intervals[which(dev.names%in%"cold.dev")])
newod<-(-8)+0:19
nrep<-1200

remove.cols<-which(names(abund.newDat) %in% c("cold.dev","on.dev"))
abund.newDat3<-cbind(abund.newDat[,-remove.cols], cold.dev=rep(newcold,60), on.dev=rep(newod,each=60))
abund.newDat3$pred <- predict(best.model, abund.newDat3, allow.new.levels =T)

lims.ab<-naba.1 %>%
  group_by(code) %>%
  summarize(minyr=min(year, na.rm=T),maxyr=max(year, na.rm=T),
            minpc1=min(pc1.dev, na.rm=T),maxpc1=max(pc1.dev, na.rm=T),
            minpc2=min(pc2.dev, na.rm=T), maxpc2=max(pc2.dev, na.rm=T),
            mincold=min(cold.dev, na.rm=T), maxcold=max(cold.dev, na.rm=T),
            minondev=min(on.dev, na.rm=T)-1, maxondev=max(on.dev, na.rm=T)+1)

abund.newDat.F3<-inner_join(abund.newDat3,lims.ab) %>%
  mutate(f1=ifelse(on.dev>=minondev,ifelse(on.dev<=maxondev,1,0),0)+ifelse(pc2.dev>=minpc2,ifelse(pc2.dev<=maxpc2,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel = factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel, as.numeric(code)))


##for fig 
library(viridis)
F.abund.3<-ggplot(data=abund.newDat.F3, aes(x=on.dev, y=pc2.dev, fill=pred)) + 
  geom_tile() +  labs(x="Onset deviation (wks)", y="PC2 deviation") + 
  scale_fill_viridis(name="Log abundance") + 
  facet_wrap(~codelabel)
F.abund.3

save(F.abund.3,file="output/abund.fig.3.png")







newwarmearly<-round(c( ((min(naba.1$pc1.dev, na.rm=T))*5):((max(naba.1$pc1.dev, na.rm=T))*5)/5),2)
newab.py<-mean(naba.1$logab.py, na.rm=T)
newpc2<-round(c( ((min(naba.1$pc2.dev, na.rm=T))*10):((max(naba.1$pc2.dev, na.rm=T))*10)/10),2)
newondev<-round(c( ((round(min(naba.1$on.dev, na.rm=T)))):((round(max(naba.1$on.dev, na.rm=T))))),2)
newdurdev<-round(c( ((round(min(naba.1$dur.dev, na.rm=T)))):((round(max(naba.1$dur.dev, na.rm=T))))),2)

nrep<-length(newcode)*length(newwarmearly)*length(newpc2)
newDatpc <- data.frame(cell = rep(1,nrep), CountID=rep(1,nrep), 
                     abslag=rep(0,nrep), dur.dev=rep(0,nrep),
                     MonthF=as.factor(rep(6,nrep)),ObsDay=rep(1,nrep),
                     cold.dev=rep(0,nrep),year=rep(10,nrep),
                     logab.py=rep(newab.py,nrep), 
                     on.dev=rep(0, nrep),
                     lag.dev=rep(mean(naba.1$lag.dev,na.rm=T),nrep),
                     pc1.dev=rep(newwarmearly, 3*length(newpc2)),
                     pc2.dev=rep(newpc2,each=length(newwarmearly),3),
                     code=rep(newcode, each=length(newwarmearly)*length(newpc2)),
                     codelabel=rep(codelabel, each=length(newwarmearly)*length(newpc2)))

newDatpc$pred <- predict(ab.yr.2, newDatpc,allow.new.levels =T)
library(viridis)

newDat1<-inner_join(newDatpc,lims.pc1.ondev) %>%
  mutate(f1=ifelse(pc1.dev>=minPC1,ifelse(pc1.dev<=maxPC1,1,0),0)+ifelse(on.dev>=minond,ifelse(on.dev<=maxond,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel=factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel,as.numeric(code)))


ggplot(data=newDat1, aes(x=pc2.dev, y=pc1.dev, fill=pred)) + 
  geom_tile() +  labs(x="PC2 deviation", y="PC1 deviation") + 
  scale_fill_viridis(name="log abundance") + 
  facet_wrap(~codelabel)


##### on.dev and dur.dev
nrep<-length(newcode)*length(newondev)*length(newdurdev)

newDatpheno <- data.frame(cell = rep(1,nrep), CountID=rep(1,nrep), 
                          lag.dev=rep(mean(naba.1$lag.dev,na.rm=T),nrep),
                          abslag=rep(0,nrep), 
                       MonthF=as.factor(rep(6,nrep)),ObsDay=rep(1,nrep),
                       cold.dev=rep(0,nrep),year=rep(10,nrep),
                       logab.py=rep(newab.py,nrep), 
                       on.dev=rep(newondev, 3*length(newdurdev)),
                       dur.dev=rep(newdurdev,each=length(newondev),3),
                       pc1.dev=rep(0, nrep),
                       pc2.dev=rep(0, nrep),
                       code=rep(newcode, each=length(newondev)*length(newdurdev)),
                       codelabel=rep(codelabel, each=length(newondev)*length(newdurdev)))

newDatpheno$pred <- predict(ab.yr.2, newDatpheno,allow.new.levels =T)

newDat2<-inner_join(newDatpheno,lims.pc1.ondev) %>%
  mutate(f1=ifelse(dur.dev>=mindurd,ifelse(dur.dev<=maxdurd,1,0),0)+ifelse(on.dev>=minond,ifelse(on.dev<=maxond,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel=factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel,as.numeric(code)))


ggplot(data=newDat2, aes(x=on.dev, y=dur.dev, fill=pred)) + 
  geom_tile() +  labs(x="Onset deviation", y="Duration deviation") + 
  scale_fill_viridis(name="log abundance") + 
  facet_wrap(~codelabel)






##### on.dev and dur.dev
newcold<-round(c( ((round(min(naba.1$cold.dev, na.rm=T)))):((round(max(naba.1$cold.dev, na.rm=T))))),2)

newyear<-c(0:17)
nrep<-length(newcode)*length(newcold)*length(newyear)

newDat3 <- data.frame(cell = rep(1,nrep), CountID=rep(1,nrep), 
                          abslag=rep(0,nrep), 
                          MonthF=as.factor(rep(6,nrep)),ObsDay=rep(1,nrep),
                          logab.py=rep(newab.py,nrep), 
                          on.dev=rep(0, nrep),
                          dur.dev=rep(0,nrep),
                          pc1.dev=rep(0, nrep),
                          pc2.dev=rep(0, nrep),
                          cold.dev=rep(newcold,3*length(newyear)),
                          year=rep(newyear,each=length(newcold),3),
                          code=rep(newcode, each=length(newcold)*length(newyear)),
                          codelabel=rep(codelabel, each=length(newcold)*length(newyear)))

newDat3$pred <- predict(ab.final1, newDat3,allow.new.levels =T)

newDat3b<-inner_join(newDat3,lims.pc1.ondev) %>%
  mutate(f1=ifelse(cold.dev>=mincold,ifelse(cold.dev<=maxcold,1,0),0)+ifelse(year>=minyr,ifelse(year<=maxyr,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel=factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel,as.numeric(code)),year=year+2000)


ggplot(data=newDat3b, aes(x=year, y=cold.dev, fill=pred)) + 
  geom_tile() +  labs(x="Year", y="Winter cold deviation") + 
  scale_fill_viridis(name="log abundance") + 
  facet_wrap(~codelabel)






#######################3333



newcode<-c("RE","RL","RP")
codelabel<-c("egg","caterpillar","pupa")
newwarmearly<-round(c( ((min(naba.1$warmearly, na.rm=T))*2):((max(naba.1$warmearly, na.rm=T))*2)/2),2)
newab.py<-mean(naba.1$logab.py, na.rm=T)
newwarmlateopen<-round(c( ((min(naba.1$warmlateopen, na.rm=T))*2):((max(naba.1$warmlateopen, na.rm=T))*2)/2),2)
newondev<-0
nrep<-length(newcode)*length(newwarmearly)*length(newwarmlateopen)
newDat <- data.frame(cell = rep(1,nrep), CountID=rep(1,nrep), 
                     abslag=rep(0,nrep), dur.dev=rep(0,nrep),
                     MonthF=as.factor(rep(7,nrep)),ObsDay=rep(1,nrep),
                     FR.dev=rep(0,nrep),year=rep(10,nrep),
                     logab.py=rep(newab.py,nrep), warmlateopen=rep(newwarmlateopen, each=length(length(newwarmearly)),length(newcode)),
                     on.dev=rep(newondev, nrep),
                     warmearly=rep(newwarmearly, 3*length(newwarmlateopen)),
                     code=rep(newcode, each=nrep/3),
                     codelabel=rep(codelabel, each=nrep/3))

newDat$pred <- predict(ab.final, newDat,allow.new.levels =T)
library(viridis)

newDat3<-inner_join(newDat,lims.pc1.ondev) %>%
  mutate(f1=ifelse(warmearly>=minPC1,ifelse(warmearly<=maxPC1,1,0),0)+ifelse(on.dev>=minond,ifelse(on.dev<=maxond,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel=factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel,as.numeric(code)))


ggplot(data=newDat3, aes(x=warmlateopen, y=warmearly, fill=pred)) + 
  geom_tile() +  labs(x="GDD + open canopy lag", y="GDD + greenup advance") + 
  scale_fill_viridis(name="log abundance") + 
  facet_wrap(~codelabel)


############ ONSET DEVIATION & YEAR

newcode<-c("RE","RL","RP")
codelabel<-c("egg","caterpillar","pupa")
newwarmearly<-mean(naba.1$warmearly, na.rm=T)
newab.py<-mean(naba.1$logab.py, na.rm=T)
newwarmlateopen<-mean(naba.1$warmlateopen, na.rm=T)
newondev<-round(c( ((round(min(naba.1$on.dev, na.rm=T)))):((round(max(naba.1$on.dev, na.rm=T))))),2)
newyear<-c(0:17)
nrep<-length(newcode)*length(newyear)*length(newondev)
newDat <- data.frame(cell = rep(1,nrep), CountID=rep(1,nrep), 
                     abslag=rep(0,nrep), dur.dev=rep(0,nrep),
                     MonthF=as.factor(rep(7,nrep)),ObsDay=rep(1,nrep),
                     FR.dev=rep(0,nrep),year=rep(newyear,nrep/length(newyear)),
                     logab.py=rep(newab.py,nrep), warmlateopen=rep(newwarmlateopen,nrep),
                     on.dev=rep(newondev, each=length(newyear),3),
                     warmearly=rep(newwarmearly,nrep),
                     code=rep(newcode, each=length(newyear)*length(newondev)),
                     codelabel=rep(codelabel, each=length(newyear)*length(newondev)))

newDat$pred <- predict(ab.final, newDat,allow.new.levels =T)
library(viridis)

newDat4<-inner_join(newDat,lims.pc1.ondev) %>%
  mutate(f1=ifelse(warmearly>=minPC1,ifelse(warmearly<=maxPC1,1,0),0)+ifelse(on.dev>=minond,ifelse(on.dev<=maxond,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel=factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel,as.numeric(code)))


ggplot(data=newDat4, aes(x=on.dev, y=year, fill=pred)) + 
  geom_tile() +  labs(x="Onset delay from baseline", y="Year") + 
  scale_fill_viridis(name="log abundance") + 
  facet_wrap(~codelabel)
















##################



newwarmearly<-round(c( ((min(naba.1$warmearly, na.rm=T))*2):((max(naba.1$warmearly, na.rm=T))*2)/2),2)
newab.py<-mean(naba.1$logab.py, na.rm=T)
newwarmlateopen<-round(c( ((min(naba.1$warmlateopen, na.rm=T))*2):((max(naba.1$warmlateopen, na.rm=T))*2)/2),2)
newondev<-0
nrep<-length(newcode)*length(newwarmearly)*length(newwarmlateopen)
newDat <- data.frame(cell = rep(1,nrep), CountID=rep(1,nrep), 
                     abslag=rep(0,nrep), dur.dev=rep(0,nrep),
                     MonthF=as.factor(rep(7,nrep)),ObsDay=rep(1,nrep),
                     FR.dev=rep(0,nrep),year=rep(10,nrep),
                     logab.py=rep(newab.py,nrep), warmlateopen=rep(newwarmlateopen, each=length(length(newwarmearly)),length(newcode)),
                     on.dev=rep(newondev, nrep),
                     warmearly=rep(newwarmearly, 3*length(newwarmlateopen)),
                     code=rep(newcode, each=nrep/3),
                     codelabel=rep(codelabel, each=nrep/3))

newDat$pred <- predict(ab.final, newDat,allow.new.levels =T)
library(viridis)

newDat3<-inner_join(newDat,lims.pc1.ondev) %>%
  mutate(f1=ifelse(warmearly>=minPC1,ifelse(warmearly<=maxPC1,1,0),0)+ifelse(on.dev>=minond,ifelse(on.dev<=maxond,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel=factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel,as.numeric(code)))


ggplot(data=newDat3, aes(x=warmlateopen, y=warmearly, fill=pred)) + 
  geom_tile() +  labs(x="GDD + open canopy lag", y="GDD + greenup advance") + 
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



#### FIG: winter and onset deviation
#pc1 x year
newcold<-dev.minvals[which(dev.names%in%"cold.dev")]+0:19*(dev.intervals[which(dev.names%in%"cold.dev")])
newod<-(-8)+0:19
nrep<-1200

remove.cols<-which(names(abund.newDat) %in% c("cold.dev","on.dev"))
abund.newDat3<-cbind(abund.newDat[,-remove.cols], cold.dev=rep(newcold,60), on.dev=rep(newod,each=60))
abund.newDat3$pred <- predict(best.model, abund.newDat3, allow.new.levels =T)

lims.ab<-naba.1 %>%
  group_by(code) %>%
  summarize(minyr=min(year, na.rm=T),maxyr=max(year, na.rm=T),
            minpc1=min(pc1.dev, na.rm=T),maxpc1=max(pc1.dev, na.rm=T),
            minpc2=min(pc2.dev, na.rm=T), maxpc2=max(pc2.dev, na.rm=T),
            mincold=min(cold.dev, na.rm=T), maxcold=max(cold.dev, na.rm=T),
            minondev=min(on.dev, na.rm=T)-1, maxondev=max(on.dev, na.rm=T)+1)

abund.newDat.F3<-inner_join(abund.newDat3,lims.ab) %>%
  mutate(f1=ifelse(on.dev>=minondev,ifelse(on.dev<=maxondev,1,0),0)+ifelse(pc2.dev>=minpc2,ifelse(pc2.dev<=maxpc2,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel = factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel, as.numeric(code)))


##for fig 
library(viridis)
F.abund.3<-ggplot(data=abund.newDat.F3, aes(x=on.dev, y=cold.dev, fill=pred)) + 
  geom_tile() +  labs(x="Onset deviation (wks)", y="Winter cold deviation") + 
  scale_fill_viridis(name="Log abundance") + 
  facet_wrap(~codelabel)
F.abund.3

save(F.abund.3,file="output/abund.fig.3.png")




#######################3
abundx<-na.omit(naba.1)
ab.full<-lmer(log.abund~-1+code*(pc1.dev+logab.py+pc2.dev+on.dev+abslag+dur.dev+MonthF*ObsDay+year+cold.dev)+(1|cell) + (1|CountID:cell), data=abundx, na.action=na.fail)


xd1<-dredge(ab.full, beta = c("none"), evaluate = TRUE,
            rank = "AICc", fixed = c("code"), m.lim = NULL,
            trace = FALSE)

summary(xd1)
xd2<-xd1[xd1$weight>0.01,]
x3<-model.avg(xd2)

print(x3)

