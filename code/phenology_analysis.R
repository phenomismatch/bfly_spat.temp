#Analysis of Presence-only Phenology Metrics
#eButterfly, NABA Butterflies I've Seen, iNaturalist (research grade)
## Resident butterfly species groups defined by overwinter stage (OWS)
## Phenometrics estimated from Weibull (M Belitz)
## Species OWS traits compiled by GU Ries Lab
#E Larsen, Georgetown U, Updated 2021-09




#libraries
library(tidyverse)
library(ggplot2)
library(ggcorrplot)
library(lme4)
library(lmerTest)
library(MuMIn)
library(car)
library(sp)
library(sjPlot)
library(viridis)
theme_set(theme_sjplot())
###PHENO DATA
load("data/derived/pheno.RData")

(hexyrs<-pheno.quant %>% dplyr::select(year, cell)  %>% group_by(year, cell) %>% tally())
nrow(pheno.quant)
(dev.hexyrs<-pheno.dev %>% dplyr::select(year, cell)  %>% group_by(year, cell) %>% tally())
nrow(pheno.dev)


#Environmental data
#env.var<-read_csv("data/derived/envir.csv")
env.var<-read_csv("data/derived/envir2.csv") %>%
  mutate(summer.gdd=summer.gdd/100)

## Rescale summer gdd and calculate devisions
env.dev<-env.var %>%
  mutate(include=ifelse(year>2015,1,NA),baselinepc1=pc1*include, baselinepc2=pc2*include, basecold=colddays*include, basewarm=warmdays*include, basesummer=summer.gdd*include, baselag=gr_mn_lag*include) %>%
  group_by(cell) %>%
  mutate( pc1.dev=pc1-mean(baselinepc1, na.rm=T),pc2.dev=pc2-mean(baselinepc2, na.rm=T),
          cold.dev=colddays-mean(basecold, na.rm=T), warm.dev=warmdays-mean(basewarm, na.rm=T),summer.dev=summer.gdd-mean(basesummer, na.rm=T), lag.dev=gr_mn_lag-mean(baselag, na.rm=T)) %>%
  dplyr::select(year, cell, cell_lat, cell_lng, cellyr,pc1.dev, pc2.dev, open.lag=gr_mn_lag, lag.dev,cold.dev,warm.dev,summer.dev)

#env.dev<-env.var %>%
#  group_by(cell) %>%
#  mutate( pc1.dev=warmearly-mean(warmearly, na.rm=T), pc2.dev=warmlateopen-mean(warmlateopen, na.rm=T))


pheno.input.dev<-merge(pheno.dev, env.dev, by=c("year","cell"), all.x=T) %>%
  filter(code %in% c("RE","RL","RP"), year<2018) %>%
  mutate(year=year-2000, onset.ci.wt=1/(onset.ci+1),med.ci.wt=1/(median.ci+1),dur.ci.wt=1/(dur.ci+1))

ggplot(data=pheno.input.dev, aes(x=year, y=onset.dev, color=code)) + geom_jitter()
ggplot(data=pheno.input.dev, aes(x=year, y=median.dev, color=code)) + geom_jitter()
ggplot(data=pheno.input.dev, aes(x=year, y=dur.dev, color=code)) + geom_jitter()

meanon<-pheno.quant %>% group_by(year, code) %>% summarize(meanon=mean(q5, na.rm=T))
ggplot(data=pheno.input.dev, aes(x=year, y=onset.dev, color=code)) + geom_smooth(method="lm") + labs(x="Year", y="Onset deviation")
ggplot(data=pheno.quant, aes(x=year, y=q5, color=code)) + geom_smooth(method="lm") + labs(x="Year", y="Onset")
ggplot(data=pheno.quant, aes(x=year, y=q5, color=code)) + geom_smooth(method="lm") + labs(x="Year", y="Onset") + geom_line(data=meanon, aes(x=year, y=meanon))

pheno.input<-merge(pheno.quant, env.var, by=c("year","cell")) %>%
  filter(code %in% c("RE","RL","RP"), year<2018) %>%
  mutate(year=year-2000, q5.ci.wt=1/(q5_ci+1),q50.ci.wt=1/(q50_ci+1),qdur.ci.wt=1/(qdur_ci+1))


mean.on<-pheno.input %>% group_by(code) %>% summarize(meanon=mean(q5, na.rm=T),meandur=mean(qdur, na.rm=T), nmet=n())
vardata<-pheno.input %>% group_by(code) %>% summarize(sdmed=sd(q50, na.rm=T), sddur=sd(qdur, na.rm=T))

save(pheno.input, pheno.input.dev, file="data/derived/model.input.RData")
save(env.var, env.dev, file="data/derived/env.input.RData")

### BASE for visualizations

newcode<-c("RE","RL","RP")
codelabel<-c("BOE","BOL","BOP")
nrep<-1200

newDat <- data.frame(cell = rep(507,nrep), 
                     uniqObsDays=rep(median(pheno.input$uniqObsDays, na.rm=T),nrep),
                     warmdays=rep(median(pheno.input$warmdays, na.rm=T),nrep),
                     gr_mn_lag=rep(median(pheno.input$gr_mn_lag, na.rm=T),nrep),
                     pc1=rep(median(pheno.input$pc1, na.rm=T),nrep),
                     pc2=rep(median(pheno.input$pc2, na.rm=T),nrep),
                     year=rep(10, nrep),
                     code=rep(newcode, each=20,20),
                     codelabel=rep(codelabel, each=20,20)) %>%
  mutate(codelabel = factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel, as.numeric(code)))

maxvals<-apply(na.omit(pheno.input[,c(32,34:40)]),2,FUN=max)
minvals<-apply(na.omit(pheno.input[,c(32,34:40)]),2,FUN=min)
intervals<-round((maxvals-minvals)/20,3)

dev.newDat <- data.frame(cell = rep(507,nrep), 
                     uniqObsDays=rep(median(pheno.input.dev$uniqObsDays, na.rm=T),nrep),
                     warm.dev=rep(0,nrep),
                     lag.dev=rep(0,nrep),
                     pc1.dev=rep(0,nrep),
                     pc2.dev=rep(0,nrep),
                     year=rep(10, nrep),
                     code=rep(newcode, each=20,20),
                     codelabel=rep(codelabel, each=20,20)) %>%
  mutate(codelabel = factor(codelabel),code=factor(code)) #%>%
 #mutate(codelabel=fct_reorder(codelabel, as.numeric(code)))

dev.maxvals<-apply(na.omit(pheno.input.dev[,c(1,4,8,11,16:22)]),2,FUN=max)
dev.minvals<-apply(na.omit(pheno.input.dev[,c(1,4,8,11,16:22)]),2,FUN=min)
dev.intervals<-round((dev.maxvals-dev.minvals)/20,3)







################################################################################
#### PHENOLOGICAL MODELS
################################################################################

### Analysis
#### WEIGHTED LMER MODELS
library(MASS)

################################################################################
#### PHENOLOGICAL MODELS: ONSET
################################################################################
onset.full<-lmer(q5~-1+code+code*(warmdays+pc1+pc2+ year + gr_mn_lag + uniqObsDays) + (1|cell), data=pheno.input, weights=q5.ci.wt)

extractAIC(onset.full)
(onset.step <- step(onset.full))
onset.best <- get_model(onset.step) #stepAIC(ab.yr.full)
summary(onset.best)
onset.final<-lmer(q5~-1+code+code:warmdays +code:pc1 + code:year + code:gr_mn_lag + uniqObsDays + (1|cell), data=pheno.input, weights=q5.ci.wt)
r.squaredGLMM(onset.final)
summary(onset.final)
extractAIC(onset.final)
onset.pc1r2<-lmer(q5~-1+code+code:warmdays +code:year + code:gr_mn_lag + uniqObsDays + (1|cell), data=pheno.input, weights=q5.ci.wt)
r.squaredGLMM(onset.final)[1]-r.squaredGLMM(onset.pc1r2)[1]
onset.yrr2<-lmer(q5~-1+code+code:warmdays +code:pc1 + code:gr_mn_lag + uniqObsDays + (1|cell), data=pheno.input, weights=q5.ci.wt)
r.squaredGLMM(onset.final)[1]-r.squaredGLMM(onset.yrr2)[1]
onset.openlagr2<-lmer(q5~-1+code+code:warmdays +code:year + code:pc1 + uniqObsDays + (1|cell), data=pheno.input, weights=q5.ci.wt)
r.squaredGLMM(onset.final)[1]-r.squaredGLMM(onset.openlagr2)[1]
onset.wwdr2<-lmer(q5~-1+code+code:pc1 +code:year + code:gr_mn_lag + uniqObsDays + (1|cell), data=pheno.input, weights=q5.ci.wt)
r.squaredGLMM(onset.final)[1]-r.squaredGLMM(onset.wwdr2)[1]

### 

#require(nlme)
#fm <- lmer(q5 ~code + code:[], random = ~ 1|Subject,data=Orthodont,method='ML')
#Terms <- terms(fm)
#todrop <-  1:length(attr(Terms,'term.labels'))
#subs <- unlist(sapply(todrop,function(p)
#  combn(todrop,p,simplify=F)),recursive =F)
#fm.subList <- lapply(subs[-length(subs)],function(s,...){
#  newf<- formula(drop.terms(terms(fm),s,keep.response = TRUE))
#  update(fm,newf)
#})
#names(fm.subList) <- sapply(fm.subList, function(x) paste('fm',attr(
#  terms(x),'term.labels'),sep='.'))
#sort(sapply(fm.subList,BIC))   

#best.subsets<-regsubsets(q5 ~ ., data = model.input[c(1:12),],nbest=1,nvmax=NULL, method="exhaustive")
#summary(tx.best.subset)
#test1<-summary(tx.best.subset)
#vars<-names(model.input[c(1:12),])[-1]

#full.aic<-lm(avg10~., data=model.input)


onset.vif<-lmer(q5~code+warmdays + pc1 + gr_mn_lag + year + uniqObsDays + (1|cell), data=pheno.input, weights=q5.ci.wt)
vif(onset.vif)

plot_model(onset.final, type = "eff", terms = c("year", "code"), title="Onset~Year")
plot_model(onset.final, type = "eff", terms = c("warmdays", "code"), title="Onset~warm winter days")
plot_model(onset.final, type = "eff", terms = c("pc1", "code"), title="Onset~PC1")
plot_model(onset.final, type = "eff", terms = c("gr_mn_lag", "code"), title="Onset~open canopy lag")

onset.output<-as_tibble(summary(onset.final)$coefficients) %>%
  dplyr::mutate(param=row.names(summary(onset.final)$coefficients), sig=ifelse(`Pr(>|t|)`>0.05,0,1), r2m=round(r.squaredGLMM(onset.final)[1],2),r2c=round(r.squaredGLMM(onset.final)[2],2))
write.csv(onset.output, file="output/onset.model0921.csv")
#write.csv(summary(onset.final)$coefficients, file="output/onset0921.csv")


######################################################
##ONSET VISUALIZATION

#SPATIAL PRED: WARM WINTER DAYS, OPEN CANOPY LAG, PC1, YEAR
new.warmdays<-minvals["warmdays"]+(0:19)*intervals["warmdays"]
new.lags<-minvals["gr_mn_lag"]+(0:19)*intervals["gr_mn_lag"]
new.pc1<-minvals["pc1"]+(0:19)*intervals["pc1"]
new.yr<-c(0:19)

remove.cols<-which(names(newDat) %in% c("warmdays","pc1"))
newdat.clim<-cbind(newDat[,-remove.cols], warmdays=rep(new.warmdays, 60), pc1=rep(new.pc1, each=60))
newdat.clim$pred<-predict(onset.final, newdat.clim,allow.new.levels =T)

lims.2<-pheno.input %>%
  group_by(code) %>%
  summarize(minwd=min(warmdays, na.rm=T),maxwd=max(warmdays, na.rm=T),
            minpc1=min(pc1, na.rm=T),maxpc1=max(pc1, na.rm=T))

newDat2<-inner_join(newdat.clim,lims.2) %>%
  mutate(f1=ifelse(warmdays>=minwd,ifelse(warmdays<=maxwd,1,0),0)+ifelse(pc1>=minpc1,ifelse(pc1<=maxpc1,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel = factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel, as.numeric(code)))


F.onset.spat1<-ggplot(data=newDat2, aes(x=warmdays, y=pc1, fill=pred)) + 
  geom_tile() +  labs(y="PC1", x="Winter Warm Days") + 
  scale_fill_viridis(name="Onset") + 
  facet_wrap(~codelabel)
F.onset.spat1
save(F.onset.spat1,file="output/onset.fig.1.png")


## pc1 x year

remove.cols<-which(names(newDat) %in% c("pc1","year","warmdays","gr_mn_lag"))
newdat.clim<-cbind(newDat[,-remove.cols], year=rep(new.yr, 60), pc1=rep(new.pc1, each=60), warmdays=rep(75,nrep),gr_mn_lag=rep(40,nrep))
newdat.clim$pred<-predict(onset.final, newdat.clim,allow.new.levels =T)

lims.on2<-pheno.input %>%
  group_by(code) %>%
  summarize(minyr=min(year, na.rm=T),maxyr=max(year, na.rm=T),
            minpc1=min(pc1, na.rm=T),maxpc1=max(pc1, na.rm=T))

newDat.on2<-inner_join(newdat.clim,lims.on2) %>%
  mutate(f1=ifelse(year>=minyr,ifelse(year<=maxyr,1,0),0)+ifelse(pc1>=minpc1,ifelse(pc1<=maxpc1,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel = factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel, as.numeric(code)))


F.onset.spat2<-ggplot(data=newDat.on2, aes(x=year, y=pc1, fill=pred)) + 
  geom_tile() +  labs(y="PC1", x="Year") + 
  #scale_fill_viridis(name="Onset") + 
  scale_fill_gradient2() +  
  facet_wrap(~codelabel)
F.onset.spat2
save(F.onset.spat2,file="output/onset.pc1.year.fig.png")


################################################################################
#### PHENOLOGICAL MODELS: ONSET DEVIATION
################################################################################

onset.dev.full<-lmer(onset.dev~-1+code+warm.dev + pc1.dev + pc2.dev + lag.dev + year+ uniqObsDays + code:warm.dev + code:lag.dev + code:pc1.dev + code:pc2.dev + code:year  + code:uniqObsDays + (1|cell), data=pheno.input.dev, weights=onset.ci.wt)
extractAIC(onset.dev.full)
r.squaredGLMM(onset.dev.full)     
(t2<-as_tibble(bind_cols(Parameters=row.names(summary(onset.dev.full)$coefficients),summary(onset.dev.full)$coefficients)) %>% arrange(abs(`t value`)))

onset.dev.best<-get_model(step(onset.dev.full))

(t2<-as_tibble(bind_cols(Parameters=row.names(summary(t1)$coefficients),summary(t1)$coefficients)) %>% arrange(abs(`t value`)))
summary(onset.dev.best)
summary(onset.dev.best)$call

onset.dev.final<-lmer(onset.dev~ -1 + code + uniqObsDays + code:pc1.dev + code:pc2.dev + code:year +(1|cell), data=pheno.input.dev, weights=onset.ci.wt)
summary(onset.dev.final)
extractAIC(onset.dev.final)
r.squaredGLMM(onset.dev.final)     
onset.dev.pc1<-lmer(onset.dev~ -1 + code + uniqObsDays +code:pc2.dev  + code:year +(1|cell), data=pheno.input.dev, weights=onset.ci.wt)
r.squaredGLMM(onset.dev.final)[1]-r.squaredGLMM(onset.dev.pc1)[1]
onset.dev.pc2<-lmer(onset.dev~ -1 + code + uniqObsDays + code:pc1.dev + code:year +(1|cell), data=pheno.input.dev, weights=onset.ci.wt)
r.squaredGLMM(onset.dev.final)[1]-r.squaredGLMM(onset.dev.pc2)[1]
onset.dev.yr<-lmer(onset.dev~ -1 + code + uniqObsDays + code:pc1.dev + code:pc2.dev +(1|cell), data=pheno.input.dev, weights=onset.ci.wt)
r.squaredGLMM(onset.dev.final)[1]-r.squaredGLMM(onset.dev.yr)[1]
onset.dev.eff<-lmer(onset.dev~ -1 + code +  code:pc1.dev + code:pc2.dev +  code:year +(1|cell), data=pheno.input.dev, weights=onset.ci.wt)
r.squaredGLMM(onset.dev.final)[1]-r.squaredGLMM(onset.dev.eff)[1]

test.vif<-lm(onset.dev~code+year+pc1.dev+pc2.dev+uniqObsDays, data=pheno.input.dev)
vif(test.vif)

plot_model(onset.dev.final, type = "eff", terms = c("year", "code"), title="OnDev~Year")
plot_model(onset.dev.final, type = "eff", terms = c("pc2.dev", "code"), title="OnDev~PC2 dev")
plot_model(onset.dev.final, type = "eff", terms = c("pc1.dev", "code"), title="OnDev~PC1 dev")

onsetdev.output<-as_tibble(summary(onset.dev.final)$coefficients) %>%
  dplyr::mutate(param=row.names(summary(onset.dev.final)$coefficients), sig=ifelse(`Pr(>|t|)`>0.05,0,1), r2m=round(r.squaredGLMM(onset.dev.final)[1],2),r2c=round(r.squaredGLMM(onset.dev.final)[2],2))
write.csv(onsetdev.output, file="output/onset.dev.model0921.csv")

##################################
#ONSET DEV PLOT: PC1, PC2, Year are interesting
#predict
##
new.pc1<-dev.minvals["pc1.dev"]+(0:19)*dev.intervals["pc1.dev"]
new.pc2<-dev.minvals["pc2.dev"]+(0:19)*dev.intervals["pc2.dev"]
#new.pc1<-minvals["pc1"]+(0:19)*intervals["pc1"]
new.yr<-c(0:19)

remove.cols<-which(names(dev.newDat) %in% c("pc1.dev","pc2.dev"))
newdat.od<-cbind(dev.newDat[,-remove.cols], pc1.dev=rep(new.pc1, 60), pc2.dev=rep(new.pc2, each=60))
newdat.od$pred <- predict(onset.dev.final, newdat.od, allow.new.levels =T)


lims.dev<-pheno.input.dev %>%
  group_by(code) %>%
  summarize(minyr=min(year, na.rm=T),maxyr=max(year, na.rm=T),
            minpc1=min(pc1.dev, na.rm=T),maxpc1=max(pc1.dev, na.rm=T),
            minpc2=min(pc2.dev, na.rm=T), maxpc2=max(pc2.dev, na.rm=T))

newDat.dev<-inner_join(newdat.od,lims.dev) %>%
  mutate(f1=ifelse(pc2.dev>=minpc2,ifelse(pc2.dev<=maxpc2,1,0),0)+ifelse(pc1.dev>=minpc1,ifelse(pc1.dev<=maxpc1,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel = factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel, as.numeric(code)))


##for fig 
F.ondev.1<-ggplot(data=newDat.dev, aes(x=pc2.dev, y=pc1.dev, fill=pred)) + 
  geom_tile() +  labs(x="PC2 Deviation", y="PC1 Deviation", fill="Onset Deviation") + 
  scale_fill_gradient2() +  
  facet_wrap(~codelabel)
F.ondev.1
save(F.ondev.1,file="output/onset.dev.fig.1.png")


#next fig
remove.cols<-which(names(dev.newDat) %in% c("pc1.dev","year"))
newdat.od2<-cbind(dev.newDat[,-remove.cols], pc1.dev=rep(new.pc1, 60), year=rep(new.yr, each=60))
newdat.od2$pred <- predict(onset.dev.final, newdat.od2, allow.new.levels =T)


F.ondev.2<-ggplot(data=filter(newdat.od2, year<18), aes(x=year+2000, y=pc1.dev, fill=pred)) + 
  geom_tile() +  labs(y="PC1 Deviation", x="Year", title="Onset deviation") + 
  #scale_fill_gradient2(name="Predicted dev.") +  
  scale_fill_viridis(name="Predicted dev.") + 
  facet_wrap(~codelabel,nrow=3) + theme(legend.position="bottom")
F.ondev.2
#save(F.ondev.2,file="output/ondev.pc1.year.fig.png")

F.onset.spat2<-ggplot(data=newDat.on2, aes(x=year+2000, y=pc1, fill=pred)) + 
  geom_tile() +  labs(y="PC1", x="Year", title="Adult onset") + 
  #scale_fill_gradient2(name="Prediction") +
  #scale_fill_viridis(name="Prediction") + 
  facet_wrap(~codelabel, nrow=3) + theme(legend.position="bottom")
F.onset.spat2


library(gridExtra)
(figonset<-grid.arrange(F.onset.spat2,F.ondev.2, nrow=1))
#ggsave(figonset, width=6,height=7, units="in", file="output/fig.onsetdiverg.pdf")
#ggsave(figonset, width=6,height=7, units="in", file="output/fig.onsetdiverg.png")
ggsave(figonset, width=6,height=7, units="in", file="output/fig.onset.pdf")
ggsave(figonset, width=6,height=7, units="in", file="output/fig.onset.png")

(pheno.dev.power<-pheno.input.dev %>% group_by(code, cell) %>% summarize(n=1) %>% group_by(code) %>% tally())


###########################################################################
### PHENOLOGICAL MODELS: MEDIAN
###########################################################################
#ph2<-pheno.input
#pheno.input<-ph2 %>% mutate(summer.gdd=summer.gdd/500)
med.full<-lmer(q50~-1+code+code*(warmdays+pc1+pc2+summer.gdd+year + gr_mn_lag + uniqObsDays) + (1|cell), data=pheno.input, weights=q50.ci.wt)
(med.step <- step(med.full))
med.best <- get_model(med.step) #stepAIC(ab.yr.full)
summary(med.best)
med.final<-lmer(q50~-1+code+warmdays +pc1 + gr_mn_lag + code:summer.gdd +  (1|cell), data=pheno.input, weights=q50.ci.wt)
r.squaredGLMM(med.final)
summary(med.final)
extractAIC(med.final)

med.vif<-lmer(q50~code+warmdays + pc1 + gr_mn_lag + summer.gdd + (1|cell), data=pheno.input, weights=q50.ci.wt)
vif(med.vif)

plot_model(med.final, type = "eff", terms = c("warmdays", "code"), title="med~warm winter days")
plot_model(med.final, type = "eff", terms = c("pc1", "code"), title="med~PC1")
plot_model(med.final, type = "eff", terms = c("gr_mn_lag", "code"), title="med~open canopy lag")
plot_model(med.final, type = "eff", terms = c("summer.gdd", "code"), title="med~summer.gdd")

med.output<-as_tibble(summary(med.final)$coefficients) %>%
  dplyr::mutate(param=row.names(summary(med.final)$coefficients), sig=ifelse(`Pr(>|t|)`<0.05,1,0), r2m=round(r.squaredGLMM(med.final)[1],2),r2c=round(r.squaredGLMM(med.final)[2],2))
write.csv(med.output, file="output/med.model0921.csv")
#write.csv(summary(med.final)$coefficients, file="output/med0921.csv")


#SPATIAL PRED: WARM WINTER DAYS, SUMMER GDD
newcode<-c("RE","RL","RP")
codelabel<-c("egg","caterpillar","pupa")
newffd<-round(c( ((round(min(pheno.input$warmdays, na.rm=T)))):((round(max(pheno.input$warmdays, na.rm=T))))),2)
newgdd<-round(c( ((round(min(pheno.input$summer.gdd, na.rm=T)*50))):((round(max(pheno.input$summer.gdd, na.rm=T)*50)))),2)/50
newyr<-c(10)
newlag<-20
nrep<-length(newcode)*length(newgdd)*length(newffd)
newDatm <- data.frame(cell = rep(622,nrep), 
                       uniqObsDays=rep(median(pheno.input$uniqObsDays, na.rm=T),nrep),
                       gr_mn_lag=rep(median(pheno.input$gr_mn_lag, na.rm=T),nrep),
                       pc1=rep(median(pheno.input$pc1, na.rm=T),nrep),
                       pc2=rep(median(pheno.input$pc2, na.rm=T),nrep),
                       year=rep(newyr,nrep),
                       warmdays=rep(newffd,length(newgdd)*length(newcode)), ##
                       summer.gdd=rep(newgdd,each=length(newffd),length(newcode)), ##
                       code=rep(newcode, each=length(newgdd)*length(newffd)),
                       codelabel=rep(codelabel, each=length(newgdd)*length(newffd)))

newDatm$pred <- predict(med.final, newDatm,allow.new.levels =T)

lims.2<-pheno.input %>%
  mutate(codelabel=ifelse(code=="RE","egg",ifelse(code=="RL","caterpillar","pupa"))) %>% 
  group_by(code, codelabel) %>%
  summarize(miny=floor(min(summer.gdd, na.rm=T)*5)/5,maxy=ceiling(max(summer.gdd, na.rm=T)*5)/5,
            minx=floor(min(warmdays, na.rm=T)*5)/5,maxx=ceiling(max(warmdays, na.rm=T)*5)/5)

newDatm2<-inner_join(newDatm,lims.2) %>%
  mutate(f1=ifelse(warmdays>=minx,ifelse(warmdays<=maxx,1,0),0)+ifelse(summer.gdd>=miny,ifelse(summer.gdd<=maxy,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel = factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel, as.numeric(code)), summer.gdd=summer.gdd*500)


(F.med.spat1<-ggplot(data=newDatm2, aes(y=warmdays, x=summer.gdd, fill=pred)) + 
  geom_tile() +  labs(x="Summer GDD", y="Winter Warm Days") + 
  scale_fill_viridis(name="med") + 
  facet_wrap(~codelabel) )

save(F.med.spat1,file="output/med.fig.1.png")

###########################################################################
### PHENOLOGICAL MODELS: MEDIAN DEVIATION
###########################################################################


med.dev.full<-lmer(median.dev~-1+code+warm.dev + summer.dev + pc1.dev + pc2.dev + lag.dev + year+ uniqObsDays + code:warm.dev + code:summer.dev +code:lag.dev + code:pc1.dev + code:pc2.dev + code:year + uniqObsDays + code:uniqObsDays + (1|cell), data=pheno.input.dev, weights=med.ci.wt)
extractAIC(med.dev.full)
r.squaredGLMM(med.dev.full)     
(t2<-as_tibble(bind_cols(Parameters=row.names(summary(med.dev.full)$coefficients),summary(med.dev.full)$coefficients)) %>% arrange(abs(`t value`)))

t1<-get_model(step(med.dev.full))

(t2<-as_tibble(bind_cols(Parameters=row.names(summary(t1)$coefficients),summary(t1)$coefficients)) %>% arrange(abs(`t value`)))
summary(t1)$call
med.dev.final<-lmer(median.dev~ -1 + code + pc1.dev + warm.dev + code:lag.dev +(1|cell), data=pheno.input.dev, weights=med.ci.wt)
summary(med.dev.final)
extractAIC(med.dev.final)
r.squaredGLMM(med.dev.final)     
test.vif<-lm(median.dev~code+pc1.dev+warm.dev+lag.dev, data=pheno.input.dev)
vif(test.vif)

plot_model(med.dev.final, type = "eff", terms = c("warm.dev", "code"), title="MedDev~Lag winter warm days")
plot_model(med.dev.final, type = "eff", terms = c("lag.dev", "code"), title="MedDev~Lag dev")
plot_model(med.dev.final, type = "eff", terms = c("pc1.dev", "code"), title="MedDev~PC1 dev")

meddev.output<-as_tibble(summary(med.dev.final)$coefficients) %>%
  dplyr::mutate(param=row.names(summary(med.dev.final)$coefficients), sig=ifelse(`Pr(>|t|)`<0.05,1,0), r2m=round(r.squaredGLMM(med.dev.final)[1],2),r2c=round(r.squaredGLMM(med.dev.final)[2],2))
write.csv(meddev.output, file="output/med.dev.model0921.csv")




###########################################################################
### PHENOLOGICAL MODELS: DURATION
###########################################################################
#ph2<-pheno.input
#pheno.input<-ph2 %>% mutate(summer.gdd=summer.gdd/500)
dur.full<-lmer(qdur~-1+code+code*(warmdays+pc1+pc2+summer.gdd+year + gr_mn_lag + uniqObsDays) + (1|cell), data=pheno.input, weights=qdur.ci.wt)
(dur.step <- step(dur.full))
dur.best <- get_model(dur.step) #stepAIC(ab.yr.full)
summary(dur.best)
dur.final<-lmer(qdur~-1+code+pc2+code:warmdays + code:gr_mn_lag + code:summer.gdd + code:year +code:uniqObsDays+ (1|cell), data=pheno.input, weights=qdur.ci.wt)
r.squaredGLMM(dur.final)
summary(dur.final)
extractAIC(dur.final)

dur.vif<-lmer(q50~code+warmdays + pc2 + uniqObsDays + warmdays + year + gr_mn_lag + summer.gdd + (1|cell), data=pheno.input, weights=q50.ci.wt)
vif(dur.vif)

plot_model(dur.final, type = "eff", terms = c("warmdays", "code"), title="dur~warm winter days")
plot_model(dur.final, type = "eff", terms = c("year", "code"), title="dur~year")
plot_model(dur.final, type = "eff", terms = c("gr_mn_lag", "code"), title="dur~open canopy lag")

plot_model(dur.final, type = "eff", terms = c("summer.gdd", "code"), title="dur~summer.gdd")

dur.output<-as_tibble(summary(dur.final)$coefficients) %>%
  dplyr::mutate(param=row.names(summary(dur.final)$coefficients), sig=ifelse(`Pr(>|t|)`<0.05,1,0), r2m=round(r.squaredGLMM(dur.final)[1],2),r2c=round(r.squaredGLMM(dur.final)[2],2))
write.csv(dur.output, file="output/dur.model0921.csv")
#write.csv(summary(dur.final)$coefficients, file="output/dur0921.csv")



###########################################################################
### PHENOLOGICAL MODELS: DURATION DEVIATION
###########################################################################


dur.dev.full<-lmer(dur.dev~-1+code+ code*(warm.dev + summer.dev + pc1.dev + pc2.dev + lag.dev + year+ uniqObsDays) + (1|cell),
                   data=pheno.input.dev, weights=dur.ci.wt)
extractAIC(dur.dev.full)
r.squaredGLMM(dur.dev.full)     
(t2<-as_tibble(bind_cols(Parameters=row.names(summary(dur.dev.full)$coefficients),summary(dur.dev.full)$coefficients)) %>% arrange(abs(`t value`)))

dur.best<-get_model(step(dur.dev.full))

(t2<-as_tibble(bind_cols(Parameters=row.names(summary(dur.best)$coefficients),summary(dur.best)$coefficients)) %>% arrange(abs(`t value`)))
summary(dur.best)$call
dur.dev.final<-lmer(dur.dev~ -1 + code + warm.dev + code:summer.dev + code:pc1.dev + code:pc2.dev + code:lag.dev + code:year + code:uniqObsDays +  (1|cell), data=pheno.input.dev, weights=dur.ci.wt)
summary(dur.dev.final)
extractAIC(dur.dev.final)
r.squaredGLMM(dur.dev.final)     
test.vif<-lm(dur.dev~code+pc1.dev+pc2.dev+warm.dev+lag.dev+summer.dev+year+uniqObsDays, data=pheno.input.dev)
vif(test.vif)

plot_model(dur.dev.final, type = "eff", terms = c("warm.dev", "code"), title="durDev~Lag winter warm days")
plot_model(dur.dev.final, type = "eff", terms = c("lag.dev", "code"), title="durDev~Lag dev")
plot_model(dur.dev.final, type = "eff", terms = c("pc1.dev", "code"), title="durDev~PC1 dev")

durdev.output<-as_tibble(summary(dur.dev.final)$coefficients) %>%
  dplyr::mutate(param=row.names(summary(dur.dev.final)$coefficients), sig=ifelse(`Pr(>|t|)`<0.05,1,0), r2m=round(r.squaredGLMM(dur.dev.final)[1],2),r2c=round(r.squaredGLMM(dur.dev.final)[2],2))
write.csv(durdev.output, file="output/dur.dev.model0921.csv")




########################################################
##  DATA FIGS
#######################################################
ggplot(data=pheno.input, aes(x=year, y=qdur,color=code)) + geom_point()

ggplot(data=pheno.input, aes(y=qdur,x=code)) + geom_boxplot()

ggplot(data=pheno.input, aes(y=q5,x=code)) + geom_boxplot()
pheno.input<-pheno.input %>% mutate(dec=as.factor(ifelse(year>9,"yes","no")),
                                    codelabels=ifelse(code=="RE","Egg",ifelse(code=="RL","Caterpillar","Pupa")))

ggplot(data=filter(pheno.input, year<10), aes(y=q5,x=qdur, color=cell_lat, shape=dec, size=uniqObsDays)) + geom_point() + 
  scale_shape_manual(name=">2009", values=c(1,15)) + 
  labs(x="Duration", y="Adult Onset") + 
  scale_size_continuous(guide="none") + 
  facet_wrap(~fct_reorder(codelabels,as.numeric(as.factor(code))))

### 
Sys.time()

























##MEDIAN DEVIATION

med.dev.full<-lmer(median.dev~-1+code*(FFD.dev + pc1.dev + pc2.dev + summer.dev + year)+ uniqObsDays+ (1|cell), data=pheno.input.dev, weights=med.ci.wt)
r.squaredGLMM(med.dev.full)     
summary(med.dev.full)
med.dev.step <- step(med.dev.full)
med.dev.best <- get_model(med.dev.step) #stepAIC(ab.yr.full)
summary(med.dev.best)
r.squaredGLMM(med.dev.best)     
plot_model(med.dev.best)
med.dev.final<-lmer(median.dev~-1+code + pc1.dev + code:pc2.dev + summer.dev +  (1|cell), data=pheno.input.dev, weights=med.ci.wt)
med.dev.vif<-lmer(median.dev~code + pc1.dev + pc2.dev + summer.dev + (1|cell), data=pheno.input.dev, weights=med.ci.wt)
vif(med.dev.vif)
plot_model(med.dev.final, type = "eff", terms = c("pc2.dev", "code"), title="Onset~PC2")


#med
med.full<-lmer(q50~code*(FFD + warmearly + warmlateopen + summer.dev  + year) + (1|cell), data=pheno.input, weights=q50.ci.wt)
(med.step <- step(med.full))
med.best <- get_model(med.step) #stepAIC(ab.yr.full)
summary(med.best)
med.vif<-lmer(q50 ~ code + FFD + warmlateopen + summer.dev + (1|cell), data=pheno.input, weights=q50.ci.wt)
vif(med.vif)
r.squaredGLMM(med.best)     
plot_model(med.best, type = "eff", terms = c("warmlateopen","summer.dev"), title="Median")
plot_model(med.best, type = "eff", terms = c("FFD","code"), title="Median")


##DURATION DEVIATION
dur.dev.full<-lmer(dur.dev~-1+code*(FFD.dev + pc1.dev + pc2.dev + summer.dev + year) + uniqObsDays + (1|cell), data=pheno.input.dev, weights=dur.ci.wt)
r.squaredGLMM(dur.dev.full)     
summary(dur.dev.full)
dur.dev.step <- step(dur.dev.full)
dur.dev.best <- get_model(dur.dev.step) #stepAIC(ab.yr.full)
summary(dur.dev.best)
dur.dev.final<-lmer(dur.dev~-1 + code + pc1.dev + code:pc2.dev + code:year + code:summer.dev + (1|cell), data=pheno.input.dev, weights=dur.ci.wt)
r.squaredGLMM(dur.dev.final)     
summary(dur.dev.final)
extractAIC(dur.dev.final)     
dur.dev.vif<-lmer(dur.dev~code + pc1.dev + pc2.dev + year + summer.dev + (1|cell), data=pheno.input.dev, weights=dur.ci.wt)
vif(dur.dev.vif)
plot_model(dur.dev.final, type = "eff", terms = c("pc2.dev", "code"), title="Duration deviation~ pc2 deviation")
plot_model(dur.dev.final, type = "eff", terms = c("year", "code"), title="Duration deviation~ year")
plot_model(dur.dev.final, type = "eff", terms = c("summer.dev", "code"), title="Duration deviation~ summer deviation")

#DURATION
dur.full<-lmer(qdur~-1+code*(FFD + warmearly + warmlateopen + summer.dev  + year) + (1|cell), data=pheno.input, weights=qdur.ci.wt)
(dur.step <- step(dur.full))
dur.best <- get_model(dur.step) #stepAIC(ab.yr.full)
summary(dur.best)
r.squaredGLMM(dur.best)     
dur.vif<-lmer(qdur~code + FFD + warmearly + warmlateopen + year + (1|cell), data=pheno.input, weights=qdur.ci.wt)
vif(dur.vif)
dur.final<-lmer(qdur~-1+code:FFD + code:warmearly + warmlateopen + code:year + (1|cell), data=pheno.input, weights=qdur.ci.wt)
plot_model(dur.final, type="eff", terms=c("year","code"))
plot_model(dur.final, type="eff", terms=c("warmearly","code"))
plot_model(dur.final, type="eff", terms=c("FFD","code"))




#all plot_models
models<-list(onset.best, onset.dev.best, med.best, med.dev.best, dur.best, dur.dev.best)

plot.panel<-function(model.i) {
  coefs<-summary(model.i)$coefficients
  min.row<-which(coefs[,1]==min(coefs[,1]))
  ylab<-coefs[min.row,1]-2*(coefs[min.row,2])
  plot_model(model.i, show.values=T, show.p=T) + 
    theme_bw() + geom_hline(yintercept=0) + 
    #scale_color_manual(values=thiscolor[c(2,4)]) + 
    labs(subtitle = paste("marginal r-sq.=",round(r.squaredGLMM(model.i),2)))
}

plot.panel(onset.best)
fixedlabels<-lapply(models, label.sig)

plots<-lapply(models, plot.panel)
pdf("output/pheno.models1.pdf")
grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]], nrow=3)
dev.off()

## Getting interaction coefficients
onset.best<-lmer(onset~-1+code + warmearly + code:year + (1|cell), data=pheno.input)
onset.final<-get_model(step(onset.best))
r.squaredGLMM(onset.final)     

ondev.best<-lmer(dev_onset~-1+code + warmearly + code:year + (1|cell), data=pheno.input.dev)
ondev.final<-get_model(step(ondev.best))
r.squaredGLMM(ondev.final)     

med.best<-lmer(fiftieth~-1+code + warmlateopen + summer.dev + (1|cell), data=pheno.input)
med.final<-get_model(step(med.best))
r.squaredGLMM(med.final)     

meddev.best<-lmer(dev_median~-1+code +warmlateopen + summer.dev + (1|cell), data=pheno.input.dev)
meddev.final<-get_model(step(meddev.best))
r.squaredGLMM(meddev.final)     

dur.best<-lmer(duration~-1+code + warmearly + code:warmlateopen + code:FFD + year + (1|cell), data=pheno.input)
dur.final<-get_model(step(dur.best))
r.squaredGLMM(dur.final)     

durdev.best<-lmer(dur.dev~-1+code + warmearly + year + (1|cell), data=pheno.input.dev)
durdev.final<-get_model(step(durdev.best))
r.squaredGLMM(durdev.final)     


finalmodels<-list(onset.final, ondev.final, med.final, meddev.final, dur.final, durdev.final)

plots<-lapply(finalmodels, plot.panel)
plots[[3]]<-plots[[3]] + scale_color_manual(values=thiscolor[c(4,2)])
plots[[5]]<-plots[[5]] + scale_color_manual(values=thiscolor[c(4,2)])
pdf("output/pheno.models2.pdf")
grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]], nrow=3)
dev.off()



### PLOT PHENOLOGY RESULTS

temp1<-predict(onset.dev.best) #, data=pheno.input.dev)
pheno1<-pheno.input.dev %>%
  dplyr::select(code,warmearly,year,dev_onset) 
pheno1<-na.omit(pheno1) %>% dplyr::mutate(fit=temp1)
(plot.ondev1<-ggplot(data=pheno1, aes(x=year, y=dev_onset, color=code)) + 
    geom_jitter() + xlim(0,20) + 
    geom_smooth(method="lm")) + theme_classic() + scale_color_viridis_d() + labs(x="Year", y="Onset Deviation")


(plot.onset1<-ggplot(data=filter(pheno.input), aes(x=year, y=onset, color=code)) + 
    geom_smooth(method="lm")) + theme_classic() + scale_color_viridis_d() + labs(x="Year", y="Onset")
(plot.onset1<-ggplot(data=filter(pheno.input.dev), aes(x=year, y=dev_onset, color=code)) + 
    geom_smooth(method="lm")) + theme_classic() + scale_color_viridis_d() + labs(x="Year", y="Onset Deviation")

(plot.onset1<-ggplot(data=filter(pheno.input, code=="RL"), aes(x=cell_lat, y=spring.dev, color=onset)) + 
  geom_point()) + theme_classic() + scale_color_viridis() + labs(x="Hex Latitude", y="Spring GDD deviation", title="Onset (RL)")
(plot.onset1<-ggplot(data=filter(pheno.input, code=="RE"), aes(x=cell_lat, y=spring.dev, color=onset)) + 
    geom_point()) + theme_classic() + scale_color_viridis() + labs(x="Hex Latitude", y="Spring GDD deviation", title="Onset (RE)")

(plot.onset1<-ggplot(data=filter(pheno.input), aes(x=cell_lat, y=spring.dev, color=onset)) + 
    geom_point()) + theme_classic() + scale_color_viridis() + labs(x="Hex Latitude", y="Spring GDD deviation", title="Onset") + facet_wrap(~code)



(plot.abund1<-ggplot(data=filter(pheno.input), aes(x=year, y=onset, color=code)) + 
    geom_smooth(method="lm")) + theme_classic() + scale_color_viridis_d() + labs(x="Year", y="Onset")




(plot.abund1<-ggplot(data=filter(pheno.input.dev, year>09), aes(x=year, y=onset, color=cell, size=obsDays)) + geom_point() +  
    geom_smooth(method="lm")) + theme_classic() + scale_color_viridis() + labs(x="Year", y="Onset deviation") + facet_wrap(~code)

(plot.abund1<-ggplot(data=filter(pheno.input.dev, code=="RE"), aes(x=year, y=onset, color=cell, size=obsDays)) + geom_point() +  
    geom_smooth(method="lm")) + theme_classic() + scale_color_viridis() + labs(x="Year", y="Onset deviation")

(plot.abund1<-ggplot(data=filter(pheno.input.dev, code=="RE"), aes(x=year, y=onset, color=cell, size=obsDays)) + geom_point() +  
    geom_smooth(method="lm")) + theme_classic() + scale_color_viridis() + labs(x="Year", y="Onset deviation")






#### ALT MODELS

ONSET

onset.full<-lmer(q5~-1+code+code*(warmdays+pc1+pc2+ year + gr_mn_lag + uniqObsDays) + (1|cell), data=pheno.input, weights=q5.ci.wt)

onset.nocode<-lmer(q5~-1+code+(warmdays+pc1+pc2+ year + gr_mn_lag + uniqObsDays) + (1|cell), data=pheno.input, weights=q5.ci.wt)

Terms<-c("warmdays","pc1","pc2","gr_mn_lag","year","uniqObsDays")

model.o<-lmer(q5~-1+code + (1|cell), data=pheno.input, weights=q5.ci.wt)
aic1<-extractAIC(model.o)[2]
ll1<-logLik(model.o, method="lme")
for (i in 1:length(Terms)) {
  thisformula<-formula(paste("q5 ~ -1 + code + ",Terms[i]," + (1|cell)"),sep=" ")
  onset.1<-lmer(thisformula, data=pheno.input, weights=q5.ci.wt)
  if(logLik(onset.1, method="lme")>ll1) { model.o<-onset.1; ll1<-logLik(onset.1, method="lme")}
}
for (i in 1:length(Terms)) {
  thisformula<-formula(paste("q5 ~ -1 + code + code:",Terms[i]," + (1|cell)"),sep="")
  onset.1<-lmer(thisformula, data=pheno.input, weights=q5.ci.wt)
  if(logLik(onset.1, method="lme")>ll1) { model.o<-onset.1; ll1<-logLik(onset.1, method="lme")}
}
for (i in c(1,3:length(Terms))) {
  thisformula<-formula(paste("q5 ~ -1 + code + code:pc1 + ",Terms[i]," + (1|cell)"),sep="")
  onset.1<-lmer(thisformula, data=pheno.input, weights=q5.ci.wt)
  if(logLik(onset.1, method="lme")>ll1) { model.o<-onset.1; ll1<-logLik(onset.1, method="lme")}
}
for (i in c(1,3:5)) {
  thisformula<-formula(paste("q5 ~ -1 + code + code:pc1 + uniqObsDays + code:",Terms[i]," + (1|cell)"),sep="")
  onset.1<-lmer(thisformula, data=pheno.input, weights=q5.ci.wt)
  if(logLik(onset.1, method="lme")>ll1) { model.o<-onset.1; ll1<-logLik(onset.1, method="lme")}
}

summary(model.o)

onset.a1<-lmer(q5~-1+code+code*(warmdays+pc1+pc2+ gr_mn_lag) + uniqObsDays + (1|cell), data=pheno.input, weights=q5.ci.wt)
onset.a2<-lmer(q5~-1+code+code*(warmdays+pc1+pc2+ year) + uniqObsDays + (1|cell), data=pheno.input, weights=q5.ci.wt)
onset.a3<-lmer(q5~-1+code+code*(warmdays+pc1+pc2+ year + gr_mn_lag) + (1|cell), data=pheno.input, weights=q5.ci.wt)
onset.a4<-lmer(q5~-1+code+code*(pc1+pc2+ year + gr_mn_lag) + uniqObsDays + (1|cell), data=pheno.input, weights=q5.ci.wt)
onset.a5<-lmer(q5~-1+code+code*(warmdays+pc2+ year + gr_mn_lag) + uniqObsDays + (1|cell), data=pheno.input, weights=q5.ci.wt)

onset.b1<-lmer(q5~-1+code+code*(pc1+pc2+ gr_mn_lag) + uniqObsDays + (1|cell), data=pheno.input, weights=q5.ci.wt)
onset.b2<-lmer(q5~-1+code+code*(pc1+pc2+ year) + uniqObsDays + (1|cell), data=pheno.input, weights=q5.ci.wt)
onset.b3<-lmer(q5~-1+code+code*(pc1+pc2+ year + gr_mn_lag) + (1|cell), data=pheno.input, weights=q5.ci.wt)
onset.b4<-lmer(q5~-1+code+code*(pc1+pc2+ year + gr_mn_lag) + uniqObsDays + (1|cell), data=pheno.input, weights=q5.ci.wt)
onset.b5<-lmer(q5~-1+code+code*(warmdays+pc2+ year + gr_mn_lag) + uniqObsDays + (1|cell), data=pheno.input, weights=q5.ci.wt)





################
pheno.input1<-na.omit(pheno.input)
onset.full<-lmer(q5~-1+code+code*(warmdays+pc1+pc2+ year + gr_mn_lag) + uniqObsDays + (1|cell), data=pheno.input1, weights=q5.ci.wt, na.action=na.fail)

x1<-dredge(onset.full, beta = c("none"), evaluate = TRUE,
       rank = "AICc", fixed = c("-1","code","uniqObsDays"), m.lim = NULL,
       trace = FALSE)

summary(x1)
x2<-x1[x1$weight>0.01,]
model.avg(x2)$coefficients[1,]

##
pheno.input.dev1<-na.omit(pheno.input.dev)
onsetdev.full<-lmer(onset.dev~-1+code+code*(warm.dev+pc1.dev+pc2.dev+ year + lag.dev) + uniqObsDays + (1|cell), data=pheno.input.dev1, weights=onset.ci.wt, na.action=na.fail)

xd1<-dredge(onsetdev.full, beta = c("none"), evaluate = TRUE,
           rank = "AICc", fixed = c("code","uniqObsDays"), m.lim = NULL,
           trace = FALSE)

summary(xd1)
xd2<-xd1[xd1$weight>0.01,]
model.avg(xd2)$coefficients[1,]



ggplot(data=pheno.input, aes(x=year, y=pc2)) + geom_point()
