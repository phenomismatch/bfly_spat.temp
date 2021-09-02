#Analysis of Presence-only Phenology Metrics
#eButterfly, NABA Butterflies I've Seen, iNaturalist (research grade)
## Resident butterfly species groups defined by overwinter stage (OWS)
## Phenometrics estimated from Weibull (M Belitz)
## Species OWS traits compiled by GU Ries Lab
#E Larsen, Georgetown U, Updated 2021-08


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
env.var<-read_csv("data/derived/envir.csv")


env.dev<-env.var %>%
  mutate(include=ifelse(year>2015,1,NA),baselinepc1=warmearly*include, baselinepc2=warmlateopen*include) %>%
  group_by(cell) %>%
  mutate( pc1.dev=warmearly-mean(baselinepc1, na.rm=T),pc2.dev=warmlateopen-mean(baselinepc2, na.rm=T))

#env.dev<-env.var %>%
#  group_by(cell) %>%
#  mutate( pc1.dev=warmearly-mean(warmearly, na.rm=T), pc2.dev=warmlateopen-mean(warmlateopen, na.rm=T))


pheno.input.dev<-merge(pheno.dev, env.dev, by=c("year","cell")) %>%
  filter(code %in% c("RE","RL","RP")) %>%
  mutate(year=year-2000, onset.ci.wt=1/(onset.ci+1),med.ci.wt=1/(median.ci+1),dur.ci.wt=1/(dur.ci+1))

ggplot(data=pheno.input.dev, aes(x=year, y=onset.dev, color=code)) + geom_jitter()
ggplot(data=pheno.input.dev, aes(x=year, y=median.dev, color=code)) + geom_jitter()
ggplot(data=pheno.input.dev, aes(x=year, y=dur.dev, color=code)) + geom_jitter()


pheno.input<-merge(pheno.quant, env.var, by=c("year","cell")) %>%
  filter(code %in% c("RE","RL","RP")) %>%
  mutate(year=year-2000, q5.ci.wt=1/(q5_ci+1),q50.ci.wt=1/(q50_ci+1),qdur.ci.wt=1/(qdur_ci+1))


### Analysis
library(MASS)

#### WEIGHTED LMER MODEL

onset.dev.full<-lmer(onset.dev~-1+code+FFD.dev + pc1.dev + pc2.dev + year+ uniqObsDays + code:FFD.dev + code:pc1.dev + code:pc2.dev + code:year + (1|cell), data=pheno.input.dev, weights=onset.ci.wt)
extractAIC(onset.dev.full)
r.squaredGLMM(onset.dev.full)     
t2<-as_tibble(bind_cols(Parameters=row.names(summary(onset.dev.full)$coefficients),summary(onset.dev.full)$coefficients)) %>% arrange(abs(`t value`))

t1<-get_model(step(onset.dev.full))

t2<-as_tibble(bind_cols(Parameters=row.names(summary(t1)$coefficients),summary(t1)$coefficients)) %>% arrange(abs(`t value`))

onset.dev.final<-lmer(onset.dev~ -1 + code + uniqObsDays + pc1.dev + code:FFD.dev  + code:year +(1|cell), data=pheno.input.dev, weights=onset.ci.wt)
summary(onset.dev.final)
extractAIC(onset.dev.final)
r.squaredGLMM(onset.dev.final)     
test.vif<-lm(onset.dev~code+year+pc1.dev+FFD.dev+uniqObsDays, data=pheno.input.dev)
vif(test.vif)

plot_model(onset.dev.final, type="est", title="Onset deviation: param. estimates")
plot_model(onset.dev.final, type = "eff", terms = c("pc1.dev", "code"), title="Onset deviation~Warm early deviation")

write.csv(summary(onset.dev.final)$coefficients, file="output/onsetdev.csv")

#predict
newcode<-c("RE","RL","RP")
codelabel<-c("egg","caterpillar","pupa")
newpc1.dev<-round(c( ((min(pheno.input.dev$pc1.dev, na.rm=T))*2):((max(pheno.input.dev$pc1.dev, na.rm=T))*2)/2),2)
newffd.dev<-round(c( ((round(min(pheno.input.dev$FFD.dev, na.rm=T)))):((round(max(pheno.input.dev$FFD.dev, na.rm=T))))),2)
newyr.dev<-c(0:20)
nrep<-length(newcode)*length(newyr.dev)*length(newffd.dev)
newDat <- data.frame(cell = rep(622,nrep), uniqObsDays=rep(median(pheno.input.dev$uniqObsDays, na.rm=T),nrep),
                     FFD.dev=rep(newffd.dev,length(newyr.dev)*length(newcode)),
                     pc1.dev<-rep(median(pheno.input.dev$pc1.dev, na.rm=T),nrep),
                     year=rep(newyr.dev,each=length(newffd.dev),length(newcode)),
                     code=rep(newcode, each=length(newyr.dev)*length(newffd.dev)),
                     codelabel=rep(codelabel, each=length(newyr.dev)*length(newffd.dev)))

newDat$pred <- predict(onset.dev.final, newDat,allow.new.levels =T)

lims.1<-pheno.input.dev %>%
  mutate(codelabel=ifelse(code=="RE","egg",ifelse(code=="RL","caterpillar","pupa"))) %>% 
  group_by(code, codelabel) %>%
  summarize(miny=floor(min(year, na.rm=T)*10)/10,maxy=ceiling(max(year, na.rm=T)*10)/10,
            minx=floor(min(FFD.dev, na.rm=T)*10)/10,maxx=ceiling(max(FFD.dev, na.rm=T)*10)/10)

newDat1<-inner_join(newDat,lims.1) %>%
  mutate(f1=ifelse(year>=miny,ifelse(year<=maxy,1,0),0)+ifelse(FFD.dev>=minx,ifelse(FFD.dev<=maxx,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel = fct_reorder(codelabel, code))

ggplot(data=newDat1, aes(y=year, x=FFD.dev, fill=pred)) + 
  geom_tile() +  labs(y="Year", x="Winter Warm Days above baseline") + 
  scale_fill_viridis(name="Adult onset deviation") + 
  facet_wrap(~codelabel)


#ONSET
onset.full<-lmer(q5~-1+code+code*(FFD + warmearly + warmlateopen + year + uniqObsDays) + (1|cell), data=pheno.input, weights=q5.ci.wt)
(onset.step <- step(onset.full))
onset.best <- get_model(onset.step) #stepAIC(ab.yr.full)
summary(onset.best)
onset.final<-lmer(q5~-1+code+code:FFD +warmearly + code:warmlateopen + code:year + uniqObsDays + (1|cell), data=pheno.input, weights=q5.ci.wt)
r.squaredGLMM(onset.final)

onset.vif<-lmer(q5~code+FFD + warmearly + warmlateopen + year + uniqObsDays + (1|cell), data=pheno.input, weights=q5.ci.wt)
vif(onset.vif)

plot_model(onset.final, type = "eff", terms = c("year", "code"), title="Onset~Year")
plot_model(onset.final, type = "eff", terms = c("FFD", "code"), title="Onset~FFD")
plot_model(onset.final, type = "eff", terms = c("warmearly", "code"), title="Onset~PC1")

plot_model(onset.final, type = "eff", terms = c("warmlateopen", "code"), title="Onset~PC2")

write.csv(summary(onset.final)$coefficients, file="output/onset.csv")


#predict
newcode<-c("RE","RL","RP")
codelabel<-c("egg","caterpillar","pupa")
newffd<-round(c( ((round(min(pheno.input$FFD, na.rm=T)))):((round(max(pheno.input$FFD, na.rm=T))))),2)
newyr<-c(0:20)
nrep<-length(newcode)*length(newyr)*length(newffd)
newDatST <- data.frame(cell = rep(622,nrep), uniqObsDays=rep(median(pheno.input$uniqObsDays, na.rm=T),nrep),
                     FFD=rep(newffd,length(newyr)*length(newcode)),
                     warmearly=rep(median(pheno.input$warmearly, na.rm=T),nrep),
                     warmlateopen=rep(median(pheno.input$warmlateopen, na.rm=T),nrep),
                     year=rep(newyr,each=length(newffd),length(newcode)),
                     code=rep(newcode, each=length(newyr)*length(newffd)),
                     codelabel=rep(codelabel, each=length(newyr)*length(newffd)))

newDatST$pred <- predict(onset.final, newDatST,allow.new.levels =T)

lims.2<-pheno.input %>%
  mutate(codelabel=ifelse(code=="RE","egg",ifelse(code=="RL","caterpillar","pupa"))) %>% 
  group_by(code, codelabel) %>%
  summarize(miny=floor(min(year, na.rm=T)*10)/10,maxy=ceiling(max(year, na.rm=T)*10)/10,
            minx=floor(min(FFD.dev, na.rm=T)*10)/10,maxx=ceiling(max(FFD.dev, na.rm=T)*10)/10)

newDat2<-inner_join(newDatST,lims.2) %>%
  mutate(f1=ifelse(year>=miny,ifelse(year<=maxy,1,0),0)+ifelse(FFD>=minx,ifelse(FFD<=maxx,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel = factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel, as.numeric(code)))


ggplot(data=newDat2, aes(y=year, x=FFD, fill=pred)) + 
  geom_tile() +  labs(y="Year", x="Winter Warm Days") + 
  scale_fill_viridis(name="Onset") + 
  facet_wrap(~codelabel)


newcode<-c("RE","RL","RP")
codelabel<-c("egg","caterpillar","pupa")
newffd<-100
newyr<-c(10)
warmearly<-round(c( ((floor(min(pheno.input$warmearly, na.rm=T)))):((ceiling(max(pheno.input$warmearly, na.rm=T))))),2)
warmlateopen<-round(c( ((floor(min(pheno.input$warmlateopen, na.rm=T)))):((ceiling(max(pheno.input$warmlateopen, na.rm=T))))),2)
nrep<-length(newcode)*length(warmearly)*length(warmlateopen)
newDat <- data.frame(cell = rep(622,nrep), uniqObsDays=rep(median(pheno.input$uniqObsDays, na.rm=T),nrep),
                     FFD=rep(newffd,nrep),
                     warmearly=rep(warmearly,length(newcode)*length(warmlateopen)),
                     warmlateopen=rep(warmlateopen,each=length(warmearly),length(newcode)),
                     year=rep(newyr,nrep),
                     code=rep(newcode, each=length(warmearly)*length(warmlateopen)),
                     codelabel=rep(codelabel, each=length(warmearly)*length(warmlateopen)))

newDat$pred <- predict(onset.final, newDat,allow.new.levels =T)


ggplot(data=newDat, aes(y=warmearly, x=warmlateopen, fill=pred)) + 
  geom_tile() +  labs(y="PC1", x="PC2") + 
  scale_fill_viridis(name="Onset") + 
  facet_wrap(~codelabel)





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






#### JUST sampling effort
effort.lme<-lmer(onset.dev~code+uniqObsDays+(1|cell), data=pheno.input.dev)
effort.lm<-lm(onset.dev~code+uniqObsDays, data=pheno.input.dev)
summary(effort.lme)
summary(effort.lm)
