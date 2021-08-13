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
  group_by(cell) %>%
  mutate( pc1.dev=warmearly-mean(warmearly, na.rm=T), pc2.dev=warmlateopen-mean(warmlateopen, na.rm=T))


pheno.input.dev<-merge(pheno.dev, env.dev, by=c("year","cell")) %>%
  filter(code %in% c("RE","RL","RP")) %>%
  mutate(year=year-2000, q5.ci.wt=1/(q5_ci+1),q50.ci.wt=1/(q50_ci+1),qdur.ci.wt=1/(qdur_ci+1))

ggplot(data=pheno.input.dev, aes(x=year, y=q5_dev, color=code)) + geom_jitter()
ggplot(data=pheno.input.dev, aes(x=year, y=q50_dev, color=code)) + geom_jitter()
ggplot(data=pheno.input.dev, aes(x=year, y=qdur_dev, color=code)) + geom_jitter()


pheno.input<-merge(pheno.quant, env.var, by=c("year","cell")) %>%
  filter(code %in% c("RE","RL","RP")) %>%
  mutate(year=year-2000, q5.ci.wt=1/(q5_ci+1),q50.ci.wt=1/(q50_ci+1),qdur.ci.wt=1/(qdur_ci+1))


### Analysis
library(MASS)

#### WEIGHTED LMER MODEL

onset.dev.full<-lmer(q5_dev~-1+code+FFD.dev + pc1.dev + pc2.dev + year+ uniqObsDays + code:FFD.dev + code:pc1.dev + code:pc2.dev + code:year + (1|cell), data=pheno.input.dev, weights=q5.ci.wt)
extractAIC(onset.dev.full)
r.squaredGLMM(onset.dev.full)     

t1<-get_model(step(onset.dev.full))
summary(t1)
onset.dev.final<-lmer(q5_dev~ -1 + code + year + code:pc1.dev + (1|cell), data=pheno.input.dev, weights=q5.ci.wt)
summary(onset.dev.final)
extractAIC(onset.dev.final)
r.squaredGLMM(onset.dev.final)     
test.vif<-lm(q5_dev~code+year+pc1.dev, data=pheno.input.dev)
vif(test.vif)

plot_model(onset.dev.final, type="est", title="Onset deviation: param. estimates")
plot_model(onset.dev.final, type = "eff", terms = c("pc1.dev", "code"), title="Onset deviation~Warm early deviation")

write.csv(summary(onset.dev.final)$coefficients, file="output/onsetdev.csv")

#ONSET
onset.full<-lmer(q5~-1+code+code*(FFD + warmearly + warmlateopen + year + uniqObsDays) + (1|cell), data=pheno.input, weights=q5.ci.wt)
(onset.step <- step(onset.full))
onset.best <- get_model(onset.step) #stepAIC(ab.yr.full)
summary(onset.best)
onset.final<-lmer(q5~-1+code+code:FFD + warmearly + warmlateopen + code:year + uniqObsDays + (1|cell), data=pheno.input, weights=q5.ci.wt)
r.squaredGLMM(onset.final)

onset.vif<-lmer(q5~code+FFD + warmearly + warmlateopen + year + uniqObsDays + (1|cell), data=pheno.input, weights=q5.ci.wt)
vif(onset.vif)

plot_model(onset.final, type = "eff", terms = c("year", "code"), title="Onset~Year")
plot_model(onset.final, type = "eff", terms = c("FFD", "code"), title="Onset~FFD")
plot_model(onset.final, type = "eff", terms = c("warmearly", "code"), title="Onset~PC1")

plot_model(onset.final, type = "eff", terms = c("warmlateopen", "code"), title="Onset~PC2")

write.csv(summary(onset.final)$coefficients, file="output/onset.csv")

##MEDIAN DEVIATION

med.dev.full<-lmer(q50_dev~-1+code*(FFD.dev + pc1.dev + pc2.dev + summer.dev + year)+ uniqObsDays+ (1|cell), data=pheno.input.dev, weights=q50.ci.wt)
r.squaredGLMM(med.dev.full)     
summary(med.dev.full)
med.dev.step <- step(med.dev.full)
med.dev.best <- get_model(med.dev.step) #stepAIC(ab.yr.full)
summary(med.dev.best)
r.squaredGLMM(med.dev.best)     
plot_model(med.dev.best)
med.dev.final<-lmer(q50_dev~-1+code + pc2.dev + summer.dev + year + uniqObsDays+ (1|cell), data=pheno.input.dev, weights=q50.ci.wt)
med.dev.vif<-lmer(q50_dev~code + pc2.dev + summer.dev + year + uniqObsDays+ (1|cell), data=pheno.input.dev, weights=q50.ci.wt)
vif(med.dev.vif)
plot_model(med.dev.vif, type = "eff", terms = c("pc2.dev", "code"), title="Onset~PC2")


#med
med.full<-lmer(q50~code*(FFD + warmearly + warmlateopen + summer.dev  + year) + (1|cell), data=pheno.input, weights=q50.ci.wt)
(med.step <- step(med.full))
med.best <- get_model(med.step) #stepAIC(ab.yr.full)
summary(med.best)
med.vif<-lmer(q50 ~ warmlateopen + summer.dev + (1|cell), data=pheno.input, weights=q50.ci.wt)
vif(med.vif)
r.squaredGLMM(med.best)     
plot_model(med.best, type = "eff", terms = c("warmlateopen"), title="Median")
plot_model(med.best, type = "eff", terms = c("summer.dev"), title="Median")


##DURATION DEVIATION
dur.dev.full<-lmer(qdur_dev~-1+code*(FFD.dev + pc1.dev + pc2.dev + summer.dev + year) + uniqObsDays + (1|cell), data=pheno.input.dev, weights=qdur.ci.wt)
r.squaredGLMM(dur.dev.full)     
summary(dur.dev.full)
dur.dev.step <- step(dur.dev.full)
dur.dev.best <- get_model(dur.dev.step) #stepAIC(ab.yr.full)
summary(dur.dev.best)
dur.dev.final<-lmer(qdur_dev~-1 + code + pc1.dev + code:pc2.dev + code:year + uniqObsDays + (1|cell), data=pheno.input.dev, weights=qdur.ci.wt)
r.squaredGLMM(dur.dev.final)     
summary(dur.dev.final)
extractAIC(dur.dev.final)     
dur.dev.vif<-lmer(qdur_dev~code + pc1.dev + pc2.dev + year + uniqObsDays + (1|cell), data=pheno.input.dev, weights=qdur.ci.wt)
vif(dur.dev.vif)
plot_model(dur.dev.final, type = "eff", terms = c("pc2.dev", "code"), title="Duration deviation~ pc2 deviation")
plot_model(dur.dev.final, type = "eff", terms = c("year", "code"), title="Duration deviation~ pc2 deviation")
plot_model(dur.dev.final, type = "eff", terms = c("year", "code"), title="Duration deviation~ pc2 deviation")

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
effort.lme<-lmer(dev_onset~code+obsDays+(1|cell), data=pheno.input.dev)
effort.lm<-lm(dev_onset~code+obsDays, data=pheno.input.dev)
summary(effort.lme)
summary(effort.lm)
