#Analysis of Presence-only Phenology Metrics
#eButterfly, NABA Butterflies I've Seen, iNaturalist (research grade)
## Resident butterfly species groups defined by overwinter stage (OWS)
## Phenometrics estimated from Weibull (M Belitz)
## Species OWS traits compiled by GU Ries Lab
#E Larsen, Georgetown U, Updated 2022-01

#libraries
library(lme4)
library(lmerTest)
library(MASS)
library(MuMIn)
library(car)
library(sp)
library(sjPlot)
library(viridis)
library(tidyverse)
library(ggeffects)
library(ggpubr)
library(gridExtra)
theme_set(theme_sjplot())

rundat<-Sys.Date()

###Input Files
pheno.data<-"data/derived/phenoDev.RData"
pheno.data.st<-"data/derived/pheno.RData"
env.dev.csv<-"data/derived/envDevs.csv"

##Parameters
mycolors<-viridis_pal()(8)[c(1,4,7)]

## Output files
onset.dev.output.csv<-"output/ondev.params.csv"
onset.dev.model.file<-"output/ondev.model.RData"

#Load phenology data
load(pheno.data)

(dev.hexyrs<-pheno.dev %>% dplyr::select(year, cell)  %>% group_by(year, cell) %>% tally())
nrow(pheno.dev)

#Load environmental data
env.dev<-read_csv(env.dev.csv)

pheno.input.dev<-merge(pheno.dev, env.dev, by=c("year","cell"), all.x=T) %>%
  filter(code %in% c("RE","RL","RP"), year<2018) %>%
  mutate(year=year-2000, onset.ci.wt=1/(onset.ci+1),med.ci.wt=1/(median.ci+1),dur.ci.wt=1/(dur.ci+1))

save(pheno.input.dev, file="data/derived/pheno.dev.modelinput.csv")

ggplot(data=pheno.input.dev, aes(x=year, y=onset.dev, color=code)) + geom_jitter() + scale_color_manual(values=mycolors)
ggplot(data=pheno.input.dev, aes(x=year, y=median.dev, color=code)) + geom_jitter() + scale_color_manual(values=mycolors)
ggplot(data=pheno.input.dev, aes(x=year, y=dur.dev, color=code)) + geom_jitter() + scale_color_manual(values=mycolors)

###########################################################################
# Phenology model: Annual deviance of onset -------------------------------

#weighted LMER with interactions of variables with overwintering code
#random effect of cell
onset.dev.full<-lmer(onset.dev~-1+code+warm.dev + dev.pc1 + dev.pc2 + open.lag.dev + year+ uniqObsDays + code:warm.dev + code:open.lag.dev + code:dev.pc1 + code:dev.pc2 + code:year  + code:uniqObsDays + (1|cell), data=pheno.input.dev, weights=onset.ci.wt)
extractAIC(onset.dev.full)
r.squaredGLMM(onset.dev.full)     
(t1<-as_tibble(bind_cols(Parameters=row.names(summary(onset.dev.full)$coefficients),summary(onset.dev.full)$coefficients)) %>% arrange(abs(`t value`)))

#stepwise model selection
onset.dev.best<-get_model(step(onset.dev.full))

(t2<-as_tibble(bind_cols(Parameters=row.names(summary(onset.dev.best)$coefficients),summary(onset.dev.best)$coefficients)) %>% arrange(abs(`t value`)))
summary(onset.dev.best)
summary(onset.dev.best)$call

#same model but structured for interpretability
#include either effect across all overwintering groups, or with group-specific model fitting
onset.dev.final<-lmer(onset.dev~ -1 + code + uniqObsDays + dev.pc1 + code:dev.pc2 + code:year +(1|cell), data=pheno.input.dev, weights=onset.ci.wt)
summary(onset.dev.final)
extractAIC(onset.dev.final)
r.squaredGLMM(onset.dev.final)     

#check for variable collinearity
test.vif<-lm(onset.dev~code+year+dev.pc1+dev.pc2+uniqObsDays, data=pheno.input.dev)
vif(test.vif)

#look at partial r2
onset.dev.pc1<-lmer(onset.dev~ -1 + code + uniqObsDays +code:dev.pc2  + code:year +(1|cell), data=pheno.input.dev, weights=onset.ci.wt)
print(paste("partial r2 for pc1:",round(r.squaredGLMM(onset.dev.final)[1]-r.squaredGLMM(onset.dev.pc1)[1],3)))
onset.dev.pc2<-lmer(onset.dev~ -1 + code + uniqObsDays + dev.pc1 + code:year +(1|cell), data=pheno.input.dev, weights=onset.ci.wt)
print(paste("partial r2 for pc2:",round(r.squaredGLMM(onset.dev.final)[1]-r.squaredGLMM(onset.dev.pc2)[1],3)))
onset.dev.yr<-lmer(onset.dev~ -1 + code + uniqObsDays + dev.pc1 + code:dev.pc2 +(1|cell), data=pheno.input.dev, weights=onset.ci.wt)
print(paste("partial r2 for year:",round(r.squaredGLMM(onset.dev.final)[1]-r.squaredGLMM(onset.dev.yr)[1],3)))
onset.dev.eff<-lmer(onset.dev~ -1 + code +  dev.pc1 + code:dev.pc2 +  code:year +(1|cell), data=pheno.input.dev, weights=onset.ci.wt)
print(paste("partial r2 for data density:",round(r.squaredGLMM(onset.dev.final)[1]-r.squaredGLMM(onset.dev.eff)[1],3)))

#save results and parameter table
(onset.params<-as_tibble(summary(onset.dev.final)$coefficients) %>%
  mutate(Estimate=round(Estimate, 4), StdError=round(`Std. Error`,4),parameter=row.names(summary(onset.dev.final)$coefficients)))
write_csv(onset.params, file=onset.dev.output.csv)
save(onset.dev.final, file=onset.dev.model.file)

# Figure 2: Phenology (annual deviation) model results ---------------------------------------

(plotyr<-plot_model(onset.dev.final, type = "eff", terms = c("year", "code"), alpha=0.5, size=2, colors=mycolors, se=TRUE, title="A. Onset Deviation ~ Year") + 
 theme_classic() + 
    scale_x_continuous(breaks = c(5, 10, 15),label = c("2005", "2010", "2015")) + 
  labs(x="Year", y="Adult onset deviation", color="Group", fill="Group") + 
    theme(legend.position="none"))

(plot.pc1<-plot_model(onset.dev.final, type = "eff", terms = c("dev.pc1", "code"), alpha=0.5, size=2, colors=mycolors, se=TRUE, title="B. Onset Deviation ~ PC1") + 
    theme_classic() + 
    labs(x="PC1 --> early greenup, warm spring", y="Adult onset deviation", color="Group", fill="Group") +
    theme(legend.position="none"))

(plot.pc2<-plot_model(onset.dev.final, type = "eff", terms = c("dev.pc2", "code"), alpha=0.5, size=2 , colors=mycolors, se=TRUE, title="C. Onset Deviation ~ PC2") + 
    theme_classic()  +
    labs(x="PC2 --> greenup later than expected", y="Adult onset deviation", color="Group", fill="Group")+ 
    theme(legend.position="none")) 

(legend.2<-ggplot(data=pheno.input.dev, aes(x=dev.pc1, y=onset.dev, color=code, fill=code)) + 
  geom_smooth(method="lm", alpha=0.5) + scale_color_manual(values=mycolors, aesthetics=c("color","fill"), labels=c("BOE","BOL","BOP")) +
   theme(legend.position="right",legend.key.size = unit(1, 'cm'),legend.title = element_text(size=14)) + labs(color="Group", fill="Group") )

l2<-get_legend(legend.2)

(figonset<-grid.arrange(plotyr, plot.pc1, plot.pc2,l2, nrow=1, widths=c(1,1,1,.5)))
ggsave(figonset, width=12,height=4, units="in", file=paste0("output/figures/FIG2pheno.png"))


##############



###########################################################################
# Phenology model: spatiotemporal onset -------------------------------
#Load environmental data
env<-read_csv("data/derived/spatemp_env.csv")
load(pheno.data.st)
pheno.quant<- pheno.quant %>%
  mutate(qdur=q95-q5, qdur_low=q95_low-q5_high, qdur_high=q95_high-q5_low,
       onset.ci=q5_high-q5_low,median.ci=q50_high-q50_low,
       q95_ci=q95_high-q95_low,dur.ci=qdur_high-qdur_low)


pheno.input<-merge(pheno.quant, env, by=c("year","cell"), all.x=T) %>%
  filter(code %in% c("RE","RL","RP"), year<2018) %>%
  rename(onset=q5, median=q50, duration=qdur) %>%
  mutate(Year=year, year=year-2000, open.lag=gr_mn_open-forest.greenup,
         onset.ci.wt=1/(onset.ci+1),med.ci.wt=1/(median.ci+1),dur.ci.wt=1/(dur.ci+1))

save(pheno.input, file="data/derived/pheno.ST.input.csv")

#weighted LMER with interactions of variables with overwintering code
#random effect of cell
onset.full<-lmer(onset~-1+code+warmdays + ST.PC1 + ST.PC2 + open.lag + year+ uniqObsDays + code:warmdays + code:open.lag + code:ST.PC1 + code:ST.PC2 + code:year  + code:uniqObsDays + (1|cell), data=pheno.input, weights=onset.ci.wt)
extractAIC(onset.full)
r.squaredGLMM(onset.full)     
(t1<-as_tibble(bind_cols(Parameters=row.names(summary(onset.full)$coefficients),summary(onset.full)$coefficients)) %>% arrange(abs(`t value`)))

#stepwise model selection
onset.best<-get_model(step(onset.full))

(t2<-as_tibble(bind_cols(Parameters=row.names(summary(onset.best)$coefficients),summary(onset.best)$coefficients)) %>% arrange(abs(`t value`)))
summary(onset.best)
summary(onset.best)$call

#same model but structured for interpretability
#include either effect across all overwintering groups, or with group-specific model fitting
onset.st.final<-lmer(onset~ -1 + code + code:ST.PC1 + code:warmdays + code:year + uniqObsDays + (1|cell), data=pheno.input, weights=onset.ci.wt)
summary(onset.st.final)
extractAIC(onset.st.final)
r.squaredGLMM(onset.st.final)     
write.csv(as.data.frame(summary(onset.st.final)$coefficients), file="output/onset.st.param.csv")

#check for variable collinearity
test.vif<-lm(onset~code + ST.PC1 + warmdays + year + uniqObsDays, data=pheno.input)
vif(test.vif)


# Supplemental: spatiotemporal onset model results ---------------------------------------

(plotyr<-plot_model(onset.st.final, type = "eff", terms = c("year", "code"), alpha=0.5, size=2, colors=mycolors, se=TRUE, title="A. Onset Deviation ~ Year") + 
   theme_classic() + 
   scale_x_continuous(breaks = c(5, 10, 15),label = c("2005", "2010", "2015")) + 
   labs(x="Year", y="Adult onset", color="Group", fill="Group") + 
   theme(legend.position="none"))

(plot.pc1<-plot_model(onset.st.final, type = "eff", terms = c("dev.pc1", "code"), alpha=0.5, size=2, colors=mycolors, se=TRUE, title="B. Onset Deviation ~ PC1") + 
    theme_classic() + 
    labs(x="PC1 --> early greenup, warm spring", y="Adult onset", color="Group", fill="Group") +
    theme(legend.position="none"))

(plot.pc2<-plot_model(onset.st.final, type = "eff", terms = c("dev.pc2", "code"), alpha=0.5, size=2 , colors=mycolors, se=TRUE, title="C. Onset Deviation ~ PC2") + 
    theme_classic()  +
    labs(x="PC2 --> greenup later than expected", y="Adult onset", color="Group", fill="Group")+ 
    theme(legend.position="none")) 

(legend.2<-ggplot(data=pheno.input.dev, aes(x=dev.pc1, y=onset.dev, color=code, fill=code)) + 
    geom_smooth(method="lm", alpha=0.5) + scale_color_manual(values=mycolors, aesthetics=c("color","fill"), labels=c("BOE","BOL","BOP")) +
    theme(legend.position="right",legend.key.size = unit(1, 'cm'),legend.title = element_text(size=14)) + labs(color="Group", fill="Group") )

l2<-get_legend(legend.2)

(figonset<-grid.arrange(plotyr, plot.pc1, plot.pc2,l2, nrow=1, widths=c(1,1,1,.5)))
ggsave(figonset, width=12,height=4, units="in", file=paste0("output/figures/supp.onset.fig.png"))


##############



###########################################################################
### PHENOLOGICAL MODELS: DURATION
###########################################################################
#ph2<-pheno.input
#pheno.input<-ph2 %>% mutate(summer.gdd=summer.gdd/500)
dur.full<-lmer(duration~-1+code+code*(warmdays + ST.PC1 + ST.PC2 + open.lag + summer.gdd + year+ uniqObsDays) + (1|cell), data=pheno.input, weights=dur.ci.wt)
extractAIC(dur.full)
r.squaredGLMM(dur.full)     
(t1<-as_tibble(bind_cols(Parameters=row.names(summary(dur.full)$coefficients),summary(dur.full)$coefficients)) %>% arrange(abs(`t value`)))

#stepwise model selection
dur.best<-get_model(step(dur.full))

(t2<-as_tibble(bind_cols(Parameters=row.names(summary(dur.best)$coefficients),summary(dur.best)$coefficients)) %>% arrange(abs(`t value`)))
summary(dur.best)
summary(dur.best)$call

#same model but structured for interpretability
#include either effect across all overwintering groups, or with group-specific model fitting
dur.st.final<-lmer(duration ~ -1 + code + ST.PC1 + code:warmdays + code:open.lag + code:summer.gdd + code:year + code:uniqObsDays + (1|cell), data=pheno.input, weights=dur.ci.wt)
summary(dur.st.final)
extractAIC(dur.st.final)
r.squaredGLMM(dur.st.final)     
write.csv(as.data.frame(summary(dur.st.final)$coefficients), file="output/dur.st.param.csv")


###########################################################################
### PHENOLOGICAL MODELS: DURATION DEVIATION
###########################################################################


dur.dev.full<-lmer(dur.dev~-1+code+ code*(warm.dev + summer.gdd.dev + dev.pc1 + dev.pc2 + open.lag.dev + year+ uniqObsDays) + (1|cell),
                   data=pheno.input.dev, weights=dur.ci.wt)
extractAIC(dur.dev.full)
r.squaredGLMM(dur.dev.full)     
(t2<-as_tibble(bind_cols(Parameters=row.names(summary(dur.dev.full)$coefficients),summary(dur.dev.full)$coefficients)) %>% arrange(abs(`t value`)))

dur.best<-get_model(step(dur.dev.full))

(t2<-as_tibble(bind_cols(Parameters=row.names(summary(dur.best)$coefficients),summary(dur.best)$coefficients)) %>% arrange(abs(`t value`)))
summary(dur.best)$call
dur.dev.final<-lmer(dur.dev~ -1 + code + warm.dev + code:summer.gdd.dev + dev.pc1 + code:dev.pc2 + code:year + code:uniqObsDays +  (1|cell), data=pheno.input.dev, weights=dur.ci.wt)
summary(dur.dev.final)
extractAIC(dur.dev.final)
r.squaredGLMM(dur.dev.final)     
test.vif<-lm(dur.dev~code+warm.dev+dev.pc1+dev.pc2+summer.gdd.dev+year+uniqObsDays, data=pheno.input.dev)
vif(test.vif)

plot_model(dur.dev.final, type = "eff", terms = c("warm.dev", "code"), title="durDev~Lag winter warm days")
plot_model(dur.dev.final, type = "eff", terms = c("summer.gdd.dev", "code"), title="durDev~Lag dev")
plot_model(dur.dev.final, type = "eff", terms = c("dev.pc2", "code"), title="durDev~PC1 dev")

durdev.output<-as_tibble(summary(dur.dev.final)$coefficients) %>%
  dplyr::mutate(param=row.names(summary(dur.dev.final)$coefficients), sig=ifelse(`Pr(>|t|)`<0.05,1,0), r2m=round(r.squaredGLMM(dur.dev.final)[1],2),r2c=round(r.squaredGLMM(dur.dev.final)[2],2))
write.csv(durdev.output, file=paste0("output/dur.dev.model",rundat,".csv"))


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
onsetdev.full<-lmer(onset.dev~-1+code+code*(warm.dev+dev.pc1+dev.pc2+ year + lag.dev) + uniqObsDays + (1|cell), data=pheno.input.dev1, weights=onset.ci.wt, na.action=na.fail)

xd1<-dredge(onsetdev.full, beta = c("none"), evaluate = TRUE,
           rank = "AICc", fixed = c("code","uniqObsDays"), m.lim = NULL,
           trace = FALSE)

summary(xd1)
xd2<-xd1[xd1$weight>0.01,]
model.avg(xd2)$coefficients[1,]



ggplot(data=pheno.input, aes(x=year, y=pc2)) + geom_point()
