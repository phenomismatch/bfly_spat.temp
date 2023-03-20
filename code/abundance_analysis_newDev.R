#Modeling patterns in butterfly abundance
#Elise Larsen, Georgetown U, Updated 2012-08
#2022 Sept - Oct : adding precipitation

#libraries
library(MASS)
library(lme4)
library(lmerTest)
library(r2glmm)
library(sjPlot)
library(MuMIn)
library(car) #(for VIF function)
library(tidyverse)
library(viridis)
library(ggpubr)
library(gridExtra)
library(png)
library(magick)
library(cowplot)


#Parameters
ows.colors<-viridis_pal()(8)[c(1,4,7)]
survey.months<-c(6:8)
theme_set(theme_sjplot())

#Input files
abundance.file<-"data/derived/naba_OWS_abundances.csv"
study.cells.file<-"data/spatial.domain.RData"
pheno.data<-"data/derived/pheno.input23.RData"
env.dev.csv<-"data/derived/envDevs.csv"
load("data/derived/precipSum.RData")

#Output files
abund.model.input<-"data/abund.input23.RData"
abundance.model.fit<-"output/abundance.finalmodel23.csv"
#Load, filter, & format abundance data
load(study.cells.file)

abund<-read_csv(abundance.file) %>%
  dplyr::select(cell, year=ObsYear, doy, CountID, SurveyID, ObsMonth, code=group, abund.bph, log.abund, SR) %>%
  filter(ObsMonth %in% survey.months, cell %in% STUDYCELLS) %>%
  mutate(code=ifelse(code=="RE","EO",ifelse(code=="RL","LO","PO"))) %>%
  group_by(cell, year, CountID, code) %>%
  mutate(x=1, seq=cumsum(x), nsurveys=max(seq))

(a1<-abund %>% group_by(code) %>% summarize(meanbph=mean(abund.bph, na.rm=T), sdbph=sd(abund.bph, na.rm=T)))

##Data visualization - base for FIGURE 1 panel (created in figures.R file)
a2<-abund %>% group_by(year, code) %>% summarize(ct=log(mean(abund.bph, na.rm=T)))
ggplot(data=a2, aes(x=year, y=ct, color=code)) + geom_line() + geom_smooth(method="lm") + scale_color_manual(values=ows.colors)


#previous year abundance
abund.py<-abund %>%
  mutate(year=year+1) %>%
  filter(year < 2018) %>%
  dplyr::select(cell, year, CountID, code, doy.py=doy, abund.py=abund.bph, logab.py=log.abund, seq.py=seq, nsurv.py=nsurveys) 

#Add previous year abundance to abundance table
abund1<-merge(abund, abund.py, by=c("year", "cell","CountID", "code"), all.x=T) %>%
  mutate(keep=ifelse(seq==seq.py | nsurv.py==1,1,0)) %>% filter(keep==1) %>%
  dplyr::select(year:SR, abund.py, logab.py)

summary(abund1)

#Add environmental data
env.dev<-read_csv(env.dev.csv)

##LOAD PRECIP SUM
precip.s1<-precip.sum %>% rename(year=yr.1) %>%
  pivot_wider(id_cols=c(year, cell), names_from=season,values_from=pcp.dev.std)
env.dev<-merge(env.dev, precip.s1, all.x=T, by=intersect(names(env.dev), names(precip.s1)))
summary(env.dev)
abund.env<-na.omit(merge(abund1, env.dev, by=c("year", "cell"), all.x=T))

##Filter and Add phenology data
###PHENO DATA
load(pheno.data)
pheno.dev<-pheno.dev %>% 
  mutate(code=ifelse(code=="RE","EO",ifelse(code=="RL","LO","PO"))) 

abund.pheno.1<-merge(abund.env, pheno.dev, by=c("year","cell","code"), all.x=T) %>%
  mutate(Year=ifelse(year>1900,year-2000,year), year=ifelse(year<1900,year+2000,year), abs.surveylag=abs(q50-doy))

summary(abund.pheno.1)

naba.1<-(abund.pheno.1) %>% 
  mutate(code=as.factor(code),daylag=doy-q50) %>%
  rename(summer.dev=summer.gdd.dev, on.dev=onset.dev, dur.dev=dur.dev, lag.dev=open.lag.dev, cold.dev=cold.dev,abslag=abs.surveylag) %>%
  mutate(MonthF=as.factor(ObsMonth), ObsDay=as.numeric(format(as.Date(doy, origin = paste((year+1999),"12-31",sep="-")),format="%d")))


save(naba.1, file="data/derived/abund.model.input23.RData")


#Visualization: basis for Figure 1D
vis.1<-naba.1 %>% group_by(code, year) %>% summarize(abundance=log(mean(abund.bph, na.omit=T)))
ggplot(data=vis.1, aes(x=year, y=abundance, color=code)) + geom_line() + geom_smooth(method="lm") + scale_color_manual(values=ows.colors)
ab.yr<-lmer(log.abund~-1+code+code:year+(1|cell) + (1|CountID:cell), data=naba.1)
summary(ab.yr)


# Abundance model ---------------------------------------------------------
#weighted LMER with interactions of variables with overwintering code
#random effect of cell
#Additional model selection further below
abundance.fullmodel<-lmer(log.abund~-1+code*(dev.pc1+dev.pc2+logab.py+lag.dev+on.dev+abslag+dur.dev+MonthF*ObsDay+Year+cold.dev+warm.dev+summer.dev+spring.p+summer.p)+(1|cell) + (1|CountID:cell), data=na.omit(naba.1))
summary(abundance.fullmodel)
abundance.best.1<-get_model(step(abundance.fullmodel))
summary(abundance.best.1)$call

#same model but with parameters for either all overwintering groups together, or interactive effect, for interpretability
naba.2<-naba.1 %>% 
  mutate(od.wk=on.dev/7, sumGDD=summer.dev/100, 
         sampl.lag=abslag/7, dd.wk=dur.dev/7)

best1b<-lmer(log.abund ~ -1 + code + logab.py + od.wk + sumGDD +
             (1 | cell) + (1 | CountID:cell) + 
  code:dev.pc1 + code:dev.pc2 + code:lag.dev + 
  code:sampl.lag + code:dd.wk + code:MonthF:ObsDay + 
  code:Year + code:cold.dev + code:spring.p, data = na.omit(naba.2))

abund.best<-get_model(step(best1b))
extractAIC(best1b)
extractAIC(abund.best)


r.squaredGLMM(abund.best)
bestmodel<-abund.best
summary(abund.best)

abund.vif<-lmer(log.abund ~ code + logab.py + od.wk + sumGDD +
                  dev.pc1 + dev.pc2 + lag.dev + 
                  sampl.lag + dd.wk + MonthF:ObsDay + 
                  Year + cold.dev + spring.p + (1 | cell) + (1 | CountID:cell), data = na.omit(naba.2))

vif(abund.vif)

write.csv(summary(bestmodel)$coefficients, file=abundance.model.fit)

# FIGURE 3: Model results for abundance -----------------------------------
(plot.pc1<-plot_model(bestmodel, type = "eff", terms = c("dev.pc1", "code"), colors=ows.colors, alpha=0.5, se=TRUE, title="Spring GDD + forest greenup") + 
      theme_classic() + ylim(0,4.3) + 
      labs(x="PC1  ---> warm spring \n + early greenup (days)", y="", color="Group", fill="Group") +
      theme(legend.position="none"))

(plot.pc2<-plot_model(bestmodel, type = "eff", terms = c("dev.pc2", "code"), colors=ows.colors, alpha=0.5,se=TRUE, title="GDD-Greenup decoupling") + 
      theme_classic() +  ylim(0,4.3) + 
      labs(x="PC2 --> greenup later \n than expected (days)",y="", color="Group", fill="Group") +
      theme(legend.position="none"))
  
(plot.lag<-plot_model(bestmodel, type = "eff", terms = c("lag.dev", "code"), colors=ows.colors, alpha=0.5, se=TRUE, title="Open canopy greenup lag") + 
      theme_classic() + ylim(0,4.3) + 
      labs(x="Open canopy greenup \n lag (days)", y="", color="Group", fill="Group") +
      theme(legend.position="none"))
  
(plot.dur<-plot_model(bestmodel, type = "eff", terms = c("dd.wk", "code"), colors=ows.colors, alpha=0.5, se=TRUE, title="Adult flight duration") + 
      theme_classic() + ylim(0,4.3) +  
      labs(x="\u25B3 Flight duration (weeks)", y="Log adult abundance",  color="Group", fill="Group") +
      theme(legend.position="none"))
  
(plot.yr<-plot_model(bestmodel, type = "eff", terms = c("Year", "code"), colors=ows.colors, alpha=0.5, se=TRUE, title="Year") + 
      theme_classic() + ylim(0,4.3) + 
      labs(x="Year", y="Log adult abundance", color="Group", fill="Group") +
      theme(legend.position="none")) + scale_x_continuous(breaks=c(5,10,15), labels=c("2005","2010","2015"))
(plot.rain<-plot_model(bestmodel, type = "eff", terms = c("spring.p", "code"), colors=ows.colors, alpha=0.5, se=TRUE, title="Spring precipitation") + 
    theme_classic() +  ylim(0,4.3) + 
    labs(x="\u25B3 spring precipitation", y="", color="Group", fill="Group") +
    theme(legend.position="none"))
(plot.cold<-plot_model(bestmodel, type = "eff", terms = c("cold.dev", "code"), colors=ows.colors, alpha=0.5, se=TRUE, title="Winter cold") + 
    theme_classic() +  ylim(0,4.3) +  
    labs(x="\u25B3 # winter days \n below freezing", y="", color="Group", fill="Group") +
    theme(legend.position="none"))


(legend.2<-ggplot(data=naba.1, aes(x=dev.pc1, y=on.dev, color=code, fill=code)) + 
      geom_smooth(method="lm", alpha=0.5) + scale_color_manual(values=ows.colors, aesthetics=c("color","fill"), labels=c("EO","LO","PO")) +
      theme(legend.position="right",legend.key.size = unit(1, 'cm'),legend.title = element_text(size=14)) + labs(color="Group", fill="Group") )
l2<-get_legend(legend.2)

#Create a panel legend with shapes for each group
x11<-data.frame(group=c("EO","LO","PO"),x1=c(1,1,1),y1=c(1,2,3),pt=c(1.5,2.5,3.5))

#overwinter group images
BOEraster = readPNG("data/images/boe_art.png")
BOLraster = readPNG("data/images/bol_art.png")
BOPraster = readPNG("data/images/bop_art.png")
x1<-c(1,1,1)
y1<-c(1:3)

rev.col<-rev(ows.colors)
(legend3<- ggplot(data=x11, aes(x=x1, y=y1)) +  
    geom_rect(xmin=x1, xmax=x1+1,ymin=y1,ymax=y1+1, fill=rev.col, alpha=0.5) +
    geom_segment(x=x1,xend=x1+1,y=y1+0.5,yend=y1+0.5, color=rev.col, alpha=1) +
    ylim(-.5,6.5) + xlim(0,5) +  theme_map(12) + labs(title="Adult Abundance") +
    #annotate("text", x = 2, y = 5.5, label = c("Adult Abundance"), size=6) + 
    annotate("text", x = 1, y = 5, label = c("Overwinter stage"), size=5, hjust="left") + 
    annotate("text", x = rep(2.6,3), y = c(1:3)+.5, label = c("PO", "LO", "EO"), size=4) + 
    theme(plot.background = element_rect(color="white", fill = "white")) +
    annotation_raster(BOPraster, xmin = 3.1, xmax=4.5,
                      ymin = 1.2, ymax = 1.8, interpolate=T) +
    annotation_raster(BOLraster, xmin = 3.1, xmax = 4.3, 
                      ymin = 2.2, ymax = 2.8, interpolate = T) + 
  annotation_raster(BOEraster, xmin = 3.1, xmax = 4.5, 
                    ymin = 3.2, ymax = 3.8, interpolate = T) 
)
  
margin = theme(plot.margin = unit(c(0.1,0,0.5,0), "cm")) #top, right, bottom, left
p.a1<-legend3
p.a2<-plot.yr+margin
p.a3<-plot.dur+margin + labs(y="") + xlim(-10,10)
p.a4<-plot.rain+margin
p.b1<-plot.cold + labs(y="Log adult abundance") +xlim(-30,40)   +margin
p.b2<-plot.pc1 +xlim(-30,40) +margin
p.b3<-plot.pc2 +xlim(-30,40) +margin
p.b4<-plot.lag +xlim(-30,40) +margin
(figabund<-grid.arrange(p.a1,p.a2,p.a3,p.a4,p.b1,p.b2,p.b3,p.b4, nrow=2))
ggsave(figabund, width=10,height=4.5, units="in", file=paste0("output/figures/abund.dev",rundat,".png"))


legend<-legend3 #+labs(title="Adult Abundance")
(figabund<-grid.arrange(fig3a, fig3b, legend3, plot.pc2, plot.dur,plot.yr,nrow=2))
(figabundrain<-grid.arrange(fig3a, fig3b, plot.pc2,legend3, plot.dur,plot.yr,plot.rainq1, plot.rainq3,nrow=2))
ggsave(figabund, width=10,height=6, units="in", file=paste0("output/figures/abund.dev3",rundat,".png"))
ggsave(figabundrain, width=12,height=4, units="in", file=paste0("output/figures/abund.rain",rundat,".png"))

(figabund<-grid.arrange(legend,plot.yr, plot.cold, plot.rain,plot.dur,plot.pc1,  plot.pc2, plot.lag,  nrow=2))
ggsave(figabund, width=12,height=6, units="in", file=paste0("output/figures/abund.dev",rundat,".png"))

                        
                      


# Model Selection ---------------------------------------------------------


  ##REGSUBSETS MODELS
  b<-regsubsets(log.abund~-1+code+code:dev.pc1+code:dev.pc2+code:open.lag.dev+code:on.dev+code:abs.surveylag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+code:logab.py+code:Year+(1|cell), data=naba.1,nbest=1)
  b<-regsubsets(log.abund~-1+dev.pc1+dev.pc2+open.lag.dev+on.dev+abs.surveylag+dur.dev+MonthF + MonthF:ObsDay+cold.dev+logab.py+Year, data=filter(naba.1,code=="RE"),nbest=1)
  
  
  
  
  ab.yr.1d<-lmer(log.abund~-1+code+code:dev.pc1+code:dev.pc2+code:on.dev+abs.surveylag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+logab.py+code:Year+(1|cell) + (1|CountID:cell), data=naba.1)
  summary(ab.yr.1d)
  extractAIC(ab.yr.1d)
  
  ab.yr.1e<-lmer(log.abund~-1+code+code:dev.pc1+code:dev.pc2+code:open.lag.dev+abs.surveylag+code:dur.dev+code:MonthF+code:cold.dev+logab.py+code:Year+(1|cell) + (1|CountID:cell), data=naba.1)
  summary(ab.yr.1e)
  extractAIC(ab.yr.1e)
  
  
  
  
    ab.yr.1b<-lmer(log.abund~-1+code+code:logab.py+(1|cell) + (1|CountID:cell), data=naba.1)
  extractAIC(ab.yr.1b)
  
  
    best.model<-ab.yr.1
  ab.yr.noint1<-lmer(log.abund~-1+code+dev.pc1+code:logab.py+code:dev.pc2+code:on.dev+code:abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+code:year+(1|cell) + (1|CountID:cell), data=naba.1)
  extractAIC(ab.yr.noint1)
  if(extractAIC(ab.yr.noint1)[2]<extractAIC(best.model)[2]) {best.model<-ab.yr.noint1}
  
  ab.yr.nointlag<-lmer(log.abund~-1+code+code:dev.pc1+code:logab.py+code:dev.pc2+code:on.dev+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+code:year+(1|cell) + (1|CountID:cell), data=naba.1)
  extractAIC(ab.yr.nointlag)
  if(extractAIC(ab.yr.nointlag)[2]<extractAIC(best.model)[2]) {best.model<-ab.yr.nointlag}
  extractAIC(best.model)
  
  ab.yr.nointyr<-lmer(log.abund~-1+code+code:dev.pc1+code:logab.py+code:dev.pc2+code:on.dev+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+year+(1|cell) + (1|CountID:cell), data=naba.1)
  extractAIC(ab.yr.nointyr)
  if(extractAIC(ab.yr.nointyr)[2]<extractAIC(best.model)[2]) {best.model<-ab.yr.nointyr}
  extractAIC(best.model)
 
  ab.yr.noyr<-lmer(log.abund~-1+code+code:dev.pc1+code:logab.py+code:dev.pc2+code:on.dev+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+(1|cell) + (1|CountID:cell), data=naba.1)
  extractAIC(ab.yr.noyr)
  if(extractAIC(ab.yr.noyr)[2]<extractAIC(best.model)[2]) {best.model<-ab.yr.noyr}
  extractAIC(best.model)

  ab.yr.nopyint<-lmer(log.abund~-1+code+code:dev.pc1+logab.py+code:dev.pc2+code:on.dev+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+(1|cell) + (1|CountID:cell), data=naba.1)
  extractAIC(ab.yr.nopyint)
  if(extractAIC(ab.yr.nopyint)[2]<extractAIC(best.model)[2]) {best.model<-ab.yr.nopyint}
  extractAIC(best.model)
  
  
  
  
  
  (yrpr2<-r.squaredGLMM(best.model)[1]-r.squaredGLMM(ab.yr.noyr)[1])
  ab.yr.nopc1<-lmer(log.abund~-1+code+code:logab.py+code:dev.pc2+code:on.dev+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+year+(1|cell) + (1|CountID:cell), data=naba.1)
  (pc1pr2<-r.squaredGLMM(best.model)[1]-r.squaredGLMM(ab.yr.nopc1)[1])
  ab.yr.nopc2<-lmer(log.abund~-1+code+code:logab.py+code:dev.pc1+code:on.dev+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+year+(1|cell) + (1|CountID:cell), data=naba.1)
  (pc2pr2<-r.squaredGLMM(best.model)[1]-r.squaredGLMM(ab.yr.nopc2)[1])
  ab.yr.noon<-lmer(log.abund~-1+code+code:logab.py+code:dev.pc1+code:dev.pc2+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+year+(1|cell) + (1|CountID:cell), data=naba.1)
  (onpr2<-r.squaredGLMM(best.model)[1]-r.squaredGLMM(ab.yr.noon)[1])
  ab.yr.py<-lmer(log.abund~-1+code+code:on.dev+code:dev.pc1+code:dev.pc2+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+year+(1|cell) + (1|CountID:cell), data=naba.1)
  (pypr2<-r.squaredGLMM(best.model)[1]-r.squaredGLMM(ab.yr.py)[1])
  ab.yr.dd<-lmer(log.abund~-1+code+code:on.dev+code:dev.pc1+code:dev.pc2+abslag+code:logab.py+code:MonthF + code:MonthF:ObsDay+code:cold.dev+year+(1|cell) + (1|CountID:cell), data=naba.1)
  (ddpr2<-r.squaredGLMM(best.model)[1]-r.squaredGLMM(ab.yr.dd)[1])
  ab.yr.cd<-lmer(log.abund~-1+code+code:on.dev+code:dev.pc1+code:dev.pc2+abslag+code:logab.py+code:MonthF + code:MonthF:ObsDay+code:dur.dev+year+(1|cell) + (1|CountID:cell), data=naba.1)
  (cdpr2<-r.squaredGLMM(best.model)[1]-r.squaredGLMM(ab.yr.cd)[1])
  
    ab.vif<-lmer(log.abund~logab.py+code+dev.pc1+dev.pc2+on.dev+abslag+dur.dev+ObsDay:MonthF+cold.dev+year+(1|cell) + (1|CountID:cell), data=naba.1)
    vif(ab.vif)
    
    lm1<-lm(log.abund ~ -1 + code + code:dev.pc1 + code:dev.pc2 + code:open.lag.dev +  
              abs.surveylag + code:dur.dev + code:MonthF + code:MonthF:ObsDay +  
              code:cold.dev + logab.py + code:Year, data=naba.1)
    
best.model<-best1b

    
abund.output<-as_tibble(summary(best.model)$coefficients) %>%
  dplyr::mutate(param=row.names(summary(best.model)$coefficients), sig=ifelse(`Pr(>|t|)`<0.05,1,0), r2m=round(r.squaredGLMM(best.model)[1],2),r2c=round(r.squaredGLMM(best.model)[2],2))
write.csv(abund.output, file="output/abund.model1027.csv")

yr.only<-lmer(log.abund ~ -1+code+code:year + code:logab.py+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay + (1|cell)+(1|CountID:cell), data=naba.1)
summary(yr.only)              


plot_model(best.model, type="eff", terms=c("ObsDay","code","MonthF"))

plot_model(best.model, type="eff", terms=c("on.dev","code"))

plot_model(best.model, type="eff", terms=c("dur.dev","code"))
plot_model(best.model, type="eff", terms=c("cold.dev","code"))
plot_model(best.model, type="eff", terms=c("dev.pc1","code"))
plot_model(best.model, type="eff", terms=c("dev.pc2","code"))
plot_model(best.model, type="eff", terms=c("logab.py","code"))
plot_model(best.model, type="eff", terms=c("year","code"))
plot_model(best.model, type="eff", terms=c("on.dev","code","warmearly"))
plot_model(best.model, type="eff", terms=c("dev.pc1","on.dev","code"))

}

names(naba.1)

#### PREDICT FOR VISUALIZATIONS



lims.pc1.ondev<-naba.1 %>%
#  select(code,logab.py,cell, CountID,warmearly,on.dev,warmlateopen,MonthF,year,ObsDay,FR.dev) %>%
  mutate(codelabel=ifelse(code=="RE","BOE",ifelse(code=="RL","BOL","BOP"))) %>% 
#         logab.py=median(logab.py, na.rm=T),cell=median(cell, na.rm=T),
#         CountID=median(CountID, na.rm=T),) %>%
  group_by(code, codelabel) %>%
  summarize(minPC1=floor(min(dev.pc1, na.rm=T)*10)/10,maxPC1=ceiling(max(dev.pc1, na.rm=T)*10)/10,
            minPC2=floor(min(dev.pc2, na.rm=T)*10)/10,maxPC2=ceiling(max(dev.pc2, na.rm=T)*10)/10,
            minond=floor(min(on.dev, na.rm=T)*10)/10,maxond=ceiling(max(on.dev, na.rm=T)*10)/10,
            mindurd=floor(min(dur.dev, na.rm=T)*10)/10,maxdurd=ceiling(max(dur.dev, na.rm=T)*10)/10,
            minyr=floor(min(year, na.rm=T)*10)/10,maxyr=ceiling(max(year, na.rm=T)*10)/10,
            mincold=floor(min(cold.dev, na.rm=T)*10)/10,maxcold=ceiling(max(cold.dev, na.rm=T)*10)/10  )

  
  

newcode<-c("RE","RL","RP")
codelabel<-c("EO","LO","PO")
nrep<-1200

abund.newDat <- data.frame(cell = rep(507,nrep), CountID=rep(497,nrep),
                         logab.py=rep(median(naba.1$logab.py, na.rm=T),nrep),
                         cold.dev=rep(0,nrep),
                         abslag=rep(10,nrep),
                         dev.pc1=rep(0,nrep),
                         dev.pc2=rep(0,nrep),
                         on.dev=rep(0,nrep), dur.dev=rep(0,nrep),
                         MonthF=as.factor(rep(6,nrep)),ObsDay=rep(30,nrep),
                         year=rep(10, nrep),
                         code=rep(newcode, each=20,20),
                         codelabel=rep(codelabel, each=20,20)) %>%
  mutate(codelabel = factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel, as.numeric(code)))

dev.maxvals<-(as.numeric(apply(na.omit(naba.1[,c(2,12,18:20,22:26,28,34,40:41)]),2,FUN=max)))
dev.minvals<-(as.numeric(apply(na.omit(naba.1[,c(2,12,18:20,22:26,28,34,40:41)]),2,FUN=min)))
dev.intervals<-round((dev.maxvals-dev.minvals)/20,3)
dev.names<-names(naba.1)[c(2,12,18:20,22:26,28,34,40:41)]

#pc1 x year
newpc1<- dev.minvals[which(dev.names=="dev.pc1")]+0:19*dev.intervals[which(dev.names=="dev.pc1")]
newpc2<- dev.minvals[which(dev.names=="dev.pc2")]+0:19*dev.intervals[which(dev.names=="dev.pc2")]


best.model<-ab.yr.1c


newyr<-0:19
nrep<-1200
minpc1<-naba.1 %>% group_by(year) %>% summarize(minpc1=min(dev.pc1, na.rm=T),maxpc1=max(dev.pc1, na.rm=T))
remove.cols<-which(names(abund.newDat) %in% c("dev.pc1","year"))
abund.newDat1<-cbind(abund.newDat[,-remove.cols], dev.pc1=rep(newpc1,60), year=rep(newyr,each=60))
abund.newDat1$pred <- predict(best.model, abund.newDat1, allow.new.levels =T)

lims.ab<-naba.1 %>%
  group_by(code) %>%
  summarize(minyr=min(year, na.rm=T),maxyr=max(year, na.rm=T),
            minpc1=min(dev.pc1, na.rm=T),maxpc1=max(dev.pc1, na.rm=T),
            minpc2=min(dev.pc2, na.rm=T), maxpc2=max(dev.pc2, na.rm=T))

abund.newDat.F1<-inner_join(abund.newDat1,lims.ab) %>%
  mutate(f1=ifelse(year>=minyr,ifelse(year<=maxyr,1,0),0)+ifelse(dev.pc1>=minpc1,ifelse(dev.pc1<=maxpc1,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel = factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel, as.numeric(code)))


##for fig 
pc1.yr<-naba.1 %>% filter(year>0) %>%group_by(year) %>% 
  summarize(minpc1=min(dev.pc1, na.rm=T),maxpc1=max(dev.pc1, na.rm=T),
            minpc2=min(dev.pc2, na.rm=T),maxpc2=max(dev.pc2, na.rm=T))

newpc1<-dev.minvals[which(dev.names%in%"dev.pc1")]+0:19*(dev.intervals[which(dev.names%in%"dev.pc1")])
newcold<-dev.minvals[which(dev.names%in%"cold.dev")]+0:19*(dev.intervals[which(dev.names%in%"cold.dev")])
nrep<-1200

remove.cols<-which(names(abund.newDat) %in% c("dev.pc1","cold.dev"))
abund.newDat2<-cbind(abund.newDat[,-remove.cols], dev.pc1=rep(newpc1,60), cold.dev=rep(newcold,each=60))
abund.newDat2$pred <- predict(best.model, abund.newDat2, allow.new.levels =T)

lims.ab<-naba.1 %>%
  group_by(code) %>%
  summarize(minyr=min(year, na.rm=T),maxyr=max(year, na.rm=T),
            minpc1=min(dev.pc1, na.rm=T),maxpc1=max(dev.pc1, na.rm=T),
            minpc2=min(dev.pc2, na.rm=T), maxpc2=max(dev.pc2, na.rm=T),
            mincold=min(cold.dev, na.rm=T), maxcold=max(cold.dev, na.rm=T),
            minondev=min(on.dev, na.rm=T), maxondev=max(on.dev, na.rm=T))

abund.newDat.F2<-inner_join(abund.newDat2,lims.ab) %>%
  mutate(f1=ifelse(cold.dev>=mincold,ifelse(cold.dev<=maxcold,1,0),0)+ifelse(dev.pc1>=minpc1,ifelse(dev.pc1<=maxpc1,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel = factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel, as.numeric(code)))




library(viridis)
library(metR)
n11<-naba.1 %>% filter(year<18) %>% mutate(Year=year+2000) %>% group_by(cell, year, Year) %>% summarize (dev.pc1=mean(dev.pc1, na.rm=T))
F.abund0<-ggplot(data=n11, aes(x=Year, y=dev.pc1, group=Year)) + geom_boxplot() + labs(y=" ")                 

ab.col<-c("wheat1","coral4") ## coral4 ##viridis_pal(option="B")(10)[c(10,6)]
and1<-abund.newDat.F1 %>% mutate(pred1=exp(pred))
  
 
F.abund.1<-ggplot(data=and1, aes(x=year+2000, y=dev.pc1, fill=pred1)) + 
  geom_tile() +  labs(x="Year", y="Later greenup + colder <----- AD.PC1 -----> Early greenup + warmer") + 
  scale_fill_gradient(name="Butterflies per hour", low=ab.col[1], high=ab.col[2], limits=c(0,27)) + 
  #scale_fill_viridis(name="Log abundance", option="A", begin=0.2,end=1, limits=c(0.3,3.6)) + #min(abund.newDat.F1$pred),max(abund.newDat.F1$pred))) +
  geom_contour(aes(x=year+2000, y=dev.pc1, z=pred1), color="black", binwidth=2) + theme_classic() + 
  geom_text_contour(aes(z=pred1), nudge_y=-.2,binwidth=2) + theme(legend.position="bottom") + 
  #  geom_line(data=pc1.yr,  inherit.aes = FALSE,aes(x=year+2000, y=maxpc2), color="black") + 
#  geom_line(data=pc1.yr,  inherit.aes = FALSE,aes(x=year+2000, y=minpc2), color="white") + 
#  geom_point(data=minpc1, aes(x=year, y=minpc1), shape=24, fill=NA, color="white") + 
#  geom_point(data=minpc1, aes(x=year, y=maxpc1), shape=25,fill=NA, color="white") + 
  facet_wrap(~codelabel, nrow=3)
F.abund.1

#save(F.abund.1,file="output/abund.fig.1027.png")

#### FIG: pc1 and cold winter days
#pc1 x year

##for fig 
#F.abund.2<-ggplot(data=abund.newDat.F2, aes(x=cold.dev, y=dev.pc1, fill=pred)) + 
#  geom_tile() +  labs(x="Cold winter days", y="PC1 deviation") + 
#  #geom_point(x=2,y=0.5,shape=24, color="white", size=1.5) + 
#  scale_fill_gradient(name="Log abundance", low=ab.col[1], high=ab.col[2], limits=c(0.3,3.6)) + 
  #scale_fill_viridis(name="Log abundance", option="A", begin=0.2,end=1,limits=c(0.3,3.6)) + #
#  facet_wrap(~codelabel)
#F.abund.2

#save(F.abund.2,file="output/abund.fig.2.png")


#### FIG: winter and onset deviation
#pc1 x year
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
            minpc1=min(dev.pc1, na.rm=T),maxpc1=max(dev.pc1, na.rm=T),
            minpc2=min(dev.pc2, na.rm=T), maxpc2=max(dev.pc2, na.rm=T),
            mincold=min(cold.dev, na.rm=T), maxcold=max(cold.dev, na.rm=T),
            minondev=min(on.dev, na.rm=T)-1, maxondev=max(on.dev, na.rm=T)+1)

abund.newDat.F3<-inner_join(abund.newDat3,lims.ab) %>%
  mutate(f1=ifelse(on.dev>=minondev,ifelse(on.dev<=maxondev,1,0),0)+ifelse(dev.pc2>=minpc2,ifelse(dev.pc2<=maxpc2,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel = factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel, as.numeric(code)))

abund.newDat.F3<-abund.newDat.F3 %>% mutate(pred1=exp(pred))

##for fig 
library(viridis)
F.abund.3<-ggplot(data=abund.newDat.F3, aes(x=on.dev, y=cold.dev, fill=pred1)) + 
  geom_tile() +  labs(x="Onset deviation (wks)", y="Winter cold deviation") + 
  scale_fill_viridis(name="Butterflies per hour") + 
  facet_wrap(~codelabel)
F.abund.3


##for fig 
abund.newDat.F3<-abund.newDat.F3 %>% mutate(pred1=round(exp(pred),3), ondevday=on.dev*7)

#Our transformation function
scaleFUN <- function(x) round(x,2)

F.abund.3<-ggplot(data=filter(abund.newDat.F3), aes(y=ondevday, x=cold.dev, fill=pred1)) + 
  geom_tile() + 
  scale_fill_gradient(name="Butterflies per hour", low=ab.col[1], high=ab.col[2], limits=c(0,27)) + 
  #scale_fill_viridis(name="Log abundance", option="A", begin=0.2,end=1,limits=c(0,4)) + #limits = c(min(abund.newDat.F1$pred),max(abund.newDat.F1$pred))) +
  scale_x_continuous(name="Winter cold days deviation", limits=c(-4.5,6.5), breaks=c(-4,0,4)) + 
  scale_y_continuous(name="Adult onset deviation (days)", limits=c(-60,80), breaks=c(-50,0,50)) + 
  geom_contour2(aes(y=ondevday, x=cold.dev, z=pred1),breaks=MakeBreaks(binwidth = NULL, bins = 4, exclude = NULL),
                global.breaks = FALSE, color="black") + #,label.placer = label_placer_n(n = 3)) + 
  geom_text_contour(aes(y=ondevday, x=cold.dev, z=pred1,label=scaleFUN(..level..)),
                     global.breaks = FALSE, color="black") + theme_classic() + 
  theme(legend.position="bottom") +   
  #geom_text_contour(aes(z=pred1),breaks=MakeBreaks(binwidth = NULL, bins = 5, exclude = NULL),global.breaks = FALSE, nudge_y=-.2) + 
  facet_wrap(~codelabel, nrow=3)
F.abund.3

ows.colors<-viridis_pal()(8)[c(1,5,7)] #c('#440154','#39568C','#FDE725')  #viridis(20, option="D")[c(1,6,20)]

Fabund30<-ggplot(data=(naba.1), aes(x=on.dev, y=cold.dev, color=code)) + 
  geom_contour2(aes(x = on.dev, y = cold.dev), alpha = 0.8) + 
  #geom_point(aes(shape=code)) + 
  scale_color_manual(values=mycolors, aesthetics = c("colour","fill")) + 
  #scale_shape_manual(values=c(3,8,21)) + 
  theme(legend.position="none", axis.title=element_blank()) 
Fabund30

#from https://stackoverflow.com/questions/46068074/double-box-plots-in-ggplot2
plot.x <- ggplot(naba.1) + geom_boxplot(aes(code, on.dev))
plot.y <- ggplot(naba.1) + geom_boxplot(aes(code, cold.dev))
plot.x <- layer_data(plot.x)[,1:6]
plot.y <- layer_data(plot.y)[,1:6]
colnames(plot.x) <- paste0("x.", gsub("y", "", colnames(plot.x)))
colnames(plot.y) <- paste0("y.", gsub("y", "", colnames(plot.y)))
df <- cbind(plot.x, plot.y); rm(plot.x, plot.y)
df$category <- sort(unique(naba.1$code))

df.outliers <- df %>%
  dplyr::select(category, x.middle, x.outliers, y.middle, y.outliers) %>%
  data.table::data.table()
df.outliers <- df.outliers[, list(x.outliers = unlist(x.outliers), y.outliers = unlist(y.outliers)), 
                           by = list(category, x.middle, y.middle)]

#
Fabund31<-ggplot(df, aes(fill = category, color = category)) +
  
  # 2D box defined by the Q1 & Q3 values in each dimension, with outline
  geom_rect(aes(xmin = x.lower, xmax = x.upper, ymin = y.lower, ymax = y.upper), alpha = 0.3) +
  geom_rect(aes(xmin = x.lower, xmax = x.upper, ymin = y.lower, ymax = y.upper), 
            color = "black", fill = NA) +
  
  # whiskers for x-axis dimension with ends
  geom_segment(aes(x = x.min, y = y.middle, xend = x.max, yend = y.middle)) + #whiskers
  geom_segment(aes(x = x.min, y = y.lower, xend = x.min, yend = y.upper)) + #lower end
  geom_segment(aes(x = x.max, y = y.lower, xend = x.max, yend = y.upper)) + #upper end
  
  # whiskers for y-axis dimension with ends
  geom_segment(aes(x = x.middle, y = y.min, xend = x.middle, yend = y.max)) + #whiskers
  geom_segment(aes(x = x.lower, y = y.min, xend = x.upper, yend = y.min)) + #lower end
  geom_segment(aes(x = x.lower, y = y.max, xend = x.upper, yend = y.max)) + #upper end
  scale_color_manual(values=mycolors, aesthetics=c("colour","fill"), guide="none" ) +
  # outliers 
  #geom_point(data = df.outliers, aes(x = x.outliers, y = y.middle), size = 3, shape = 1) + # x-direction
  #geom_point(data = df.outliers, aes(x = x.middle, y = y.outliers), size = 3, shape = 1) + # y-direction
  theme(legend.position="none", axis.title = element_blank()) +
  xlab(" ") + ylab(" ") +
  #coord_cartesian(xlim = c(-10,10), ylim = c(-10,10)) +
  theme_classic() 


Fabund31
###



F.abund.3a<-ggplot(data=filter(abund.newDat.F3, code=="RE"), aes(x=on.dev, y=cold.dev, fill=pred)) + 
  geom_tile() +  labs(y=" ") + 
  scale_fill_gradient(name="Log abundance", low=ab.col[1], high=ab.col[2], limits=c(0.3,3.6)) + 
  #scale_fill_viridis(name="Log abundance", option="A", begin=0.2,end=1,limits=c(0,4)) + #limits = c(min(abund.newDat.F1$pred),max(abund.newDat.F1$pred))) +
  theme(legend.position="none", axis.title.x=element_blank()) +  
  xlim(-10,12) + ylim(-4,8) + 
  geom_line(data=filter(lims.obs, code=="RE"),  inherit.aes = FALSE,aes(x=ond1, y=maxwc), color="white") + 
  geom_line(data=filter(lims.obs, code=="RE"),  inherit.aes = FALSE,aes(x=ond1, y=minwc), color="white") + 
  geom_contour(aes(x=on.dev, y=cold.dev, z=pred), color="black", binwidth=.2) + 
  geom_text_contour(aes(z=pred), nudge_y=-.2,binwidth=.2) + 
  facet_grid(~codelabel)
F.abund.3a
F.abund.3b<-ggplot(data=filter(abund.newDat.F3, code=="RL"), aes(x=on.dev, y=cold.dev, fill=pred)) + 
  geom_tile() +  labs( y="Std. # days below 0 C (winter)") + 
  scale_fill_gradient(name="Log abundance", low=ab.col[1], high=ab.col[2], limits=c(0.3,3.6)) + 
  #scale_fill_viridis(name="Log abundance", option="A", begin=0.2,end=1,limits=c(0,4)) + #limits = c(min(abund.newDat.F1$pred),max(abund.newDat.F1$pred))) +
  theme(legend.position="none", axis.title.x=element_blank()) + 
  xlim(-10,12) + ylim(-4,8) + 
  geom_line(data=filter(lims.obs, code=="RL"),  inherit.aes = FALSE,aes(x=ond1, y=maxwc), color="white") + 
  geom_line(data=filter(lims.obs, code=="RL"),  inherit.aes = FALSE,aes(x=ond1, y=minwc), color="white") + 
  geom_contour(aes(x=on.dev, y=cold.dev, z=pred), color="black", binwidth=.1) + 
  geom_text_contour(aes(z=pred), nudge_y=-.2,binwidth=.1) + 
  facet_grid(~codelabel)
F.abund.3b
F.abund.3c<-ggplot(data=filter(abund.newDat.F3, code=="RP"), aes(x=on.dev, y=cold.dev, fill=pred)) + 
  geom_tile() +  labs(x="Onset deviation (wks) ") + 
  scale_fill_gradient(name="Log abundance", low=ab.col[1], high=ab.col[2], limits=c(0.3,3.6)) + 
  #scale_fill_viridis(name="Log abundance", option="A", begin=0.2,end=1,limits=c(0,4)) + #limits = c(min(abund.newDat.F1$pred),max(abund.newDat.F1$pred))) +
  xlim(-10,12) + ylim(-4,8) + theme(axis.title.y = element_blank(),legend.position="bottom") + 
  geom_line(data=filter(lims.obs, code=="RP" & ond1),  inherit.aes = FALSE,aes(x=ond1, y=minwc), color="white") + 
  geom_line(data=filter(lims.obs, code=="RP"),  inherit.aes = FALSE,aes(x=ond1, y=maxwc), color="white") +   
  geom_contour(aes(x=on.dev, y=cold.dev, z=pred), color="black", binwidth=.1) + 
  geom_text_contour(aes(z=pred), nudge_y=-.2,binwidth=.1) + 
  facet_wrap(~codelabel)
F.abund.3c


library(gridExtra)
(figabund<-grid.arrange(F.abund.3,F.abund.1, nrow=1, ncol=2, heights=c(7.5,2))) 
ggsave(figabund,file="output/figures/fig3.10.png", height=10, width=7)

(figabund<-grid.arrange(F.abund.1,F.abund.3, nrow=1, ncol=2)) 
ggsave(figabund,file="output/figures/fig3.11.png", height=10, width=7)

                        

newwarmearly<-round(c( ((min(naba.1$dev.pc1, na.rm=T))*5):((max(naba.1$dev.pc1, na.rm=T))*5)/5),2)
newab.py<-mean(naba.1$logab.py, na.rm=T)
newpc2<-round(c( ((min(naba.1$dev.pc2, na.rm=T))*10):((max(naba.1$dev.pc2, na.rm=T))*10)/10),2)
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
                     dev.pc1=rep(newwarmearly, 3*length(newpc2)),
                     dev.pc2=rep(newpc2,each=length(newwarmearly),3),
                     code=rep(newcode, each=length(newwarmearly)*length(newpc2)),
                     codelabel=rep(codelabel, each=length(newwarmearly)*length(newpc2)))

newDatpc$pred <- predict(ab.yr.2, newDatpc,allow.new.levels =T)
library(viridis)

newDat1<-inner_join(newDatpc,lims.pc1.ondev) %>%
  mutate(f1=ifelse(dev.pc1>=minPC1,ifelse(dev.pc1<=maxPC1,1,0),0)+ifelse(on.dev>=minond,ifelse(on.dev<=maxond,1,0),0)) %>%
  filter(f1==2) %>% mutate(codelabel=factor(codelabel),code=factor(code)) %>%
  mutate(codelabel=fct_reorder(codelabel,as.numeric(code)))


ggplot(data=newDat1, aes(x=dev.pc2, y=dev.pc1, fill=pred)) + 
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
                       dev.pc1=rep(0, nrep),
                       dev.pc2=rep(0, nrep),
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
                          dev.pc1=rep(0, nrep),
                          dev.pc2=rep(0, nrep),
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



#######################3
abundx<-na.omit(naba.1)
ab.full<-lmer(log.abund~-1+code*(dev.pc1+logab.py+dev.pc2+on.dev+abslag+dur.dev+MonthF*ObsDay+year+cold.dev)+(1|cell) + (1|CountID:cell), data=abundx, na.action=na.fail)


xd1<-dredge(ab.full, beta = c("none"), evaluate = TRUE,
            rank = "AICc", fixed = c("code"), m.lim = NULL,
            trace = FALSE)

summary(xd1)
xd2<-xd1[xd1$weight>0.01,]
x3<-model.avg(xd2)

print(x3)

