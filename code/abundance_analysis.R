#Modeling patterns in butterfly abundance
#Elise Larsen, Georgetown U
#Run in R 4.1.2

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
abund.model.input<-"data/abund.input.RData"
abundance.model.fit<-"output/abundance.finalmodel.csv"
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


#add previous year abundance to input data
abund.py<-abund %>%
  mutate(year=year+1) %>%
  filter(year < 2018) %>%
  dplyr::select(cell, year, CountID, code, doy.py=doy, abund.py=abund.bph, logab.py=log.abund, seq.py=seq, nsurv.py=nsurveys) 

#Add previous year abundance to abundance table
abund1<-merge(abund, abund.py, by=c("year", "cell","CountID", "code"), all.x=T) %>%
  mutate(keep=ifelse(seq==seq.py | nsurv.py==1,1,0)) %>% filter(keep==1) %>%
  dplyr::select(year:SR, abund.py, logab.py)

summary(abund1)

#Add environmental deviations
env.dev<-read_csv(env.dev.csv)

#Add precipitation data
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


save(naba.1, file="data/derived/abund.model.input.RData")


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

### OUTPUT VISUALIZATION

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


# Model Selection Exploration---------------------------------------------------------
  ##REGSUBSETS MODELS
  b<-regsubsets(log.abund~-1+code+code:dev.pc1+code:dev.pc2+code:open.lag.dev+code:on.dev+code:abs.surveylag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+code:logab.py+code:Year+(1|cell), data=naba.1,nbest=1)

  ab.yr.1d<-lmer(log.abund~-1+code+code:dev.pc1+code:dev.pc2+code:on.dev+abs.surveylag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+logab.py+code:Year+(1|cell) + (1|CountID:cell), data=naba.1)
  summary(ab.yr.1d)
  extractAIC(ab.yr.1d)
  
  ab.yr.1e<-lmer(log.abund~-1+code+code:dev.pc1+code:dev.pc2+code:open.lag.dev+abs.surveylag+code:dur.dev+code:MonthF+code:cold.dev+logab.py+code:Year+(1|cell) + (1|CountID:cell), data=naba.1)
  summary(ab.yr.1e)
  extractAIC(ab.yr.1e)
  
  ab.yr.1b<-lmer(log.abund~-1+code+code:logab.py+(1|cell) + (1|CountID:cell), data=naba.1)
  extractAIC(ab.yr.1b)
  
  ab.yr.noint1<-lmer(log.abund~-1+code+dev.pc1+code:logab.py+code:dev.pc2+code:on.dev+code:abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+code:year+(1|cell) + (1|CountID:cell), data=naba.1)
  extractAIC(ab.yr.noint1)
  if(extractAIC(ab.yr.noint1)[2]<extractAIC(best.model)[2]) {best.model<-ab.yr.noint1}
  
  ab.yr.nointlag<-lmer(log.abund~-1+code+code:dev.pc1+code:logab.py+code:dev.pc2+code:on.dev+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+code:year+(1|cell) + (1|CountID:cell), data=naba.1)
  if(extractAIC(ab.yr.nointlag)[2]<extractAIC(best.model)[2]) {best.model<-ab.yr.nointlag}
  extractAIC(best.model)
  
  ab.yr.nointyr<-lmer(log.abund~-1+code+code:dev.pc1+code:logab.py+code:dev.pc2+code:on.dev+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+year+(1|cell) + (1|CountID:cell), data=naba.1)
  if(extractAIC(ab.yr.nointyr)[2]<extractAIC(best.model)[2]) {best.model<-ab.yr.nointyr}
  extractAIC(best.model)
 
  ab.yr.noyr<-lmer(log.abund~-1+code+code:dev.pc1+code:logab.py+code:dev.pc2+code:on.dev+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+(1|cell) + (1|CountID:cell), data=naba.1)
  if(extractAIC(ab.yr.noyr)[2]<extractAIC(best.model)[2]) {best.model<-ab.yr.noyr}
  extractAIC(best.model)

  ab.yr.nopyint<-lmer(log.abund~-1+code+code:dev.pc1+logab.py+code:dev.pc2+code:on.dev+abslag+code:dur.dev+code:MonthF + code:MonthF:ObsDay+code:cold.dev+(1|cell) + (1|CountID:cell), data=naba.1)
  if(extractAIC(ab.yr.nopyint)[2]<extractAIC(best.model)[2]) {best.model<-ab.yr.nopyint}
  extractAIC(best.model)

#output table  
abund.output<-as_tibble(summary(best.model)$coefficients) %>%
  dplyr::mutate(param=row.names(summary(best.model)$coefficients), sig=ifelse(`Pr(>|t|)`<0.05,1,0), r2m=round(r.squaredGLMM(best.model)[1],2),r2c=round(r.squaredGLMM(best.model)[2],2))
write.csv(abund.output, file="output/abund.model.csv")

