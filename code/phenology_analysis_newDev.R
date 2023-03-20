#Analysis of Presence-only Phenology Metrics
#eButterfly, NABA Butterflies I've Seen, iNaturalist (research grade)
## Resident butterfly species groups defined by overwinter stage (OWS)
## Phenometrics estimated from Weibull (M Belitz)
## Species OWS traits compiled by GU Ries Lab
#E Larsen, Georgetown U, Updated 2022-01
#2022 Sept - Oct : adding precipitation
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
ows.colors<-mycolors
## Output files
onset.dev.output.csv<-"output/ondev.params23.csv"
onset.dev.model.file<-"output/ondev.model23.RData"

#Load phenology data
load(pheno.data)

(dev.hexyrs<-pheno.dev %>% dplyr::select(year, cell)  %>% group_by(year, cell) %>% tally())
nrow(pheno.dev)

#Load environmental data
load("data/derived/precipSum.RData")
env.dev<-read_csv(env.dev.csv)
precip.s1<-precip.sum %>% rename(year=yr.1) %>%
  #mutate(qtr=ifelse(season>1,ifelse(season>2,ifelse(season==4,"fourth","third"),"second"),"first")) %>%
  pivot_wider(id_cols=c(year, cell), names_from=season,values_from=pcp.dev.std)
env.dev<-merge(env.dev, precip.s1, all.x=T, by=intersect(names(env.dev), names(precip.s1)))
save(env.dev, file="data/derived/pheno.input23.RData")

pheno.input.dev<-merge(pheno.dev, env.dev, by=c("year","cell"), all.x=T) %>%
  filter(code %in% c("RE","RL","RP"), year<2018) %>%
  mutate(code=ifelse(code=="RE","EO",ifelse(code=="RL","LO","PO"))) %>%
  mutate(year=year-2000, onset.ci.wt=1/(onset.ci+1),med.ci.wt=1/(median.ci+1),dur.ci.wt=1/(dur.ci+1))

save(pheno.input.dev, file="data/derived/pheno.dev.modelinput23.RData")


ggplot(data=pheno.input.dev, aes(x=year, y=onset.dev, color=code)) + geom_jitter() + scale_color_manual(values=mycolors)
ggplot(data=pheno.input.dev, aes(x=year, y=median.dev, color=code)) + geom_jitter() + scale_color_manual(values=mycolors)
ggplot(data=pheno.input.dev, aes(x=year, y=dur.dev, color=code)) + geom_jitter() + scale_color_manual(values=mycolors)

###########################################################################
# Phenology model: Annual deviance of onset -------------------------------

#weighted LMER with interactions of variables with overwintering code
#random effect of cell
onset.dev.full<-lmer(onset.dev~-1+code+warm.dev + dev.pc1 + dev.pc2 + open.lag.dev + year+ uniqObsDays + code:warm.dev + code:open.lag.dev + code:dev.pc1 + code:dev.pc2 + code:year  + code:uniqObsDays + code:(winter.p + spring.p) + (1|cell), data=pheno.input.dev, weights=onset.ci.wt)
onset.dev.full<-lmer(onset.dev~-1+code*(warm.dev + dev.pc1 + dev.pc2 + open.lag.dev + year+ uniqObsDays + winter.p + spring.p) + (1|cell), data=pheno.input.dev, weights=onset.ci.wt)

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
onset.dev.final<-lmer(onset.dev~ -1 + code + warm.dev + dev.pc1 + uniqObsDays +
                       code:dev.pc2 + code:year + code:winter.p + code:spring.p + (1|cell), 
                      data=pheno.input.dev, weights=onset.ci.wt)
onset.dev.best<-get_model(step(onset.dev.final))
summary(onset.dev.best)
extractAIC(onset.dev.final)
extractAIC(onset.dev.best)
r.squaredGLMM(onset.dev.best)     

#check for variable collinearity
test.vif<-lm(onset.dev~code+year+dev.pc1+dev.pc2+year+warm.dev+winter.p + spring.p+uniqObsDays, data=pheno.input.dev)
vif(test.vif)

#look at partial r2
onset.dev.pc1<-lmer(onset.dev~ -1 + code + uniqObsDays + warm.dev + code:dev.pc2 + code:year + code:winter.p + code:spring.p + (1|cell), data=pheno.input.dev, weights=onset.ci.wt)
print(paste("partial r2 for pc1:",round(r.squaredGLMM(onset.dev.final)[1]-r.squaredGLMM(onset.dev.pc1)[1],3)))
onset.dev.pc2<-lmer(onset.dev~ -1 + code + uniqObsDays + dev.pc1 + warm.dev + code:year + code:winter.p + code:spring.p + (1|cell), data=pheno.input.dev, weights=onset.ci.wt)
print(paste("partial r2 for pc2:",round(r.squaredGLMM(onset.dev.final)[1]-r.squaredGLMM(onset.dev.pc2)[1],3)))
onset.dev.yr<-lmer(onset.dev~ -1 + code + uniqObsDays + dev.pc1 + warm.dev+code:dev.pc2 + code:winter.p + code:spring.p + (1|cell), data=pheno.input.dev, weights=onset.ci.wt)
print(paste("partial r2 for year:",round(r.squaredGLMM(onset.dev.final)[1]-r.squaredGLMM(onset.dev.yr)[1],3)))
onset.dev.eff<-lmer(onset.dev~ -1 + code  + dev.pc1 + warm.dev+code:dev.pc2 + code:year + code:winter.p + code:spring.p + (1|cell), data=pheno.input.dev, weights=onset.ci.wt)
print(paste("partial r2 for data density:",round(r.squaredGLMM(onset.dev.final)[1]-r.squaredGLMM(onset.dev.eff)[1],3)))
onset.dev.wp<-lmer(onset.dev~ -1 + code + uniqObsDays + dev.pc1 + warm.dev+code:dev.pc2 + code:year + code:spring.p + (1|cell), data=pheno.input.dev, weights=onset.ci.wt)
print(paste("partial r2 for winter precip:",round(r.squaredGLMM(onset.dev.final)[1]-r.squaredGLMM(onset.dev.wp)[1],3)))
onset.dev.sp<-lmer(onset.dev~ -1 + code + uniqObsDays + dev.pc1 + warm.dev+code:dev.pc2 + code:year + code:winter.p + (1|cell), data=pheno.input.dev, weights=onset.ci.wt)
print(paste("partial r2 for spring precip:",round(r.squaredGLMM(onset.dev.final)[1]-r.squaredGLMM(onset.dev.sp)[1],3)))

#save results and parameter table
(onset.params<-as_tibble(summary(onset.dev.best)$coefficients) %>%
  mutate(Estimate=round(Estimate, 4), StdError=round(`Std. Error`,4),parameter=row.names(summary(onset.dev.final)$coefficients)))
write_csv(onset.params, file=onset.dev.output.csv)
#save(onset.dev.best, file=onset.dev.model.file)
plot_model(onset.dev.best,show.values = TRUE, value.offset = .5)


# Figure 2: Phenology (annual deviation) model results ---------------------------------------

(plotyr<-plot_model(onset.dev.final, type = "eff", terms = c("year", "code"), alpha=0.5, size=2, colors=mycolors, se=TRUE, title="Year") + 
 theme_classic() + ylim(-25,25) + 
    scale_x_continuous(breaks = c(5, 10, 15),label = c("2005", "2010", "2015")) + 
  labs(x="Year", y="\u25B3 Adult onset", color="Group", fill="Group") + 
    theme(legend.position="none"))

(plot.pc1<-plot_model(onset.dev.final, type = "eff", terms = c("dev.pc1", "code"), alpha=0.5, size=2, colors=mycolors, se=TRUE, title="Spring GDD + forest greenup") + 
    theme_classic() + ylim(-25,25) + 
    labs(x="PC1 --> early greenup (days), warm spring", y="\u25B3 Adult onset", color="Group", fill="Group") +
    theme(legend.position="none"))

(plot.pc2<-plot_model(onset.dev.final, type = "eff", terms = c("dev.pc2", "code"), alpha=0.5, size=2 , colors=mycolors, se=TRUE, title="GDD-Greenup decoupling") + 
    theme_classic()  + ylim(-25,25) + 
    labs(x="PC2 --> greenup later than expected (days)", y="", color="Group", fill="Group")+ 
    theme(legend.position="none")) 

(plot.warm<-plot_model(onset.dev.final, type = "eff", terms = c("warm.dev", "code"), alpha=0.5, size=2 , colors=mycolors, se=TRUE, title="Winter warm days") + 
    theme_classic()  + ylim(-25,25) + 
    labs(x="\u25B3 Winter days above 0C", y="", color="Group", fill="Group")+ 
    theme(legend.position="none")) 

(plot.winterp<-plot_model(onset.dev.final, type = "eff", terms = c("winter.p", "code"), alpha=0.5, size=2 , colors=mycolors, se=TRUE, title="Winter precipitation") + 
    theme_classic()  + ylim(-25,25) + 
    labs(x="\u25B3 Winter precipitation", y="", color="Group", fill="Group")+ 
    theme(legend.position="none")) 

(plot.springp<-plot_model(onset.dev.final, type = "eff", terms = c("spring.p", "code"), alpha=0.5, size=2 , colors=mycolors, se=TRUE, title="Spring precipitation") + 
    theme_classic()  +ylim(-25,25) + 
    labs(x="\u25B3 Spring precipitation", y="", color="Group", fill="Group")+ 
    theme(legend.position="none")) 

(legend.2<-ggplot(data=pheno.input.dev, aes(x=dev.pc1, y=onset.dev, color=code, fill=code)) + 
  geom_smooth(method="lm", alpha=0.5) + scale_color_manual(values=mycolors, aesthetics=c("color","fill"), labels=c("EO","LO","PO")) +
   theme(legend.position="right",legend.key.size = unit(1, 'cm'),legend.title = element_text(size=14)) + labs(color="Group", fill="Group") )

l2<-get_legend(legend.2)


(legend.2<-ggplot(data=pheno.input.dev, aes(x=dev.pc1, y=onset.dev, color=code, fill=code)) + 
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
    ylim(-.5,6.5) + xlim(0,8) +  theme_map(12) + labs(title="\u25B3 Adult flight onset") +
    annotate("text", x = 2, y = 4.5, label = c("Overwinter stage"), size=5) + 
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
pp.a1<-legend3


(figonset<-grid.arrange(pp.a1, plot.pc1, plot.pc2, plotyr,plot.winterp,plot.springp, nrow=2))
#ggsave(figonset, width=12,height=4, units="in", file=paste0("output/figures/onset.dev.result.23.png"))

ggsave(figonset, width=12,height=5, units="in", file=paste0("output/figures/onset.dev.result.23.png"))


(figonset<-grid.arrange(plotyr, plot.pc1, plot.pc2,l2, nrow=1, widths=c(1,1,1,.5)))
#ggsave(figonset, width=10,height=6, units="in", file=paste0("output/figures/FIG2pheno.png"))

figon.10<-grid.arrange(plotyr, plot.pc1, plot.pc2,plot.winterp,plot.springp,pp.a1, nrow=2, widths=c(1,1,1))
#ggsave(figon.10, width=10,height=4, units="in", file=paste0("output/figures/onset.rain.pdf"))
ggsave(figon.10, width=10,height=4, units="in", file=paste0("output/figures/onset.rain23.png"))
##############



###########################################################################
# Phenology model: spatiotemporal onset -------------------------------
#Load environmental data
env<-read_csv("data/derived/spatemp_env.csv")
env<-merge(env, precip.s1, all.x=T, by=intersect(names(env), names(precip.s1)))

load(pheno.data.st)
pheno.quant<- pheno.quant %>%
  mutate(qdur=q95-q5, qdur_low=q95_low-q5_high, qdur_high=q95_high-q5_low,
       onset.ci=q5_high-q5_low,median.ci=q50_high-q50_low,
       q95_ci=q95_high-q95_low,dur.ci=qdur_high-qdur_low)


pheno.input<-merge(pheno.quant, env, by=c("year","cell"), all.x=T) %>%
  filter(code %in% c("RE","RL","RP"), year<2018) %>%
  mutate(code=ifelse(code=="RE","EO",ifelse(code=="RL","LO","PO"))) %>%
  rename(onset=q5, median=q50, duration=qdur) %>%
  mutate(Year=year, year=year-2000, open.lag=gr_mn_open-forest.greenup,
         onset.ci.wt=1/(onset.ci+1),med.ci.wt=1/(median.ci+1),dur.ci.wt=1/(dur.ci+1))

save(pheno.input, file="data/derived/pheno.ST.input23.csv")

#weighted LMER with interactions of variables with overwintering code
#random effect of cell
onset.full<-lmer(onset~-1+code+warmdays + ST.PC1 + ST.PC2 + open.lag + spring.p+winter.p+year+ uniqObsDays + code:warmdays + code:open.lag + code:ST.PC1 + code:ST.PC2 + code:spring.p + code:winter.p + code:year  + code:uniqObsDays + (1|cell), data=pheno.input, weights=onset.ci.wt)
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
onset.st.final<-lmer(onset~ -1 + code + code:ST.PC1 + code:warmdays + code:year + winter.p + uniqObsDays + (1|cell), data=pheno.input, weights=onset.ci.wt)
summary(onset.st.final)
extractAIC(onset.st.final)
r.squaredGLMM(onset.st.final)     
write.csv(as.data.frame(summary(onset.st.final)$coefficients), file="output/onset.st.param23.csv")

#check for variable collinearity
test.vif<-lm(onset~code + ST.PC1 + warmdays + year + winter.p + uniqObsDays, data=pheno.input)
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

(plot.wp<-plot_model(onset.st.final, type = "eff", terms = c("winter.p", "code"), alpha=0.5, size=2 , colors=mycolors, se=TRUE, title="C. Onset Deviation ~ PC2") + 
    theme_classic()  +
    labs(x="Winter precipitation", y="Adult onset", color="Group", fill="Group")+ 
    theme(legend.position="none")) 


(legend.2<-ggplot(data=pheno.input.dev, aes(x=dev.pc1, y=onset.dev, color=code, fill=code)) + 
    geom_smooth(method="lm", alpha=0.5) + scale_color_manual(values=mycolors, aesthetics=c("color","fill"), labels=c("EO","LO","PO")) +
    theme(legend.position="right",legend.key.size = unit(1, 'cm'),legend.title = element_text(size=14)) + labs(color="Group", fill="Group") )

l2<-get_legend(legend.2)

(figonset<-grid.arrange(plotyr, plot.pc1, plot.pc2,l2, nrow=1, widths=c(1,1,1,.5)))
ggsave(figonset, width=12,height=4, units="in", file=paste0("output/figures/supp.onset.23.png"))


##############



###########################################################################
### PHENOLOGICAL MODELS: DURATION
###########################################################################
#ph2<-pheno.input
if(max(pheno.input$summer.gdd, na.rm=T)>1000) {
  pheno.input<-pheno.input %>% mutate(summer.gdd=summer.gdd/100)
}
dur.full<-lmer(duration~-1+code+code*(warmdays + ST.PC1 + ST.PC2 + winter.p + spring.p + summer.p + open.lag + summer.gdd + year)+ uniqObsDays + (1|cell), data=pheno.input, weights=dur.ci.wt)
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
dur.st.final<-lmer(duration ~ -1 + code + ST.PC2 + code:warmdays + code:open.lag + summer.p + 
                     code:spring.p + code:summer.gdd + code:year + uniqObsDays + (1|cell), 
                   data=pheno.input, weights=dur.ci.wt)
summary(dur.st.final)
extractAIC(dur.st.final)
r.squaredGLMM(dur.st.final)     
write.csv(as.data.frame(summary(dur.st.final)$coefficients), file="output/dur.st.param23.csv")


###########################################################################
### PHENOLOGICAL MODELS: DURATION DEVIATION
###########################################################################


dur.dev.full<-lmer(dur.dev~-1+code+ code*(warm.dev + summer.gdd.dev + dev.pc1 + dev.pc2  + spring.p + summer.p + open.lag.dev + year+ uniqObsDays) + (1|cell),
                   data=pheno.input.dev, weights=dur.ci.wt)
extractAIC(dur.dev.full)
r.squaredGLMM(dur.dev.full)     
(t2<-as_tibble(bind_cols(Parameters=row.names(summary(dur.dev.full)$coefficients),summary(dur.dev.full)$coefficients)) %>% arrange(abs(`t value`)))

dur.best<-get_model(step(dur.dev.full))

(t2<-as_tibble(bind_cols(Parameters=row.names(summary(dur.best)$coefficients),summary(dur.best)$coefficients)) %>% arrange(abs(`t value`)))
summary(dur.best)$call
dur.dev.final<-lmer(dur.dev~ -1 + code + warm.dev + summer.gdd.dev + code:dev.pc1 + code:dev.pc2 + code:spring.p + summer.p + code:open.lag.dev + code:year + code:uniqObsDays +  (1|cell), data=pheno.input.dev, weights=dur.ci.wt)
summary(dur.dev.final)
dur.dev.best<-get_model(step(dur.dev.final))
extractAIC(dur.dev.final)
r.squaredGLMM(dur.dev.final)     
extractAIC(dur.dev.best)
r.squaredGLMM(dur.dev.best)     
test.vif<-lm(dur.dev~code+warm.dev+dev.pc1+dev.pc2+spring.p+summer.p+open.lag.dev+summer.gdd.dev+year+uniqObsDays, data=pheno.input.dev)
vif(test.vif)
summary(dur.dev.best)

#plot_model(dur.dev.final, type = "eff", terms = c("warm.dev", "code"), title="durDev~winter warm days") + 
#  scale_fill_manual(values=ows.colors) + scale_color_manual(values=ows.colors)  
(pp1<-plot_model(dur.dev.final, type = "eff", terms = c("summer.gdd.dev", "code"), title="Summer GDD", alpha=0.5, size=2, colors=ows.colors, se=TRUE) +
  theme_classic() + ylim(-90,90) + theme(legend.position="none") + labs(x="\u25B3 Summer GDD",y="", color="Group", fill="Group"))

(pp.pc1<-plot_model(dur.dev.final, type = "eff", terms = c("dev.pc1", "code"), title="Spring GDD + forest greenup", alpha=0.5, size=2, colors=ows.colors, se=TRUE) +
    theme_classic() + ylim(-90,90) +theme(legend.position="none") + labs(x="PC2 --> warm spring, early greenup (days)",y="\u25B3 Flight period duration", color="Group", fill="Group"))
(pp2<-plot_model(dur.dev.final, type = "eff", terms = c("dev.pc2", "code"), title="GDD-Greenup decoupling", alpha=0.5, size=2, colors=ows.colors, se=TRUE) +
  theme_classic() + ylim(-90,90) +theme(legend.position="none") + labs(x="PC2 --> greenup later than expected (days)",y="", color="Group", fill="Group"))

(pp3<-plot_model(dur.dev.final, type = "eff", terms = c("spring.p", "code"), title="Spring precipitation", alpha=0.5, size=2, colors=ows.colors, se=TRUE) +
  theme_classic() +ylim(-90,90) + theme(legend.position="none") + labs(x="\u25B3 Spring precipitation",y="", color="Group", fill="Group"))
(pp4<-plot_model(dur.dev.final, type = "eff", terms = c("summer.p", "code"), title="Summer precipitation", alpha=0.5, size=2, colors=ows.colors, se=TRUE) +
  theme_classic() +ylim(-90,90) + theme(legend.position="none") + labs(x="\u25B3 summer precipitation",y="", color="Group", fill="Group"))
(pp5<-plot_model(dur.dev.final, type = "eff", terms = c("warm.dev", "code"), title="Winter warm days", alpha=0.5, size=2, colors=ows.colors, se=TRUE) +
  theme_classic() +ylim(-90,90) + theme(legend.position="none") + labs(x="\u25B3 days above 0C",y="", color="Group", fill="Group"))
(pp6<-plot_model(dur.dev.final, type = "eff", terms = c("open.lag.dev", "code"), title="Open canopy greenup lag", alpha=0.5, size=2, colors=ows.colors, se=TRUE) +
  theme_classic() + ylim(-90,90) +theme(legend.position="none") + labs(x="\u25B3 Open canopy lag",y="", color="Group", fill="Group") )
(pp7<-plot_model(dur.dev.final, type = "eff", terms = c("year", "code"), title="Year", alpha=0.5, size=2, colors=ows.colors, se=TRUE) +
  theme_classic() +ylim(-90,90) + theme(legend.position="none") + scale_x_continuous(breaks=c(5,10,15),labels=c("2005","2010","2015")) +
    labs(x="Year",y="\u25B3 Flight period duration", color="Group", fill="Group") )
(pp8<-plot_model(dur.dev.final, type = "eff", terms = c("uniqObsDays", "code"), title="N.Obs", alpha=0.5, size=2, colors=ows.colors, se=TRUE) +
  theme_classic() +ylim(-90,90) + theme(legend.position="none")+ labs(x="n. Obs Days",y="\u25B3 Flight pd duration", color="Group", fill="Group") )

(fig.dur<-grid.arrange(arrangeGrob(legend3,pp.pc1,pp2,pp6,pp7,pp3,pp4,pp1, nrow=2, heights=c(2,2))))
(fig.dur<-grid.arrange(arrangeGrob(legend3,pp3,pp7,pp.pc1,pp2,pp6, nrow=2, heights=c(2,2))))

ggsave(fig.dur, width=12,height=4, units="in", file=paste0("output/figures/dur.rain.pdf"))
ggsave(fig.dur, width=12,height=4, units="in", file=paste0("output/figures/dur.rain23.png"))

durdev.output<-as_tibble(summary(dur.dev.final)$coefficients) %>%
  dplyr::mutate(param=row.names(summary(dur.dev.final)$coefficients), sig=ifelse(`Pr(>|t|)`<0.05,1,0), r2m=round(r.squaredGLMM(dur.dev.final)[1],2),r2c=round(r.squaredGLMM(dur.dev.final)[2],2))
write.csv(durdev.output, file=paste0("output/dur.dev.model23",rundat,".csv"))


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
