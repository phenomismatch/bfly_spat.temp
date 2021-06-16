#Occurrence data phenology Metrics
## Resident butterfly species groups defined by overwinter stage (OWS)
## Phenometrics Onset (1%), Median (50%), Termination (99%), Duration (99%-1%)
## Code written by M Belitz, UFL
## Species OWS traits compiled by GU Ries Lab
# Modifications by E Larsen, Georgetown U, Updated 2021-06

library(tidyverse)
library(splitstackshape)

# read adult bfly pdfs
pdfs <- read.csv("data/derived/simpleton_pheno_pdfs-OutlierDetection_CIs.csv")

head(pdfs)

# split into list by year, overwintering code, and hexcell
pdf_list <- split(pdfs, 
                  f = list(pdfs$year, 
                           pdfs$code,
                           pdfs$HEXcell), 
                  drop = TRUE)

# write function to calculate 10% & 50% phenometrics
pheno_fun <- function(x){
  
  pdf_df <- pdf_list[[x]] %>% 
    mutate(
      cum_prob = cumsum(meanPDF),
      cum_perc = cum_prob / max(cum_prob),
      onset = DOY[which.max(cum_perc >= 0.01)],
      fiftieth = DOY[which.max(cum_perc >= 0.5)],
      term = DOY[which.max(cum_perc >= 0.99)],
      
      cum_prob_lowCI = cumsum(PDF2.5),
      cum_perc_lowCI = cum_prob_lowCI / max(cum_prob_lowCI),
      onset_highCI = DOY[which.max(cum_perc_lowCI >= 0.01)],
      fiftieth_lowCI = DOY[which.max(cum_perc_lowCI >= 0.5)],
      term_lowCI = DOY[which.max(cum_perc_lowCI >= 0.99)],
      
      cum_prob_highCI = cumsum(PDF97.5),
      cum_perc_hichCI = cum_prob_highCI / max(cum_prob_highCI),
      onset_lowCI = DOY[which.max(cum_perc_hichCI >= 0.01)],
      fiftieth_highCI = DOY[which.max(cum_perc_hichCI >= 0.5)],
      term_highCI = DOY[which.max(cum_perc_hichCI >= 0.99)],
      
      fiftieth_CI_days = fiftieth_highCI - fiftieth_lowCI,
      onset_CI_days = onset_highCI - onset_lowCI,
      term_CI_days = term_highCI - term_lowCI
    ) %>% 
    dplyr::select(year, HEXcell, code, obsDays, abundance,
                  onset, onset_lowCI, onset_highCI, 
                  fiftieth, fiftieth_lowCI, fiftieth_highCI,
                  term, term_lowCI, term_highCI,
                  onset_CI_days, fiftieth_CI_days, term_CI_days) 
  
  pdf_df2 <- pdf_df[1,]
  
  return(pdf_df2)
}

# calculate phenometrics for each grouping
pheno_list <- lapply(X = 1:length(pdf_list), FUN = pheno_fun)
# combine outputs into df
pheno_df_output <- do.call(bind_rows, pheno_list)

write.csv(x = pheno_df_output, "data/derived/adult_bfly_metrics.csv", row.names = F)


# data exploration, examining the effect of observation days and total number of 
# observations on the number of days the CI spans

#ggplot(pheno_df_output) +
#  geom_point(mapping = aes(x = abundance, y = abs(onset_CI_days), color = code)) +
#  geom_smooth(mapping = aes(x = abundance, y = abs(onset_CI_days), color = code)) +
#  geom_hline(yintercept = 0) +
#  scale_x_log10() +
#  theme_bw()
#
#ggplot(pheno_df_output) +
#  geom_point(mapping = aes(x = obsDays, y = abs(onset_CI_days), color = code)) +
#  geom_smooth(mapping = aes(x = obsDays, y = abs(onset_CI_days), color = code)) +
#  geom_hline(yintercept = 0) +
#  scale_x_log10() +
#  theme_bw()

# read in data
bfly_pheno <- pheno_df_output %>% 
  mutate(cell_code = paste(code, HEXcell, sep = "_")) # add hexcell for joining

#filter to 2016 - 2020
pheno_2016 <- filter(bfly_pheno, year >= 2016)

consq_yrs2 <- pheno_2016 %>% 
  group_by(HEXcell, code) %>% 
  summarise(years = n(), 
            mean10 = mean(onset),
            mean50 = mean(fiftieth),
            mean90 = mean(term))

tot_2016.2020 <- filter(consq_yrs2, years == 5) %>% 
  mutate(cell_code = paste(code, HEXcell, sep = "_")) # add hexcell for joining/antijoining

head(tot_2016.2020)

#### Filter to only keep cells with at least 5 years
bfly_pheno_prun <- bfly_pheno %>% 
  filter(cell_code %in% tot_2016.2020$cell_code)

bfly_pheno_mean <- left_join(bfly_pheno_prun, tot_2016.2020)

# calculate deviations from mean
bfly_pheno_mean <- bfly_pheno_mean %>% 
  mutate(dev_10 = onset - mean10,
         dev_50 = fiftieth - mean50,
         dev_90 = term - mean90)

write.csv(bfly_pheno_mean,
          "data/derived/adult_bfly_metrics_deviations.csv", row.names = F)

## do some real quick visualizations
#ggplot(bfly_pheno_mean) + 
#  geom_point(mapping = aes(x = year, y = dev_10, color = code)) +
#  geom_smooth(mapping = aes(x = year,y = dev_10, color = code), method = "lm") + 
#  facet_wrap(~HEXcell)
#
#ggplot(bfly_pheno_mean) + 
#  geom_point(mapping = aes(x = year, y = dev_50, color = code)) +
#  geom_smooth(mapping = aes(x = year,y = dev_50, color = code), method = "lm") + 
#  facet_wrap(~HEXcell)
#

# There are some major outliers. These should probably be removed before analysis.
# i.e., remove any years that have have a deviation >60 (?) days 