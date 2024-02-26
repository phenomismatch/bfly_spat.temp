# Butterfly_Abundance_and_Phenology

to do: Add Zenodo DOI

Code for analyzing drivers of aggregate butterfly phenology and abundance by overwintering stage

**Analysis published in:**

Larsen, EA, MW Belitz, GJ Di Cecco, J Glassberg, AH Hurlbert, L Ries, & RP Guralnick. Accepted 2023-12. [Overwintering strategy regulates phenological sensitivity and consequences for ecological services in a clade of temperate North American insects ](https://besjournals.onlinelibrary.wiley.com/journal/13652435) __*Functional Ecology*__


**Repository structure:**

* `data/` - Datasets relevant for project
  * `derived/`  - data that have undergone processing


  * `envir/`  - environmental greenup phenology data
    * `MidGreenup-2021-11-30-forest.rds` - input file for forest canopy greenup dates
    * `MidGreenup-2021-11-30-open.rds` - input file for open canopy greenup dates

  * `images/`  - cartoon art used in figure legends
    * `boe_art.png` - image file for representing butterflies overwintering as eggs
    * `bol_art.png` - image file for representing butterflies overwintering as larvae
    * `bop_art.png` - image file for representing butterflies overwintering as pupae

  * `maps/`  - files for spatial alignment &  aggregation
    * `hex_grid_crop.shp` - shapefiles for hexagonal grid used for data aggregation, from [Bird Phenology repo] (https://github.com/phenomismatch/Bird_Phenology/tree/master/Data/hex_grid_crop)
    * `spatial.domain.RData` - R data file for spatial domain of analysis
    
  * `naba/` - public data related to NABA circles in analysis
    * `naba_circles.csv` - spreadsheet of 777 NABA circles with mapping to hex cell IDs. Includes ID, Name, State,	Country, Lat,	Lng, hex cell
    * `naba_names.csv` - spreadsheet of 780 scientific name entries in NABA data mapped to Species IDs.  Includes ScientificName, n (# rows in raw datasheet), heirarchy (0 = above species, 1 = species, 2 = subspeices), SpeciesID (species level name)


* `code/` - Datasets relevant for project
  * `abundance_analysis.R` - code for analysis of drivers of abundance metrics
  * `cpc.process.R` - code for calculating thermal metrics from CPC data
  * `figures.R` - code for Figure 1 and supplemental figures
  * `gather_envir_vars.R` - code for preparing environmental covariates for analysis
  * `gdd_variations.R` - code for supplemental analysis of variations in GDD if using different thresholds, or considering within-hex-cell variation of daily tmin and tmax values
  * `phenology_analysis.R` - code for analyzing drivers of phenology metrics
  * `phenology_metrics.R` - code for estimating phenology metrics
  
  * `src/`  - reference code files from previous analysis
    * `degday1.R` - function for calculating single-sine estimation of degree days

* `output/` - output files from analysis (figures, tables)
  * `abundance.finalmodel.csv`  - csv spreadsheet of model parameters and statistics for final abundance model
  * `onsetDeviation.finalmodel.csv` - csv spreadsheet of model parameters and statistics for final model of adult onset phenology (start of flight period) deviation

  * `figures/`  - directory for image files of output figures
    * `Fig1.2024.png` - final Figure 1 (data summary)
    * `Fig2.2024.png` - final Figure 2 (onset deviation model)
    * `Fig3.2024.png` - final Figure 3 (abundance model)
    * `FigS31.png` - supplementary figure 3.1
    * `FigS32.png` - supplementary figure 3.2
    * `FigS33.png` - supplementary figure 3.3
