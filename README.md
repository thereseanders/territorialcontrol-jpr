<p float="left">
  <img src="https://github.com/thereseanders/territorialcontrol-jpr/blob/master/plots/hmm_col_yearly.png" width="325"/>
  <img src="https://github.com/thereseanders/territorialcontrol-jpr/blob/master/plots/hmm_nga_yearly.png" width="325"/>
</p>

# Replication code

The repository contains the data and code to replicate the analysis for: 

:page_facing_up: Anders, Therese (2020): [**Territorial control in civil wars: Theory and measurement using machine learning**](https://journals.sagepub.com/doi/abs/10.1177/0022343320959687?journalCode=jpra) In *Journal of Peace Research 57*(6): 701-714.

[:paperclip: **Online appendix**](http://therese.rbind.io/publication/anders2019_jpr/territorialcontrol_jpr_appendix_200608.pdf)


<p float="center">
  <img src="https://github.com/thereseanders/territorialcontrol-jpr/blob/master/plots/territorialcontrol.png" width="500"/>
</p>

---

This repository leverages the [`renv`](https://rstudio.github.io/renv/articles/renv.html) package for `R` package dependency management. 

*Note: The raw GTD and GED data is not included in this repository. Please follow the download instructions below to obtain the data from the original sources. The repository contains the relevant subset of events, as well as the code to reproduce the pre-cleaning steps to get from the raw data to the conflict event exposure values.*


## Pre-cleaning steps

### GTD (2018) data
Raw files `gtd_96to13_0718dist.csv` and `gtd_14to17_0718dist.csv` (downloaded from https://www.start.umd.edu/gtd/). A subset of the data for Nigeria and Colombia with attached PRIO GRID 2.0 grid cell id is shared under `./raw/gtd_attached_prioid.rds`. PRIO GRID IDs were attached to each event using file `priogrid_cell.shp` (downloaded from https://grid.prio.org/#/download). 

### GED data Global version 18.1
Downloaded from https://ucdp.uu.se/downloads/index.html#ged_global.

### ACLED
Note: some hand-coding was necessary to attribute ACLED events to the respective perpetrators and targets of violence (files available upon request). The repository contains the cleaned data for alternative aggregation schemes. 

### Steps
*(code but not raw data included)*

1. `./precleaning/hexgrid25_ged_gtd.R`: Subset data to relevant units, create hexagonal grid, and associate event data with each grid cell. Note: This file takes a while to run due to computationally intensive spatial operations. There is room for improving efficiency in future iterations of this code.
2. `./precleaning/weights_nganorth.R` and `./precleaning/weights_col.R` (as well as `./precleaning/weights_col_noassess.R` and `./precleaning/weights_col_noofficial.R` for robustness checks): Compute exposure of grid cells to conflict events using spatial and temporal weights.

## Compute HMM
Run `./est/hmm_col.R` `./est/hmm_nga.R` (as well as `./est/hmm_col_noassess.R` and `./est/hmm_col_noofficial.R` for robustness checks) to compute HMM. The analysis is powered by the [`HMM`](https://cran.r-project.org/web/packages/HMM/index.html) package.

## Graphs and tables

### Article
The code to reproduce all figures in the [article](https://journals.sagepub.com/doi/abs/10.1177/0022343320959687?journalCode=jpra#articleCitationDownloadContainer) can be found in the `./analysis/figures_tables.R` script.

### Online appendix
The code to reproduce all figures and tables in the appendix (with the exception of the Colombia validation exercise and sensitivity analysis) can be found in the `./analysis/figures_tables_appendix.R` script. 

Note that because the PRIO Grid 2.0 data is very large, only the relevant subset of the data for the analysis of the underreporting bias in the online appendix. Please recall that the PRIO GRID ids are added to the GTD data in a pre-cleaning step.

### Validation of Colombia estimates using deforestation data.
`./analysis/validation_col.R` performs the validation exercise. Deforestation are obtained from the *Instituto de Hidrologia, Meteorologiıa y Estudios Ambientales* and in geoTiff format and have been aggregated to hexagonal grid cells. Estimation of clustered standard errors based on [`clusterSEs`](https://github.com/cran/clusterSEs/blob/master/R/clusterBS.glm.R) package. Note that the estimation of clustered standard errors takes some time, so this code is shared separate from the `R` scripts to replicate the figures and tables in the article and online appendix.

### Sensitivity analysis
The script `./analysis/vary_emissions.R` performs sensitivity of the Nigeria HMM estimates to variations in the emission probabilities respectively. Note that the estimation of `N = 389` HMMs takes some time, therefore this code is shared separate from the scripts to replicate the figures in the online appendix. The code leverages the [`furrr`](https://github.com/DavisVaughan/furrr) package for multicore processing.

<img src="https://github.com/thereseanders/territorialcontrol-jpr/blob/master/plots/appendix/plot2_evvariance_rubin.png" width="1000"/> <figcaption>Sensitivity of estimates to changes in the matrix of emission probabilities (NE Nigeria 2013–2017). </figcaption>
