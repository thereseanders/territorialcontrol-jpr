############################
# Territorial Control paper
# Prio GRID data cleaning
#
# Therese Anders
############################

library(rgeos)
library(raster)
library(sf)
library(tidyverse)

##################
#
# Reading data
# 
##################


# Creating a subset of the PRIO Grid 2.0 data (for population values)
# See README.md for details
# downloaded from https://grid.prio.org/#/download
# PRIO GRID is WGS84
# priorgrid <- readr::read_csv("~/raw/PRIO-GRID\ Yearly\ Variables\ for\ 1989-2014\ -\ 2016-11-20.csv")

# Data available for years 1990, 1995, 2000, and 2005
priogrid_sub <- priorgrid %>%
  dplyr::select(priogrid = gid, 
                year,
                gwno,
                pop_gpw_sum) %>%
  mutate(case = case_when(
    gwno == 475 ~ "nga",
    gwno == 100 ~ "col"
  )) %>%
  filter(!is.na(case),
         !is.na(pop_gpw_sum)) %>%
  mutate(ln_pop = log(pop_gpw_sum))
saveRDS(priogrid_sub, here::here("prepped", "priogrid_pop_sub.rds"))

