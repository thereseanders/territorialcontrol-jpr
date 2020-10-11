##########################################################
# Territorial control
# Computing event influence using distance and time decay
# Colombia (main specification)
#
# Therese Anders
##########################################################

library(tidyverse)
library(parallel)
library(sf)

source(here::here("utils","func_weights.R"))

#######################################
# Reading in data (25km grid, colombia)
centroids <- readRDS(here::here("prepped","grids_hex_col_hex25.rds")) %>%
  filter(countrypoly == 1) %>%
  mutate(longitude = X,
         latitude = Y) 

# table(centroids$countrypoly)
events <- readRDS(here::here("prepped","events_col_noassass_hex25.rds")) %>%
  filter(year %in% seq(1997, 2017)) %>%
  mutate(longitude = as.numeric(longitude),
         latitude = as.numeric(latitude)) %>%
  st_as_sf()

events_trans <- events %>%
  mutate(time = paste(month, year, sep = "-"),
         lonr_v = as.numeric(deg2rad(longitude)),
         latr_v = as.numeric(deg2rad(latitude))) %>%
  dplyr::select(lonr_v,
                latr_v,
                time,
                type,
                year,
                month)

df_months <- data.frame(time = apply(expand.grid(seq(1,12), 
                                                 seq(1997, max(events_trans$year))), 
                                     1, paste, collapse = "-")) %>%
  mutate(index = seq(1,nrow(.)))

events_new <- events_trans %>%
  left_join(df_months, by = "time")

events_terr <- events_new %>%
  filter(type == "terrorism") %>%
  arrange(year, month)

events_comb <- events_new %>%
  filter(type == "conventional") %>%
  arrange(year, month)


centroids_new <- centroids %>%
  mutate(lonr_c = as.numeric(deg2rad(longitude)),
         latr_c = as.numeric(deg2rad(latitude))) %>%
  dplyr::select(gid,
                lonr_c,
                latr_c,
                countrypoly,
                colwest)

# spatial 7, -0.35; temporal 8, -2.5
df_s_t <- func_decay(spatial_a = 7,
                     spatial_b = -0.35,
                     temporal_a = 8, 
                     temporal_b = -2.5, 
                     countrypoly = T)

df_s_t_trunc <- df_s_t %>%
  mutate(weighted_terrorism_trunc = ifelse(weighted_terrorism < 0.05, 0, weighted_terrorism),
         weighted_combats_trunc = ifelse(weighted_combats < 0.05, 0, weighted_combats))

saveRDS(df_s_t_trunc, here::here("prepped","events_decay_col_noassass_hex25_logistic_trunc.rds"))

