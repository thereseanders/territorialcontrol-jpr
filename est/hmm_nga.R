############################
# Territorial Control paper
# HMM analysis
# Nigeria
#
# Therese Anders
############################


library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(HMM)
library(lubridate)
library(VGAM)
library(sf)

##################################
#
#
# PART A) Preprocessing
#
#
##################################

source(here::here("utils","func_hmm.R"))

### Reading in data
grid_og <- readRDS(here::here("prepped","grids_hex_nga_hex25.rds")) %>%
  dplyr::select(-X, -Y) %>%
  filter(nganorth == 1)

crs <- st_crs(grid_og)
crs

gid_sub <- unique(grid_og)

events <- readRDS(here::here("prepped","events_decay_nga_hex25_logistic_trunc.rds"))

#####################
# new frame
to <- max(events$timeindex)

time <- data.frame(timeindex = seq(1,to),
                   month = rep(seq(1,12))) %>%
  mutate(year_trans = rep(seq(1, nrow(.)/12), each = 12)) %>%
  mutate(year = (2008 + year_trans)-1) %>%
  dplyr::select(-year_trans)

dat_raw <- events %>%
  st_as_sf(coords = c("lonr_c", "latr_c"),
           crs = crs) %>%
  left_join(time, by = "timeindex") %>%
  mutate(gwno = 475) %>%
  dplyr::select(gid, 
                gwno,
                year,
                month,
                timeindex,
                nganorth,
                gtd = weighted_terrorism_trunc,
                ged = weighted_combats_trunc) %>%
  dplyr::mutate(gtd = replace(gtd, is.na(gtd), 0),
                ged = replace(ged, is.na(ged), 0)) %>%
  dplyr::mutate(gtd_floor = floor(gtd),
                ged_floor = floor(ged),
                gtd_round = round(gtd),
                ged_round = round(ged)) %>%
  st_set_geometry(NULL)

dat_mean_info <- dat_raw %>%
  gather(indicator, value, gtd:ged) %>%
  
  # Computing monthly averages
  group_by(gwno, indicator, timeindex) %>%
  
  # Monthly means
  dplyr::summarise(mean = mean(value, na.rm =T),
                   var = var(value, na.rm = T)) %>%
  gather(statistic, value, mean:var) %>%
  unite(indc_stat, indicator, statistic) %>%
  spread(indc_stat, value) %>%
  dplyr::select(gwno,
                mean_t = gtd_mean,
                mean_c = ged_mean,
                var_t = gtd_var,
                var_c = ged_var,
                timeindex) %>%
  mutate(mean_t = replace(mean_t, is.na(mean_t), 0),
         mean_c = replace(mean_c, is.na(mean_c), 0),
         var_t = replace(var_t, is.na(var_t), 0),
         var_c = replace(var_c, is.na(var_c), 0))
summary(dat_mean_info$mean_t)
summary(dat_mean_info$mean_c)


dat <- dat_raw  %>%
  left_join(dat_mean_info, by = c("gwno", "timeindex")) %>%
  left_join(grid_og) %>%
  dplyr::mutate(cdf_T = pzipois(floor(gtd), mean_t), # This rounds gtd (ged) down to the next lowest full number
                cdf_C = pzipois(floor(ged), mean_c),
                gid = as.numeric(str_replace_all(gid, "ngahex", ""))) %>%
  filter(nganorth == 1)



##################################
#
#
# PART B) Setting HMM parameters
#
#
##################################

#############################
# Setting up empirical states
ev <- c("O1", "O2", "O3", "O4")

##########################################
# Setting up values of the latent variable
## States are the possible levels of territorial control (K = 5)
states <- c("R", "DR", "D", "DG", "G")

###############################
## Setting up evidence matrices
# Evidence probabilities (translating observations into probabilities, based on theory)

ev_mat <- matrix(nrow = length(states), ncol = length(ev))
ev_mat[1,] <- c(0.6, 0.175, 0.175, 0.05)
ev_mat[2,] <- c(0.05, 0.6, 0.175, 0.175)
ev_mat[3,] <- c(0.05, 0.175, 0.6, 0.175)
ev_mat[4,] <- c(0.05, 0.175, 0.175, 0.6)
ev_mat[5,] <- c(0.6, 0.05, 0.175, 0.175)

###############################
# Setting up initial probabilities
# For Nigeria: Strong prior for government control
d_init <- c(0.025, 0.025, 0.025, 0.025, 0.9)

####### Transition matrix ########
# Modified empirical transition matrix
trans_mat <- data.frame(matrix(nrow = length(states), ncol = length(states)))
trans_mat[1,] <- c(0.25, 0.5, 0.025, 0.2, 0.025)
trans_mat[2,] <- c(0.25, 0.15, 0.075, 0.5, 0.025)
trans_mat[3,] <- c(0.05, 0.025, 0.05, 0.85, 0.025)
trans_mat[4,] <- c(0.025, 0.075, 0.15, 0.125, 0.625)
trans_mat[5,] <-c(0.05, 0.075, 0.475, 0.025, 0.375)
row.names(trans_mat) <- states
names(trans_mat) <- states



##################################
#
#
# PART C) Implementing the HMM
#
#
##################################


###################################
# Setting up overlap and trimming value

# Setting country codes
ccode <- 475 
imar <- 0.05
ixs <- 0.1

hmm_raw <- func_preprocess_hmm(df = dat,
                               countrycode = ccode,
                               numstates = length(states),
                               evidence_matrix = ev_mat,
                               transition_matrix = t(trans_mat),
                               marval = imar,
                               xsval = ixs)

hmm <- do.call(rbind, hmm_raw) %>%
  mutate(gid = paste0("ngahex", gid))

hmm_grid <-  left_join(hmm, grid_og, by = "gid") %>%
  left_join(time)
table(hmm_grid$control) 
round(table(hmm_grid$control)/length(hmm_grid$gid)*100, 2)
table(hmm_grid$nganorth)

hmm_grid$control_lab <- factor(hmm_grid$control,
                               levels = c("R",
                                          "DR",
                                          "D",
                                          "DG",
                                          "G"),
                               labels = c("Full rebel control",
                                          "Contested, closer to rebels",
                                          "Contested",
                                          "Contested, closer to government",
                                          "Full government control"))


saveRDS(hmm_grid, here::here("prepped","hmm_nga.rds"))
