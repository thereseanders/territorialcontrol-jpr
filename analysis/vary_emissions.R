############################
# Territorial Control paper
# Sensitivity analysis NGA
#
# Therese Anders
############################

# NOTE: do not run this in RStudio due to forked processing 
# issues when using the furrr package (and future).

setwd("/Users/thereseanders/Dropbox/Dissertation/TerritorialControl/Manuscript/submissions/jpr_march2019/final_apr2020/territorialcontrol-jpr")

library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(HMM)
library(lubridate)
library(VGAM)
library(sf)
library(withr)
library(xtable)
library(furrr)
library(tictoc)
library(purrr)

source("./utils/func_hmm.R")

plotdirapp <- "~/Dropbox/Dissertation/TerritorialControl/Manuscript/submissions/jpr_march2019/final_apr2020/territorialcontrol-jpr/plots/appendix"


##################################
#
#
# PART A) Preprocessing
#
#
##################################


## Loading data and functions
grid_nga <- readRDS("./prepped/grids_hex_nga_hex25.rds") %>%
  filter(nganorth == 1)
gids <- unique(grid_nga$gid)
events_nga <-
  readRDS("./prepped/events_decay_nga_hex25_logistic_trunc.rds")
summary(events_nga$weighted_terrorism)

## Prepping data
crs <- st_crs(grid_nga)

to <- max(events_nga$timeindex)

time <- data.frame(timeindex = seq(1, to),
                   month = rep(seq(1, 12))) %>%
  mutate(year_trans = rep(seq(1, nrow(.) / 12), each = 12)) %>%
  mutate(year = (2008 + year_trans) - 1) %>%
  dplyr::select(-year_trans)

dat_raw <- events_nga %>%
  st_as_sf(coords = c("lonr_c", "latr_c"),
           crs = crs) %>%
  left_join(time, by = "timeindex") %>%
  mutate(gwno = 475) %>%
  dplyr::select(gid,
                gwno,
                year,
                month,
                timeindex,
                gtd = weighted_terrorism_trunc,
                ged = weighted_combats_trunc) %>%
  dplyr::mutate(gtd = replace(gtd, is.na(gtd), 0),
                ged = replace(ged, is.na(ged), 0)) %>%
  st_set_geometry(NULL)


dat_mean_info <- dat_raw %>%
  gather(indicator, value, gtd:ged) %>%
  group_by(gwno, indicator, timeindex) %>%
  dplyr::summarise(mean = mean(value, na.rm = T),
                   var = var(value, na.rm = T)) %>%
  gather(statistic, value, mean:var) %>%
  unite(indc_stat, indicator, statistic) %>%
  spread(indc_stat, value) %>%
  dplyr::select(
    gwno,
    mean_t = gtd_mean,
    mean_c = ged_mean,
    var_t = gtd_var,
    var_c = ged_var,
    timeindex
  ) %>%
  mutate(
    mean_t = replace(mean_t, is.na(mean_t), 0),
    mean_c = replace(mean_c, is.na(mean_c), 0),
    var_t = replace(var_t, is.na(var_t), 0),
    var_c = replace(var_c, is.na(var_c), 0)
  )

dat <- dat_raw  %>%
  left_join(dat_mean_info, by = c("gwno", "timeindex")) %>%
  left_join(grid_nga) %>%
  dplyr::mutate(cdf_T = pzipois(floor(gtd), mean_t),
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
ev_mat[1, ] <- c(0.6, 0.175, 0.175, 0.05)
ev_mat[2, ] <- c(0.05, 0.6, 0.175, 0.175)
ev_mat[3, ] <- c(0.05, 0.175, 0.6, 0.175)
ev_mat[4, ] <- c(0.05, 0.175, 0.175, 0.6)
ev_mat[5, ] <- c(0.6, 0.05, 0.175, 0.175)

###############################
# Setting up initial probabilities
# For Nigeria: Strong prior for government control
d_init <- c(0.025, 0.025, 0.025, 0.025, 0.9)

####### Transition matrix ########
# Modified empirical transition matrix
trans_mat <-
  data.frame(matrix(nrow = length(states), ncol = length(states)))
trans_mat[1, ] <- c(0.25, 0.5, 0.025, 0.2, 0.025)
trans_mat[2, ] <- c(0.25, 0.15, 0.075, 0.5, 0.025)
trans_mat[3, ] <- c(0.05, 0.025, 0.05, 0.85, 0.025)
trans_mat[4, ] <- c(0.025, 0.075, 0.15, 0.125, 0.625)
trans_mat[5, ] <- c(0.05, 0.075, 0.475, 0.025, 0.375)
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

hmm_base_ls <- func_preprocess_hmm(
      df = dat,
      countrycode = ccode,
      numstates = length(states),
      evidence_matrix = ev_mat,
      transition_matrix = t(trans_mat),
      marval = imar,
      xsval = ixs
    )
hmm_base <- do.call(rbind, hmm_base_ls) %>%
  mutate(gid = paste0("ngahex", gid)) %>%
  mutate(control_num = as.numeric(case_when(
    control == "R" ~ 0,
    control == "DR" ~ .25,
    control == "D" ~ .5,
    control == "DG" ~ .75,
    control == "G" ~ 1
  )))
table(hmm_base$control)

## Randomization
n <- 30
xsam <- 500
seed <- 12345
iseeds <- withr::with_seed(seed,
                           sample(100:1e+6, xsam))
probs <- list()
vals <- list()
varv <- c()
i = 1
for (x in 1:length(iseeds)) {
  value <- withr::with_seed(iseeds[x],
                            sample(
                              c(1, 2, 3, 4),
                              size = n,
                              replace = T,
                              prob = c(.6, .175, .175, .05)
                            ))
  val <- as.numeric(table(value) / n)
  if (length(val) == 4) {
    probs[[i]] <- val
    vals[[i]] <- value
    varv[i] <- var(value)
    i = i + 1
  }
}


sd06 <- sd(map_dbl(probs, 1))
sd0175_1 <- sd(map_dbl(probs, 2))
sd0175_2 <- sd(map_dbl(probs, 3))
sd005 <- sd(map_dbl(probs, 4))

df_sd <- data.frame(prob = c(.6, .175, .175, .05),
                    sd = c(sd06, sd0175_1, sd0175_2, sd005))

####################
# Appendix Table
print(xtable(df_sd), include.rownames = FALSE)
####################

## Creating alternative evidence matrices (n = 389)
# func_ev_mat in func_hmm.R
ev_mat_list <- purrr::map(probs,
                          ~ withr::with_seed(seed, func_ev_mat(.)))

# Computing HMMs on each of these evidence matrices
# WARNING: This code takes a while to run! 
# One iteration ~ 34s; 389 iterations ~ >1.5h
# supportsMulticore()
# tic()
# plan(multiprocess)
# hmm_list <- furrr::future_map(
#   ev_mat_list,
#   ~ func_preprocess_hmm(
#     df = dat,
#     countrycode = ccode,
#     numstates = length(states),
#     evidence_matrix = .x,
#     transition_matrix = t(trans_mat),
#     marval = imar,
#     xsval = ixs
#   ), .progress = TRUE
# )
# toc()
# saveRDS(hmm_list, "./prepped/hmm_list_nga_ta_x500_n30.rds")


hmm_list <- readRDS("./prepped/hmm_list_nga_ta_x500_n30.rds")


# testing base hmm vs distribution
test_dist <- append(hmm_list, list(hmm_base_ls)) %>%
  furrr::future_map(., ~bind_rows(.)) %>%
  purrr::map(~ dplyr::mutate(., gid = paste0("ngahex", gid))) %>%
  purrr::map(~ mutate(., control_num = as.numeric(case_when(
    control == "R" ~ 0,
    control == "DR" ~ .25,
    control == "D" ~ .5,
    control == "DG" ~ .75,
    control == "G" ~ 1
  ))))


test_constellation <- purrr::map(test_dist,
                                 ~as.numeric(table(.[["control"]])))

# below removing "base" HMM
test_summary <- test_dist[-length(test_dist)] %>%
  map_dfc(., dplyr::select, c("control_num")) %>%
  as_tibble() %>%
  mutate(sd = apply(dplyr::select(., contains("control")), 1, sd, na.rm=TRUE)) %>%
  mutate(av = apply(dplyr::select(., contains("control")), 1, mean, na.rm=TRUE)) %>%
  mutate(gid = test_dist[[1]]$gid,
         timeindex = test_dist[[1]]$timeindex) %>%
  dplyr::select(gid, timeindex, sd, av) %>%
  left_join(time)

comp <- left_join(test_summary, hmm_base) %>%
  mutate(cilo = av - sd,
         cihi = av + sd) %>%
  rowwise() %>%
  mutate(test = ifelse(between(control_num, cilo, cihi), 1, 0))
table(comp$test)/nrow(comp) #

## Using Rubin's rules
# https://bookdown.org/mwheymans/bookmi/rubins-rules.html
test_rubin <- test_dist[-length(test_dist)] %>%
  map_dfc(., dplyr::select, c("control_num")) %>%
  as_tibble() %>%
  mutate(gid = test_dist[[1]]$gid,
         timeindex = test_dist[[1]]$timeindex) %>%
  left_join(time) %>%
  dplyr::select(gid, year, month, everything(), -timeindex) %>%

  # assumption: parameter estimate is annual mean of control
  pivot_longer(cols = -c(gid, year, month), names_to = "sample", values_to = "value") %>%
  group_by(gid, year, sample) %>%
  summarise(yr_mean = mean(value),
            yr_sd = sd(value),
            yr_nobs = n()) %>%
  mutate(yr_se = yr_sd/sqrt(yr_nobs)) %>%
  ungroup()


m <- length(unique(test_rubin$sample))
gid_nganorth <- unique(grid_nga$gid[grid_nga$nganorth == 1])


test_rubin_pooledmean <- test_rubin %>%
  group_by(year, gid) %>%
  summarise(v_mean = mean(yr_mean),
            v_w = (1/m)*(sum(yr_se^2))) %>%
  ungroup()

test_rubin_res <- left_join(test_rubin, test_rubin_pooledmean) %>%
  mutate(res = (v_mean - yr_mean)^2) %>%
  group_by(gid, year) %>%
  summarise(v_b = sum(res)/(m-1)) %>%
  ungroup()

test_rubin_summary <- left_join(test_rubin_pooledmean, test_rubin_res) %>%
  mutate(v_total = v_w + v_b + v_b/m) %>%
  mutate(sd_total = sqrt(v_total)) %>%
  mutate(gid_new = dplyr::group_indices(., gid))

mycols <- rev(c("#2c7bb6", #Blue
                "#abd9e9",
                "#ffffbf", # Yellow
                "#fdae61",
                "#d7191c"))

# Plotting it
p1 <- ggplot(subset(test_rubin_summary, year %in% seq(2008, 2012)),
       aes(x = gid_new,
           y = v_mean,
           color = v_mean)) +
  geom_hline(yintercept = 0, color = "grey", alpha = 0.5, size = 0.2) +
  geom_hline(yintercept = 0.25, color = "grey", alpha = 0.5, size = 0.2) +
  geom_hline(yintercept = 0.5, color = "grey", alpha = 0.5, size = 0.2) +
  geom_hline(yintercept = 0.75, color = "grey", alpha = 0.5, size = 0.2) +
  geom_hline(yintercept = 1, color = "grey", alpha = 0.5, size = 0.2) +
  geom_point(size = 0.1) +
  geom_linerange(aes(ymin = v_mean - 1.95*sd_total, ymax = v_mean + 1.95*sd_total),
                 size = 0.1) +
  facet_wrap(~year, ncol = 1) +
  scale_color_gradientn(name = "",
                        colors = mycols,
                        breaks = c(0, 0.5, 1),
                        limits = c(0, 1),
                        labels = c("Full rebel\ncontrol",
                                   "Highly\ncontested",
                                   "Full government\ncontrol")) +
  theme_light() +
  theme(legend.position = "top",
        panel.grid = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  labs(y = "Territorial control",
       x = "Grid cell ID")  +
  scale_y_continuous(breaks = c(0, 0.25, .5, .75, 1))
ggsave(paste0(plotdirapp, "/", "plot1_evvariance_rubin.png"),
       width = 12, height = 8, dpi = 400, p1)


p2 <- ggplot(subset(test_rubin_summary, year %in% seq(2013, 2017)),
       aes(x = gid_new,
           y = v_mean,
           color = v_mean)) +
  geom_hline(yintercept = 0, color = "grey", alpha = 0.5, size = 0.2) +
  geom_hline(yintercept = 0.25, color = "grey", alpha = 0.5, size = 0.2) +
  geom_hline(yintercept = 0.5, color = "grey", alpha = 0.5, size = 0.2) +
  geom_hline(yintercept = 0.75, color = "grey", alpha = 0.5, size = 0.2) +
  geom_hline(yintercept = 1, color = "grey", alpha = 0.5, size = 0.2) +
  geom_point(size = 0.1) +
  geom_linerange(aes(ymin = v_mean - 1.95*sd_total, ymax = v_mean + 1.95*sd_total),
                 size = 0.1) +
  facet_wrap(~year, ncol = 1) +
  scale_color_gradientn(name = "",
                        colors = mycols,
                        breaks = c(0, 0.5, 1),
                        limits=c(0,1),
                        labels = c("Full rebel\ncontrol",
                                   "Highly\ncontested",
                                   "Full government\ncontrol")) +
  theme_light() +
  theme(legend.position = "top",
        panel.grid = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  labs(y = "Territorial control",
       x = "Grid cell ID")  +
  scale_y_continuous(breaks = c(0, 0.25, .5, .75, 1))
ggsave(paste0(plotdirapp, "/", "plot2_evvariance_rubin.png"),
       width = 12, height = 8, dpi = 400, p2)


p_scatter <- ggplot(test_rubin_summary,
       aes(x = v_mean,
           y = sd_total,
           fill = v_mean)) +
  geom_smooth(color = "black", alpha = 0.6) +
  geom_point(size = 1,
             shape = 21,
             color = "black") +
  scale_fill_gradientn(name = "",
                       colors = mycols,
                       breaks = c(0, 0.5, 1),
                       limits=c(0,1),
                       labels = c("Full rebel\ncontrol",
                                  "Highly\ncontested",
                                  "Full government\ncontrol")) +
  theme_light() +
  theme(legend.position = "top",
        legend.key.width = unit(1.2,"cm"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  labs(x = "Pooled estimate of territorial control",
       y = "Total standard deviation")
ggsave(paste0(plotdirapp, "/", "plot2_est_sd.png"),
       width = 8, height = 5, dpi = 400, p_scatter)
