############################
# Territorial Control paper
# Graphs and figures appendix
#
# Therese Anders
############################

library(raster)
library(tidyverse)
library(sf)
library(lubridate)
library(rgeos)
library(xtable)
library(scales)
library(ggpubr)
library(ggalluvial)
library(VGAM)
library(haven)

source(here::here("utils","func_weights.R"))
source(here::here("utils","func_hmm.R"))

plotdirapp <- here::here("plots", "appendix")

####################################################
# Loading data
####################################################
nigeria_districts <- raster::getData("GADM", country = "NGA", level = 2)
events_nga <- readRDS(here::here("prepped","events_nga_hex25.rds"))
events_col <- readRDS(here::here("prepped","events_col_hex25.rds"))
events_nga_cont <- readRDS(here::here("prepped","events_decay_nga_hex25_logistic_trunc.rds"))
events_col_cont <- readRDS(here::here("prepped","events_decay_col_hex25_logistic_trunc.rds"))
hmm_col <- readRDS(here::here("prepped","hmm_col.rds"))
hmm_nga <- readRDS(here::here("prepped","hmm_nga.rds"))
grid_col <- readRDS(here::here("prepped","grids_hex_col_hex25.rds"))
grid_nga <- readRDS(here::here("prepped","grids_hex_nga_hex25.rds"))
df_priogrid <- readRDS(here::here("prepped","priogrid_pop_sub.rds"))


mycols <- rev(c("#2c7bb6", #Blue
                "#abd9e9",
                "#ffffbf", # Yellow
                "#fdae61",
                "#d7191c"))
names(mycols) <- levels(hmm_col$control_lab)

# For validation and robustness checks: ACLED territorial control estimates
files <- list.files(here::here("prepped","acled"))
acled_ls <- list()
for (f in 1:length(files)) {
  name <- files[f] %>%
    str_remove(".rds") %>%
    str_remove("acledtesting_nganorth_hex25_control_flip_")
  acled_ls[[f]] <- readRDS(paste0(here::here("prepped","acled"), "/", files[f])) %>%
    dplyr::rename(control_test = control) %>%
    mutate(type = name) %>%
    separate(type, c("temporal", "duration", "type"))
   names(acled_ls)[[f]] <- name
}
acled <- bind_rows(acled_ls)

#######################################################
# Appendix Figure 1: Temporal distribution of events
#######################################################
events_merged <- bind_rows(events_nga, events_col) %>%
  mutate(select = case_when(
    str_detect(gid, "nga") & nganorth == 1 ~ 1,
    str_detect(gid, "col") & colwest == 1 ~ 1,
    T ~ 0
    )) %>%
  filter(select == 1) %>%
  mutate(country_label = ifelse(str_detect(gid, "nga"), 
                          "Boko Haram associated events in NE Nigeria", 
                          "FARC associated events in Colombia (without Orinoco and Amazon)"))

events_overtime <- events_merged %>%
  group_by(year, month, type, country_label) %>%
  summarize(no_events = n()) %>%
  ungroup() %>%
  mutate(ym = ymd(paste0(year, "-", month), truncated = 1))

min <- min(events_overtime$ym)
max <- max(events_overtime$ym)

df_overtime <- data.frame(expand.grid(year = seq(2009, 2017),
                                      month = seq(1,12))) %>%
  mutate(ym = ymd(paste0(year, "-", month), truncated = 1)) %>%
  filter(ym >= min,
         ym <= max) %>%
  arrange(ym) %>%
  left_join(events_overtime) %>%
  mutate(no_events = ifelse(is.na(no_events), 0, no_events)) %>%
  filter(!is.na(type))

p_overtime <- ggplot(events_overtime,
                     aes(x = ym,
                         y = no_events,
                         color = type)) + 
  geom_line(alpha = 0.9) +
  facet_wrap(~country_label, ncol = 1) +
  theme_light() +
  labs(title = "Monthly event count",
       x = "",
       y = "Number of events") +
  scale_color_manual(labels = c("Conventional (GED)",
                                "Terrorism (GTD)"),
                     name = "",
                     values = c("#d7191c", #red
                                "#2c7bb6")) +
  theme(legend.position = "top")
ggsave(paste0(plotdirapp, "/", "plot_overtime_numevents.png"), 
       width = 7, height = 6, dpi = 400, p_overtime)  

#######################################################
# Appendix Figure 2: Spatial distribution of events (NGA)
#######################################################
events_nganorth_sub <- events_nga %>%
  filter(nganorth == 1) %>%
  dplyr::select(gid, type, month, year, eventid, fatalities) %>%
  st_set_geometry(NULL) %>%
  group_by(gid, year, type) %>%
  summarize(no_events = n()) %>%
  ungroup()

grid_nganorth <- grid_nga %>%
  filter(nganorth == 1)

gid_nganorth <- data.frame(expand.grid(gid = unique(grid_nganorth$gid),
                                       type = c("conventional", "terrorism"),
                                       year = seq(min(events_nganorth_sub$year), 
                                                  max(events_nganorth_sub$year)))) %>%
  left_join(events_nganorth_sub) %>%
  mutate(no_events = replace(no_events, is.na(no_events), 0)) 

events_nganorth_discrete <- gid_nganorth %>%
  group_by(gid, type) %>%
  summarise(av = mean(no_events)) %>%
  mutate(type_label = ifelse(type == "conventional",
                             "Conventional (GED)",
                             "Terrorism (GTD)")) %>%
  ungroup() %>%
  mutate(aggregation = "Discrete") %>%
  left_join(grid_nganorth) %>%
  st_as_sf()

df_time <- data.frame(expand.grid(year = seq(2008, 2017),
                                  month = seq(1,12))) %>%
  arrange(year, month) %>%
  mutate(timeindex = seq(1:nrow(.)))

events_nganorth_contin <- events_nga_cont %>%
  filter(nganorth == 1) %>%
  left_join(df_time) %>%
  filter(year != 0) %>%
  dplyr::select(gid, timeindex, weighted_terrorism_trunc, weighted_combats_trunc) %>%
  pivot_longer(c(weighted_terrorism_trunc, weighted_combats_trunc)) %>%
  group_by(gid, name) %>%
  summarise(av = mean(value)) %>%
  mutate(type_label = ifelse(name == "weighted_combats_trunc",
                             "Conventional (GED)",
                             "Terrorism (GTD)")) %>%
  ungroup() %>%
  mutate(aggregation = "Continuous") %>%
  left_join(grid_nga) %>%
  st_as_sf()

p_nganorth_contin <- ggplot(events_nganorth_contin) +
  geom_sf(aes(fill = av),
          color = "grey",
          size = 0.1) +
  theme_light() +
  scale_fill_gradient(low = "white",
                      high = "black",
                      trans = "log10",
                      na.value = "transparent",
                      name = "Average annual events",
                      limits = c(0.1, 30)) +
  facet_wrap(~type_label) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  labs(title = "Weighted assignment to grid cells",
       subtitle = "Boko Haram related events in NE Nigeria 2009-2017")
ggsave(paste0(plotdirapp, "/", "plot_nganorth_continuous.png"), 
       width = 7, height = 3, dpi = 400, p_nganorth_contin) 


p_nganorth_discrete <- ggplot(events_nganorth_discrete) +
  geom_sf(aes(fill = av),
          color = "grey",
          size = 0.1) +
  theme_light() +
  scale_fill_gradient(low = "white",
                      high = "black",
                      trans = "log10",
                      na.value = "transparent",
                      name = "Average annual events") +
  facet_wrap(~type_label) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  labs(title = "Discrete assignment to grid cells",
       subtitle = "Boko Haram related events in NE Nigeria 2009-2017")
ggsave(paste0(plotdirapp, "/", "plot_nganorth_discrete.png"), 
       width = 7, height = 3, dpi = 400, p_nganorth_discrete) 



#######################################################
# Appendix Figure 3: Spatial distribution of events (COL)
#######################################################
events_colwest_sub <- events_col %>%
  filter(colwest == 1) %>%
  dplyr::select(gid, type, month, year, eventid, fatalities) %>%
  st_set_geometry(NULL) %>%
  group_by(gid, year, type) %>%
  summarize(no_events = n()) %>%
  ungroup()

gid_colwest <- data.frame(expand.grid(gid = unique(grid_col$gid[grid_col$colwest == 1]),
                                      type = c("conventional", "terrorism"),
                                      year = seq(min(events_colwest_sub$year), 
                                                 max(events_colwest_sub$year)))) %>%
  left_join(events_colwest_sub) %>%
  mutate(no_events = replace(no_events, is.na(no_events), 0)) 


events_colwest_discrete <- gid_colwest %>%
  group_by(gid, type) %>%
  summarise(av = mean(no_events)) %>%
  mutate(type_label = ifelse(type == "conventional",
                             "Conventional (GED)",
                             "Terrorism (GTD)")) %>%
  ungroup() %>%
  mutate(aggregation = "Discrete") %>%
  left_join(grid_col) %>%
  filter(colwest == 1) %>%
  st_as_sf()


events_colwest_contin <- events_col_cont %>%
  dplyr::select(gid, timeindex, weighted_terrorism_trunc, weighted_combats_trunc) %>%
  pivot_longer(c(weighted_terrorism_trunc, weighted_combats_trunc)) %>%
  group_by(gid, name) %>%
  summarise(av = mean(value)) %>%
  mutate(type_label = ifelse(name == "weighted_combats_trunc",
                             "Conventional (GED)",
                             "Terrorism (GTD)")) %>%
  ungroup() %>%
  mutate(aggregation = "Continuous") %>%
  left_join(grid_col) %>%
  filter(colwest == 1) %>%
  st_as_sf()

p_colwest_contin <- ggplot(events_colwest_contin) +
  geom_sf(aes(fill = av),
          color = "grey",
          size = 0.1) +
  theme_light() +
  scale_fill_gradient(low = "white",
                      high = "black",
                      #trans = "log10",
                      na.value = "transparent",
                      name = "Average annual events") +
  facet_wrap(~type_label) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  labs(title = "Weighted assignment to grid cells",
       subtitle = "FARC related events in Colombia (except Orinoco and Amazon) 1997-2017") +
  coord_sf(xlim = c(-80, -70),
           ylim = c(0, 12.75))
ggsave(paste0(plotdirapp, "/", "plot_colwest_continuous.png"), 
       width = 7, height = 3, dpi = 400, p_colwest_contin) 


p_colwest_discrete <- ggplot(events_colwest_discrete) +
  geom_sf(aes(fill = av),
          color = "grey",
          size = 0.1) +
  theme_light() +
  scale_fill_gradient(low = "white",
                      high = "black",
                      #trans = "log10",
                      na.value = "transparent",
                      name = "Average annual events") +
  facet_wrap(~type_label) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  labs(title = "Discrete assignment to grid cells",
       subtitle = "FARC related events in Colombia (except Orinoco and Amazon) 1997-2017") +
  coord_sf(xlim = c(-80, -70),
           ylim = c(0, 12.75))
ggsave(paste0(plotdirapp, "/","plot_colwest_discrete.png"), 
       width = 7, height = 3, dpi = 400, p_colwest_discrete) 





#######################################################
# Appendix Figure 5: decay space & time
#######################################################
df_s <- data.frame(d = seq(1,40,1)) %>%
  mutate(w = func_weight(d, 7, -0.35),
         type = "Spatial distance (in km)")

df_t <- data.frame(d = seq(1, 12, 0.1)) %>%
  mutate(w = func_weight(d, 8, -2.5),
         type = "Temporal distance (in months)")

p_s <- ggplot(df_s,
       aes(x = d, y = w)) +
  geom_line() +
  labs(title = "Logistic decay for spatial dimension",
       x = "Distance in km",
       y = "Weight") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
ggsave(paste0(plotdirapp, "/", "decay_space.png"), 
       width = 4, height = 3, dpi = 400, p_s)

p_t <- ggplot(df_t,
       aes(x = d, y = w)) +
  geom_line() +
  labs(title = "Logistic decay for temporal dimension",
       x = "Distance in months",
       y = "Weight") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
ggsave(paste0(plotdirapp, "/",  "decay_time.png"), 
       width = 4, height = 3, dpi = 400, p_t)


#######################################################
# Appendix Figure 7: Plotting transition probabilities
#######################################################
states <- c("R", "DR", "D", "DG", "G")
trans_mat <- data.frame(matrix(nrow = length(states), ncol = length(states)))
trans_mat[1,] <- c(0.25, 0.5, 0.025, 0.2, 0.025)
trans_mat[2,] <- c(0.25, 0.15, 0.075, 0.5, 0.025)
trans_mat[3,] <- c(0.05, 0.025, 0.05, 0.85, 0.025)
trans_mat[4,] <- c(0.025, 0.075, 0.15, 0.125, 0.625)
trans_mat[5,] <-c(0.05, 0.075, 0.475, 0.025, 0.375)
row.names(trans_mat) <- states
names(trans_mat) <- states

trans_df <- as.data.frame(trans_mat) %>%
  mutate(q_tm1 = states %>% factor(levels = states)) %>%
  pivot_longer(-q_tm1, names_to = "temp") %>%
  mutate(q_t = temp %>% factor(levels = states)) 

p_trans <- ggplot(trans_df,
       aes(axis1 = q_tm1, axis2 = q_t, y = value)) +
  geom_flow(aes(fill = value), alpha = 0.9, color = "black") +
  geom_stratum(fill = "#9d827d", alpha = 0.5) +
  geom_label(stat = "stratum", infer.label = TRUE) +
  facet_wrap(~q_tm1, nrow = 1) +
  scale_fill_gradient(low = "white", high = "black",
                      name = "Transition probability") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top") +
  labs(title = "Transition probabilities based on Kalyvas (2006)",
       x = "",
       y = "") +
  scale_x_discrete(limits = c("q_{t-1}", "q_t"))
ggsave(paste0(plotdirapp, "/", "plot_trans.png"), 
       width = 7, height = 5.5, dpi = 400, p_trans)

#######################################################
# Appendix Figure 7: Plotting emission probabilities
#######################################################

ev <- c("O1", "O2", "O3", "O4")
ev_mat <- matrix(nrow = length(states), ncol = length(ev))
ev_mat[1,] <- c(0.6, 0.175, 0.175, 0.05)
ev_mat[2,] <- c(0.05, 0.6, 0.175, 0.175)
ev_mat[3,] <- c(0.05, 0.175, 0.6, 0.175)
ev_mat[4,] <- c(0.05, 0.175, 0.175, 0.6)
ev_mat[5,] <- c(0.6, 0.05, 0.175, 0.175)

ev_df <- as.data.frame(ev_mat) %>%
  mutate(q_t = states %>% factor(levels = states)) %>%
  pivot_longer(-q_t, names_to = "o_t") %>%
  mutate(o_t = str_replace_all(o_t, "V", "O") %>% factor(levels = ev))

p_ev <- ggplot(ev_df,
       aes(axis1 = q_t, axis2 = o_t, y = value)) +
  geom_flow(aes(fill = value), alpha = 0.9, color = "black") +
  geom_stratum(fill = "#9d827d", alpha = 0.5) +
  geom_label(stat = "stratum", infer.label = TRUE) +
  facet_wrap(~q_t, nrow = 1) +
  scale_fill_gradient(low = "white", high = "black",
                      name = "Emission probability") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top") +
  labs(title = "Emission probabilities",
       x = "",
       y = "") +
  scale_x_discrete(limits = c("q_t", "o_t"))
ggsave(paste0(plotdirapp, "/", "plot_ev.png"), 
       width = 7, height = 5.5, dpi = 400, p_ev)



#######################################################
# Appendix Figure 9 a) & b): ACLED validation data (6m)
#######################################################
acled_yearly <- acled %>%
  filter(temporal == "yearly") %>%
  dplyr::select(gid, year, month, control_test, type, duration) %>%
  left_join(grid_nga) %>%
  st_as_sf()

p_acled_6m_full <- ggplot(acled_yearly %>%
                            filter(type == "full", duration == "6m")) +
  geom_sf(aes(fill = control_test), size = 0.01, color = "lightgrey") +
  facet_wrap(~year, nrow = 2) +
  theme_bw() +
  scale_fill_gradientn(name = "",
                       colors = mycols,
                       breaks = c(0, 0.5, 1),
                       limits = c(0,1),
                       labels = c("Full rebel\ncontrol",
                                  "Highly\ncontested",
                                  "Full government\ncontrol")) +
  theme(legend.position = "bottom",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(5,5,2,2),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.box.background = element_blank(),
        legend.key.width = unit(2, "cm")) +
  labs(title = "") +
  coord_sf(datum = NA)
ggsave(paste0(plotdirapp, "/", "acled_6m_full.png"), 
       width = 7, height = 4, dpi = 400, p_acled_6m_full)


p_acled_6m_select <- ggplot(acled_yearly %>%
                            filter(type == "select", duration == "6m")) +
  geom_sf(aes(fill = control_test), size = 0.01, color = "lightgrey") +
  facet_wrap(~year, nrow = 2) +
  theme_bw() +
  scale_fill_gradientn(name = "",
                       colors = mycols,
                       breaks = c(0, 0.5, 1),
                       limits = c(0,1),
                       labels = c("Full rebel\ncontrol",
                                  "Highly\ncontested",
                                  "Full government\ncontrol")) +
  theme(legend.position = "bottom",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(5,5,2,2),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.box.background = element_blank(),
        legend.key.width = unit(2, "cm")) +
  labs(title = "") +
  coord_sf(datum = NA)
ggsave(paste0(plotdirapp, "/", "acled_6m_select.png"),
       width = 7, height = 4, dpi = 400, p_acled_6m_select)


#######################################################
# Appendix Figure 10 a) & b): Acled validation data (12m)
#######################################################


p_acled_12m_full <- ggplot(acled_yearly %>%
                            filter(type == "full", duration == "12m")) +
  geom_sf(aes(fill = control_test), size = 0.01, color = "lightgrey") +
  facet_wrap(~year, nrow = 2) +
  theme_bw() +
  scale_fill_gradientn(name = "",
                       colors = mycols,
                       breaks = c(0, 0.5, 1),
                       limits = c(0,1),
                       labels = c("Full rebel\ncontrol",
                                  "Highly\ncontested",
                                  "Full government\ncontrol")) +
  theme(legend.position = "bottom",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(5,5,2,2),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.box.background = element_blank(),
        legend.key.width = unit(2, "cm")) +
  labs(title = "") +
  coord_sf(datum = NA)
ggsave(paste0(plotdirapp, "/", "acled_12m_full.png"), 
       width = 7, height = 4, dpi = 400, p_acled_12m_full)


p_acled_12m_select <- ggplot(acled_yearly %>%
                              filter(type == "select", duration == "12m")) +
  geom_sf(aes(fill = control_test), size = 0.01, color = "lightgrey") +
  facet_wrap(~year, nrow = 2) +
  theme_bw() +
  scale_fill_gradientn(name = "",
                       colors = mycols,
                       breaks = c(0, 0.5, 1),
                       limits = c(0,1),
                       labels = c("Full rebel\ncontrol",
                                  "Highly\ncontested",
                                  "Full government\ncontrol")) +
  theme(legend.position = "bottom",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(5,5,2,2),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.box.background = element_blank(),
        legend.key.width = unit(2, "cm")) +
  labs(title = "") +
  coord_sf(datum = NA)
ggsave(paste0(plotdirapp, "/", "acled_12m_select.png"), 
       width = 7, height = 4, dpi = 400, p_acled_12m_select)

####################################################
# Appendix Figure 11: Correlation w ACLED validation
####################################################
hmm_nga_monthly <- hmm_nga %>%
  mutate(control_num = case_when(
    control == "R" ~ 0,
    control == "DR" ~ 0.25,
    control == "D" ~ 0.5,
    control == "DG" ~ 0.75,
    control == "G" ~ 1
  ))

# Computing correlations
merged_monthly <- acled %>%
  filter(temporal == "monthly") %>%
  dplyr::select(gid, year, month, control_test, type, duration) %>%
  left_join(hmm_nga_monthly)

# ## Monthly correlation?
ti <- unique(merged_monthly$timeindex)
cors_full_12m <- c()
cors_full_6m <- c()
cors_full_3m <- c()
cors_select_12m <- c()
cors_select_6m <- c()
cors_select_3m <- c()
for(y in 1:length(ti)){
  sub <- merged_monthly %>%
    filter(timeindex == ti[y])
  cors_full_12m[y] <- tryCatch({cor.test(sub$control_test[sub$type == "full" & sub$duration == "12m"], sub$control_num[sub$type == "full" & sub$duration == "12m"], 
                                         method = "spearman")$estimate},
                               error = function(e) {print(paste("Error on timeindex ", ti[y]))})
  cors_full_6m[y] <- tryCatch({cor.test(sub$control_test[sub$type == "full" & sub$duration == "6m"], sub$control_num[sub$type == "full" & sub$duration == "6m"], 
                                        method = "spearman")$estimate},
                              error = function(e) {print(paste("Error on timeindex ", ti[y]))})
  cors_full_3m[y] <- tryCatch({cor.test(sub$control_test[sub$type == "full" & sub$duration == "3m"], sub$control_num[sub$type == "full" & sub$duration == "3m"], 
                                        method = "spearman")$estimate},
                              error = function(e) {print(paste("Error on timeindex ", ti[y]))})
  cors_select_12m[y] <- tryCatch({cor.test(sub$control_test[sub$type == "select" & sub$duration == "12m"], sub$control_num[sub$type == "select" & sub$duration == "12m"], 
                                           method = "spearman")$estimate},
                                 error = function(e) {print(paste("Error on timeindex ", ti[y]))})
  cors_select_6m[y] <- tryCatch({cor.test(sub$control_test[sub$type == "select" & sub$duration == "6m"], sub$control_num[sub$type == "select" & sub$duration == "6m"], 
                                          method = "spearman")$estimate},
                                error = function(e) {print(paste("Error on timeindex ", ti[y]))})
  cors_select_3m[y] <- tryCatch({cor.test(sub$control_test[sub$type == "select" & sub$duration == "3m"], sub$control_num[sub$type == "select" & sub$duration == "3m"], 
                                          method = "spearman")$estimate},
                                error = function(e) {print(paste("Error on timeindex ", ti[y]))})
  
}

df_time <- merged_monthly %>%
  dplyr::select(timeindex, year, month) %>%
  unique()

df_monthly <- data.frame(timeindex = ti,
                         full_12m = cors_full_12m,
                         full_6m = cors_full_6m,
                         full_3m = cors_full_3m,
                         select_12m = cors_select_12m,
                         select_6m = cors_select_6m,
                         select_3m = cors_select_3m) %>%
  pivot_longer(-timeindex) %>%
  left_join(df_time) %>%
  mutate(timelabel = paste(month, year, sep = "-")) %>%
  separate(name, c("type", "duration")) %>%
  mutate(type_label = case_when(
    type == "full" ~ "Full",
    type == "select" ~ "Restricted"
  )) %>%
  mutate(duration_label = case_when(
    duration == "3m" ~ "3 month lag imputation",
    duration == "6m" ~ "6 month lag imputation",
    duration == "12m" ~ "12 month lag imputation"
  )) %>%
  filter(year >= 2010)

df_monthly$duration_label <- factor(df_monthly$duration_label,
                                    levels = c("3 month lag imputation",
                                               "6 month lag imputation",
                                               "12 month lag imputation"))


x_l <- unique(df_monthly$timelabel[df_monthly$month == 1])
x_b <- unique(df_monthly$timeindex[df_monthly$month == 1])


p_acled_appendix <- ggplot(df_monthly,
                           aes(x = timeindex,
                               y = value,
                               linetype = type_label,
                               color = type_label)) +
  facet_wrap(~duration_label, nrow = 1) +
  geom_line(size = 0.25) +
  labs(title = "Monthly correlation between HMM estimates and ACLED testing data",
       subtitle = "2010 - 2017",
       x = "",
       y = "Spearman correlation coefficient") +
  scale_x_continuous(labels = x_l,
                     breaks = x_b) +
  scale_color_manual(values = c("black", "darkgrey"),
                     name = "ACLED sample") +
  scale_linetype_discrete(name = "ACLED sample") +
  theme_light() +
  geom_smooth(size = 0.5) +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 6))
ggsave(paste0(plotdirapp, "/", "cor_acled_appendix.png"), 
       width = 8, height = 4, dpi = 400, p_acled_appendix)


####################################################
# Appendix Table 6: HMM vs ACLED
####################################################
summary_dens <- merged_monthly %>%
  filter(type == "full",
         duration %in% c("6m", "12m")) %>%
  dplyr::select(-geometry) %>%
  mutate(control_test_bin = case_when(
    control_test < 0.125 ~ 0,
    control_test >= 0.125 & control_test < 0.375 ~ 0.25,
    control_test >= 0.375 & control_test < 0.625 ~ 0.5,
    control_test >= 0.625 & control_test < 0.875 ~ 0.75,
    control_test >= 0.875 ~ 1
  )) %>%
  dplyr::select(-control_test) %>%
  pivot_wider(names_from = duration, values_from = control_test_bin) %>%
  dplyr::select(gid,
                HMM = control_num,
                ACLED_6m = `6m`,
                ACLED_12m = `12m`) %>%
  pivot_longer(-gid) %>%
  group_by(name, value) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  group_by(name) %>%
  mutate(perc = round(100/sum(count)*count, 2)) %>%
  dplyr::select(name, value, perc) %>%
  pivot_wider(names_from = name, values_from = perc) %>%
  dplyr::select(`Territorial control` = value,
                HMM,
                `ACLED re-binned (6m)` = ACLED_6m,
                `ACLED re-binned (12m)` = ACLED_12m)


print(xtable(summary_dens), include.rownames = FALSE)

####################################################
# Appendix Figure 12: HMM vs ACLED
####################################################

dens <- merged_monthly %>%
  filter(type == "full",
         duration == "6m") %>%
  dplyr::select(gid,
                HMM = control_num,
                ACLED = control_test) %>%
  mutate(ACLED_bin = case_when(
    ACLED < 0.125 ~ 0,
    ACLED >= 0.125 & ACLED < 0.375 ~ 0.25,
    ACLED >= 0.375 & ACLED < 0.625 ~ 0.5,
    ACLED >= 0.625 & ACLED < 0.875 ~ 0.75,
    ACLED >= 0.875 ~ 1
  )) %>%
  pivot_longer(-gid) %>%
  group_by(name, value) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(count_log = log10(count))

dens$name <- factor(dens$name,
                    levels = c("HMM", "ACLED", "ACLED_bin"),
                    labels = c("HMM", "ACLED", "ACLED re-binned"))

p_dens <- ggplot(dens,
                 aes(x = value)) +
  geom_line(stat = "density") +
  facet_wrap(~name) +
  theme_light() +
  labs(y = "Density",
       x = "Territorial control")


scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}

p_hist <- ggplot(dens,
                 aes(x = value, y = count, fill = value)) +
  geom_bar(stat = "identity", width = 0.05, color = "grey") +
  facet_wrap(~name) +
  scale_fill_gradientn(name = "",
                       colors = mycols,
                       breaks = c(0, 1),
                       labels = c("Full rebel control",
                                  "Full government control")) +
  theme_light() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.key.width = unit(1.25, "cm"))  +
  labs(x = "Territorial control",
       y = "Observations per value of control (log base 10)") +
  scale_y_log10(label = scientific)

p_all <- ggarrange(p_dens, p_hist, nrow = 2, heights = c(1, 2))
ggsave(paste0(plotdirapp, "/", "hist.png"), 
       width = 8, height = 8, dpi = 400, p_all)


####################################################
# Appendix Figure 13: HMM vs Reuters
####################################################
### Coding districts in Yobe, Borno, Adamawa states 
# Information based on:
# http://blogs.reuters.com/data-dive/2015/05/05/mapping-boko-harams-decline-in-nigeria/
# https://www.cfr.org/blog/has-tide-turned-against-boko-haram-nigeria

districts <- sort(c("Geidam", "Magumeri", "Kaga", "Jere", "Maidugur", "Damboa", "Chibok", "Biu", "Kwaya Kusar", "Bayo", "Shani",
                    "Kukawa", "Monguno", "Marte", "Konduga", "Michika",
                    "Abadam", "Mobbar", "Guzamala", "Gubio", "Nganzai", "Mafa", "Ngala", "Kala/Balge", "Dikwa", "Bama", "Gwoza", "Madagali", "Askira/U", "Gujba", "Gulani", "Hawul"))


## February 25
contested0225 <- c("Kukawa", "Monguno", "Marte", "Konduga", "Michika")
boko0225 <- c("Abadam", "Mobbar", "Guzamala", "Gubio", "Nganzai", "Mafa", "Ngala", "Kala/Balge", "Dikwa", "Bama", "Gwoza", "Madagali", "Askira/U", "Gujba", "Gulani")

## Mar 10
contested0310 <- c("Mobbar", "Abadam", "Kukawa", "Monguno", "Marte", "Dikwa", "Mafa", "Konduga")
boko0310 <- c("Guzamala", "Gubio", "Nganzai", "Ngala", "Kala/Balge", "Bama", "Gwoza", "Askira/U","Gulani")

## Mar 18
boko0318 <- c("Abadam", "Kala/Balge", "Gwoza")

## April 24 (Sambisa forest)
boko0424 <- c("Gwoza")

df_control <- tibble(district = districts) %>%
  mutate(control0225 = case_when(
    district %in% contested0225 ~ "Contested",
    district %in% boko0225 ~ "Boko Haram",
    T ~ "Government"
  )) %>%
  mutate(control0310 = case_when(
    district %in% contested0310 ~ "Contested",
    district %in% boko0310 ~ "Boko Haram",
    T ~ "Government"
  )) %>%
  mutate(control0318 = case_when(
    district %in% boko0318 ~ "Boko Haram",
    T ~ "Government"
  )) %>%
  mutate(control0424 = case_when(
    district %in% boko0424 ~ "Boko Haram",
    T ~ "Government"
  ) %>% factor(levels = c("Boko Haram", "Contested", "Government")))


all <- districts

select_sub <- nigeria_districts[nigeria_districts$NAME_2 %in% all,] %>%
  st_as_sf()
st_crs(select_sub)

select_sf <- select_sub %>%
  left_join(df_control, by = c("NAME_2" = "district")) %>%
  gather(control, actor, contains("control")) %>%
  mutate(window = case_when(
    control == "control0225" ~ "25 February",
    control == "control0310" ~ "10 March",
    control == "control0318" ~ "18 March",
    control == "control0424" ~ "24 April"
  ) %>% factor(levels = c("25 February", "10 March", "18 March", "24 April")))


grid_nganorthsmall <- grid_nga %>%
  filter(nganorthsmall == 1) %>%
  st_transform(st_crs(select_sf))#N = 152
st_crs(grid_nganorthsmall)

# Computing area per grid cell
grid_nganorthsmall_nogeom <- grid_nganorthsmall %>%
  mutate(gid_area = as.numeric(str_replace(st_area(geometry), " m^2", ""))) %>%
  as_tibble() %>%
  dplyr::select(gid, gid_area)


# select_sf is the data containing the Reuters data
# aggregating it to the grid cell level
# Note: code below should be condensed for future iterations
crs <- st_crs(grid_nganorthsmall)

new <- select_sf %>%
  mutate(select_area = as.numeric(str_replace(st_area(geometry), " m^2", ""))) %>%
  mutate(date = str_replace_all(control, "[a-z]*", "")) %>%
  dplyr::select(date,
                actor,
                adm_area = select_area) %>%
  st_intersection(grid_nganorthsmall) %>%
  mutate(clipped_area = as.numeric(str_replace_all(st_area(geometry), " m^2", ""))) %>%
  group_by(gid, date, actor) %>%
  summarize(clipped_area = sum(clipped_area)) %>%
  ungroup() %>%
  left_join(grid_nganorthsmall_nogeom, by = "gid") %>%
  mutate(prop = 1/gid_area*clipped_area,
         control_num = actor,
         control_num = replace(control_num, actor == "Boko Haram", 0),
         control_num = replace(control_num, actor == "Contested", 0.5),
         control_num = replace(control_num, actor == "Government", 1),
         control_num = as.numeric(control_num)) %>%
  group_by(gid, date) %>%
  mutate(count_section = n()) %>%
  mutate(control_test = case_when(
    count_section == 1 ~ control_num,
    count_section > 1 ~ sum(prop*control_num)
  )) %>%
  ungroup() %>%
  mutate(clipped_area = as.numeric(str_replace_all(st_area(geometry), " m^2", ""))) %>%
  st_set_geometry(NULL) %>%
  left_join(grid_nganorthsmall) %>%
  st_as_sf() %>%
  mutate(prop = 1/gid_area*clipped_area) %>%
  filter(prop > 0.33) %>%
  mutate(type = "proportion")


new_max <- select_sf %>%
  mutate(select_area = as.numeric(str_replace(st_area(geometry), " m^2", ""))) %>%
  mutate(date = str_replace_all(control, "[a-z]*", "")) %>%
  dplyr::select(date,
                actor,
                adm_area = select_area) %>%
  st_intersection(grid_nganorthsmall) %>%
  mutate(clipped_area = as.numeric(str_replace_all(st_area(geometry), " m^2", ""))) %>%
  group_by(gid, date, actor) %>%
  summarize(clipped_area = sum(clipped_area)) %>%
  ungroup() %>%
  left_join(grid_nganorthsmall_nogeom, by = "gid") %>%
  mutate(prop = 1/gid_area*clipped_area, 2) %>%
  group_by(gid, date) %>%
  arrange(-prop) %>%
  filter(prop == max(prop)) %>%
  ungroup() %>%
  mutate(clipped_area = as.numeric(str_replace_all(st_area(geometry), " m^2", ""))) %>%
  st_set_geometry(NULL) %>%
  left_join(grid_nganorthsmall) %>%
  st_as_sf() %>%
  mutate(prop = 1/gid_area*clipped_area) %>%
  filter(prop > 0.33) %>%# coded area has to make up at least 33% of the gid area to be included
  mutate(type = "max")


hmm_nganorthsmall <- hmm_nga %>%
  filter(nganorthsmall == 1) %>%
  mutate(window = case_when(
    year == 2015 & month == 2 ~ "February",
    year == 2015 & month == 3 ~ "March",
    year == 2015 & month == 4 ~ "April"
  ) %>% factor(levels = c("February", "March", "April"))) %>%
  filter(!is.na(window)) %>%
  mutate(control_num = case_when(
    control == "R" ~ 0,
    control == "DR" ~ 0.25,
    control == "D" ~ .5,
    control == "DG" ~ .75,
    control == "G" ~ 1
  )) %>%
  mutate(control_lab = case_when(
    control == "R" ~ "R",
    control == "DR" ~ "DR",
    control == "D" ~ "D",
    control == "DG" ~ "DG",
    control == "G" ~ "G"
  ) %>% factor(levels = c("R", "DR", "D", "DG", "G"))) %>%
  mutate(control_middle = case_when(
    control == "R" ~ "R",
    control == "DR" ~ "D",
    control == "D" ~ "D",
    control == "DG" ~ "D",
    control == "G" ~ "G"
  ) %>% factor(levels = c("R", "D", "G"))) %>%
  mutate(control_middle_num = case_when(
    control_middle == "R" ~ 0,
    control_middle == "D" ~ 0.5,
    control_middle == "G" ~ 1
  )) %>%
  mutate(control_extreme_num = case_when(
    control == "R" ~ 0,
    control == "DR" ~ 0,
    control == "D" ~ 0.5,
    control == "DG" ~ 1,
    control == "G" ~ 1
  ))



mycols3 <- rev(c("#2c7bb6", #Blue
                 "#ffffbf", # Yellow
                 "#d7191c"))

new_sub_max <- new_max %>%
  mutate(window = case_when(
    date == "0225" ~ "February",
    date == "0310" ~ "March",
    date == "0318" ~ "March",
    date == "0424" ~ "April"
  )) %>%
  mutate(control_test = case_when(
    actor == "Boko Haram" ~ 0,
    actor == "Contested" ~ 0.5,
    actor == "Government" ~ 1
  )) %>%
  dplyr::select(-geometry) %>%
  st_set_geometry(NULL) %>%
  dplyr::select(gid, window, date, control_test_max = control_test)

new_sub <- new %>%
  mutate(window = case_when(
    date == "0225" ~ "February",
    date == "0310" ~ "March",
    date == "0318" ~ "March",
    date == "0424" ~ "April"
  )) %>%
  dplyr::select(-geometry) %>%
  st_set_geometry(NULL) %>%
  dplyr::select(gid, window, date, control_test)

reuters_gid <- unique(new_sub$gid)

select_sf_sub <- select_sf %>%
  mutate(date = str_replace_all(control, "[a-z]", "")) %>%
  mutate(window_new = case_when(
    date == "0225" ~ "February",
    date == "0310" ~ "March1",
    date == "0318" ~ "March2",
    date == "0424" ~ "April"
  ) %>% factor(levels = c("February", "March1", "March2", "April"))) %>%
  mutate(value = case_when(
    actor == "Boko Haram" ~ 0,
    actor == "Contested" ~ 0.5,
    actor == "Government" ~ 1
  )) %>%
  mutate(name = "reuters_og") %>%
  dplyr::select(value, name, window_new, geometry)

merged_nganorthsmall <- hmm_nganorthsmall %>%
  filter(gid %in% reuters_gid) %>%
  dplyr::select(-geometry) %>%
  left_join(new_sub_max) %>%
  left_join(new_sub) %>%
  mutate(window_new = case_when(
    window == "February" ~ "February",
    window == "March" & date == "0310" ~ "March1",
    window == "March" & date == "0318" ~ "March2",
    window == "April" ~ "April"
  ) %>% factor(levels = c("February", "March1", "March2", "April"))) %>%
  dplyr::select(gid,
                window_new,
                control_test,
                control_test_max,
                control_num,
                control_middle_num,
                control_extreme_num) %>%
  pivot_longer(-c(gid, window_new)) %>%
  left_join(grid_nganorthsmall) %>%
  dplyr::select(gid, value, name, window_new, geometry) %>%
  bind_rows(select_sf_sub) %>%
  mutate(name_label = case_when(
    name == "reuters_og" ~ "Reuters\noriginal",
    name == "control_test_max" ~ "Reuters\ngridded (max)",
    name == "control_test" ~ "Reuters\ngridded (proportional)",
    name == "control_num" ~ "HMM",
    name == "control_extreme_num" ~ "HMM\n(DR = R and DG = G)",
    name == "control_middle_num" ~ "HMM\n(DR = D and DG = D)"
  ) %>% factor(levels = c("Reuters\noriginal",
                          "Reuters\ngridded (max)",
                          "Reuters\ngridded (proportional)",
                          "HMM",
                          "HMM\n(DR = R and DG = G)",
                          "HMM\n(DR = D and DG = D)"))) %>%
  st_as_sf(crs = crs) %>%
  filter(!is.na(window_new),
         !is.na(value))


merged_nganorthsmall$window_new <- factor(merged_nganorthsmall$window_new,
                                          levels = c("February",
                                                     "March1",
                                                     "March2",
                                                     "April"))



p_hmm_reuters <- ggplot() +
  geom_sf(data = merged_nganorthsmall,
          aes(fill = value),
          alpha = 1, size = 0.2,
          color = "darkgrey") +
  facet_grid(name_label ~ window_new) +
  scale_fill_gradientn(colors = mycols,
                       name = "Territorial Control",
                       breaks = c(0,0.5,1),
                       labels = c("Boko Haram", "Contested","Government")) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.key.width = unit(1.5,"cm"))
ggsave(paste0(plotdirapp, "/", "hmm_reuters.png"), 
       width = 7, height = 10, dpi = 400, p_hmm_reuters)


####################################################
# Appendix Figure 14: MSE
####################################################
mse_estimates <- merged_nganorthsmall %>%
  st_set_geometry(NULL) %>%
  filter(name %in% c("control_extreme_num",
                     "control_middle_num",
                     "control_num")) %>%
  dplyr::select(-name_label) %>%
  unique()

mse_test <- merged_nganorthsmall %>%
  st_set_geometry(NULL) %>%
  filter(str_detect(name, "test")) %>%
  dplyr::select(-name_label) %>%
  unique() %>%
  pivot_wider(names_from = name, values_from = value)

mse_merged <- left_join(mse_estimates, mse_test) %>%
  mutate(test_base = (value - control_test)^2,
         test_max = (value - control_test_max)^2) %>%
  dplyr::select(gid,
                window_new,
                hmm_name = name,
                test_base,
                test_max) %>%
  pivot_longer(cols = c(test_base, test_max)) %>%
  group_by(window_new, hmm_name, name) %>%
  summarize(av = mean(value, na.rm = T),
            sd = sd(value, na.rm = T)) %>%
  mutate(hmm_label = case_when(
    hmm_name == "control_num" ~ "HMM",
    hmm_name == "control_extreme_num" ~ "HMM (DR = R and DG = G)",
    hmm_name == "control_middle_num" ~ "HMM (DR = D and DG = D)"
  ) %>%
    factor(levels = c("HMM",
                      "HMM (DR = R and DG = G)",
                      "HMM (DR = D and DG = D)"))) %>%
  mutate(test_label = case_when(
    name == "test_base" ~ "Proportional aggregation",
    name == "test_max" ~ "Maximum aggregation"
  ) %>% factor(levels = c("Maximum aggregation",
                          "Proportional aggregation")))

p_mse <- ggplot(mse_merged,
       aes(x = window_new,
           y = av,
           fill = hmm_label)) +
  geom_bar(stat = "identity",
           position = "dodge",
           width = 0.5,
           color = "white") +
  coord_cartesian(ylim = c(0, .5)) +
  theme_light() +
  facet_wrap(~ test_label) +
  scale_fill_brewer(name = "HMM variant") +
  labs(x = "",
       y = "Average MSE",
       title = "Mean squared error")
ggsave(paste0(plotdirapp, "/", "plot_mse.png"), 
       width = 7, height = 5, dpi = 400, p_mse)

####################################################
# Appendix Figure 18: Sensitivity to m
####################################################
cdf_C_vals <- seq(0.7, 1, 0.005)
cdf_T_vals <- seq(0.7, 1, 0.005)

out <- list()
i <- 1
for (c in cdf_C_vals) {
  for(t in cdf_T_vals) {
    for(m in c(0.025, 0.05)) {
      out[[i]] <- func_test(c, t, m)
      i <- i + 1
    }
  }
}

out_df <- out %>%
  bind_rows() %>%
  mutate(value_label = case_when(
    abs == "O2" ~ "DR",
    abs == "O3" ~ "D",
    abs == "O4" ~ "DG",
  ) %>% factor(., levels = c("DR", "D", "DG"))) %>%
  mutate(mar_label = case_when(
    mar == 0.025 ~ "m = 0.025",
    mar == 0.05 ~ "m = 0.05"
  ))


## Nigeria
to_nga <- max(events_nga_cont$timeindex)

time_nga <- data.frame(timeindex = seq(1,to_nga),
                       month = rep(seq(1,12))) %>%
  mutate(year_trans = rep(seq(1, nrow(.)/12), each = 12)) %>%
  mutate(year = (2008 + year_trans)-1) %>%
  dplyr::select(-year_trans)

# Computing averages (lambda in zero-inflated Poisson)
dat_raw_nga <- events_nga_cont %>%
  st_as_sf(coords = c("lonr_c", "latr_c"),
           crs = crs) %>%
  left_join(time_nga, by = "timeindex") %>%
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


dat_mean_info_nga <- dat_raw_nga %>%
  gather(indicator, value, gtd:ged) %>%
  group_by(gwno, indicator, timeindex) %>%
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

dat_cdf_nga <- dat_raw_nga  %>%
  left_join(dat_mean_info_nga, by = c("gwno", "timeindex")) %>%
  left_join(grid_nga) %>%
  filter(nganorth == 1) %>%
  # This rounds gtd (ged) down to the next lowest full number
  dplyr::mutate(cdf_T = pzipois(floor(gtd), mean_t), 
                cdf_C = pzipois(floor(ged), mean_c),
                gid = as.numeric(str_replace_all(gid, "ngahex", ""))) %>%
  dplyr::select(gid, gwno, year, month, gtd, ged, cdf_C, cdf_T)


## Colombia
to_col <- max(events_col_cont$timeindex)

time_col <- data.frame(timeindex = seq(1,to_col),
                       month = rep(seq(1,12))) %>%
  mutate(year_trans = rep(seq(1, nrow(.)/12),each = 12)) %>%
  mutate(year = (1997 + year_trans)-1) %>%
  dplyr::select(-year_trans)


dat_raw_col <- events_col_cont %>%
  st_as_sf(coords = c("lonr_c", "latr_c"),
           crs = crs) %>%
  left_join(time_col, by = "timeindex") %>%
  mutate(gwno = 100) %>%
  dplyr::select(gid, 
                gwno,
                year,
                month,
                timeindex,
                colwest,
                gtd = weighted_terrorism_trunc,
                ged = weighted_combats_trunc) %>%
  dplyr::mutate(gtd = replace(gtd, is.na(gtd), 0),
                ged = replace(ged, is.na(ged), 0)) %>%
  dplyr::mutate(gtd_floor = floor(gtd),
                ged_floor = floor(ged),
                gtd_round = round(gtd),
                ged_round = round(ged)) %>%
  st_set_geometry(NULL)  %>%
  filter(colwest == 1)


dat_mean_info_col <- dat_raw_col %>%
  gather(indicator, value, gtd:ged) %>%
  group_by(gwno, indicator, timeindex) %>%
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

dat_cdf_col <- dat_raw_col  %>%
  left_join(dat_mean_info_col, by = c("gwno", "timeindex")) %>%
  # This rounds gtd (ged) down to the next lowest full number
  dplyr::mutate(cdf_T = pzipois(floor(gtd), mean_t), 
                cdf_C = pzipois(floor(ged), mean_c),
                gid = as.numeric(str_replace_all(gid, "colhex", ""))) %>%
  filter(year %in% seq(1994,2017)) %>%
  dplyr::select(gid, gwno, year, month, gtd, ged, cdf_C, cdf_T)


dat <- bind_rows(dat_cdf_nga, dat_cdf_col) %>%
  mutate(gwno_label = case_when(
    gwno == 100 ~ "Colombia",
    gwno == 475 ~ "Nigeria"
  ))


p_sim <- ggplot() +
  geom_tile(data = out_df,
            aes(x = cdf_C, y = cdf_T, fill = value_label)) +
  theme_light() +
  scale_fill_manual(values = rev(c("#abd9e9",
                                   "#ffffbf", # Yellow
                                   "#fdae61")),
                    name = "Emission value") +
  labs(x = "Conventional fighting (zero-inflated Poisson)",
       y = "Terrorism (zero-inflated Poisson)") +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "top") +
  coord_fixed(expand = c(0.01,0.01),
              xlim = c(0.7, 1),
              ylim = c(0.7, 1)) +
  geom_point(dat = dat,
             aes(x = cdf_C, y = cdf_T),
             size = 0.25, alpha = 0.1) +
  facet_grid(mar_label ~ gwno_label) +
  geom_rug(dat = dat,
           aes(x = cdf_C, y = cdf_T),
           color = "black")
ggsave(paste0(plotdirapp,"/",  "plot_cdf_nga_col.png"), 
       width = 8, height = 8, dpi = 400, p_sim)
round(median(dat$cdf_C[dat$gwno == 100]), 2)
round(median(dat$cdf_T[dat$gwno == 100]), 2)
round(median(dat$cdf_C[dat$gwno == 475]), 2)
round(median(dat$cdf_T[dat$gwno == 475]), 2)


######################################################
# Appendix Figure 19: Schematic to show reporting bias
######################################################

m_t_true <- 0.7
m_t_bias <- 0.35
m_c_true <- 0.3
m_c_bias <- 0.15

df_bias <- data.frame(val = seq(0,2, 1)) %>%
  mutate(t_true = pzipois(val, m_t_true),
         t_bias = pzipois(val, m_t_bias),
         c_true = pzipois(val, m_c_true),
         c_bias = pzipois(val, m_c_bias)) %>%
  pivot_longer(-val) %>%
  separate(name, c("violence", "variance"))


p_bias <- ggplot(df_bias,
       aes(x = val, 
           y = value,
           linetype = variance,
           color = violence,
           fill = violence)) +
  geom_line(alpha = 0.2) +
  geom_hline(yintercept = pzipois(1, m_t_true), color = "black") +
  geom_hline(yintercept = pzipois(0, m_c_true), color = "darkcyan") +
  geom_hline(yintercept = pzipois(1, m_t_bias), color = "black", linetype = "dashed") +
  geom_hline(yintercept = pzipois(0, m_c_bias), color = "darkcyan", linetype = "dashed") +
  annotate("rect",
           ymin = pzipois(1, m_t_true) - 0.025/2, ymax = pzipois(1, m_t_true) + 0.025/2,
           xmin = 0 - 0.02, xmax = 1 + 0.02, alpha = 0.1, fill = "black", color = "black") +
  annotate("rect",
           ymin = pzipois(1, m_t_bias) - 0.025/2, ymax = pzipois(1, m_t_bias) + 0.025/2,
           xmin = 0- 0.02, xmax = 1+ 0.02, alpha = 0.1, fill = "black", color = "black", linetype = "dashed") +
  annotate("rect",
           ymin = pzipois(0, m_c_true) - 0.025/2, ymax = pzipois(0, m_c_true) + 0.025/2,
           xmin = 0- 0.02, xmax = 0 + 0.02, alpha = 0.1, fill = "darkcyan", color = "darkcyan") +
  annotate("rect",
           ymin = pzipois(0, m_c_bias) - 0.025/2, ymax = pzipois(0, m_c_bias) + 0.025/2,
           xmin = 0- 0.02, xmax = 0 + 0.02, alpha = 0.1, fill = "darkcyan", color = "darkcyan", linetype = "dashed") +
  theme_light() +
  scale_linetype_manual(name = "Version",
                        labels = c("Biased",
                                   "True"),
                        values = c("dashed", "solid")) +
  scale_color_manual(name = "Tactics",
                     values = c("darkcyan", "black"),
                     labels = c("Conventional", "Terrorism")) +
  scale_fill_manual(name = "Tactics",
                    values = c("darkcyan", "black"),
                    labels = c("Conventional", "Terrorism")) +
  labs(x = "Observed number of events",
       y = "Probability from zero-inflated Poisson distribution") +
  scale_x_continuous(breaks = seq(0,2)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)),
         linetype = guide_legend(override.aes = list(alpha = 1)))
ggsave(paste0(plotdirapp, "/", "plot_bias.png"),
       width = 6.5, height = 4, dpi = 300, p_bias)


######################################################
# Appendix Table 7: Poisson regression results
######################################################
# See README.md for details (only subset of PRIO GRID data included in repo)
# PRIO grid pop data available for years 1990, 1995, 2000, and 2005

df_count <- bind_rows(events_nga, events_col) %>%
  st_set_geometry(NULL) %>%
  group_by(year, type, priogrid) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(gtd = ifelse(type == "terrorism", 1, 0)) %>%
  dplyr::select(-type) %>%
  arrange(priogrid, year)


pgid_nga <- unique(df_priogrid$priogrid[df_priogrid$case == "nga"])
pgid_col <- unique(df_priogrid$priogrid[df_priogrid$case == "col"])
yr_nga <- unique(df_priogrid$year[df_priogrid$case == "nga"])
yr_col <- unique(df_priogrid$year[df_priogrid$case == "col"])


## NGA (effectively only 2010)
base_nga <- expand_grid(priogrid = pgid_nga,
                        year = yr_nga,
                        gtd = c(0, 1)) %>%
  left_join(df_priogrid) %>%
  left_join(df_count) %>%
  mutate(count = replace(count, is.na(count), 0)) %>%
  filter(year >= 2008)

m1_nga <- glm(count ~ ln_pop,
              data = base_nga, family = "poisson")
m2_nga <- glm(count ~ ln_pop + factor(gtd),
              data = base_nga, family = "poisson")
m3_nga <- glm(count ~ ln_pop*factor(gtd),
              data = base_nga, family = "poisson")


## COL
base_col <- data.frame(expand_grid(priogrid = pgid_col,
                                   year = yr_col,
                                   gtd = c(0, 1))) %>%
  left_join(df_priogrid) %>%
  left_join(df_count) %>%
  mutate(count = replace(count, is.na(count), 0)) %>%
  filter(year >= 1997)


m1_col <- glm(count ~ ln_pop + factor(year),
              data = base_col, family = "poisson")
m2_col <- glm(count ~ ln_pop + factor(year) + factor(gtd),
              data = base_col, family = "poisson")
m3a_col <- glm(count ~ ln_pop*factor(gtd) + factor(year),
               data = base_col, family = "poisson")
m3b_col <- glm(count ~ ln_pop*factor(gtd),
               data = base_col, family = "poisson")
m3c_col <- glm(count ~ ln_pop*factor(gtd),
               data = subset(base_col, year == 2010), 
               family = "poisson")

stargazer::stargazer(m1_nga, m2_nga, m3_nga,
                     m1_col, m2_col, m3a_col, m3b_col, 
                     title = "Assessment of imbalance of potential underreporting bias between the GED and GTD data sets.",
                     dep.var.labels  = "Count of events per PRIO grid cell",
                     column.labels = c("Nigeria", "Colombia"),
                     column.separate = c(3, 4),
                     notes = c("Poisson regression results with standard errors in parentheses.",
                               "Nigeria estimates are for 2010 only.",
                               "Base category for year dummies  for the Colombia estimates is 2000."),
                     out = paste0(plotdirapp, "/",  "tab_underreporting.tex"),
                     font.size = "scriptsize",
                     digits = 2,
                     column.sep.width = "-6pt",
                     order = c(3, 1, 2, 4, 5),
                     covariate.labels = c("Population (ln) x Terrorism",
                                          "Population (ln)",
                                          "Terrorism",
                                          "2005",
                                          "2010"),
                     star.cutoffs = c(0.05, 0.01, 0.001),
                     label = "tab:underreporting")

# Predicted number of events
s1 <- data.frame(ln_pop = mean(base_col$ln_pop),
                 gtd = c(0,1)) %>% mutate(int = ln_pop*gtd)

predict(m3b_col, s1, type = "response", se.fit = T)


######################################################
# Appendix Figure 20: Case selection
######################################################

# Download and unzip replication data for Polo & Gleditsch 2016
# http://file.prio.no/Journals/JPR/2016/53/6/Sara%20MT%20Polo%20&%20Kristian%20Skrede%20Gleditsch.zip
# Save in folder `pologleditsch2016`
pg <- haven::read_dta("~/Downloads/pologleditsch2016/PG_replication_final.dta")
pg_sub <- pg %>%
  dplyr::select(year, ucdpgroup, ccode, country, troopratio_all) %>%
  filter(year >= 1994)


rebgovstrength_by_country <- pg_sub %>%
  group_by(ccode, country) %>% 
  dplyr::summarize(max_rebgovstrength = max(troopratio_all, na.rm = T),
                   min_rebgovstrength = min(troopratio_all, na.rm = T),
                   av_rebgovstrength = mean(troopratio_all, na.rm = T)) %>%
  pivot_longer(max_rebgovstrength:av_rebgovstrength) %>%
  separate(name, c("measure", "indicator")) %>% 
  dplyr::mutate(value = replace(value, value %in% c("NaN", "Inf", "-Inf"), NA)) %>%
  arrange(measure, value) %>%
  mutate(include05 = ifelse(value <= 0.5, 1, 0))

p_rebgovstrength <- ggplot(subset(rebgovstrength_by_country, measure != "min"), 
       aes(x = country, y = value, fill = factor(include05))) +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() +
  facet_wrap(~measure) +
  scale_fill_manual(name = "Include at 0.5 threshold",
                    labels = c("No", "Yes"),
                    values = c("black", "red")) +
  theme(legend.position = "top") +
  labs(title = "Case selection based on power asymmetry",
       subtitle = "Data on troops ratios from Polo and Gleditsch (2016)",
       x = "",
       y = "Troop ratio")
ggsave(paste0(plotdirapp, "/", "plot_rebgovstrength_cases.png"), 
       width = 8, height = 8, dpi = 400, p_rebgovstrength)

######################################################
# Appendix Figure 21: Cases threshold
######################################################


# Inclusion of cases
steps <- seq(0, 1, 0.01)
numcases_av <- c()
numcases_max <- c()
for(i in 1:length(steps)){
  
  dat_av <- rebgovstrength_by_country %>%
    filter(measure == "av",
           value >= 0.01,
           value <= steps[i])
  
  dat_max <- rebgovstrength_by_country %>%
    filter(measure == "max",
           value >= 0.01,
           value <= steps[i])
  
  numcases_av[i] <- nrow(dat_av) 
  numcases_max[i] <- nrow(dat_max) 
  
}

df_numcases <- data.frame("step" = steps,
                          "av" = numcases_av,
                          "max" = numcases_max) %>%
  pivot_longer(cols = av:max)

p_threshold <- ggplot(df_numcases, 
       aes(x = step, y = value, color = name)) +
  geom_line() +
  theme_bw() +
  labs(title = "Number of included cases based on strenght threshold",
       x = "Strength threshold",
       y = "Number of cases") +
  scale_color_manual(name = "Threshold measure",
                     labels = c("Averaged across conflict years",
                                "Maximum across conflict years"),
                     values = c("cornflowerblue", "goldenrod2")) +
  theme(legend.position = "top")
ggsave(paste0(plotdirapp, "/", "plot_rebgovstrength_threshold.png"), 
       width = 6, height = 4, dpi = 400, p_threshold)

