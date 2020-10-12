############################
# Territorial Control paper
# Graphs and figures
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

source(here::here("utils","func_weights.R"))
source(here::here("utils","func_hmm.R"))

plotdir <- here::here("plots")


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


#################################
#
#
# Main manuscript
#
#
#################################



####################################################
# Figure 1: Territorial control and tactical choice
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
table(select_sf$control)

# Joining point and aereal data
merged_sub <- st_join(events_nga, select_sub) %>%
  filter(!is.na(NAME_1)) %>%
  mutate(window = case_when(
    date_start %within% interval(ymd("2015-02-25"), ymd("2015-03-11")) ~ "25 February",
    date_start %within% interval(ymd("2015-03-10"), ymd("2015-03-24")) ~ "10 March",
    date_start %within% interval(ymd("2015-03-18"), ymd("2015-04-01")) ~ "18 March",
    date_start %within% interval(ymd("2015-04-24"), ymd("2015-05-08")) ~ "24 April"
  ) %>% factor(levels = c("25 February", "10 March", "18 March", "24 April"))) %>%
  filter(!is.na(window)) %>%
  mutate(type_label = case_when(
    type == "terrorism" ~ "Terrorist attack (GTD data)",
    type == "conventional" ~ "Conventional fighting (GED data)"
  ) %>% factor(levels = c("Conventional fighting (GED data)", "Terrorist attack (GTD data)")))

# Plotting map
p_map <- ggplot() +
  geom_sf(data = select_sf,
          aes(fill = actor),
          size = 0.2, 
          color = "grey98", 
          alpha = 0.2) +
  facet_wrap(~ window, nrow = 1) +
  geom_sf(data = merged_sub,
          aes(color = type_label, shape = type_label),
          alpha = 0.7,
          size = 1.7,
          show.legend = "point") +
  coord_sf(datum = NA) + 
  scale_fill_manual(name = "Territorial Control", 
                    values = rev(c("#2c7bb6",
                                   "goldenrod",
                                   "#d7191c")),
                    guide = guide_legend(override.aes = list(linetype = "blank", 
                                                             shape = NA))) +
  labs(title = "Territorial control and conflict events in NE Nigeria in 2015",
       subtitle = "Conflict events within two weeks of observing territorial control",
       x = "",
       y = "") +
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.key = element_blank()) +
  
  scale_color_manual(values = c("#d7191c", #red
                                "#2c7bb6"),
                     name = "Rebel tactics",
                     guide = guide_legend(override.aes = list(fill = NA,
                                                              size = 5))) +
  scale_shape_manual(values = c(19,17),
                     name = "Rebel tactics")
ggsave(paste0(plotdir, "/", "reuters_map.png"), width = 10, height = 4.5, dpi = 500, p_map)


# bw version
p_map_bw <- ggplot() +
  geom_sf(data = select_sf,
          aes(fill = actor),
          size = 0.1, 
          color = "black", 
          alpha = 0.5) +
  facet_wrap(~ window, nrow = 1) +
  geom_sf(data = merged_sub,
          aes(shape = type_label),
          alpha = 0.7,
          size = 1.7,
          show.legend = "point") +
  coord_sf(datum = NA) + 
  scale_fill_manual(values = c("black", "grey", "white"),
                    name = "Territorial Control",
                  guide = guide_legend(override.aes = list(shape = NA))) +
  labs(title = "Territorial control and conflict events in NE Nigeria in 2015",
       subtitle = "Conflict events within two weeks of observing territorial control",
       x = "",
       y = "") +
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.key = element_blank()) +
  # scale_color_manual(name = "Rebel tactics",
  #                    values = c("black", "black"),
  #                    guide = guide_legend(override.aes = list(fill = NA,
  #                                                             size = 5))) +
  scale_shape_manual(values = c(19,17),
                     name = "Rebel tactics",
                     guide = guide_legend(override.aes = list(fill = NA,
                                                              size = 5)))
ggsave(paste0(plotdir, "/", "reuters_map_bw.png"), width = 10, height = 4.5, dpi = 500, p_map_bw)



####################################################
# Figure 4: Territorial control estimates Colombia 
####################################################

mycols <- rev(c("#2c7bb6", #Blue
                "#abd9e9",
                "#ffffbf", # Yellow
                "#fdae61",
                "#d7191c"))
names(mycols) <- levels(hmm_col$control_lab)

mycols_bw <- rev(c("#f5f5f5", #light grey
                "#e0e0e0",
                "#a3a3a3", # Yellow
                "#666666",
                "#000000"))
names(mycols_bw) <- levels(hmm_col$control_lab)


# Numerical values for control levels
hmm_col_monthly <- hmm_col %>%
  mutate(control_num = case_when(
    control == "R" ~ 0,
    control == "DR" ~ 0.25,
    control == "D" ~ 0.5,
    control == "DG" ~ 0.75,
    control == "G" ~ 1
  ))
         
hmm_col_yearly <- hmm_col_monthly %>%
  group_by(gid, year) %>%
  summarise(control_mean = mean(control_num),
            control_median = median(control_num),
            control_sd = sd(control_num),
            control_any_s1s2 = ifelse(any(control %in% c("S1", "S2")), 1, 0),
            control_any_s1 = ifelse(any(control %in% c("S1")), 1, 0)) %>% 
  left_join(grid_col, by = "gid") %>%
  st_as_sf()

# Outline for Colombia
col_out <- getData('GADM', country = 'COL', level = 0)
col_out_simple <- gSimplify(col_out, tol = 0.01, topologyPreserve=TRUE) 
col_out_sf <- col_out_simple %>%
  st_as_sf()

## Plotting Figure 4
p_col_yearly <- ggplot() +
  geom_sf(data = col_out_sf,
          alpha = 0.2,
          size = 0.05) +
  geom_sf(data = subset(hmm_col_yearly, year %in% seq(2006,2017)),
          aes(fill = control_mean), 
          size = 0.01, color = "lightgrey",
          alpha = 1) +
  facet_wrap(~year, nrow = 3) +
  theme_bw() +
  scale_fill_gradientn(name = "",
                       colors = mycols,
                       breaks = c(0, 0.5, 1),
                       limits=c(0,1),
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
  labs(title = "Yearly averages of monthly estimates of territorial control")  +
  coord_sf(xlim = c(-79, -67), ylim = c(-4,12.5), datum = NA)
ggsave(paste0(plotdir, "/", "hmm_col_yearly.png"), width = 8, height = 10, dpi = 500, p_col_yearly)


## Plotting Figure 4 (print bw version)
p_col_yearly_bw <- ggplot() +
  geom_sf(data = col_out_sf,
          alpha = 0.2,
          size = 0.05) +
  geom_sf(data = subset(hmm_col_yearly, year %in% seq(2006,2017)),
          aes(fill = control_mean), 
          size = 0.01, 
          color = "black",
          alpha = 1) +
  facet_wrap(~year, nrow = 3) +
  theme_bw() +
  scale_fill_gradientn(name = "",
                       colors = mycols_bw,
                       breaks = c(0, 0.5, 1),
                       limits=c(0,1),
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
  labs(title = "Yearly averages of monthly estimates of territorial control")  +
  coord_sf(xlim = c(-79, -67), ylim = c(-4,12.5), datum = NA)
ggsave(paste0(plotdir, "/", "hmm_col_yearly_bw.png"), width = 8, height = 10, dpi = 400, p_col_yearly_bw)

####################################################
# Figure 5: Territorial control estimates Nigeria 
####################################################
hmm_nga_monthly <- hmm_nga %>%
  mutate(control_num = case_when(
    control == "R" ~ 0,
    control == "DR" ~ 0.25,
    control == "D" ~ 0.5,
    control == "DG" ~ 0.75,
    control == "G" ~ 1
  ))

##### Creating an annual data set
hmm_nga_yearly <- hmm_nga_monthly %>%
  group_by(gid, year) %>%
  summarise(control = mean(control_num)) %>% 
  left_join(grid_nga, by = "gid") %>%
  st_as_sf()

# Outline for Nigeria
nga_out <- getData('GADM', country = 'NGA', level = 0)
nga_out_simple <- gSimplify(nga_out, tol = 0.005, topologyPreserve=TRUE) 
nga_out_sf <- nga_out_simple %>%
  st_as_sf()

# Plotting figure 5
p_nga_yearly <- ggplot(subset(hmm_nga_yearly, year > 2008)) +
  geom_sf(data = nga_out_sf,
          alpha = 0.6,
          size = 0.1) +
  geom_sf(aes(fill = control), size = 0.01, color = "lightgrey") +
  facet_wrap(~year, nrow = 3) +
  theme_bw() +
  scale_fill_gradientn(name = "",
                       colors = mycols,
                       breaks = c(0, 0.5, 1),
                       limits=c(0,1),
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
  labs(title = "Yearly averages of monthly estimates of territorial control in NE Nigeria") +
  coord_sf(datum = NA)
ggsave(paste0(plotdir, "/", "hmm_nga_yearly.png"), width = 8, height = 8, dpi = 400, p_nga_yearly)


# Plotting figure 5
p_nga_yearly_bw <- ggplot(subset(hmm_nga_yearly, year > 2008)) +
  geom_sf(data = nga_out_sf,
          alpha = 0.9,
          size = 0.1) +
  geom_sf(aes(fill = control), size = 0.01, color = "black") +
  facet_wrap(~year, nrow = 3) +
  theme_bw() +
  scale_fill_gradientn(name = "",
                       colors = mycols_bw,
                       breaks = c(0, 0.5, 1),
                       limits=c(0,1),
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
  labs(title = "Yearly averages of monthly estimates of territorial control in NE Nigeria") +
  coord_sf(datum = NA)
ggsave(paste0(plotdir, "/", "hmm_nga_yearly_bw.png"), width = 8, height = 8, dpi = 500, p_nga_yearly_bw)



####################################################
# Figure 6: Acled validation
####################################################
table(hmm_nga_monthly$control)
acled <- bind_rows(acled_ls)
merged_monthly <- acled %>%
  filter(temporal == "monthly") %>%
  dplyr::select(gid, year, month, control_test, type, duration) %>%
  left_join(hmm_nga_monthly)

# Computing correlations
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


# Plotting Figure 6
sub_paper <- df_monthly %>%
  filter(duration == "6m",
         type == "full") %>%
  filter(year >= 2011)

# range in paper
round(range(sub_paper$value, na.rm = T), 2)


p_acled_paper <- ggplot(sub_paper,
       aes(x = timeindex,
           y = value)) +
  geom_smooth(size = 0.5, color = "black") +
  geom_line(size = 0.25) +
  labs(title = "Monthly correlation between HMM estimates and ACLED testing data",
       subtitle = "2010 - 2017",
       x = "",
       y = "Spearman correlation coefficient") +
  scale_x_continuous(labels = x_l,
                     breaks = x_b) +
  theme_light() +
  coord_cartesian(ylim = c(0, 0.55))
ggsave(paste0(plotdir, "/", "cor_acled_paper.png"), width = 7, height = 4, dpi = 500, p_acled_paper)
