#####################
# Territorial control
# Prepping event data
#
# Therese Anders
#####################


library(tidyverse)
library(stringr)
library(lubridate)
library(sf)
library(lattice)
library(raster)
library(rgeos)

source(here::here("utils", "func_grid.R"))

###########################
# GTD data 
###########################

# This is an sf data frame with PRIO Grid 2.0 IDs attached (only COL & NGA)
# Does not have lat/lon columns b/c geographic data contained in geometry
gtd <- readRDS(here::here("raw", "gtd_attached_prioid.rds"))

# 1) Preparing the raw data
gtd_sub <- gtd %>%
  
  # Rename columns
  dplyr::rename(year = iyear,
                month = imonth,
                day = iday) %>%
  
  # Code new variables
  dplyr::mutate(doubtterr = replace(doubtterr, doubtterr == -9, NA),
                date_start = ymd(paste(year, month, day)),
                type = "terrorism",
                date_end = ymd(resolution),
                duration = date_end - date_start + 1,
                duration = replace(duration, extended == 0, 1)) %>%
  
  ## Coding a time precision variable
  dplyr::mutate(date_prec = ifelse(approxdate == "", 1, NA),
                date_prec = replace(date_prec, approxdate != "", 2), #Few days - matched to less than week in other data sets
                date_prec = replace(date_prec, month != 0 & day == 0, 3),
                date_prec = replace(date_prec, month == 0, 5)) %>%
  
  # Creating lat/lon columns
  rowwise() %>%
  mutate(longitude = unlist(geometry)[1],
         latitude = unlist(geometry)[2]) %>%
  ungroup() %>%
  
  # Filtering out events 
  # a) not clearly terrorism
  # b) not attributable to min 2nd order admin region
  # c) Main model: removing attacks against military
  dplyr::filter(doubtterr == 0,
                specificity %in% c(1, 2, 3),
                !(targtype1 %in% c(4))) %>%

  dplyr::select(priogrid = gid,
                year,
                month,
                day,
                date_start,
                date_end,
                date_prec,
                extended,
                duration,
                country = country_txt,
                longitude,
                latitude,
                prec_tax = specificity,
                fatalities = nkill,
                perpetrator = gname,
                type,
                event_tax = attacktype1,
                desc = summary,
                doubtterr,
                targtype1,
                attacktype1_txt)

# 2) Creating alternative subsets (for robustness checks for Colombia data)
## Removing assassinations
gtd_sub_noassass <- gtd_sub %>%
  filter(attacktype1_txt != "Assassination")

## Additionally removing all official targets
gtd_sub_noofficial <- gtd_sub %>%
  filter(attacktype1_txt != "Assassination",
         !(targtype1 %in% c(2, 3, 7)))


# 3) Subsetting to select actors in Colombia & Nigeria
table(gtd_sub$perpetrator[gtd_sub$country == "Colombia"])
act_col <- c("Revolutionary Armed Forces of Colombia (FARC)",
             "Revolutionary Armed Forces of Colombia (FARC) dissidents")

gtd_select_col <- gtd_sub %>%
  filter(year >= 1997,
         year <= 2017,
         country == "Colombia",
         perpetrator %in% act_col)

gtd_select_noassass_col <- gtd_sub_noassass %>%
  filter(year >= 1997,
         year <= 2017,
         country == "Colombia",
         perpetrator %in% act_col)

gtd_select_noofficial_col <- gtd_sub_noofficial %>%
  filter(year >= 1997,
         year <= 2017,
         country == "Colombia",
         perpetrator %in% act_col)

gtd_select_nga <- gtd_sub %>%
  filter(year >= 2008,
         year <= 2017,
         country == "Nigeria",
         perpetrator %in% c("Boko Haram", 
                            "Al-Qaida in the Islamic Maghreb (AQIM)"))


###########################
# GED data 
###########################
# Reading data
# Note: raw data can be downloaded from the GED website
# ged <- read.csv("~/Dropbox/Data/GED/ged181.csv",
#                 stringsAsFactors = F)

ged_sub <- ged %>%

  dplyr::mutate(date_start = ymd(date_start),
                date_end = ymd(date_end),
                extended = ifelse(date_start != date_end, 1, 0),
                duration = date_end-date_start+1,
                type = "conventional",
                day = day(date_start),
                month = month(date_start),
                date_end = replace(date_end, date_end == date_start, NA)) %>%
  
  # Selecting relevant observations
  # a) state-based conflict. In state- based armed conflicts, at least one of the primary parties must be the government of a state.
  # b) not attributable to min 2nd order admin region
  dplyr::filter(type_of_violence == 1,
                where_prec %in% c(1, 2, 3)) %>%
  
  # Selecting relevant variables
  dplyr::select(year,
                month,
                day,
                date_start,
                date_end,
                date_prec,
                extended,
                duration,
                country,
                longitude,
                latitude,
                prec_tax = where_prec,
                fatalities = best,
                perpetrator = side_b,
                type,
                event_tax = type_of_violence,
                desc = source_headline,
                desc_add = source_article,
                priogrid = priogrid_gid)

# 2) Subsetting to Colombia and Nigeria
## Data for Colombia (FARC)
ged_select_col <- ged_sub %>%
  filter(year >= 1997, #doubtterr starts in 1997
         year <= 2017,
         country == "Colombia",
         perpetrator %in% c("FARC"))

ged_select_nga <- ged_sub %>%
  filter(year >= 2008, 
         year <= 2017,
         country == "Nigeria",
         perpetrator %in% c("Jama'atu Ahlis Sunna Lidda'awati wal-Jihad", 
                            "IS"))


# Merging
merged_col <- bind_rows(gtd_select_col, ged_select_col) %>%
  dplyr::mutate(long_old = longitude,
                lat_old = latitude,
                eventid = paste0("id", seq(1, nrow(.))))

merged_col_noassass <- bind_rows(gtd_select_noassass_col, ged_select_col) %>%
  dplyr::mutate(long_old = longitude,
                lat_old = latitude,
                eventid = paste0("id", seq(1, nrow(.))))

merged_col_noofficial <- bind_rows(gtd_select_noofficial_col, ged_select_col) %>%
  dplyr::mutate(long_old = longitude,
                lat_old = latitude,
                eventid = paste0("id", seq(1, nrow(.))))

merged_nga <- bind_rows(gtd_select_nga, ged_select_nga) %>%
  dplyr::mutate(long_old = longitude,
                lat_old = latitude,
                eventid = paste0("id", seq(1, nrow(.))))

# Now generate grid and aggregate point data
## Selecting common set of variables
names(ged_sub)
names(gtd_sub)

vars_select <- c("gid",
                 "year",
                 "month",
                 "day",
                 "date_start",
                 "date_end",
                 "date_prec",
                 "extended",
                 "duration",
                 "country",
                 "longitude",
                 "latitude",
                 "prec_tax",
                 "fatalities",
                 "perpetrator",
                 "type",
                 "event_tax",
                 "desc",
                 "desc_add",
                 "eventid",
                 "doubtterr",
                 "targtype1",
                 "priogrid")


############################
# 
#
# COLOMBIA
#
#
############################

#########################
# Reading geographic data
#########################
getData('ISO3')
####
# Colombia
set.seed(12345)
col <- getData('GADM', country = 'COL', level = 1)
crs <- proj4string(col)

col_sf <- col %>%
  st_as_sf() 

## Creating grid - rectangle
seed <- 12345
col_rect <- withr::with_seed(seed,
                             st_make_grid(col_sf,
                                          n = c(1,1))) %>%
  as("Spatial")

## Creating grid - hexagons 0.25 diameter
grids_col_sf <- withr::with_seed(seed,
                                 make_grid(col_rect, 
                                           type = "hexagonal", 
                                           cell_width = 0.25, clip = F)) %>%
  st_as_sf() %>%
  dplyr::mutate(gid = paste0("colhex", seq(1, nrow(.))))


# Intersection of hexagons with geographic area
intersect_col <- as_tibble(st_intersection(grids_col_sf, col_sf)) %>%
  st_sf() %>%
  st_set_geometry(NULL)

intersectgid_col <- unique(intersect_col$gid)

#### Subsetting regions
####
# Colombia west (excluding Orinoco and Amazonia)
r_andina <- c("CO.AN", 
              "CO.BY", 
              "CO.CL", 
              "CO.CU", 
              "CO.HU", 
              "CO.NS", 
              "CO.QD", 
              "CO.RI", 
              "CO.ST", 
              "CO.TO")
r_amazonia <- c("CO.AM", 
                "CO.CQ", 
                "CO.GN", 
                "CO.GV", 
                "CO.PU", 
                "CO.VP")
r_orinoquia <- c("CO.AR",
                 "CO.CS",
                 "CO.ME",
                 "CO.VD")
r_caribe <- c("CO.AT",
              "CO.BL",
              "CO.CE",
              "CO.CO",
              "CO.LG",
              "CO.MA",
              "CO.SA",
              "CO.SU")
r_pacifica <- c("CO.CA",
                "CO.CH",
                "CO.NA",
                "CO.VC")

col_sf_sub <- col_sf %>%
  mutate(region = case_when(
    HASC_1 %in% r_andina ~ "Andean",
    HASC_1 %in% r_amazonia ~ "Amazon",
    HASC_1 %in% r_orinoquia ~ "Orinoco",
    HASC_1 %in% r_caribe ~ "Caribbean",
    HASC_1 %in% r_pacifica ~ "Pacific"
  )) %>%
  filter(!(region %in% c("Amazon", "Orinoco")))

# Subsetting to just Western regions
# Intersection of hexagons with geographic area
intersect_colwest <- as_tibble(st_intersection(grids_col_sf, col_sf_sub)) %>%
  st_sf() %>%
  st_set_geometry(NULL)

intersectgid_colwest <- unique(intersect_colwest$gid)

##### Creating data frames

## 1)  Base data frame
# Geo-coding events
merged_col_sf <- merged_col %>%
  st_as_sf(coords = c("longitude", "latitude"),
           crs = crs)
col_cols <- names(merged_col_sf)

# Filtering events within the states of the study region
merged_col_sf_filter <- merged_col_sf %>%
  st_intersection(col_sf) %>%
  dplyr::select(one_of(col_cols))

# Merging events and hexagon grid cells
new_col <- merged_col_sf_filter %>%
  st_join(grids_col_sf) %>%
  dplyr::rename(longitude = long_old,
                latitude = lat_old) %>%
  filter(!is.na(eventid)) %>%
  dplyr::select(one_of(vars_select)) %>%
  mutate(colwest = ifelse(gid %in% intersectgid_colwest, 1, 0))

centroid_col <- grids_col_sf %>%
  st_centroid() %>%
  st_coordinates() %>%
  as_tibble() %>%
  dplyr::mutate(gid = paste0("colhex", seq(1, nrow(.))))

grids_col_sf_new <- grids_col_sf %>%
  mutate(countrypoly = ifelse(gid %in% intersectgid_col, 1, 0),
         colwest = ifelse(gid %in% intersectgid_colwest, 1, 0)) %>%
  left_join(centroid_col)

## Saving data
saveRDS(grids_col_sf_new, here::here("prepped", "grids_hex_col_hex25.rds"))
saveRDS(new_col, here::here("prepped", "events_col_hex25.rds"))


## 2)  No assassinations
# Geo-coding events
merged_col_noassass_sf <- merged_col_noassass %>%
  st_as_sf(coords = c("longitude", "latitude"),
           crs = crs)

# Filtering events within the states of the study region
merged_col_noassass_sf_filter <- merged_col_noassass_sf %>%
  st_intersection(col_sf) %>%
  dplyr::select(one_of(col_cols))

# Merging events and hexagon grid cells
new_col_noassass <- merged_col_noassass_sf_filter %>%
  st_join(grids_col_sf) %>%
  dplyr::rename(longitude = long_old,
                latitude = lat_old) %>%
  filter(!is.na(eventid)) %>%
  dplyr::select(one_of(vars_select)) %>%
  mutate(colwest = ifelse(gid %in% intersectgid_colwest, 1, 0))


## Saving data
saveRDS(new_col_noassass, here::here("prepped", "events_col_noassass_hex25.rds"))

## 2)  No official targets
# Geo-coding events
merged_col_noofficial_sf <- merged_col_noofficial %>%
  st_as_sf(coords = c("longitude", "latitude"),
           crs = crs)

# Filtering events within the states of the study region
merged_col_noofficial_sf_filter <- merged_col_noofficial_sf %>%
  st_intersection(col_sf) %>%
  dplyr::select(one_of(col_cols))

# Merging events and hexagon grid cells
new_col_noofficial <- merged_col_noofficial_sf_filter %>%
  st_join(grids_col_sf) %>%
  dplyr::rename(longitude = long_old,
                latitude = lat_old) %>%
  filter(!is.na(eventid)) %>%
  dplyr::select(one_of(vars_select)) %>%
  mutate(colwest = ifelse(gid %in% intersectgid_colwest, 1, 0))


## Saving data
saveRDS(new_col_noofficial, here::here("prepped", "events_col_noofficial_hex25.rds"))


############################
# 
#
# NIGERIA
#
#
############################

# Getting geographic data
nga <- getData('GADM', country = 'NGA', level = 2)
crs <- proj4string(nga)

nga_sf <- nga %>%
  st_as_sf() 

## Creating grid - rectangle
seed <- 12345
nga_rect <- withr::with_seed(seed,
                             st_make_grid(nga_sf,
                                          n = c(1,1))) %>%
  as("Spatial")

## Creating grid - hexagons 0.25 diameter
grids_nga_sf <- withr::with_seed(seed,
                                 make_grid(nga_rect, 
                                           type = "hexagonal", 
                                           cell_width = 0.25, clip = F)) %>%
  st_as_sf() %>%
  dplyr::mutate(gid = paste0("ngahex", seq(1, nrow(.))))


# Intersection of hexagons with geographic area
intersect_nga <- as_tibble(st_intersection(grids_nga_sf, nga_sf)) %>%
  st_sf() %>%
  st_set_geometry(NULL)

intersectgid_nga <- unique(intersect_nga$gid)

# Subsetting regions
## A) nganorth
states <- c("Adamawa",
            "Bauchi",
            "Benue",
            "Borno",
            "Gombe",
            "Jigawa",
            "Kaduna",
            "Kano",
            "Katsina",
            "Nassarawa",
            "Niger",
            "Plateau",
            "Taraba",
            "Yobe",
            "Federal Capital Territory")

nganorth_sf_sub <- nga_sf %>%
  filter(nga$NAME_1 %in% states)

#nganorth <- nga[nga$NAME_1 %in% states,]

# Intersection of hexagons with geographic area
intersect_nganorth <- as_tibble(st_intersection(grids_nga_sf, nganorth_sf_sub)) %>%
  st_sf() %>%
  st_set_geometry(NULL)

intersectgid_nganorth <- unique(intersect_nganorth$gid)

# b) nganorthsmall
states_control <-
  c(
    "Geidam",
    "Magumeri",
    "Kaga",
    "Jere",
    "Maidugur",
    "Damboa",
    "Chibok",
    "Biu",
    "Kwaya Kusar",
    "Bayo",
    "Shani",
    "Kukawa",
    "Monguno",
    "Marte",
    "Konduga",
    "Michika",
    "Abadam",
    "Mobbar",
    "Guzamala",
    "Gubio",
    "Nganzai",
    "Mafa",
    "Ngala",
    "Kala/Balge",
    "Dikwa",
    "Bama",
    "Gwoza",
    "Madagali",
    "Askira/U",
    "Gujba",
    "Gulani"
  )

nganorthsmall_sf_sub <- nga_sf %>%
  filter(nga$NAME_2 %in% states_control)

# nganorthsmall <- nga[nga$NAME_2 %in% states_control, ]

intersect_nganorthsmall <- as_tibble(st_intersection(grids_nga_sf, nganorthsmall_sf_sub)) %>%
  st_sf() %>%
  st_set_geometry(NULL)

intersectgid_nganorthsmall <- unique(intersect_nganorthsmall$gid)

## 1)  Base data frame
# Geo-coding events
merged_nga_sf <- merged_nga %>%
  st_as_sf(coords = c("longitude", "latitude"),
           crs = crs)
nga_cols <- names(merged_nga_sf)

# Filtering events within the states of the study region
merged_nga_sf_filter <- merged_nga_sf %>%
  st_intersection(nga_sf) %>%
  dplyr::select(one_of(nga_cols))

# Merging events and hexagon grid cells
new_nga <- merged_nga_sf_filter %>%
  st_join(grids_nga_sf) %>%
  dplyr::rename(longitude = long_old,
                latitude = lat_old) %>%
  filter(!is.na(eventid)) %>%
  dplyr::select(one_of(vars_select)) %>%
  mutate(nganorth = ifelse(gid %in% intersectgid_nganorth, 1, 0),
         nganorthsmall = ifelse(gid %in% intersectgid_nganorthsmall, 1, 0))

centroid_nga <- grids_nga_sf %>%
  st_centroid() %>%
  st_coordinates() %>%
  as_tibble() %>%
  dplyr::mutate(gid = paste0("ngahex", seq(1, nrow(.))))

grids_nga_sf_new <- grids_nga_sf %>%
  mutate(countrypoly = ifelse(gid %in% intersectgid_nga, 1, 0),
         nganorth = ifelse(gid %in% intersectgid_nganorth, 1, 0),
         nganorthsmall = ifelse(gid %in% intersectgid_nganorthsmall, 1, 0)) %>%
  left_join(centroid_nga)

## Saving data
saveRDS(grids_nga_sf_new, here::here("prepped", "grids_hex_nga_hex25.rds"))
saveRDS(new_nga, here::here("prepped", "events_nga_hex25.rds"))
