## Spatial weights
# Computing distance 
# https://www.r-bloggers.com/great-circle-distance-calculations-in-r/
# Convert degrees to radians
deg2rad <- function(deg) return(deg*pi/180)

# Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# Haversine formula (hf)
gcd.hf <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}
func_weight <- function(dist, a, b){
  return(1/(1 + exp(-(a + b*dist))))
}

func_decay <- function(spatial_a, 
                       spatial_b, 
                       temporal_a,
                       temporal_b,
                       countrypoly = T){
  
  ## Spatial weights
  value_spatial_terr <- mclapply(1:nrow(centroids_new), function(c){
    lon <- as.numeric(centroids_new$lonr_c[c])
    lat <- as.numeric(centroids_new$latr_c[c])
    mclapply(1:nrow(events_terr), function(v){
      return(gcd.hf(lon,
                    lat,
                    as.numeric(events_terr$lonr_v[v]),
                    as.numeric(events_terr$latr_v[v])))
    })
  })
  
  spatial_mat_terr <- matrix(unlist(value_spatial_terr), nrow = nrow(centroids_new), byrow = T)
  spatial_mat_terr <- func_weight(spatial_mat_terr, spatial_a, spatial_b)
  
  
  value_spatial_comb <- mclapply(1:nrow(centroids_new), function(c){
    lon <- centroids_new$lonr_c[c]
    lat <- centroids_new$latr_c[c]
    mclapply(1:nrow(events_comb), function(v){
      return(gcd.hf(lon,
                    lat,
                    as.numeric(events_comb$lonr_v[v]),
                    as.numeric(events_comb$latr_v[v])))
    })
  })
  
  spatial_mat_comb <- matrix(unlist(value_spatial_comb), nrow = nrow(centroids_new), byrow = T)
  spatial_mat_comb <- func_weight(spatial_mat_comb, spatial_a, spatial_b)
  
  
  ## Temporal weights
  value_temporal_terr <- mclapply(1:nrow(df_months), function(m){
    mclapply(1:nrow(events_terr), function(v){
      if(events_terr$index[v] > df_months$index[m]){
        return(NA)
      } else {
        return(df_months$index[m]-events_terr$index[v])
      }
    })
  })
  
  temporal_mat_terr <- matrix(unlist(value_temporal_terr), nrow = nrow(df_months), byrow = T)
  temporal_mat_terr <- func_weight(temporal_mat_terr, temporal_a, temporal_b)
  
  value_temporal_comb <- mclapply(1:nrow(df_months), function(m){
    mclapply(1:nrow(events_comb), function(v){
      if(events_comb$index[v] > df_months$index[m]){
        return(NA)
      } else {
        return(df_months$index[m]-events_comb$index[v])
      }
    })
  })
  
  temporal_mat_comb <- matrix(unlist(value_temporal_comb), nrow = nrow(df_months), byrow = T)
  temporal_mat_comb <- func_weight(temporal_mat_comb, temporal_a, temporal_b)
  
  ## Weights per centroid
  # This could be vectorized
  # terrorism
  vals_terr <- list()
  for(c in 1:nrow(spatial_mat_terr)){
    row <- spatial_mat_terr[c,]
    vals <- c()
    for(t in 1:nrow(temporal_mat_terr)){
      vals[t] <- sum(temporal_mat_terr[t,]*row, na.rm = T)
    }
    vals_terr[[c]] <- vals #influence for each time period for each centroi
  }
  
  weights_terr <- matrix(unlist(vals_terr), nrow = nrow(spatial_mat_terr), byrow = T) %>%
    as.data.frame()
  
  df_weights_terr <- bind_cols(centroids_new, weights_terr) %>%
    st_set_geometry(NULL) %>%
    pivot_longer(
      names_to = "timeindex",
      values_to = "weighted_terrorism",
      cols = starts_with("V")) %>%
    dplyr::mutate(timeindex = as.numeric(str_replace_all(timeindex, "V", "")))

  
  # combats
  vals_comb <- list()
  for(c in 1:nrow(spatial_mat_comb)){
    row <- spatial_mat_comb[c,]
    vals <- c()
    for(t in 1:nrow(temporal_mat_comb)){
      vals[t] <- sum(temporal_mat_comb[t,]*row, na.rm = T)
    }
    vals_comb[[c]] <- vals #influence for each time period for each centroi
  }
  
  weights_comb <- matrix(unlist(vals_comb), nrow = nrow(spatial_mat_comb), byrow = T) %>%
    as.data.frame()
  
  df_weights_comb <- bind_cols(centroids_new, weights_comb) %>%
    st_set_geometry(NULL) %>%
    pivot_longer(
      names_to = "timeindex",
      values_to = "weighted_combats",
      cols = starts_with("V")) %>%
    dplyr::mutate(timeindex = as.numeric(str_replace_all(timeindex, "V", "")))
  
  
  ### merging it
  merged <- full_join(df_weights_terr, df_weights_comb) 
  
  if(countrypoly == T){
    merged <- merged %>%
      dplyr::mutate(weighted_terrorism = replace(weighted_terrorism, countrypoly == 0, 0),
                    weighted_combats = replace(weighted_combats, countrypoly == 0, 0))
  }
  
  return(merged)
  
}

