####### Evidence matrix setup ########

func_ev <- function(t, c, mar, xs, cdf_T = 0, cdf_C = 0){
  
    # Preprocessing - Tolerance
    if(t <= xs){
      t <- 0
    }
    
    if(c <= xs){
      c <- 0
    }

    if((t == 0) & (c == 0)){
      return("O1")
    } else if (abs(cdf_T - cdf_C) <= mar){
      return("O3")
    } else if (cdf_C < cdf_T){
      return("O4")
    } else {
      return("O2")
  }
}

## Function to select right evidence node probabilities
ev_select <- function(input){
  output <- as.vector(ev_possibilities[row.names(ev_possibilities) == input,])
  return(output)
}



###############################################################################
# Functions for Pre-Processing the data and running the HMM for each grid cell
###############################################################################
func_preprocess_hmm <- function(df, 
                                countrycode, 
                                numstates, 
                                marval, 
                                xsval,
                                evidence_matrix,
                                transition_matrix){
  
  # Selecting the country by gwno code
  df_gwno <- df[df$gwno == countrycode,]
  
  # Selecting the number of unique grid cells
  gid <- c(unique(df_gwno$gid))
  
  # Setting up data frame
  frames <- list()
  
  # Looping through each grid cell
  for(i in 1:length(gid)){
    
    df_gwno_gid <- df_gwno[df_gwno$gid == gid[i],]
    
    # Selecting the number of time steps
    time <- sort(unique(df_gwno_gid$timeindex))
    
    ev_obs <- c()
    
    for(y in 1:length(time)){
      
      terr = df_gwno_gid$gtd[df_gwno_gid$timeindex == time[y]]
      conf = df_gwno_gid$ged[df_gwno_gid$timeindex == time[y]]
      
      cdfC <- df_gwno_gid$cdf_C[df_gwno_gid$timeindex == time[y]]
      cdfT <- df_gwno_gid$cdf_T[df_gwno_gid$timeindex == time[y]]
      
      ev_val <- func_ev(t = terr, c = conf, mar = marval, xs = xsval, cdf_T = cdfT, cdf_C = cdfC)
      ev_obs[y] <- ev_val
      
    } # Close loop populating evidence matrix for specific gid
    
    ## Implementing the HMM
    
    hmm <- initHMM(States = states, # List of states
                   Symbols = ev, # List of potential observations
                   startProbs = d_init, # Initial probabilities
                   transProbs = transition_matrix, # Transition matrix
                   emissionProbs = evidence_matrix) # Emission matrix
    
    
    
    # Compute the viterbi algorithm (Decoding)
    out_viterbi <- viterbi(hmm, ev_obs) #most probable path of states for a sequence of observations
    
    gid_hmm <- data.frame("gid" = as.numeric(gid[i]),
                          "timeindex" = as.numeric(time),
                          "control" = as.character(out_viterbi))
    
    frames[[i]] <- gid_hmm
    
  } # Close grid cell loop
  
  return(frames)
  
}



# Function for testing the sensitivity of HMM to variations in m (demo)
func_test <- function(cdf_T, cdf_C, mar, t = NULL, c = NULL, xs = NULL) {
  
  if (!is.null(t) & !is.null(c) & !is.null(xs)) {
    if(t <= xs){
      t <- 0
    }
    
    if(c <= xs){
      c <- 0
    }
  } else {
    t <- Inf
    c <- Inf
    xs <- Inf
  }
  
  if((t == 0) & (c == 0)){
    abs <- "O1"
  } else if (abs(cdf_T - cdf_C) <= mar){
    abs <- "O3"
  } else if (cdf_C < cdf_T){
    abs <- "O4"
  } else {
    abs <- "O2"
  }
  
  return(tibble(t = t,
                c = c,
                xs = xs,
                cdf_C = cdf_C,
                cdf_T = cdf_T,
                abs = abs,
                mar = mar))
  
}


# Function to test the sensitivity of the results to variations in emission matrix
func_ev_mat <- function(r){
  
  ev <- c("O1", "O2", "O3", "O4")
  states <- c("R", "DR", "D", "DG", "G")
  
  order <- sample(c(2,3), 2, replace = F)
  
  ev_mat_draw <- matrix(nrow = length(states), ncol = length(ev))
  ev_mat_draw[1,] <- c(r[[1]], r[[order[[1]]]], r[[order[[2]]]], r[[4]])
  ev_mat_draw[2,] <- c(r[[4]], r[[1]], r[[order[[1]]]], r[[order[[2]]]])
  ev_mat_draw[3,] <- c(r[[4]], r[[order[[1]]]], r[[1]], r[[order[[2]]]])
  ev_mat_draw[4,] <- c(r[[4]], r[[order[[1]]]], r[[order[[2]]]], r[[1]])
  ev_mat_draw[5,] <- c(r[[1]], r[[4]], r[[order[[1]]]], r[[order[[2]]]])
  
  return(ev_mat_draw)
}