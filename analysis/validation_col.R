###################################
# Territorial Control paper
# HMM analysis
# Validation of Colombia estimates
#
# Therese Anders
###################################

library(tidyverse)
library(sf)
library(sandwich)
library(clusterSEs)
library(stargazer)

plotdir <- here::here("plots")
plotdirapp <- here::here("plots", "appendix")

############################
# Reading data
############################

## Deforestation
deforest_raw <- readRDS(here::here("prepped","deforestation_1213_1314_1415_1516.rds")) %>%
  dplyr::filter(colwest == 1) %>%
  st_as_sf()

## Grid data
grid_og <- readRDS(here::here("prepped","grids_hex_col_hex25.rds")) %>%
  dplyr::select(-X, -Y) %>%
  dplyr::filter(colwest == 1)
crs <- st_crs(grid_og)
gid_hmm <- unique(grid_og$gid)

hmm_grid <- readRDS(here::here("prepped","hmm_col.rds")) %>%
  dplyr::rename(base = control)

hmm_grid_noassass <- readRDS(here::here("prepped","hmm_col_noassess.rds")) %>%
  dplyr::select(gid,
                timeindex,
                noassass = control)

hmm_grid_noofficial <- readRDS(here::here("prepped","hmm_col_noofficial.rds")) %>%
  dplyr::select(gid,
                timeindex,
                noofficial = control)

hmm_all <- left_join(hmm_grid, hmm_grid_noassass) %>%
  left_join(hmm_grid_noofficial) %>%
  pivot_longer(cols = c(base, noassass, noofficial),
               names_to = "name",
               values_to = "control") %>%
  mutate(control_num = case_when(
    control == "R" ~ 0,
    control == "DR" ~ 0.25,
    control == "D" ~ 0.5,
    control == "DG" ~ 0.75,
    control == "G" ~ 1
  )) 

hmm_all_wide <- hmm_all %>%
  pivot_wider(values_from = c(control_num, control), names_from = name)

# merging
merged_deforest <- hmm_all %>%
  dplyr::group_by(gid, year, name) %>%
  dplyr::summarise(control_mean = mean(control_num),
                   control_median = median(control_num),
                   control_sd = sd(control_num)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(year <= 2016) %>%
  dplyr::arrange(gid, year) %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(control_mean_lag = dplyr::lag(control_mean),
         control_median_lag = dplyr::lag(control_median)) %>%
  dplyr::mutate(diff_control_mean = control_mean - control_mean_lag,
         diff_control_median = control_median - control_median_lag) %>%
  dplyr::mutate(diff_control_mean_lag = dplyr::lag(diff_control_mean)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(subset(deforest_raw, gid %in% gid_hmm)) %>%
  dplyr::select(-geometry, -X, -Y) %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(count_deforest_lag = dplyr::lag(count_deforest),
         count_unstable_lag = dplyr::lag(count_unstable),
         count_stable_lag = dplyr::lag(count_stable),
         prop_deforest_lag = dplyr::lag(prop_deforest),
         any_deforest = ifelse(count_deforest >= 1, 1, 0)) %>%
  dplyr::mutate(any_deforest_lag = dplyr::lag(any_deforest)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(year %in% seq(2013, 2016)) %>%
  dplyr::mutate(treat = ifelse(year == 2016, 1, 0))


#########################
# Summary statistics
#########################

tab_sub <- merged_deforest %>%
  filter(name == "base") %>%
  dplyr::select(control_mean,
                diff_control_mean,
                treat,
                any_deforest) %>%
  na.omit()

stargazer(as.data.frame(tab_sub),
          font.size = "small",
          digits = 4,
          covariate.labels = c("$\\text{Control}_{i,t}$",
                               "$\\Delta\\text{Control}_{i,t}$",
                               "$\\text{Peace}_{t}$",
                               "$\\text{Deforestation}_{i,t}$"),
          out = paste0(plotdirapp, "/", "tab_summary_deforest.tex"),
          title = "Summary statistics for the logistic regression model of deforestation in Colombia on changes in territorial control as a result of the 2016 peace agreement. The unit of analysis for territorial control is annual averages of monthly-level estimates for 0.25 degree hexagonal grid cells.")


##########################
# Full Model on difference
##########################

## base model
merged_deforest_base <- merged_deforest %>%
  filter(name == "base")

m0_base <- glm(any_deforest ~ diff_control_mean + treat,
               data = merged_deforest_base,
               family = binomial(link='logit')) 
summary(m0_base)

m1_base <- glm(any_deforest ~ diff_control_mean + treat + diff_control_mean*treat,
               data = merged_deforest_base,
               family = binomial(link='logit'))
summary(m1_base) 

m2_base <- glm(any_deforest ~ diff_control_mean + treat + diff_control_mean*treat + any_deforest_lag,
               data = merged_deforest_base,
               family = binomial(link='logit'))
summary(m2_base)

## noassass model
merged_deforest_noassass <- merged_deforest %>%
  filter(name == "noassass")

m0_noassass <- glm(any_deforest ~ diff_control_mean + treat,
                   data = merged_deforest_noassass,
                   family = binomial(link='logit')) 
summary(m0_noassass)

m1_noassass <- glm(any_deforest ~ diff_control_mean + treat + diff_control_mean*treat,
                   data = merged_deforest_noassass,
                   family = binomial(link='logit')) 
summary(m1_noassass) 

m2_noassass <- glm(any_deforest ~ diff_control_mean + treat + diff_control_mean*treat + any_deforest_lag,
                   data = merged_deforest_noassass,
                   family = binomial(link='logit')) 
summary(m2_noassass)


## noofficial model
merged_deforest_noofficial <- merged_deforest %>%
  filter(name == "noofficial")

m0_noofficial <- glm(any_deforest ~ diff_control_mean + treat,
                     data = merged_deforest_noofficial,
                     family = binomial(link='logit')) 
summary(m0_noofficial)

m1_noofficial <- glm(any_deforest ~ diff_control_mean + treat + diff_control_mean*treat,
                     data = merged_deforest_noofficial,
                     family = binomial(link='logit')) 
summary(m1_noofficial)

m2_noofficial <- glm(any_deforest ~ diff_control_mean + treat + diff_control_mean*treat + any_deforest_lag,
                     data = merged_deforest_noofficial,
                     family = binomial(link='logit'))
summary(m2_noofficial)


## Cluster robust SEs
# Trying it manually
# https://github.com/cran/clusterSEs/blob/master/R/clusterBS.glm.R
cluster.bs.glm<-function(mod, dat, cluster, ci.level = 0.95, boot.reps = 1000, stratify = FALSE, 
                         cluster.se = TRUE, report = TRUE, prog.bar = TRUE, output.replicates = FALSE){
  
  form <- mod$formula                                                 # what is the formula of this model?  
  variables <- all.vars(form)                                         # what variables are in this model?
  clust.name <- all.vars(cluster)                                     # what is the name of the clustering variable?
  used.idx <- which(rownames(dat) %in% rownames(mod$model))           # what were the actively used observations in the model?
  dat <- dat[used.idx,]                                               # keep only active observations (drop the missing)
  clust <- as.vector(unlist(dat[[clust.name]]))                       # store cluster index in convenient vector
  G<-length(unique(clust))                                            # how many clusters are in this model?
  ind.variables.full <- names(coefficients(mod))                      # what independent variables are in this model?
  ind.variables <- rownames(summary(mod)$coefficients)                # what non-dropped independent variables in this model?
  
  
  # load in a function to create clustered standard errors
  # by Mahmood Arai: http://thetarzan.wordpress.com/2011/06/11/clustered-standard-errors-in-r/
  cl   <- function(dat, fm, cluster){
    #require(sandwich, quietly = TRUE)
    #require(lmtest, quietly = TRUE)
    M <- length(unique(cluster))
    N <- length(cluster)
    K <- fm$rank
    dfc <- (M/(M-1))
    uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
    vcovCL <- dfc*sandwich(fm, meat.=crossprod(uj)/N)
    coeftest(fm, vcovCL) }
  
  if(cluster.se == T){
    
    se.clust <- cl(dat, mod, clust)[ind.variables,2]               # retrieve the clustered SEs
    beta.mod <- coefficients(mod)[ind.variables]                   # retrieve the estimated coefficients
    w <- beta.mod / se.clust                                       # calculate the t-test statistic
    
  }else{
    
    se.beta <- summary(mod)$coefficients[ind.variables,2]          # retrieve the vanilla SEs
    beta.mod <- coefficients(mod)[ind.variables]                   # retrieve the estimated coefficients
    w <- beta.mod / se.beta                                        # calculate the t-test statistic
    
  }
  
  w.store <- matrix(data=NA, nrow=boot.reps, ncol=length(ind.variables))  # store bootstrapped test statistics
  
  # keep track of the beta bootstrap replicates for possible output
  rep.store <- matrix(data=NA, nrow=boot.reps, ncol=length(beta.mod))
  colnames(rep.store) <- ind.variables
  
  if(prog.bar==TRUE){pb <- txtProgressBar(min = 0, max = boot.reps, initial = 0, style = 3)}
  for(i in 1:boot.reps){
    
    if(prog.bar==TRUE){setTxtProgressBar(pb, value=i)}
    
    boot.sel <- sample(1:G, size=G, replace=T)                            # randomly select clusters
    
    # pick the observations corresponding to the randomly selected clusters
    boot.ind <- c()                                                       # where the selected obs will be stored
    boot.clust <- c()                                                     # create new cluster index for the bootstrap data
    
    for(k in 1:G){
      
      obs.sel <- which(clust == unique(clust)[boot.sel[k]])               # which observations are in the sampled cluster?
      if(stratify==T){
        
        obs.samp <- sample(obs.sel, size = length(obs.sel), replace=T)    # sample randomly from the selected cluster
        boot.ind <- c(boot.ind, obs.samp)                                 # append the selected obs index to existing index
        
      }else{
        
        boot.ind <- c(boot.ind, obs.sel)                                  # append the selected obs index to existing index
        
      }
      boot.clust <- c(boot.clust, rep(k, length(obs.sel)))                # store the new bootstrap cluster index
      
    }
    
    boot.dat <- dat[boot.ind,]                                            # create the bootstrapped data
    
    # run a model on the bootstrap replicate data
    boot.mod <- suppressWarnings(tryCatch(glm(form, data = boot.dat, family = mod$family), 
                                          error = function(e){return(NULL)}))                                    
    
    if(is.null(boot.mod) == FALSE ){
      if(boot.mod$converged == 0){boot.mod <- NULL}                    # judge GLM as failure if convergence not achieved
    }
    fail <- is.null(boot.mod)                                          # determine whether the GLM process created an error
    
    if(fail==0){                                                     # proceed if the GLM model was not in error
      
      if(cluster.se == T){
        
        se.boot <- tryCatch(cl(boot.dat, boot.mod, boot.clust)[ind.variables,2],
                            error = function(e){return(NA)}, 
                            warning = function(w){return(NA)})                              # retrieve the bootstrap clustered SE
        beta.boot <- tryCatch(coefficients(boot.mod)[ind.variables],
                              error = function(e){return(NA)}, 
                              warning = function(w){return(NA)})                            # store the bootstrap beta coefficient
        w.store[i,] <- (beta.boot-beta.mod) / se.boot                              # store the bootstrap test statistic
        
        rep.store[i,] <- beta.boot                                                 # store the bootstrap beta for output
        
      }else{
        
        se.boot <- tryCatch(summary(boot.mod)$coefficients[ind.variables,2],
                            error = function(e){return(NA)}, 
                            warning = function(w){return(NA)})                               # retrieve the bootstrap vanilla SE
        beta.boot <- tryCatch(coefficients(boot.mod)[ind.variables],
                              error = function(e){return(NA)}, 
                              warning = function(w){return(NA)})                             # retrieve the bootstrap beta coefficient
        w.store[i,] <- (beta.boot-beta.mod) / se.boot                               # calculate the t-test statistic
        
        rep.store[i,] <- beta.boot                                                  # store the bootstrap beta for output
        
      }
      
    }else{
      w.store[i,] <- NA                                                  # if model didn't converge, store NA as a result 
      rep.store[i,] <- NA
    }
    
  }
  if(prog.bar==TRUE){close(pb)}
  
  num.fail <- length(attr(na.omit(w.store), "na.action"))         # count the number of times something went wrong
  w.store <- na.omit(w.store)                                     # drop the erroneous bootstrap replicates
  
  
  comp.fun<-function(vec2, vec1){as.numeric(vec1>vec2)}                              # a simple function comparing v1 to v2
  p.store.s <- t(apply(X = abs(w.store), FUN=comp.fun, MARGIN = 1, vec1 = abs(w)))   # compare the BS test stats to orig. result
  p.store <- 1 - ( colSums(p.store.s) / dim(w.store)[1] )                            # calculate the cluster bootstrap p-value
  
  # compute critical t-statistics for CIs
  crit.t <- apply(X=abs(w.store), MARGIN=2, FUN=quantile, probs=ci.level )
  if(cluster.se == TRUE){
    ci.lo <- beta.mod - crit.t*se.clust
    ci.hi <- beta.mod + crit.t*se.clust
  }else{
    ci.lo <- beta.mod - crit.t*se.beta
    ci.hi <- beta.mod + crit.t*se.beta
  }
  
  out.se <- cbind(ind.variables, se.clust)
  
  print.ci <- cbind(ind.variables, ci.lo, ci.hi)
  print.ci <- rbind(c("variable name", "CI lower", "CI higher"), print.ci)
  
  out.ci <- cbind(ci.lo, ci.hi)
  rownames(out.ci) <- ind.variables
  colnames(out.ci) <- c("CI lower", "CI higher")
  
  out <- matrix(p.store, ncol=1)
  colnames(out) <- c("clustered bootstrap p-value")
  rownames(out) <- ind.variables
  out.p <- cbind(ind.variables, out)
  out.p <- rbind(c("variable name", "cluster bootstrap p-value"), out.p)
  
  
  printmat <- function(m){
    write.table(format(m, justify="right"), row.names=F, col.names=F, quote=F, sep = "   ")
  }
  
  if(report==T){
    
    if(num.fail!=0){
      cat("\n", "\n", "\n", "****", "Warning: ", num.fail, " out of ", boot.reps, " bootstrap replicate models failed to estimate.", "****", "\n", sep="")
    }
    
    cat("\n", "Cluster Bootstrap p-values: ", "\n", "\n")
    printmat(out.p)
    
    cat("\n", "Confidence Intervals (derived from bootstrapped t-statistics): ", "\n", "\n")
    printmat(print.ci)
    
    cat("\n", "Cluster standard errors (derived from bootstrapped t-statistics): ", "\n", "\n")
    printmat(out.se)
    
    if(length(ind.variables) < length(ind.variables.full)){
      cat("\n", "\n", "****", "Note: ", length(ind.variables.full) - length(ind.variables), " variables were unidentified in the model and are not reported.", "****", "\n", sep="")
      cat("Variables not reported:", "\n", sep="")
      cat(ind.variables.full[!ind.variables.full %in% ind.variables], sep=", ")
      cat("\n", "\n")
    }
    
  }
  
  out.list<-list()
  out.list[["p.values"]]<-out
  out.list[["ci"]] <- out.ci
  out.list[["se"]] <- out.se
  if(output.replicates == TRUE){out.list[["replicates"]] <- rep.store}
  return(invisible(out.list))
  
}

## base
clust.bs.p0_base <- cluster.bs.glm(m0_base, merged_deforest_base, ~gid, ci.level = 0.95, boot.reps = 1000, stratify = FALSE, 
                                   cluster.se = TRUE, report = TRUE, prog.bar = TRUE, output.replicates = FALSE)
clust.bs.p1_base <- cluster.bs.glm(m1_base, merged_deforest_base, ~gid, ci.level = 0.95, boot.reps = 1000, stratify = FALSE, 
                                   cluster.se = TRUE, report = TRUE, prog.bar = TRUE, output.replicates = FALSE)
clust.bs.p2_base <- cluster.bs.glm(m2_base, merged_deforest_base, ~gid, ci.level = 0.95, boot.reps = 1000, stratify = FALSE, 
                                   cluster.se = TRUE, report = TRUE, prog.bar = TRUE, output.replicates = FALSE)

se_clust.bs.p0_base <- as.numeric(clust.bs.p0_base$se[,2])
se_clust.bs.p1_base <- as.numeric(clust.bs.p1_base$se[,2])
se_clust.bs.p2_base <- as.numeric(clust.bs.p2_base$se[,2])

## noassass
clust.bs.p0_noassass <- cluster.bs.glm(m0_noassass, merged_deforest_noassass, ~gid, ci.level = 0.95, boot.reps = 1000, stratify = FALSE, 
                                       cluster.se = TRUE, report = TRUE, prog.bar = TRUE, output.replicates = FALSE)
clust.bs.p1_noassass <- cluster.bs.glm(m1_noassass, merged_deforest_noassass, ~gid, ci.level = 0.95, boot.reps = 1000, stratify = FALSE, 
                                       cluster.se = TRUE, report = TRUE, prog.bar = TRUE, output.replicates = FALSE)
clust.bs.p2_noassass <- cluster.bs.glm(m2_noassass, merged_deforest_noassass, ~gid, ci.level = 0.95, boot.reps = 1000, stratify = FALSE, 
                                       cluster.se = TRUE, report = TRUE, prog.bar = TRUE, output.replicates = FALSE)

se_clust.bs.p0_noassass <- as.numeric(clust.bs.p0_noassass$se[,2])
se_clust.bs.p1_noassass <- as.numeric(clust.bs.p1_noassass$se[,2])
se_clust.bs.p2_noassass <- as.numeric(clust.bs.p2_noassass$se[,2])

## noofficial
clust.bs.p0_noofficial <- cluster.bs.glm(m0_noofficial, merged_deforest_noofficial, ~gid, ci.level = 0.95, boot.reps = 1000, stratify = FALSE, 
                                         cluster.se = TRUE, report = TRUE, prog.bar = TRUE, output.replicates = FALSE)
clust.bs.p1_noofficial <- cluster.bs.glm(m1_noofficial, merged_deforest_noofficial, ~gid, ci.level = 0.95, boot.reps = 1000, stratify = FALSE, 
                                         cluster.se = TRUE, report = TRUE, prog.bar = TRUE, output.replicates = FALSE)
clust.bs.p2_noofficial <- cluster.bs.glm(m2_noofficial, merged_deforest_noofficial, ~gid, ci.level = 0.95, boot.reps = 1000, stratify = FALSE, 
                                         cluster.se = TRUE, report = TRUE, prog.bar = TRUE, output.replicates = FALSE)

se_clust.bs.p0_noofficial <- as.numeric(clust.bs.p0_noofficial$se[,2])
se_clust.bs.p1_noofficial <- as.numeric(clust.bs.p1_noofficial$se[,2])
se_clust.bs.p2_noofficial <- as.numeric(clust.bs.p2_noofficial$se[,2])


library(stargazer)
stargazer(m0_base, m1_base, m2_base,
          se=list(se_clust.bs.p0_base,
                  se_clust.bs.p1_base,
                  se_clust.bs.p2_base),
          covariate.labels = c("$\\Delta\\text{Control}_{i,t} \\times \\text{Peace}_{t}$",
                               "$\\Delta\\text{Control}_{i,t}$",
                               "$\\text{Peace}_{t}$",
                               "$\\text{Deforestation}_{i,t-1}$"),
          order = c(4, 1, 2, 3),
          font.size = "small",
          digits = 2,
          omit.stat = c("f","ser","rsq"),
          column.sep.width = "-6pt",
          dep.var.labels.include = F,
          star.cutoffs = c(0.05, 0.01, 0.001),
          no.space = TRUE,
          model.names = FALSE,
          dep.var.caption  = "$\\text{Deforestation}_{i,t}$",
          title = "Relationship between rebel territorial control and deforestation in Colombia.",
          notes = c("Logistic regression coefficients with", 
                    "bootstrapped clustered standard",
                    "errors by grid cell in parentheses."),
          label = "tab:deforest",
          out = paste0(plotdir, "/",  "tab_deforest.tex"))


### Robustness check
stargazer(m0_base, m1_base, m2_base,
          m0_noassass, m1_noassass, m2_noassass,
          m0_noofficial, m1_noofficial, m2_noofficial,
          se = list(se_clust.bs.p0_base, se_clust.bs.p1_base, se_clust.bs.p2_base,
                    se_clust.bs.p0_noassass, se_clust.bs.p1_noassass, se_clust.bs.p2_noassass,
                    se_clust.bs.p0_noofficial, se_clust.bs.p1_noofficial, se_clust.bs.p2_noofficial),
          covariate.labels = c("$\\Delta\\text{Control}_{i,t} \\times \\text{Peace}_{t}$",
                               "$\\Delta\\text{Control}_{i,t}$",
                               "$\\text{Peace}_{t}$",
                               "$\\text{Deforestation}_{i,t-1}$"),
          column.labels = c("Base sample", "Exclude assassination", "Exclude government targets"),
          column.separate = c(3, 3, 3),
          order = c(4, 1, 2, 3),
          font.size = "footnotesize",
          digits = 2,
          omit.stat = c("f","ser","rsq"),
          column.sep.width = "-6pt",
          dep.var.labels.include = F,
          star.cutoffs = c(0.05, 0.01, 0.001),
          no.space = TRUE,
          model.names = FALSE,
          dep.var.caption  = "$\\text{Deforestation}_{i,t}$",
          title = "Relationship between rebel territorial control and deforestation in Colombia.",
          notes = c("Logistic regression coefficients with bootstrapped clustered standard errors by grid cell", 
                    " in parentheses."),
          label = "tab:deforest_robust",
          out = paste0(plotdirapp, "/",  "tab_deforest_robustness.tex"))



