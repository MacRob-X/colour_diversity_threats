# Functions for examining individual colour metric - threat relationships
# Statistics and plotting


# Test for holistic differences between focal and non-focal distributions using e.g. Wasserstein's
# distance for a single metric axis and threat
# This function uses the inbuilt bootstrapping capabilities of twosamples to generate a p-value
test_metric_threat_twosamples <- function(
    metric_name, 
    threat_colour_long, 
    threat_type, 
    focal_sex = "both_sexes", 
    n_boots = 1000, 
    stat_test = c("cvm", "dts", "wass", "ks", "kuiper", "ad"), 
    threatened_spp_only = FALSE
){
  
  # display focal combination
  message(paste0("Working on ", focal_sex, ", ", threat_type, ", ", metric_name))
  
  # Set twosamples stats test to use
  stat_func <- switch(
    stat_test,
    "cvm" = twosamples::cvm_test,
    "dts" = twosamples::dts_test,
    "wass" = twosamples::wass_test,
    "ks"  = twosamples::ks_test,
    "kuiper" = twosamples::kuiper_test,
    "ad"  = twosamples::ad_test,
    stop("Error: Invalid test type specified. Choose 'cvm', 'wass', 'dts', 'ks', 'kuiper' or 'ad'.")
  )
  
  # filter to threatened species only
  # NOTE that this includes species which are threatened by 
  # non-significant drivers of extinction, as 
  # determined by Stewart et al 2025 Nat Ecol Evol (see Supplementary_Dataset.xlsx for details)
  if(threatened_spp_only == TRUE){
    analysis_dat <- threat_colour_long |> 
      filter(
        iucn_cat %in% c("CR", "EN", "VU"), 
        metric == metric_name
      ) 
  } else {
    # if not using threatened species only, just restrict to metric axis of interest
    analysis_dat <- threat_colour_long |> 
      filter(
        metric == metric_name
      ) 
  }
  
  
  # Filter to sex of interest if specific sex requested
  if(focal_sex != "both_sexes"){
    analysis_dat <- analysis_dat |> 
      filter(
        sex == focal_sex
      )
  }
  
  
  # Extract distribution of focal threat species and non-focal-threat species
  focal_distrib <- analysis_dat[analysis_dat[, threat_type] == 1, ]
  non_focal_distrib <- analysis_dat[analysis_dat[, threat_type] == 0, ] 
  
  # extract metric values
  focal_distrib <- focal_distrib[["metric_value"]]
  non_focal_distrib <- non_focal_distrib[["metric_value"]]
  
  
  # Calculate test statistic with count-based p-value
  test_stat <- stat_func(
    a = focal_distrib,
    b = non_focal_distrib,
    nboots = n_boots,
    keep.boots = F
  )
  test_stat <- as.data.frame(t(as.data.frame(test_stat)))
  test_stat$metric <- metric_name
  test_stat$sex <- focal_sex
  
  return(test_stat)
  
}

# Permutation test for mean shift and standard deviation difference
test_metric_meanshift_varshift <- function(
    combo_number,
    sig_combos_df,
    filtered_colour_data, # this must be filtered to the focal metric only before calling, otherwise it will time out
    n_boots = 1000, 
    threatened_spp_only = FALSE,
    focal_sex = "both_sexes"
){
  
  # define focal combination
  focal_combo <- sig_combos_df[combo_number, ]
  focal_metric <- focal_combo[["metric"]]
  focal_threat <- focal_combo[["extinction_driver"]]
  
  # determine if the dataset has been restricted to the focal metric only
  # if not, throw a warning as it will take ages to run
  metrics_pres <- unique(filtered_colour_data$metric)
  if(length(metrics_pres) > 1){
    error("Dataset not filtered to focal metric only. Filter to focal metric before running function.")
  } else if(metrics_pres != focal_metric){
    error("metric in focal combo does not match metric in data.")
  }
  
  # display focal combination
  # message(paste0("Working on ", focal_threat, ", ", focal_metric))
  
  # discard unused columns in data
  filtered_colour_data <- filtered_colour_data[, c("jetz_species", "sex", "iucn_cat", focal_threat, "metric", "metric_value")]
  
  # filter to threatened species only
  # NOTE that this includes species which are threatened by 
  # non-significant drivers of extinction, as 
  # determined by Stewart et al 2025 Nat Ecol Evol (see Supplementary_Dataset.xlsx for details)
  if(threatened_spp_only == TRUE){
    filtered_colour_data <- filtered_colour_data |> 
      filter(
        iucn_cat %in% c("CR", "EN", "VU")
      ) 
  }
  
  
  # Filter to sex of interest if specific sex requested
  if(focal_sex != "both_sexes"){
    filtered_colour_data <- filtered_colour_data |> 
      filter(
        sex == focal_sex
      )
  }
  
  
  # Extract distribution of focal threat species and non-focal-threat species
  focal_distrib <- filtered_colour_data[filtered_colour_data[, focal_threat] == 1, ]
  non_focal_distrib <- filtered_colour_data[filtered_colour_data[, focal_threat] == 0, ] 
  
  
  # Calculate test statistics (mean shift and log variance [F] ratio)
  # +ve mean shift means the mean metric value for the focal distribution is higher than for the non-focal
  obs_mean_shift <- mean(focal_distrib$metric_value) -  mean(non_focal_distrib$metric_value)
  # a log ratio >0 means the focal distribution has higher variance than the non-focal, <1 means
  # the focal distirbution has lower variance than the non-focal
  # If the focal distrib has higher variance than non-focal, the focal extinction driver tends
  # to act on phenotypes on the extremes of this metric (more or less, - it's a little more 
  # nuanced)
  obs_var_ratio <- log(var(focal_distrib$metric_value) / var(non_focal_distrib$metric_value))
  
  # randomly reassign focal/non-focal extinction driver to generate null distributions of metric value/
  # extinction driver combinations
  # Should this sample random male and female pairs rather than random rows?
  null_distribs <- lapply(1:n_boots, 
                          function(i){
                            # get unique species in full sample
                            unique_spp <- unique(filtered_colour_data$jetz_species)
                            # get number of species in focal distribution
                            sample_size <- length(focal_distrib$metric_value)
                            if(focal_sex == "both_sexes"){
                              sample_size <- sample_size / 2
                            }
                            
                            # sample from unique species
                            sample_spp <- sample(
                              unique_spp,
                              size = sample_size,
                              replace = FALSE
                            )
                            
                            sample_focal <- filtered_colour_data[filtered_colour_data$jetz_species %in% sample_spp, ]
                            sample_focal[, focal_threat] <- 1
                            sample_non_focal <- filtered_colour_data[!(filtered_colour_data$jetz_species %in% sample_spp), ]
                            sample_non_focal[, focal_threat] <- 0
                            
                            sample_distrib <- rbind(sample_focal, sample_non_focal)
                            
                            return(sample_distrib)
                          }
  )
  
  # Get null distribution of mean shifts
  null_mean_shifts <- unlist(
    lapply(
      null_distribs,
      function(sample){
        focal_sample <- sample[sample[[focal_threat]] == 1, "metric_value"]
        non_focal_sample <- sample[sample[[focal_threat]] == 0, "metric_value"]
        null_mean_shift <- mean(focal_sample$metric_value) - mean(non_focal_sample$metric_value)
        return(null_mean_shift)
      }
    )
  )
  
  # Calculate SES for mean shifts
  null_ms_mean <- mean(null_mean_shifts)
  null_ms_sd <- sd(null_mean_shifts)
  mean_shift_es <- obs_mean_shift - null_ms_mean
  mean_shift_ses <- mean_shift_es / null_ms_sd
  
  # Calculate p-value for t statistics
  p_val_ms <- (length(which(abs(null_mean_shifts) >= abs(obs_mean_shift))) + 1) / (n_boots + 1)
  
  # Get null distribution of log variance ratios
  null_var_ratios <- unlist(
    lapply(
      null_distribs,
      function(sample){
        focal_sample <- sample[sample[[focal_threat]] == 1, "metric_value"]
        non_focal_sample <- sample[sample[[focal_threat]] == 0, "metric_value"]
        null_var_ratio <- log(var(focal_sample$metric_value) / var(non_focal_sample$metric_value))
        return(null_var_ratio)
      }
    )
  )
  
  # Calculate SES for variance ratios
  null_vr_mean <- mean(null_var_ratios)
  null_vr_sd <- sd(null_var_ratios)
  var_ratio_es <- obs_var_ratio - null_vr_mean
  var_ratio_ses <- var_ratio_es / null_vr_sd
  
  # Calculate p-value for t statistics
  p_val_vr <- (length(which(abs(null_var_ratios) >= abs(obs_var_ratio))) + 1) / (n_boots + 1)
  
  # Calculate standard deviation of focal and non-focal combo (to double check variance
  # inequality directionality)
  sd_focal <- sd(focal_distrib$metric_value)
  sd_non_focal <- sd(non_focal_distrib$metric_value)
  delta_sd = sd_focal - sd_non_focal
  
  
  # define results dataframe
  res <- data.frame(
    metric = focal_metric, 
    ex_driver = focal_threat,
    sex = focal_sex,
    mean_shift_obs = obs_mean_shift, 
    mean_shift_null_mean = null_ms_mean, 
    mean_shift_null_sd = null_ms_sd, 
    mean_shift_es = mean_shift_es,
    mean_shift_ses = mean_shift_ses,
    mean_shift_p = p_val_ms,
    log_var_ratio_obs = obs_var_ratio, 
    log_var_ratio_null_mean = null_vr_mean, 
    log_var_ratio_sd = null_vr_sd, 
    log_var_ratio_es = var_ratio_es,
    log_var_ratio_ses = var_ratio_ses,
    log_var_ratio_p = p_val_vr,
    sd_focal = sd_focal,
    sd_non_focal = sd_non_focal,
    delta_sd = delta_sd)
  
  return(res)
}