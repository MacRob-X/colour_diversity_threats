# Functions for examining individual PC axis - threat relationships
# Statistics and plotting

# function to get proportional density plot of species threatened by a particular threat
# proportional to density of all species
prop_dens_2d <- function(threat_col_data, focal_threat, x_axis, y_axis, threatened_spp_only = FALSE, n_bins = 200, return_params = FALSE){
  
  # make copy of full data, to extract axes limits later
  full_dat <- threat_col_data
  
  # If making proportional to all other THREATENED species, get rows of all  
  # species for which we have threat data and which are threatened (i.e. not LC or NT)
  if(threatened_spp_only == TRUE){
    threat_col_data <- threat_col_data[which(threat_col_data$iucn_cat %in% c("CR", "EN", "VU")), ]
  }
  
  # add row ID column for subsetting
  threat_col_data$row_id <- 1:nrow(threat_col_data)
  
  # Get rows corresponding to species threatened by focal threat
  focal_threat_rows <- threat_col_data[which(threat_col_data$ex_driver == focal_threat), ][["row_id"]]
  # Get rows corresponding to species which are NOT threatened by the focal threat (i.e.,
  # species which are not in the focal threat rows)
  non_focal_threat_rows <- threat_col_data[which(!(threat_col_data$jetz_species %in% threat_col_data[focal_threat_rows, ]$jetz_species)), ][["row_id"]]
  # need to remove duplicates - make sure each species is represented only once per sex, rather
  # than once per sex per extinction driver
  non_focal_dat <- threat_col_data[threat_col_data$row_id %in% non_focal_threat_rows, ]
  non_focal_dat <- dplyr::distinct(non_focal_dat, jetz_species, sex, .keep_all = TRUE)
  non_focal_threat_rows <- non_focal_dat$row_id
  
  # Get x and y limits - use the entire range of the x and y axes
  x_lim <- range(full_dat[[x_axis]])
  y_lim <- range(full_dat[[y_axis]])
  
  # Get 2D density of non-focal rows
  non_focal_2dkd <- MASS::kde2d(
    x = threat_col_data[non_focal_threat_rows, x_axis], 
    y = threat_col_data[non_focal_threat_rows, y_axis],
    n = n_bins,
    lims = c(x_lim, y_lim)
  )
  # Get 2D density of focal rows
  focal_2dkd <- MASS::kde2d(
    x = threat_col_data[focal_threat_rows, x_axis], 
    y = threat_col_data[focal_threat_rows, y_axis],
    n = n_bins,
    lims = c(x_lim, y_lim)
  )
  # calculate proportional 2D density (i.e., 2D density of focal species - 2D density of 
  # non-focal species)
  prop_2dkd <- list(
    x = focal_2dkd$x,
    y = focal_2dkd$y,
    z = focal_2dkd$z - non_focal_2dkd$z
  )
  
  if(return_params == FALSE){
    return(prop_2dkd)
  } else {
    params = list(
      focal_threat = focal_threat,
      x_axis = x_axis,
      y_axis = y_axis,
      threatened_spp_only = threatened_spp_only,
      n_bins = n_bins
    )
    res <- list(
      prop_2dkd = prop_2dkd,
      params = params
    )
    return(res)
  }
  
  
}

# Plot a single proportional density plot
plot_prop_dens_2d <- function(prop_dens_obj, focal_threat = NULL, x_axis = NULL, y_axis = NULL, neg_colour = "darkblue", pos_colour = "darkred", plot_contours = FALSE, title = FALSE, max_density_val = NULL){
  
  # Check if only a MASS::kde2d object is provided - if so, need to either provide 
  # focal threat and axes data or set to defaults
  if(length(prop_dens_obj) == 3){
    
    # Check if focal threat name and axes provided
    if(is.null(focal_threat)){
      warning("Focal threat not provided. Labelling as 'Focal threat'.")
      focal_threat <- "not_provided"
    }
    if(is.null(x_axis)){
      warning("X-axis not provided. Labelling as 'X-axis'.")
      x_axis <- "not_provided"
    }
    if(is.null(y_axis)){
      warning("Y-axis not provided. Labelling as 'Y-axis'.")
      y_axis <- "not_provided"
    }
    
    prop_2dkd <- prop_dens_obj
    
  } else if(length(prop_dens_obj) == 2){
    
    focal_threat <- prop_dens_obj$params$focal_threat
    x_axis <- prop_dens_obj$params$x_axis
    y_axis <- prop_dens_obj$params$y_axis
    
    prop_2dkd <- prop_dens_obj$prop_2dkd
    
  }
  
  # Set up colour palette
  # Check if max density value provided (if plotting multiple plots with the same 
  # colour scale)
  if(is.null(max_density_val)){
    
    # get min and max density values from proportional density
    den_min <- min(prop_2dkd$z)
    den_max <- max(prop_2dkd$z)
    
    col_max <- max(abs(den_min), abs(den_max))
  } else {
    col_max = max_density_val
  }
  
  # create symmetrical colour palette
  n_cols <- 100
  col_min <- -(col_max)
  col_breaks <- seq(col_min, col_max, length.out = n_cols + 1)
  col_pal_neg <- colorRampPalette(c(neg_colour, "white"))(n_cols / 2)
  col_pal_pos <- colorRampPalette(c("white", pos_colour))(n_cols / 2)
  col_pal <- c(col_pal_neg, col_pal_pos)
  
  # Set up axes labels
  if(x_axis == "not_provided"){
    x_lab <- "X-axis"
  } else {
    x_lab <- x_axis
  }
  if(y_axis == "not_provided"){
    y_lab <- "Y-axis"
  } else {
    y_lab <- y_axis
  }
  
  # Plot as heatmap
  image(
    prop_2dkd,
    col = col_pal,
    breaks = col_breaks,
    xlab = x_lab, ylab = y_lab, 
  )
  
  if(title == TRUE){
    
    # Set plot title
    if(focal_threat == "clim_chan"){
      plot_title <- "Climate change"
    } else if(focal_threat == "hab_loss"){
      plot_title <- "Habitat loss"
    } else if(focal_threat == "pollut"){
      plot_title <- "Pollution"
    } else if(focal_threat == "hunt_col"){
      plot_title <- "Hunting & collection"
    } else if(focal_threat == "no_sig_driver"){
      plot_title <- "No significant extinction drivers"
    } else if(focal_threat == "acc_mort"){
      plot_title <- "Accidental mortality"
    } else if(focal_threat == "invas_spec"){
      plot_title <- "Invasive species"
    } else if(focal_threat == "other"){
      plot_title <- "Other threats"
    }
    
    # Add plot title
    title(plot_title)
    
  }
  
  if(plot_contours == TRUE){
    
    contour(
      prop_2dkd,
      add = TRUE,
      col = "darkgrey"
    )
    
  }
  
  
  
  
}

# Calculate SES or test statistic of twosamples distributional difference test
# (e.g. Wasserstein's distance) for a single PC axis and threat
# OBSOLETE - NOT RUN
test_pc_threat <- function(pc_axis, threat_colour_long, threat, n = 1000, seed = NULL, return_distrib = TRUE, stat_test = c("cvm", "dts", "wass", "ks", "kuiper", "ad")){
  
  # Set twosamples stats test to use
  stat_func <- switch(
    stat_test,
    "cvm" = twosamples::cvm_stat,
    "dts" = twosamples::dts_stat,
    "wass" = twosamples::wass_stat,
    "ks"  = twosamples::ks_stat,
    "kuiper" = twosamples::kuiper_stat,
    "ad"  = twosamples::ad_stat,
    stop("Error: Invalid test type specified. Choose 'cvm', 'wass', 'dts', 'ks', 'kuiper' or 'ad'.")
  )
  
  # Extract distribution of focal threat species and non-focal-threat species
  focal_distrib <- threat_colour_long |> 
    filter(
      iucn_cat != "LC", 
      PC == pc_axis, 
      ex_driver == threat
    ) |> 
    pull(PC_value)
  non_focal_distrib <- threat_colour_long |> 
    filter(
      iucn_cat != "LC", 
      PC == pc_axis, 
      ex_driver != threat
    ) |> 
    pull(PC_value)
  # NOTE that this excludes species for which ex_driver is NA - this includes species
  # which truly have no threats (i.e. LC species) AND threatened species for which there is
  # no threat data
  # This is the correct thing to do because we cannot say for sure that these threatened species
  # are not threatened by the focal threat
  # It includes species which are threatened by non-significant drivers of extinction, as 
  # determined by Stewart et al 2025 Nat Ecol Evol (see Supplementary_Dataset.xlsx for details)
  
  # Extract distribution of all threatened species (focal and non-focal)
  # We will use this to generate our null distributions
  all_distrib <- threat_colour_long |> 
    filter(
      iucn_cat != "LC",
      PC == pc_axis
    ) |> 
    pull(PC_value)
  
  # Check that lengths of distributions match
  distrib_match <- function(focal_distrib, non_focal_distrib, all_distrib){
    assertthat::are_equal(
      length(focal_distrib) + length(non_focal_distrib), 
      length(all_distrib)
    )
  }
  assertthat::on_failure(distrib_match) <- function(call, env){
    "Lengths of non-focal and focal distributions do not equal length of all-threat distribution. \nCheck if species threats are correctly coded (e.g. that threatened species with no listed threats do not have NA ex_driver value)."
  }
  assertthat::assert_that(distrib_match(focal_distrib, non_focal_distrib, all_distrib))
  
  # get lengths of focal and non-focal distributions
  n_focal <- length(focal_distrib)
  n_non_focal <- length(non_focal_distrib)
  n_all <- n_focal + n_non_focal
  
  # calculate observed CVM value
  test_stat_obs <- stat_func(
    a = focal_distrib,
    b = non_focal_distrib,
  )
  
  # set seed (if specified)
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  
  # Generate set of null distribution from non-focal-
  null_test_stats <- unlist(
    lapply(
      1:n, 
      function(i){
        all_sample <- sample(
          all_distrib,
          size = n_all,
          replace = FALSE
        )
        focal_sample <- all_sample[1:n_focal]
        non_focal_sample <- all_sample[(n_focal + 1):n_all]
        
        # Calculate test statistic for each random distribution, vs the original non-focal-threat distribution
        null_test_stat <- stat_func(
          a = focal_sample,
          b = non_focal_sample,
        )
        return(null_test_stat)
      }
    )
  )
  
  # Check if distribution of null test statistics is normal
  norm_dist <- shapiro.test(null_test_stats)$p.value >= 0.05
  if(norm_dist == FALSE){
    
    warning("Null test statistic distribution is non-normal. Using median-based Standardized \nEffect Size.")
    
    # Use median-based standardized effect size if non-normal distribution
    # with median and median absolute deviation - robust to non-normal data
    # Calculate SES (obs - median(null)) / MAD(null)
    # SES > 1.96 indicates deviation from difference between distributions expected by chance
    random_median <- median(null_test_stats)
    random_mad <- mad(null_test_stats)
    es <- test_stat_obs - random_median
    ses <- es / random_mad
  } else {
    
    # Calculate SES (obs - mean(null)) / sd(null)
    # SES > 1.96 indicates deviation from difference between distributions expected by chance
    random_mean <- mean(null_test_stats)
    random_sd <- sd(null_test_stats)
    es <- test_stat_obs - random_mean
    ses <- es / random_sd
    
  }
  
  
  # calculate p-value with Laplace smoothing
  # method: rank/counting: proportion of null values equal to or larger than observed CVM statistic
  # This makes no assumptions about the shape of the null distribution, unlike standard
  # p-value calculations (which assume normal distribution). CVM stat distribution is likely
  # a chi-square distribution, not normal
  # This is the method used under the hood in twosamples::cvm_test
  p_val <- (length(which(null_test_stats >= test_stat_obs)) + 1) / (n + 1)
  
  # bind results together for output
  if(norm_dist == FALSE){
    res <- list(
      extinction_driver = threat,
      ES = es,
      SES = ses,
      test_statistic = stat_test,
      test_stat_observed = test_stat_obs,
      test_stat_null_median = random_median,
      test_stat_null_mad = random_mad,
      p_value = p_val
    )
  } else {
    
    res <- list(
      extinction_driver = threat,
      ES = es,
      SES = ses,
      test_statistic = stat_test,
      test_stat_observed = test_stat_obs,
      test_stat_null_mean = random_mean,
      test_stat_null_sd = random_sd,
      p_value = p_val
    )
    
  }
  
  # bind to distributions, if desired
  if(return_distrib == TRUE){
    res$null_test_stats <- null_test_stats
  }
  
  return(res)
  
}

# Null distribution and observed test value plotting function
plot_cvm_distrib <- function(cvm_res, bins = 30){
  hist(cvm_res$cvm_null, breaks = bins)
  abline(v = cvm_res$cvm_observed, col = "red")
}

# Test for holistic differences between focal and non-focal distributions using e.g. Wasserstein's
# distance for a single PC axis and threat
# This function uses the inbuilt bootstrapping capabilities of twosamples to generate a p-value
test_pc_threat_twosamples <- function(
    pc_axis, 
    threat_colour_long, 
    threat_type, 
    focal_sex = "both_sexes", 
    n_boots = 1000, 
    stat_test = c("cvm", "dts", "wass", "ks", "kuiper", "ad"), 
    threatened_spp_only = FALSE){
  
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
  # NOTE that this excludes species for which ex_driver is NA - this includes species
  # which truly have no threats (i.e. LC species) AND threatened species for which there is
  # no threat data
  # This is the correct thing to do because we cannot say for sure that these threatened species
  # are not threatened by the focal threat
  # It includes species which are threatened by non-significant drivers of extinction, as 
  # determined by Stewart et al 2025 Nat Ecol Evol (see Supplementary_Dataset.xlsx for details)
  if(threatened_spp_only == TRUE){
    analysis_dat <- threat_colour_long |> 
      filter(
        iucn_cat %in% c("CR", "EN", "VU"), 
        PC == pc_axis
      ) |> 
      filter(
        !is.na(ex_driver)
      )
  } else {
    # if not using threatened species only, still filter out threatened species for which there is
    # no threat data
    # Again, we cannot say for sure that these threatened species are not threatened by the 
    # focal threat
    analysis_dat <- threat_colour_long |> 
      filter(
        PC == pc_axis
      ) |> 
      filter(
        !is.na(ex_driver) | iucn_cat == "LC" # anything that is not listed as NT should have an associated threat since we know they're threatened - so these are rows with missing threat data 
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
  focal_rows <- analysis_dat |> 
    mutate(
      row_num = row_number()
    ) |> 
    filter(
      ex_driver == threat_type
    ) |> 
    pull(
      row_num
    )
  focal_distrib <- analysis_dat[focal_rows, ]
  non_focal_distrib <- analysis_dat[-focal_rows, ] |> 
    filter(
      !(jetz_species %in% unique(focal_distrib$jetz_species))
    )
  # get a single row per species/sex - otherwise species with more than one threat will be counted
  # multiple times
  # only necessary for non-focal species as focal species will by definition only have one (focal) 
  # threat
  non_focal_distrib <- non_focal_distrib |> 
    distinct(
      jetz_species, sex, .keep_all = TRUE
    )
  
  # extract pc values
  focal_distrib <- focal_distrib[["PC_value"]]
  non_focal_distrib <- non_focal_distrib[["PC_value"]]
  

  # Calculate test statistic with count-based p-value
  test_stat <- stat_func(
    a = focal_distrib,
    b = non_focal_distrib,
    nboots = n_boots,
    keep.boots = F
  )
  test_stat <- as.data.frame(t(as.data.frame(test_stat)))
  test_stat$PC <- pc_axis
  test_stat$sex <- focal_sex
  
  return(test_stat)
  
}
