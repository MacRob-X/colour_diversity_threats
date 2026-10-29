# Compare mean colour pattern diversity for species threatened by different IUCN threat groups
# Robert MacDonald

# Clear environment
rm(list=ls())

# Load libraries ----
library(dplyr)
library(ggplot2)
library(ggridges)

# Load custom functions ----
source(
  here::here(
    "02_scripts", "R", "06_pc_axes_threat_functions.R"
  )
)
source(
  here::here(
    "02_scripts", "R", "08_image_feature_metrics_threats_functions.R"
  )
)

## EDITABLE CODE ##
# Use latest IUCN assessment data or use most recent assessment data pre- specified cutoff year?
latest <- TRUE
# If not using latest assessment data, specify a cutoff year. Set to NULL if using latest.
cutoff_year <- NULL
# Clade to focus on ("Aves", "Neognaths", "Neoaves", "Passeriformes")
clade <- "Aves"
# Colour space to use
space <- "lab"
# Exclude past threats?
exclude_past_threats <- FALSE
# Remove extinct (EX) species?
remove_extinct <- TRUE

# Load data ----

# Load threat data with centroid distances (binary variables version)
if(remove_extinct){
  ext_par <- "remove_extinct_"
} else {
  ext_par <- "keep_extinct_"
}
if(exclude_past_threats) {
  past_threat_par <- "exclude_past_threats_"
} else {
  past_threat_par <- NULL
}
if(latest == TRUE){
  filename <- paste0(ext_par, past_threat_par, "centroid_jetz_threat_matrix_latest_2026-01-07.csv")
} else if(latest == FALSE){
  filename <- paste0(ext_par, past_threat_par, "centroid_jetz_threat_matrix_", cutoff_year, "_cutoff_year.csv")
}
threat_data <- read.csv(
  here::here(
    "03_output_data", "03a_threat_centroid",
    filename
  )
)


# load handcrafted colour metric data
# currently this requires running 07_recolorize_feature_metrics.R prior to current script!
metrics_filename <- if(remove_extinct){
  paste("handcrafted_colour_metrics", "remove_EX.csv", sep = "_")
} else if(!remove_extinct){
  paste("handcrafted_colour_metrics", "keep_EX.csv", sep = "_")
}
metrics <- read.csv(
  here::here(
    "03_output_data", "07_recolorize_image_feature_metrics",
    metrics_filename
  ),
  row.names = 1
)

# Data preparation ----


# Exclude Data Deficient (DD) species (already done at taxonomy matching stage,
# duplicating here for completeness)
threat_data <- threat_data %>%
  filter(
    !(iucn_cat %in% c("DD"))
  )



# get threat data with handcrafted colour metrics
# remove PCs from colour metrics before attaching to threat data (no need to duplicate)
metric_names <- c("mean_lum", "max_lum", "countershading", "colourfulness", paste0("p_", 1:15))
metrics <- metrics[, which(colnames(metrics) %in% metric_names)]
metrics$species <- sapply(strsplit(rownames(metrics), split = "-"), "[", 1)
metrics$sex <- sapply(strsplit(rownames(metrics), split = "-"), "[", 2)

# attach handcrafted colour metrics to threat data - we don't need PC values here so remove first
threat_metrics <- threat_data %>% 
  select(
    -starts_with("PC")
  ) %>% 
  inner_join(
    metrics,
    by = join_by("jetz_species" == "species", "sex")
  )

# pivot longer on metric types
threat_metrics_long <- threat_metrics %>% 
  tidyr::pivot_longer(
    cols = all_of(metric_names),
    names_to = "metric",
    values_to = "metric_value"
  ) 



# Handcrafted colour metric statistics ----

## Wasserstein distance permutation test

threat_types <- c(hab_loss = "hab_loss", hunt_col = "hunt_col", invas_spec = "invas_spec", clim_chan = "clim_chan", pollut = "pollut", acc_mort = "acc_mort")

sexes <- c("both_sexes")
threatened_spp_only_par <- FALSE
res_twosamples <- lapply(
  sexes,
  function(sex){
    
    sex_res <- lapply(
      
      threat_types,
      
      function(threat_type){
        
        metric_stats <- lapply(
          metric_names,
          function(metric){
            # filter to focal metric here to save overhead
            dat <- threat_metrics_long[threat_metrics_long$metric == metric, ]
            stat <- test_metric_threat_twosamples(metric, 
                                              threat_colour_long = dat, 
                                              threat_type = threat_type,
                                              focal_sex = sex,
                                              n_boots = 1000,
                                              stat_test = "wass",
                                              threatened_spp_only = threatened_spp_only_par)
            return(stat)
          }
        )
        
        threat_res <- do.call(rbind, metric_stats)
      }
      
    )
    sex_res <- data.table::rbindlist(sex_res, idcol = "extinction_driver")
  }
)
# format results
res_twosamples <- do.call(rbind, res_twosamples)
colnames(res_twosamples) <- c("extinction_driver", "test_stat", "p_value", "metric", "sex")
res_twosamples <- res_twosamples %>% 
  select(metric, extinction_driver, sex, test_stat, p_value)

# Bonferroni correction for p-values (within sex)
res_twosamples <- res_twosamples %>% 
  group_by(
    sex
  ) %>% 
  mutate(
    p_adjusted = p.adjust(p_value, method = "bonferroni")
  )
# NOTE that this is not necessary - I only want to identify which threat types might differentially
# threaten which PC axes - better to keep the net wide here as the mean shift/variance tests
# will add an additional filter
# It's also not robust, since there is a hard lower limit for p-values of 0.0005
# based on the laplace smoothing in the count-based p-value calculation

# Save results as CSV
# set filename based on parameters
if(threatened_spp_only_par == TRUE){
  wass_filename <- "metrics_ex_driver_wasserstein_distance_threatened_spp_only.csv"
} else if(threatened_spp_only_par == FALSE){
  wass_filename <- "metrics_ex_driver_wasserstein_distance_all_spp.csv"
}

write.csv(
  res_twosamples,
  file = here::here(
    "03_output_data", "08_image_feature_metrics_threats",
    wass_filename
  )
)


## Mean shift and variance shift permutation tests ----

# Keep only the axes for which p-value (from twosamples results) < 0.05
# Use this because SES is unreliable for skewed distributions of null statistics
sig_metric_exdrive <- res_twosamples %>% 
  filter(
    p_value < 0.05
  )
sig_combos <- sig_metric_exdrive %>% 
  select(
    metric, extinction_driver, sex
  ) %>% 
  mutate(
    combo = paste(metric, extinction_driver, sex, sep = "_")
  )

# Test for SPECIFIC distributional differences in each of these metric/threat combinations
# --> Mean shift (via difference in means, with permutation test)
# --> Shift in variance (via difference in standard deviation, via permutation test)

# 2026-10-07 New method for detecting mean shifts and variance differences
# test function
# first filter data to only the focal PC (essential to save overhead, otherwise function
# times out if you call it in an lapply because the dataset is so big)
focal_combo_num <- 9
focal_pc <- sig_combos[[focal_combo_num, "PC"]]
test_data <- threat_pc_long[threat_pc_long[["PC"]] == focal_pc, ]
debug(test_meanshift_varshift)
test_res <- test_meanshift_varshift(
  combo_number = focal_combo_num,
  sig_combos_df = sig_combos,
  filtered_colour_data = test_data,
  n_boots = 1000,
  threatened_spp_only = FALSE,
  focal_sex = "both_sexes"
)

# now run across all combos
spec_test_res <- lapply(
  1:nrow(sig_combos),
  function(focal_combo_num){
    # get focal PC
    focal_metric <- sig_combos[[focal_combo_num, "metric"]]
    # get focal threat
    focal_threat <- sig_combos[[focal_combo_num, "extinction_driver"]]
    # print which combo working on
    message("Focal combo ", focal_combo_num, " of ", nrow(sig_combos), " (", focal_metric, ", ", focal_threat, ")")
    # filter data to focal PC only
    focal_dat <- threat_metrics_long[threat_metrics_long[["metric"]] == focal_metric, ]
    # run testing function for mean shift and log variance ratio
    focal_res <- test_metric_meanshift_varshift(
      combo_number = focal_combo_num,
      sig_combos_df = sig_combos,
      filtered_colour_data = focal_dat,
      n_boots = 1000,
      threatened_spp_only = FALSE,
      focal_sex = "both_sexes"
    )
  }
)
spec_test_res <- do.call(rbind, spec_test_res)


# Save results as CSV
# set filename based on parameters
if(threatened_spp_only_par == TRUE){
  spec_filename <- "metric_ex_driver_meanshift-logvar_threatened_spp_only.csv"
} else if(threatened_spp_only_par == FALSE){
  spec_filename <- "metric_ex_driver_meanshift-logvar_all_spp.csv"
}

write.csv(
  spec_test_res,
  file = here::here(
    "03_output_data", "08_image_feature_metrics_threats",
    spec_filename
  ), 
  row.names = FALSE
)

## Plot results ----

# reload results for plotting
threatened_spp_only_par <- FALSE
# set filename based on parameters
if(threatened_spp_only_par == TRUE){
  spec_filename <- "metric_ex_driver_meanshift-logvar_threatened_spp_only.csv"
} else if(threatened_spp_only_par == FALSE){
  spec_filename <- "metric_ex_driver_meanshift-logvar_all_spp.csv"
}
spec_test_res <- read.csv(
  file = here::here(
    "03_output_data", "08_image_feature_metrics_threats",
    spec_filename
  )
)


# plot significant meanshift drivers - horizontal bars
# first the non-colour proportion ones
mean_shift_plot <- spec_test_res %>% 
  # filter(
  #   mean_shift_p < 0.05
  # ) %>% 
  mutate(
    signif = ifelse(mean_shift_p < 0.05, "y", "n")
  ) %>% 
  filter(metric %in% c("mean_lum", "max_lum", "countershading", "colourfulness")) %>% 
  ggplot(aes(y = ex_driver, x = mean_shift_ses, fill = signif)) + 
  geom_col() + 
  geom_vline(xintercept = 0) + 
  facet_wrap(~ metric, dir = "v") + 
  labs(y = "Extinction driver", x = paste0("Effect size (mean shift)")) + 
  scale_y_discrete(
    labels = c(
      "pollut" = "Pollution",
      "invas_spec" = "Invasive species",
      "hunt_col" = "Hunting & Collection",
      "hab_loss" =  "Habitat Loss",
      "clim_chan" =  "Climate Change",
      "acc_mort" = "Accidental Mortality"
    )
  ) + 
  scale_fill_discrete(palette = c("grey70", "grey25")) + 
  theme_bw() + 
  theme(legend.position = "none")
mean_shift_plot

# and significant variance inequality drivers -  horizontal bars
log_var_ratio_plot <- spec_test_res %>% 
  # filter(
  #   log_var_ratio_p < 0.05
  # ) %>% 
  mutate(
    signif = ifelse(log_var_ratio_p < 0.05, "y", "n"),
    sd_diff = ifelse(delta_sd > 0, "pos", "neg")
  ) %>% 
  filter(metric %in% c("mean_lum", "max_lum", "countershading", "colourfulness")) %>% 
  ggplot(aes(y = ex_driver, x = log_var_ratio_ses, fill = signif)) + 
  geom_col(linewidth = 1) + 
  geom_vline(xintercept = 0) + 
  facet_wrap(~ metric) + 
  labs(y = "Extinction driver", x = "Effect size (variance inequality)") + 
  scale_y_discrete(
    labels = c(
      "pollut" = "Pollution",
      "invas_spec" = "Invasive species",
      "hunt_col" = "Hunting & Collection",
      "hab_loss" =  "Habitat Loss",
      "clim_chan" =  "Climate Change",
      "acc_mort" = "Accidental Mortality"
    )
  ) +
  scale_fill_discrete(palette = c(adjustcolor("lightgrey"), "grey25")) + 
  #  scale_fill_discrete(palette = c(adjustcolor("darkblue", alpha = 0.8), adjustcolor("darkred", alpha = 0.8))) +
  theme_bw() + 
  theme(legend.position = "none")
log_var_ratio_plot

# save plots
if(threatened_spp_only_par == TRUE){
  ms_plot_filename <- "metric_ex_driver_meanshift_threatened_spp_only.svg"
  vi_plot_filename <- "metric_ex_driver_logvarratio_threatened_spp_only.svg"
} else if(threatened_spp_only_par == FALSE){
  ms_plot_filename <- "metric_ex_driver_meanshift_all_spp.svg"
  vi_plot_filename <- "metric_ex_driver_logvarratio_all_spp.svg"
}
ggsave(
  filename = ms_plot_filename,
  plot = mean_shift_plot, 
  device = "svg", 
  path = here::here(
    "04_output_plots", "08_image_feature_metrics_threats", "02_distributional_differences"
  ), 
  width = 180, height = 140, units = "mm"
)
ggsave(
  filename = vi_plot_filename,
  plot = log_var_ratio_plot, 
  device = "svg", 
  path = here::here(
    "04_output_plots", "08_image_feature_metrics_threats", "02_distributional_differences"
  ), 
  width = 180, height = 140, units = "mm"
)
