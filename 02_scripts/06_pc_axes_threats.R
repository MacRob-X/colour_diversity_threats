# Compare mean colour pattern diversity for species threatened by different IUCN threat groups
# Robert MacDonald

# Clear environment
rm(list=ls())

# Load libraries ----
library(dplyr)
library(ggplot2)
library(ggridges)

## EDITABLE CODE ##
# Use latest IUCN assessment data or use most recent assessment data pre- specified cutoff year?
latest <- TRUE
# If not using latest assessment data, specify a cutoff year. Set to NULL if using latest.
cutoff_year <- NULL
# Clade to focus on ("Aves", "Neognaths", "Neoaves", "Passeriformes")
clade <- "Aves"

# Load data ----

# Load threat data (Jetz taxonomy version)
if(latest == TRUE){
  jetz_threat_filename <- paste0("jetz_threat_matrix_latest_2026-01-07.csv")
} else if(latest == FALSE){
  jetz_threat_filename <- paste("jetz_threat_matrix", cutoff_year, "cutoff_year.csv", sep = "_")
}
threat_matrix <- read.csv(
  file = here::here(
    "03_output_data", jetz_threat_filename
  )
)

# load colour pattern space (created in Chapter 1 - patch-pipeline)
colspace_path <- paste0("G:/My Drive/patch-pipeline/2_Patches/3_OutputData/", clade, "/2_PCA_ColourPattern_spaces/1_Raw_PCA/", clade, ".matchedsex.patches.250716.PCAcolspaces.rds")
colour_space <- readRDS(colspace_path)[["lab"]][["x"]]

# Data preparation ----


# Exclude Data Deficient (DD), Extinct (EX) and Extinct in the Wild (EW) species
threat_matrix <- threat_matrix |> 
  filter(
    !(iucn_cat %in% c("DD", "EX", "EW"))
  )

# add second-order threat codes to threat matrix (derived from third-order codes)
threat_matrix$second_ord_code <- stringr::str_extract(threat_matrix$code, "[^_]*_[^_]*")


# assign second-order IUCN threat types to grouped 'driver of extinction' categories
# same system as Stewart et al 2025 Nat Ecol Evol - grouping is provided in Supplementary Dataset 1
# of that paper
# There are many threats that aren't assigned to one of these groups - this is because these
# threats were non-significant in predicting IUCN threat level in Stewart et al 2025
# We assign these as 'FLAG' in case we want to do anything with them later

# Groups
# Accidental mortality and disturbance
acc_mort_codes <- c("4_2", "5_4", "6_3")
# Climate change and severe weather
clim_chan_codes <- c("11_1", "11_4")
# Habitat loss and degradation
hab_loss_codes <- c("1_2", "1_3", "2_1", "2_2", "2_3", "5_3", "7_1", "7_2")
# Hunting and collecting
hunt_col_codes <- "5_1"
# Invasive species and disease
invas_spec_codes <- c("8_1", "8_2")
# Other ["Threats that affected ten or fewer species were grouped with other threats"]
other_codes <- c("10_1", "10_2", "10_3", "12_1")
# Pollution
pollut_codes <- "9_3"


# First, let's make a version with a binary variable for each extinction driver
threat_matrix_bin <- threat_matrix |> 
  mutate(
    acc_mort = ifelse(second_ord_code %in% acc_mort_codes, 1, 0),
    clim_chan = ifelse(second_ord_code %in% clim_chan_codes, 1, 0),
    hab_loss = ifelse(second_ord_code %in% hab_loss_codes, 1, 0),
    hunt_col = ifelse(second_ord_code %in% hunt_col_codes, 1, 0),
    invas_spec = ifelse(second_ord_code %in% invas_spec_codes, 1, 0),
    other = ifelse(second_ord_code %in% other_codes, 1, 0),
    pollut = ifelse(second_ord_code %in% pollut_codes, 1, 0),
  ) |> 
  select(
    jetz_species, iucn_cat, acc_mort, clim_chan, hab_loss, hunt_col, invas_spec, other, pollut 
  ) |> 
  distinct() |> 
  group_by(jetz_species, iucn_cat) %>%
  # Apply the 'max' function across all remaining columns
  # This ensures if a threat is 1 in ANY row, it becomes 1 in the result
  summarise(across(everything(), max), .groups = 'drop')

# can make a long version of this
threat_bin_long <- threat_matrix_bin |> 
  tidyr::pivot_longer(
    c(acc_mort, clim_chan, hab_loss, hunt_col, invas_spec, other, pollut),
    names_to = "threat_type",
    values_to = "threat_present"
  )

# Alternate version with a categorical variable instead of binary
threat_matrix <- threat_matrix |> 
  mutate(
    ex_driver = second_ord_code
  ) |> 
  mutate( # surely there's a more elegant way to do this
    ex_driver = ifelse(
      ex_driver %in% acc_mort_codes,
      "acc_mort",
      ifelse(
        ex_driver %in% clim_chan_codes,
        "clim_chan",
        ifelse(
          ex_driver %in% hab_loss_codes,
          "hab_loss",
          ifelse(
            ex_driver %in% hunt_col_codes,
            "hunt_col",
            ifelse(
              ex_driver %in% invas_spec_codes,
              "invas_spec",
              ifelse(
                ex_driver %in% other_codes,
                "other",
                ifelse(
                  ex_driver %in% pollut_codes,
                  "pollut",
                  ifelse(
                    is.na(ex_driver),
                    NA,
                    "FLAG"
                  )
                )
              )
            )
          )
        )
      )
    )
  )


# For now, let's make all the flagged threats (i.e. those which are not significant predictors
# of extinction risk) NA, as we're not interested them
# threat_matrix <- threat_matrix |> 
#   mutate(
#     ex_driver = ifelse(
#       ex_driver == "FLAG",
#       NA,
#       ex_driver
#     )
#   )
# NOTE 18/02/2026
# I don't actually want to do the above snippet, as I need to treat these species differently to
# species for which we truly have no threat data (which will be NA)
# Instead I will mark these as "no_sig_driver" so that I can include them in downstream 
# analysis with confidence
threat_matrix <- threat_matrix |>
  mutate(
    ex_driver = ifelse(
      ex_driver == "FLAG",
      "no_sig_driver",
      ex_driver
    )
  )




# Get threat data with PC axes
colour_space_sppsex <- data.frame(species = sapply(strsplit(rownames(colour_space), split = "-"), "[", 1),
                                  sex = sapply(strsplit(rownames(colour_space), split = "-"), "[", 2), 
                                  colour_space)
threat_colour <- threat_matrix |> 
  inner_join(colour_space_sppsex, by = join_by("jetz_species" == "species"))
# this throws a warning but it's just because we have male and female data together
# - it's not a problem

# remove duplicates based on second-order code
threat_colour_clean <- threat_colour |> 
  distinct(second_ord_code, jetz_species, sex, .keep_all = TRUE)


# Pivot longer, so we can facet by PC axis
threat_colour_long <- threat_colour_clean |> 
  tidyr::pivot_longer(
    cols = starts_with("PC"),
    names_to = "PC",
    values_to = "PC_value"
  ) 


# Plotting ----

# boxplot of PC values by threat type - exclude LC species
# Compare to all threatened species
threat_colour_long |> 
  bind_rows(
    threat_colour_long |> 
      mutate(
        ex_driver = "all"
      )
  ) |> 
  filter(
    iucn_cat != "LC",
    !(ex_driver %in% c(NA, "other")),
    PC %in% paste0("PC", 2)
  ) |> 
  mutate(
    ex_driver = factor(ex_driver, levels = c("all", "hab_loss", "hunt_col", "clim_chan", "invas_spec", "acc_mort", "pollut")),
    sex = factor(sex, levels = c("M", "F"))
  ) |> 
  ggplot(aes(x = ex_driver, y = PC_value, fill = ex_driver)) + 
  geom_boxplot(outliers = F) + 
  facet_grid(~ sex)

# Density plot of all species (excluding LC) vs species threatened by each threat type,
# faceted by PC axis (first 7 PCs only)
threat_colour_long |> 
  bind_rows(
    threat_colour_long |> 
      filter(
        ex_driver != "hunt_col"
      ) |> 
      mutate(
        ex_driver = "all"
      ) |> 
      distinct()
  ) |> 
  filter(
    iucn_cat != "LC",
    !(ex_driver %in% c(NA, "other")),
    PC %in% paste0("PC", 1)
  ) |> 
  mutate(
    ex_driver = factor(ex_driver, levels = c("all", "hab_loss", "hunt_col", "clim_chan", "invas_spec", "acc_mort", "pollut")),
    sex = factor(sex, levels = c("M", "F"))
  ) |> 
  filter(
    ex_driver %in% c("all", "hunt_col")
  ) |> 
  ggplot(aes(x = PC_value, fill = ex_driver)) + 
  geom_density(alpha = 0.4) + 
  geom_segment(x = 0, xend = 0, y = 0, yend = 0.04, lwd = 0.3, colour = "grey30") + # Vertical line at 0
  facet_grid(rows = vars(PC), cols = vars(sex)) + 
  theme_minimal()

# function to do the above
plot_pc_distrib <- function(pc_axis, threat_type, long_data, absolute = FALSE) {
  
  if(absolute == TRUE){
    long_data$PC_value <- abs(long_data$PC_value)
  }
  
  p <- long_data |> 
    bind_rows(
      long_data |> 
        filter(
          ex_driver != threat_type
        ) |> 
        mutate(
          ex_driver = "all"
        ) |> distinct()
    ) |> 
    filter(
      iucn_cat != "LC",
      !(ex_driver %in% c(NA, "other")),
      PC == pc_axis
    ) |> 
    mutate(
      ex_driver = factor(ex_driver, levels = c("all", "hab_loss", "hunt_col", "clim_chan", "invas_spec", "acc_mort", "pollut")),
      sex = factor(sex, levels = c("M", "F"))
    ) |> 
    filter(
      ex_driver %in% c("all", threat_type)
    ) |> 
    ggplot(aes(x = PC_value, fill = ex_driver)) + 
    geom_density(alpha = 0.4) + 
    scale_fill_discrete(name = "Threats", labels = c("All threats", "Specific threat")) + 
    labs(title = names(threat_types)[threat_types == threat_type]) + 
    geom_segment(x = 0, xend = 0, y = 0, yend = 0.04, lwd = 0.3, colour = "grey30") + # Vertical line at 0
    facet_grid(cols = vars(sex)) + 
    theme_minimal() + 
    theme(
      #      legend.position = "none",
      plot.title = element_text(size = 10, vjust = -4, hjust = -0.1),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    )
  
  return(p)
}

threat_types <- c(hab_loss = "hab_loss", hunt_col = "hunt_col", invas_spec = "invas_spec", clim_chan = "clim_chan", pollut = "pollut", acc_mort = "acc_mort")
pc_axis <- "PC2"
pc_plots <- lapply(
  threat_types,
  plot_pc_distrib,
  pc_axis = pc_axis,
  long_data = threat_colour_long,
  absolute = FALSE
)
p_pc1_allthreats <- ggpubr::ggarrange(plotlist = pc_plots, common.legend = TRUE, ncol = 2, nrow = 3)
ggpubr::annotate_figure(
  p_pc1_allthreats,
  bottom = ggpubr::text_grob(pc_axis),
  left = ggpubr::text_grob("Density", rot = 90)
)


# Same plot, but with absolute PC values - this is effectively distance to centroid of each PC
abs_pc_plots <- lapply(
  threat_types,
  plot_pc_distrib,
  pc_axis = pc_axis,
  long_data = threat_colour_long,
  absolute = TRUE
)
abs_allthreats <- ggpubr::ggarrange(plotlist = abs_pc_plots, common.legend = TRUE, ncol = 2, nrow = 3)
ggpubr::annotate_figure(
  abs_allthreats,
  bottom = ggpubr::text_grob(pc_axis),
  left = ggpubr::text_grob("Density", rot = 90)
)


# STATISTICS ----

# Perform a Cramer-von Mises test to check for difference in distribution of PC1 values between 
# species threatened by hunting and collection and species not threatened by hunting and collection
# Exclude LC species

pc_axis <- "PC1"
focal_driver <- "hunt_col"

focal_distrib <- threat_colour_long |> 
  filter(
    iucn_cat != "LC", 
    PC == pc_axis, 
    ex_driver == focal_driver
    ) |> 
  pull(PC_value)
non_focal_distrib <- threat_colour_long |> 
  filter(
    iucn_cat != "LC", 
    PC == pc_axis, 
    ex_driver != focal_driver
    ) |> 
  pull(PC_value)

cvm <- twosamples::cvm_test(
  a = focal_distrib,
  b = non_focal_distrib,
)
cvm
plot(cvm)

# Calculate Standard Effect Size of CVM test statistic by randomly sampling from the
# non-hunting&collection distribution 
# Note that I think this is essentially what twosamples::cvm_test does under the hood anyway
# first calculate CVM test statistic for the observed distribution
cvm_obs <- twosamples::cvm_stat(
  a = hc_distrib,
  b = non_hc_distrib,
)
random_samples <- lapply(1:2000, 
                         function(i){
                           sample(
                             non_hc_distrib,
                             size = length(hc_distrib),
                             replace = FALSE
                           )
                         }
                         ) 
cvm_stats_random <- unlist(
  lapply(
    random_samples,
    function(sample){
      cvm <- twosamples::cvm_stat(
        a = sample,
        b = non_hc_distrib,
      )
      return(cvm)
    }
  )
)

# Now calculate SES by calculating (obs - mean(null)) / sd(null)
# SES > 1.96 indicates deviation from difference between distributions expected by chance
ses <- (cvm_obs - mean(cvm_stats_random)) / sd(cvm_stats_random)
# plot histogram of random sample distribution with abserved value as red dotted line
hist(cvm_stats_random, breaks = 50)
abline(v = cvm_obs, col = "red")
# compare with the bootstrap distribution of twosamples::cvm_test
plot(cvm)
# it looks like my SES code is doing essentially the same thing as twosamples::cvm_test
# So, no need to manually do the random sampling myself - I can let the premade function do it


# Set it up to loop over the different PC axes, for a single threat
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
  random_samples <- lapply(1:n, 
                           function(i){
                             sample(
                               all_distrib,
                               size = length(focal_distrib),
                               replace = TRUE
                             )
                           }
  )
  
  # Calculate CVM for each random distribution, vs the original non-focal-threat distribution
  null_test_stats <- unlist(
    lapply(
      random_samples,
      function(sample){
        cvm <- stat_func(
          a = sample,
          b = non_focal_distrib,
        )
        return(cvm)
      }
    )
  )
  
  # Calculate SES (obs - mean(null)) / sd(null)
  # SES > 1.96 indicates deviation from difference between distributions expected by chance
  random_mean <- mean(null_test_stats)
  random_sd <- sd(null_test_stats)
  ses <- (test_stat_obs - random_mean) / random_sd
  # calculate p-value with Laplace smoothing
  # method: rank/counting: proportion of null values equal to or larger than observed CVM statistic
  # This makes no assumptions about the shape of the null distribution, unlike standard
  # p-value calculations (which assume normal distribution). CVM stat distribution is likely
  # a chi-square distribution, not normal
  # This is the method used under the hood in twosamples::cvm_test
  p_val <- (length(which(null_test_stats >= test_stat_obs)) + 1) / (n + 1)
  
  # bind results together for output
  res <- list(
    extinction_driver = threat,
    ses = ses,
    test_statistic = stat_test,
    test_stat_observed = test_stat_obs,
    test_stat_null_mean = random_mean,
    test_stat_null_sd = random_sd,
    p_value = p_val
  )
  
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

# Now apply across each extinction driver and PC
# could parallelise this but it only takes ~30s for 7 axes
# define threat types and focal PCs first 
threat_types <- c(hab_loss = "hab_loss", hunt_col = "hunt_col", invas_spec = "invas_spec", clim_chan = "clim_chan", pollut = "pollut", acc_mort = "acc_mort")
pcs <- paste0("PC", 1:3)

# Run across PCs/extinction drivers
cvm_res <- pbapply::pblapply(
  
  threat_types,
  
  function(threat_type){
    
    cvm <- pbapply::pblapply(
      pcs, 
      test_pc_threat,
      threat_colour_long = threat_colour_long,
      threat = threat_type,
      n = 1000,
      seed = 42,
      return_distrib = FALSE,
      stat_test = "dts"
    )
    names(cvm) <- pcs
    cvm <- data.table::rbindlist(cvm, idcol = "PC")
    # Bonferroni correction for p-values
    # cvm$p_adjusted <- p.adjust(cvm$p_value, method = "bonferroni")
    return(cvm)
    
  }
  
)
cvm_res <- do.call(rbind, cvm_res)
# Bonferroni correction for p-values
cvm_res$p_adjusted <- p.adjust(cvm_res$p_value, method = "bonferroni")

# inspect results
cvm_res |> 
  # filter(
  #   abs(ses) > 1.96
  # ) |> 
  arrange(
    p_adjusted
  )


### ALTERNATIVE HIERARCHICAL MODELLING APPROACH ----

# use the clean (non-long data)

# Add an extra column called 'focal' threat - this needs to be specified each time the model is run
# Here we'll use hunting and collection as the focal variable
# Remove duplicate rows and also filter out non-threatened (LC) species
focal_threat <- "hunt_col"
focal_pc <- "PC3"
mod_data <- threat_colour_clean |> 
  mutate(
    focal_threat = as.factor(ifelse(ex_driver == focal_threat, 1, 0))
  ) |> 
  distinct(
    jetz_species, focal_threat,
    .keep_all = TRUE
  ) |> 
  filter(
    iucn_cat != "LC"
  )

# set up model
form <- formula(get(focal_pc) ~ focal_threat)
mod <- lm(form, data = mod_data)
summary(mod)
# I think this is essentailly just testing for a difference in means between the two groups - 
# not holistic difference in distribution
# CONCLUSION: continue with my previous methodology (random resampling combined with CVM test)