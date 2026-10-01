# Examine how colour pattern diversity is distributed among Extinction Drivers
# 2026-09-29
# Robert MacDonald

# Clear environment
rm(list=ls())

# Load custom functions
source(
  here::here(
    "02_scripts", "R", 
    "03b_threat_centroid_SES_functions.R"
  )
)

# Load libraries ----
library(dplyr)
library(ggplot2)
library(parallel)
library(dispRity)

## EDITABLE CODE ##
## Parameters of data to load
# Use latest IUCN assessment data or use most recent assessment data pre- specified cutoff year?
latest <- TRUE
# If not using latest assessment data, specify a cutoff year. Set to NULL if using latest.
cutoff_year <- NULL
# Clade to focus on ("Aves", "Neognaths", "Neoaves", "Passeriformes")
clade <- "Aves"
# Exclude past threats?
exclude_past_threats <- FALSE
# Remove extinct (EX) species?
# Remove extinct species?
remove_extinct <- TRUE
## Parameters for current analysis
# Choose number of null simulations
n_sims <- 1000
# select sex ("M", "F", "All")
sex <- "All"
# select metric ("centr-dist", "nn-k", "nn-count")
# note that nn-k is EXTREMELY slow to run - needs parallelisation (but will probably still
# be too slow to run)
metric <- "centr-dist"

# Load data ----

# Load threat data with centroid distances
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

# Prepare data ----

# create copy of trait data to use in analyses
trait <- threat_data %>% 
  rename(
    species = jetz_species
  )

# Subset to chosen sex
if(sex != "All"){
  trait <- trait[trait[, "sex"] == sex, ]
  species <- species[species[, "sex"] == sex, ]
  unique_species <- NULL # so the parallelisation cluster export doesn't throw an error
} else {
  # get unique species (to subset species/sex pairs later)
  unique_species <- unique(trait[["species"]])
}

# set extinction drivers
ext_drivers <- c("clim_chan", "hab_loss", "pollut", "hunt_col", "acc_mort", "invas_spec")

# Initialise data frame to write results to
results <- as.data.frame(matrix(NA, nrow = 1, ncol = 9))
colnames(results) <- c("ex_driver", "species_richness", "obs_mean", "null_mean", "null_sd", "null_se", "es", "p_value" "ses")
results$ex_driver <- c("all_species")

# Get PC matrix (with species/sex as rownames)
# needs to be matrix for input to dispRity
traits <- trait %>% 
  select(
    starts_with("PC")
  ) %>% 
  as.matrix()
rownames(traits) <- paste(trait$species, trait$sex, sep = "-")  


# get collections of species/sex with those threatened by each extinction driver trimmed OUT
# these will be used to subset the trait data
trimmed_species <- list(
  all_species = paste(trait[, "species"], trait[, "sex"], sep = "-"),
  clim_chan = paste(trait[trait$clim_chan == 0, "species"], trait[trait$clim_chan == 0, "sex"], sep = "-"),
  hab_loss = paste(trait[trait$hab_loss == 0, "species"], trait[trait$hab_loss == 0, "sex"], sep = "-"),
  pollut = paste(trait[trait$pollut == 0, "species"], trait[trait$pollut == 0, "sex"], sep = "-"),
  hunt_col = paste(trait[trait$hunt_col == 0, "species"], trait[trait$hunt_col == 0, "sex"], sep = "-"),
  acc_mort = paste(trait[trait$acc_mort == 0, "species"], trait[trait$acc_mort == 0, "sex"], sep = "-"),
  invas_spec = paste(trait[trait$invas_spec == 0, "species"], trait[trait$invas_spec == 0, "sex"], sep = "-")
)

# set dispRity metric
if(metric == "centr-dist"){
  metric_get <- "centroids"
} else if (metric == "nn-k"){
  metric_get <- "mean.nn.dist"
}else if (metric == "nn-count"){
  metric_get <- "count.neighbours"
}

## Run analysis ----

# Get colour diversity (mean distance to centroid) and species richness of all species
# extract traits for species in this community
trait_vec <- traits[trimmed_species[["all_species"]], ]
allspp_sr <- get_spec_rich(trait_vec, sex)
results[results$ex_driver == "all_species", "species_richness"] <- allspp_sr
allspp_centr_dist <- mean(dispRity(trait_vec, metric = get(metric_get))$disparity[[1]][[1]])
results[results$ex_driver == "all_species", "obs_mean"] <- allspp_centr_dist



# First set up a non-parallelised lapply to apply across the extinction drivers
threat_results <- lapply(
  ext_drivers,
  function(focal_driver){
    
    # extract traits for species in the focal community
    trait_vec <- traits[trimmed_species[[focal_driver]], ]
    spec_rich <- get_spec_rich(trait_vec, sex)
    obs_mean <- mean(dispRity(trait_vec, metric = get(metric_get))$disparity[[1]][[1]])
    
    # Use parLapply to parallelise the simulated distributions
    
    # Set up a cluster using the number of cores (2 less than total number of laptop cores)
    no_cores <- parallel::detectCores() - 2
    cl <- parallel::makeCluster(no_cores)
    
    # export necessary objects to the cluster
    parallel::clusterExport(cl, 
                            c(
                              "trimmed_species", 
                              "traits", 
                              "trait_vec", 
                              "unique_species", 
                              "sex",
                              "get_random_species", 
                              "metric_get"
                              ),
                            envir = environment()
                          )
    parallel::clusterEvalQ(cl, library(dispRity))
    
    sims <- unlist(
      parallel::parLapply(
        cl = cl,
        1:n_sims, 
        function(j) {
          
          # get random species to remove
          random_sample <- get_random_species(trait_vec, unique_species, sex)
          
          random_trait_vec <- traits[random_sample, ]
          
          # Calculate the disparity value for the current simulation (column j)
          disparity_value <- mean(dispRity(random_trait_vec, metric = get(metric_get))$disparity[[1]][[1]])
          
          return(disparity_value)
        }
      )
    )
    
    parallel::stopCluster(cl)
    
    # Calculate null mean, SD, SE, ES, p-value (with Laplace smoothing), SES
    null_mean <- mean(sims)
    null_sd <- sd(sims)
    null_se <- null_sd / sqrt(n_sims)
    es <- obs_mean - null_mean
    p_value <- (length(which(abs(allspp_centr_dist - sims) >= abs(allspp_centr_dist - obs_mean))) + 1) / (n_sims + 1)
    ses <- es / null_sd
    
    to_return <- c(
      focal_driver,
      spec_rich,
      obs_mean,
      null_mean,
      null_sd,
      null_se,
      es,
      p_value,
      ses
    )
    
    return(to_return)
    
  }
)

threat_results <- do.call(rbind, threat_results)
colnames(threat_results) <- colnames(results)

results <- rbind(results, threat_results)

# convert relevant columns to numeric
results <- results %>% 
  mutate(
    across(species_richness:ses, as.numeric)
  )

# Save as CSV
if(latest == TRUE){
  res_filename <- paste0(ext_par, past_threat_par, "_centroid_threat_reduction_SES_latest_2026-01-07.csv")
} else if(latest == FALSE){
  res_filename <- paste(ext_par, past_threat_par, "_centroid_threat_reduction_SES", cutoff_year, "cutoff_year.csv", sep = "_")
}
write.csv(
  results,
  here::here(
    "03_output_data", "03b_threat_centroid_SES",
    res_filename
  ),
  row.names = FALSE
)

## Plotting ----

# Plot results
# plot significant meanshift drivers - horizontal bars
results_plot <- results %>% 
  filter(
    ex_driver != "all_species"
  ) %>% 
  mutate(
    signif = ifelse(abs(ses) > 2, "y", "n")
  ) %>% 
  ggplot(aes(y = ex_driver, x = ses, fill = signif)) + 
  geom_col() + 
  geom_vline(xintercept = 0) + 
  geom_vline(xintercept = -2, linetype = "dashed") + 
  labs(y = "Extinction driver", x = paste0("Standardised Effect Size")) + 
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
results_plot
