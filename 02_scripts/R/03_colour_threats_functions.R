# Functions for 03_colour_threats.R, 03a_threat_centroid.R and 03b_threat_centroid_SES.R scripts

# Remove past threats from threat matrix
remove_past_threats <- function(threat_mat){
  
  # filter to only past threat rows (removing NAs)
  past_threat_rows <- threat_mat %>% 
    filter(
      timing == "Past, Unlikely to Return"
    )
  
  # and only ongoing/future etc rows (including NAs)
  current_threat_rows <- threat_mat %>% 
    filter(
      !(timing %in% "Past, Unlikely to Return")
    )
  
  # get associated species
  past_threat_spp <- unique(past_threat_rows$jetz_species)
  current_threat_spp <- unique(current_threat_rows$jetz_species)
  
  # check which species have both past and current threats
  both_threats_spp <- intersect(past_threat_spp, current_threat_spp)
  
  # check which species are ONLY in the past threat rows
  # these have no current threats and can be coded as having no threats
  past_threats_only_spp <- past_threat_spp[which(!past_threat_spp %in% current_threat_spp)]
  
  # code as having no current threats (i.e. make everything NA and add 'no_current_threats' in notes)

  # threat matrix of past threat only species
  threat_mat_pto <- threat_mat[threat_mat$jetz_species %in% past_threats_only_spp, ]
  # set these threats to NA
  cols_to_na <- c("scope", "timing", "internationalTrade", "score","severity","ancestry","virus","ias","text", "code")
  threat_mat_pto[threat_mat_pto$jetz_species %in% past_threats_only_spp, cols_to_na] <- rep(NA, times = length(cols_to_na))
  threat_mat_pto[threat_mat_pto$jetz_species %in% past_threats_only_spp, "notes"] <- "past_threats_only"
  # we only need one row per species for these species, regardless of how many past threats they had
  # so let's remove any duplicate rows
  threat_mat_pto <- threat_mat_pto[!duplicated(threat_mat_pto), ]
  
  
  # For the species which have both past and current threats, I can simply remove the past threats rows
  # from the matrix, leaving the current threat rows, and rbind to the past threat only rows
  new_threat_mat <- rbind(current_threat_rows, threat_mat_pto)
  
  # check the number of rows adds up
  assertthat::assert_that(
    nrow(new_threat_mat) == nrow(threat_mat) - nrow(past_threat_rows) + nrow(threat_mat_pto),
    msg = "Number of rows inconsistent - there is some error in this function."
      
  )
  
  return(new_threat_mat)
  
}


# Assign second order threat types to extinction drivers
# second order threat types must be in a column named 'second_ord_code'
assign_ext_drivers <- function(threat_mat){
  
  # Check for second order code column
  assertthat::assert_that(
    "second_ord_code" %in% colnames(threat_mat),
    msg = "No column named 'second_ord_code' in threat matrix"
  )
  
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
  
  threat_mat <- threat_mat |> 
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
  
  # for species with no threats, add this in the ex_driver column
  threat_mat <- threat_mat %>% 
    mutate(
      ex_driver = ifelse(!is.na(notes) & notes == "no_threats", "no_threats", ex_driver)
    )
  
  # same for species with only past threats - assign these as having no threats
  threat_mat <- threat_mat %>% 
    mutate(
      ex_driver = ifelse(!is.na(notes) & notes == "past_threats_only", "no_threats", ex_driver)
    )
  
  return(threat_mat)
  
}

# Assign second order threat types to binary extinction drivers
# second order threat types must be in a column named 'second_ord_code'
assign_binary_ext_drivers <- function(threat_mat){
  
  # Check for second order code column
  assertthat::assert_that(
    "second_ord_code" %in% colnames(threat_mat),
    msg = "No column named 'second_ord_code' in threat matrix"
  )
  
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
  # All significant codes
  all_sig_codes <- c(
    acc_mort_codes, 
    clim_chan_codes, 
    hab_loss_codes, 
    hunt_col_codes,
    invas_spec_codes,
    other_codes,
    pollut_codes
    )
  
  threat_mat <- threat_mat |> 
    mutate( 
      acc_mort = ifelse(second_ord_code %in% acc_mort_codes, 1, 0),
      clim_chan = ifelse(second_ord_code %in% clim_chan_codes, 1, 0),
      hab_loss = ifelse(second_ord_code %in% hab_loss_codes, 1, 0),
      hunt_col = ifelse(second_ord_code %in% hunt_col_codes, 1, 0),
      invas_spec = ifelse(second_ord_code %in% invas_spec_codes, 1, 0),
      other_threats = ifelse(second_ord_code %in% other_codes, 1, 0),
      pollut = ifelse(second_ord_code %in% pollut_codes, 1, 0),
      non_sig_threat = ifelse(!is.na(second_ord_code) & !(second_ord_code %in% all_sig_codes), 1, 0),
      threat_data_missing = ifelse(is.na(second_ord_code) & notes != "no_threats", 1, 0)
    )
  
  return(threat_mat)
  
}

# Function to check autocorrelation in first non-zero lag of MCMCglmm fixed effects
# From MacDonald et al 2024 - Primate coloration and colour vision: a comparative approach (Supporting Information)
# https://doi.org/10.1093/biolinnean/blad089
autocorrMCMCglmm <- function(mod, var) {
  var2 <- var*var
  vals <- data.frame(coda::autocorr(mod$Sol[,1:var]))
  count <- 0
  for(i in 1:var2){
    if(vals[2, i] >= 0.1){
      return(paste0("Uh oh, problem in column <", colnames(vals[i]), ">; Lag = ", vals[2, i]))
      count <- count + 1
    }
  }
  if(count == 0){
    return("All good!")
  }
}

############### Function allowing estimation of lambda ####################
# Adapted from MacDonald et al 2024 - Primate coloration and colour vision: a comparative approach (Supporting Information)
# https://doi.org/10.1093/biolinnean/blad089
# Based on calculation from de Villemereuil and Nakagawa 2014 (from Garamszegi 2014 Modern phylogenetic comparative methods - Ch 11.2.1)
# Essentially this computes how much of the variance is explained by the random effects, which in this case is the phylogeny
lamCalc <- function(mod, phylo, n_levels_fixed){
  
  # extract names of random effects from model
  random_effects <- colnames(mod$VCV)
  
  # get phylogeny random effect name
  phylo <- as.character(phylo)
  
  # check phylogeny random effect name is actually in the model
  phylo_in_mod <- function(phylo, random_effects){
    assertthat::assert_that(phylo %in% random_effects)
  }
  assertthat::on_failure(phylo_in_mod) <- function(call, env){
    paste0(deparse(call$phylo), " is not in the random effects of the provided model.")
  }
  assertthat::assert_that(phylo_in_mod(phylo, random_effects))
  
  # Calculate Lambda 
  # ratio of phylogenetic variance to total variance [phylogenetic + other random effects + residual + fixed effects]
  # This specific formula is taken from Villemereuil 2024 - Estimation of a biological trait heritability using the animal model and MCMCglmm
  # Available at https://devillemereuil.legtux.org/downloads/
  
  # first calculate fixed effects variance
  compute_varpred <- function(beta, design_matrix) {
    var(as.vector(design_matrix %*% beta))
  }
  X <- mod[["X"]]
  var_fixed <- apply(mod[["Sol"]][, 1:(n_levels_fixed + 1)], 1, compute_varpred, design_matrix = X)
  
  # random effects variance
  var_random <- rowSums(mod$VCV)
  
  # Calculate lambda
  lambda <- mod$VCV[, phylo] / var_fixed + var_random
  
  return(lambda)
}

# Run MCMCglmm model for a single tree i in a distribution N of trees
run_itt <- function(
    i, # tree number
    phy, # tree distribution
    data_mcmcglmm,
    mod_itt, mod_thin, mod_burnin,
    n_samples_tree,
    prior
    ){
  
  # select the ith tree
  tree <- phy[[i]]
  
  animalA <- inverseA(tree)$Ainv
  
  mod <- MCMCglmm(
    log(centr_dists) ~ ex_driver + sex,   # log transform to pull in right skew
    random = ~ jetz_species + PhyloName,
    ginverse = list(PhyloName = animalA),
    prior = g_prior,
    data = data_mcmcglmm,
    rcov = ~ units,
    family = "gaussian",
    nitt = mod_itt,
    thin = mod_thin,
    burnin = mod_burnin,
    pl = TRUE,
    pr = TRUE,
    verbose = FALSE
  )
  
  # return the samples per tree
  mod_res <- list(
    VCV = mod$VCV[1:n_samples_tree, ], # [VCV is posterior distrib of covariance matrices]
    Sol = mod$Sol[1:n_samples_tree, ], # [Sol is posterior distrib of MME solutions (???) - includes fixed effects]
    Liab = mod$Liab[1:n_samples_tree, ] # [Liab is posterior distrib of latent variables])
  )
  return(mod_res)
  
}
