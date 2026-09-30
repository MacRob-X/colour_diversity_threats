# Functions for 03_colour_threats.R script


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
