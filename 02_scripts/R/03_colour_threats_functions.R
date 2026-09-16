# Functions for 03_colour_threats.R script


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
