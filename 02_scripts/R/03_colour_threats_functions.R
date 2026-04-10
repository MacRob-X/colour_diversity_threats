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


# Run MCMCglmm model for a single tree i in a distribution N of trees
run_itt <- function(i){
  
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