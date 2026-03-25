# Simulate colour pattern diversity loss for given species richness loss
# Use n_iter > 1 if want more than one sample per given species richness loss
simulate_cd_loss <- function(n_sample, spp_list, centr_dists, sex, n_iter = 1, threatened_spp_only = F, iucn_dat = NULL){
  
  # if restricting sampling to threatened species only, check that IUCN data is available
  # and get threatened species list
  if(threatened_spp_only == TRUE){
    if(is.null(iucn_dat)){
      stop("Provide IUCN category data to restrict smapling to threatened species only")
    } else{
      threatened_spp_list <- iucn_dat[iucn_dat$iucn_cat %in% c("CR", "EN", "VU"), "species_birdtree"]
      # check random sample size is smaller than number of threatened species
      n_threatened_spp <- length(threatened_spp_list)
      if(n_sample >= n_threatened_spp){
        stop(paste0("Sample size larger than number of threatened species. Please select a sample size smaller than ", n_threatened_spp))
      }
    }
  }
  
  # trim centroid distance data to sex of interest and species of interest
  centr_dists <- centr_dists[centr_dists$sex == sex & centr_dists$species %in% spp_list, ]
  
  # get baseline diversity (mean distance to centroid of all species)
  baseline_cd <- mean(centr_dists$centr_dists)
  
  # simulate one round of random extinction at a time
  random_exts <- lapply(
    1:n_iter,
    function(n_itt){
      
      # select random sample of species - restrict to threatened species, if requested
      if(threatened_spp_only == FALSE){
        random_spp <- sample(centr_dists$species, size = n_sample, replace = FALSE)
      } else{
        random_spp <- sample(threatened_spp_list, size = n_sample, replace = FALSE)
      }
      
      
      # trim centroid data to random sample
      it_dat <- centr_dists[centr_dists$species %in% random_spp, ]
      
      # get centroid distance distribution
      it_dat <- it_dat$centr_dists
      
      # get species richness loss, colour pattern diversity and colour pattern diversity loss
      sample_stat <- c(n_sample, mean(it_dat), baseline_cd - mean(it_dat))
      
      return(sample_stat)
      
    }
  )
  
  # bind into a single matrix
  random_exts <- do.call(rbind, random_exts)
  
  # name columns
  colnames(random_exts) <- c("sr_loss_abs", "mean_centr_dist", "mean_cd_loss_abs")
  
  # return result
  return(random_exts)
  
}
