# Functions to go with 03b_threat_centroid_SES.R

# get species richness from trait vector
# this is functionalised to hide the if clause which would make the code difficult to read
get_spec_rich <- function(trait_vec, sex){
  
  if(sex == "All"){
    spec_rich <- nrow(trait_vec) / 2  # divide by two because this contains males and females
  } else {
    spec_rich <- nrow(trait_vec)
  }
  
  return(spec_rich)
  
}


# get random species sample
get_random_species <- function(trait_vec, unique_species, sex, replace_par = FALSE){
  
  # if using males and females, need to get paired sex samples
  if(sex == "All"){
    spp_sample <- sample(x = unique_species, size = nrow(trait_vec)/2, replace = replace_par)
    coms <- c(paste0(spp_sample, "-M"), paste0(spp_sample, "-F"))
  } else {
    coms <- sample(x = all, size = nrow(trait_vec), replace = FALSE)
  }
  
}
