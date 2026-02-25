# Functions for examining individual PC axis - threat relationships
# Statistics and plotting

# function to get proportional density plot of species threatened by a particular threat
# proportional to density of all species
prop_dens_2d <- function(threat_col_data, focal_threat, x_axis, y_axis, threatened_spp_only = FALSE, n_bins = 200, return_params = FALSE){
  
  # Get rows corresponding to species threatened by focal threat
  focal_threat_rows <- which(threat_col_data$ex_driver == focal_threat)
  # If making proportional to ALL other species, get all other rows
  if(threatened_spp_only == FALSE){
    non_focal_threat_rows <- (1:nrow(threat_col_data))[-focal_threat_rows]
  } else {
    # If making proportional to all other THREATENED species, get rows of all other 
    # species for which we have threat data and which are threatened (i.e. not LC)
    non_focal_threat_rows <- which(threat_col_data$ex_driver != focal_threat & threat_col_data$iucn_cat != "LC")
    # this will automatically exclude NA rows
  }
  
  # Get 2D density of non-focal rows
  non_focal_2dkd <- MASS::kde2d(
    x = threat_umap_clean[non_focal_threat_rows, x_axis], 
    y = threat_umap_clean[non_focal_threat_rows, y_axis],
    n = n_bins
  )
  # Get 2D density of focal rows
  focal_2dkd <- MASS::kde2d(
    x = threat_umap_clean[focal_threat_rows, x_axis], 
    y = threat_umap_clean[focal_threat_rows, y_axis],
    n = n_bins
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
