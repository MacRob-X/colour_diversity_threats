# Functions to obtain detailed information on threats to each avian species from IUCN Red List
# Companion to script 02_scripts/01_get_threat_info.R


#' Create skeleton dataframe for single species threat data populated with NAs
#' We will use this as a skeleton for species which have no threat data to grab.
#'
#' @param binomial_name Binomial name of species WITH THREATS in form 'Genus_species'. As of 2025-12-18 "Aburria_aburri" is a good one to use.
#' @param assessment_year Year of assessment for species to get threat skeleton from. 2023 is good for Aburria_aburri
#' @param api_key Key for IUCN API
#'
#' @returns Threat skeleton in the form of a dataframe of the same structure as a threatened species but populated with NAs
#' @export
#'
#' @examples
make_threat_skeleton <- function(binomial_name, assessment_year, api_key){
  
  # Split binomial into genus and species
  genus_name <- strsplit(binomial_name, split = "_")[[1]][1]
  species_name <- strsplit(binomial_name, split = "_")[[1]][2]
  
  # Get list of assessments for species
  am_assess <- rredlist::rl_species(genus = genus_name, species = species_name, key = api_key)$assessments
  # Get assessment ID for specified year and species
  spec_assess_id <- am_assess[am_assess$year_published == as.character(assessment_year), "assessment_id"]
  # Get threats for specified assessment and species - restrict to first row only
  threats_skeleton <- rredlist::rl_assessment(id = spec_assess_id, key = api_key)$threats[1, ]
  # Populate with NAs
  threats_skeleton[1, ] = NA
  
  return(threats_skeleton)
  
}


# 
#' Create or append CSV of threat data
#' 
#' If no CSV file exists at specified writepath, creates a new CSV and adds species information in the first line. If a CSV already exists, appends to existing CSV.
#'
#' @param iucn_species List of all species for which threat data is being obtained
#' @param binomial_name Binomial name of species of interest in form 'Genus_species'
#' @param csv_writepath Path to write/append species information to CSV
#' @param threats_df Generated within 'get_threat_data' function - dataframe of threat data obtained from IUCN Red List and organised into the correct format
#'
#' @returns N/A - doesn't return anything, just writes or appends to CSV file
#' @export
#'
#' @examples 
handle_csv <- function(iucn_species, binomial_name, csv_writepath, threats_df){
  
  # if first species in list, initialise CSV to save into
  if(which(iucn_species %in% binomial_name) == 1){
    
    # check if file already exists - if it does, ask user if they want to overwrite
    # cancel operation if they don't
    if(file.exists(csv_writepath)){
      message("Warning: CSV file of this name already exists.")
      response <- readline(prompt = "Do you want to overwrite it? (y/n): ")
      clean_response <- tolower(substr(response, 1, 1)) # handle Y, y, Yes, yes, YES etc
      if(clean_response == "y"){
        message("Overwriting CSV.")
      } else {
        stop("Operation cancelled by user")
      }
    }
    
    # initialise csv
    message("Creating CSV.")
    write.csv(threats_df, file = csv_writepath)
    
  } else {
    
    # otherwise append to existing CSV
    message("Appending to CSV.")
    write.table(threats_df, file = csv_writepath, append = TRUE, col.names = FALSE, sep = ",")
    
  }
  
}


# 
#' Deal with species that have no IUCN data 
#' Some species (e.g. Cacicus leucoramphus) have no IUCN data to retrieve, probably because they are no longer a valid taxon but were in previous IUCN versions. This function returns threat information with NA in each field and 'no_data' in notes.
#'
#' @param binomial_name Binomial name of species of interest in form 'Genus_species'
#' @param threats_skeleton Skeleton structure of threat data, populated with NAs
#' @param csv_writepath Path to write/append species information to CSV. If none specified, does not write to CSV but returns dataframe of information.
#'
#' @returns Dataframe of threat information for species in the correct format. Also writes/appends to CSV if a path is specified
#' @export
#'
#' @examples
handle_missing_species <- function(binomial_name, threats_skeleton, csv_writepath = NULL){
  
  warning("No IUCN data for ", binomial_name, ". Returning NA row with 'no_data' in notes.")
  
  # assign threat skeleton
  threat_data <- threats_skeleton
  
  # remove "description" column as we don't need it and it contains commas, which messes with
  # CSV writing
  threat_data$description <- NULL
  
  # initialise df for data (following IUCN threat classification scheme)
  # https://www.iucnredlist.org/resources/threat-classification-scheme
  assessment_year <- NA
  iucn_cat <- NA
  notes <- "no_data"
  threats_df <- cbind(binomial_name, threat_data, assessment_year, iucn_cat, notes)
  # save CSV if required
  if(!is.null(csv_writepath)){
    handle_csv(iucn_species = iucn_species, binomial_name = binomial_name, csv_writepath = csv_writepath, threats_df = threats_df)
  }
  # 2 second delay makes API work better - recommended by IUCN
  Sys.sleep(2)
  return(threats_df)
  
}



# Function to deal with species which have no 'latest' data (e.g. Calliphlox_evelynae)
# 
#' Deal with species which have no 'latest' data
#' 
#' Some species have IUCN specific threat data, but with no assessments marked as 'latest'. These are probably mostly invalid taxa, e.g. Calliphlox_evelynae was moved to Nesophlox_evelynae and subsequently split into two taxa, N. evelynae and N. lyrura). This function retrieves the most recent data from before the specified cutoff year (inclusive), but specifies in the notes 'no_latest_data' to differentiate these from other species.
#'
#' @param binomial_name Binomial name of species of interest in form 'Genus_species'
#' @param species_history Species assessment history obtained from rredlist::rl_species()$assessments
#' @param year_cutoff Cutoff year later than which threat information should not be retrieved (e.g., if a year of 2023 is specified, no assessments from 2024 or later will be considered)
#' @param csv_writepath Path to write/append species information to CSV. If none specified, does not write to CSV but returns dataframe of information.
#'
#' @returns Dataframe of threat information for species in the correct format. Also writes/appends to CSV if a path is specified
#' @export
#'
#' @examples
handle_no_latest_data <- function(binomial_name, species_history, year_cutoff, csv_writepath = NULL){
  
  warning("No 'latest' IUCN data for ", binomial_name, ". Returning most recent data with 'no_latest_data' in notes.")
  
  # use the data from the latest pre-cutoff assessment
  # first check if there is data from before the cutoff year
  if(any(species_history$year_published <= year_cutoff)){
    # should already be in reverse date order but let's order just in case
    ordered_history <- species_history[order(species_history$year_published, decreasing = T), ]
    # get assessment id for relevant year
    assess_id <- ordered_history[which(ordered_history$year_published <= year_cutoff)[1], "assessment_id"]
    # get data
    species_data <- rredlist::rl_assessment(id = assess_id, key = api)
    
    
    # create and populate threats df
    
    # extract threat data
    threat_data <- species_data$threats
    # if there are no threats then we use the threats skeleton - this will then populate the threat matrix with NA
    # for this species
    if(length(threat_data) == 0){
      threat_data <- threats_skeleton
    }
    # remove "description" column as we don't need it and it contains commas, which messes with
    # CSV writing
    threat_data$description <- NULL
    
    notes <- "no_latest_data"
    assessment_year <- species_data$year_published
    iucn_cat <- species_data$red_list_category$code
    threats_df <- cbind(binomial_name, threat_data, assessment_year, iucn_cat, notes)
    
  } else {
    # if there isn't, use the skeleton but with "no_precutoff_data" in notes column
    warning("No IUCN data pre-cutoff year for ", binomial_name, ". Returning NA row with 'no_precutoff_data' in notes.")
    threat_data <- threats_skeleton
    # remove "description" column as we don't need it and it contains commas, which messes with
    # CSV writing
    threat_data$description <- NULL
    notes <- "no_precutoff_data"
    assessment_year <- NA
    iucn_cat <- NA
    threats_df <- cbind(binomial_name, threat_data, assessment_year, iucn_cat, notes)
    
  }
  
  # save CSV if required
  if(!is.null(csv_writepath)){
    handle_csv(iucn_species = iucn_species, binomial_name = binomial_name, csv_writepath = csv_writepath, threats_df = threats_df)
  }
  # 2 second delay makes API work better - recommended by IUCN
  Sys.sleep(2)
  return(threats_df)
  
}


#' Get specific threat data for a single species from IUCN Red List API
#'
#' @param binomial_name Binomial name of species of interest in form 'Genus_species'
#' @param threats_skeleton Skeleton structure of threat data, populated with NAs
#' @param latest Logical: if T, function will grab threat data from the latest IUCN assessment of the species. If F, function will grab threat data from the most recent assessment up to and including the cutoff year. E.g., with a cutoff year of 2023, if a species' has assessments from 2016, 2022, and 2025, latest = TRUE will grab threat data from the 2025 assessment, while latest = FALSE will grab threat data from the 2022 assessment. 
#' @param year_cutoff Cutoff year for non-latest data - see 'latest' param.
#' @param csv_writepath Path to write/append species information to CSV. If NULL, does not write to CSV but returns dataframe of information.
#'
#' @returns Single-row dataframe of threat data for species of interest.
#' @export
#'
#' @examples
get_threat_data <- function(binomial_name, threats_skeleton, latest = TRUE, year_cutoff = NULL, csv_writepath = NULL){
  
  # print species processing
  cat("\r", paste("Processing ", which(iucn_species %in% binomial_name), " of ", length(iucn_species), ": ",  binomial_name, sep = ""))
  
  
  # set genus and species names
  genus_name <- strsplit(binomial_name, split = "_")[[1]][1]
  species_name <- strsplit(binomial_name, split = "_")[[1]][2]
  
  
  # if using latest threat data
  if(latest == TRUE){
    # get species data
    species_data <- try(rredlist::rl_species_latest(genus = genus_name, species = species_name, key = api))
    notes <- NA
    if(class(species_data) == "try-error"){ # check if there is actually any IUCN data for the species - if not, return threats skeleton with explanatory note
      threats_df <- handle_missing_species(binomial_name, threats_skeleton, csv_writepath)
      return(threats_df)
    }
  } else if(latest == FALSE){
    # if not using latest data
    
    # get list of all assessments for species
    species_history <- try(rredlist::rl_species(genus = genus_name, species = species_name, key = api))
    if(class(species_history) == "try-error"){ # check if there is actually any IUCN data for the species - if not, return threats skeleton with explanatory note
      threats_df <- handle_missing_species(binomial_name, threats_skeleton, csv_writepath)
      return(threats_df)
    } else {
      species_history <- species_history$assessments
    }
    # check if the latest assessment is from pre-cutoff year
    # note some species seem to have multiple marked as 'latest' - use the first of these
    # note also some species (e.g. Calliphlox_evelynae) have no records marked as 'latest' - 
    # these are invalid taxa (e.g. Calliphlox_evelynae was moved to Nesophlox_evelynae and subsequently
    # split into two taxa, N. evelynae and N. lyrura). We will keep these species in but put in the notes
    # that they are now reassigned to other taxa
    if(!(any(species_history$latest == TRUE))){
      threats_df <- handle_no_latest_data(binomial_name, species_history, year_cutoff, csv_writepath)
      return(threats_df)
    }
    if(species_history[species_history$latest == TRUE, "year_published"][1] <= year_cutoff){
      # if it is, just use the latest assessment
      species_data <- rredlist::rl_species_latest(genus = genus_name, species = species_name, key = api)
      notes <- NA
    } else {
      # if it isn't, use the data from the first pre-cutoff assessment
      # first check if there is data from before the cutoff year
      if(any(species_history$year_published <= year_cutoff)){
        # should already be in reverse date order but let's order just in case
        ordered_history <- species_history[order(species_history$year_published, decreasing = T), ]
        # get assessment id for relevant year
        assess_id <- ordered_history[which(ordered_history$year_published <= year_cutoff)[1], "assessment_id"]
        # get data
        species_data <- rredlist::rl_assessment(id = assess_id, key = api)
        notes <- NA
      } else {
        # if there isn't, use the skeleton but with "no_precutoff_data" in notes column
        warning("No IUCN data pre-cutoff year for ", binomial_name, ". Returning NA row with 'no_precutoff_data' in notes.")
        threat_data <- threats_skeleton
        # remove "description" column as we don't need it and it contains commas, which messes with
        # CSV writing
        threat_data$description <- NULL
        notes <- "no_precutoff_data"
        assessment_year <- NA
        iucn_cat <- NA
        threats_df <- cbind(binomial_name, threat_data, assessment_year, iucn_cat, notes)
        # save CSV if required
        if(!is.null(csv_writepath)){
          handle_csv(iucn_species = iucn_species, binomial_name = binomial_name, csv_writepath = csv_writepath, threats_df = threats_df)
        }
        # 2 second delay makes API work better - recommended by IUCN
        Sys.sleep(2)
        return(threats_df)
      }
      
    }
  }
  
  # extract threat data
  threat_data <- species_data$threats
  # if there are no threats then we use the threats skeleton - this will then populate the threat matrix with NA
  # for this species
  if(length(threat_data) == 0){
    threat_data <- threats_skeleton
    notes <- "no_threats"
  }
  # remove "description" column as we don't need it and it contains commas, which messes with
  # CSV writing
  threat_data$description <- NULL

  # initialise df for data (following IUCN threat classification scheme)
  # https://www.iucnredlist.org/resources/threat-classification-scheme
  assessment_year <- species_data$year_published
  iucn_cat <- species_data$red_list_category$code
  threats_df <- cbind(binomial_name, threat_data, assessment_year, iucn_cat, notes)

  # 2 second delay makes API work better - recommended by IUCN
  Sys.sleep(2)
  
  # save CSV if required
  if(!is.null(csv_writepath)){
    
    handle_csv(iucn_species = iucn_species, binomial_name = binomial_name, csv_writepath = csv_writepath, threats_df = threats_df)
    
  }
  
  return(threats_df)
  
}
