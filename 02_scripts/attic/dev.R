# dev for Chapter 2 - colour loss mechanisms/threats to colour diversity

# 12/12/2025 -----
# Obtain threat info from IUCN redlist
# Bastardised/updated code from Stewart et al 2025

# set iucn API key
api <- "zNXM3miBBknHAGv5mKZraLuhPg2sUTs5c2Ay"

# get full list of birds in IUCN data
#iucn_aves = rredlist::rl_class(class = "aves", key = api)

# save, if necessary
# saveRDS(
#   iucn_aves,
#   file = here::here(
#     "3_output_data", "iucn_data_class_aves.RDS"
#     
#   )
# )

# read in again, if necessary
iucn_aves <- readRDS(
  here::here(
  "3_output_data", "iucn_data_class_aves.RDS"
  )
)

# get data
aves_data <- iucn_aves$assessments

# get unique species names with underscore
iucn_species <- unique(gsub(" ", "_", aves_data$taxon_scientific_name))
# keep only full species
iucn_species <- iucn_species[which(stringr::str_count(iucn_species, "_") == 1)]

# get IUCN threat classification scheme
iucn_threat_types <- rredlist::rl_threats(key = api)$threats
# restrict to second-level classification (after Stewart et al 2025)
second_ord_threats <- stringr::str_extract(iucn_threat_types$code, "[^_]*_[^_]*")
second_ord_threats <- unique(second_ord_threats[!is.na(second_ord_threats)])

# function to create or append CSV of threat data
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
    write.csv(threats_df, file = csv_writepath)
    
  } else {
    
    # otherwise append to existing CSV
    message("Appending to CSV.")
    write.table(threats_df, file = csv_writepath, append = TRUE, col.names = FALSE, sep = ",")
    
  }
  
}

# get specific threat data for all species
# note that this gets the latest data - this may be different to the 2023-1 IUCN Red List version
# I used in Chapter 1 for some species
# There seem to be a lot of species that were reassessed in 2024
# Probably the easiest way is to use the latest IUCN Red List data and rerun Chapter 1 global 
# diversity loss analysis with the newest IUCN version - this should show that the results are 
# near-identical and it's fine to use the newest version of the threats as well
get_threat_data <- function(binomial_name, latest = TRUE, year_cutoff = NULL, csv_writepath = NULL){
  
  # print species processing
  cat("\r", paste("Processing ", which(iucn_species %in% binomial_name), " of ", length(iucn_species), ": ",  binomial_name, sep = ""))
  
  
  # set genus and species names
  genus_name <- strsplit(binomial_name, split = "_")[[1]][1]
  species_name <- strsplit(binomial_name, split = "_")[[1]][2]
  
  # if using latest threat data
  if(latest == TRUE){
    # get species data
    species_data <- rredlist::rl_species_latest(genus = genus_name, species = species_name, key = api)
  } else if(latest == FALSE){
    # if not using latest data
    
    # get list of all assessments for species
    species_history <- rredlist::rl_species(genus = genus_name, species = species_name, key = api)$assessments
    # check if the latest assessment is from pre-cutoff year
    # note some species seem to have multiple marked as 'latest' - use the first of these
    if(species_history[species_history$latest == TRUE, "year_published"][1] <= year_cutoff){
      # if it is, just use the latest assessment
      species_data <- rredlist::rl_species_latest(genus = genus_name, species = species_name, key = api)
    } else {
      # if it isn't, use the data from the first pre-cutoff assessment
      # first check if there is data from before the cutoff year
      if(any(species_history$year_published <= year_cutoff)){
        assess_id <- species_history[which(species_history$year_published <= year_cutoff)[1], "assessment_id"]
        species_data <- rredlist::rl_assessment(id = assess_id, key = api)
      } else {
        # if there isn't, return NAs
        message("No IUCN data pre-cutoff year for ", binomial_name, ". Returning row of NAs.")
        threats_df <-  matrix(NA, nrow = 1, ncol = length(second_ord_threats) + 3)
        colnames(threats_df) <- c("species", "iucn_cat", "assessment_year", second_ord_threats)
        threats_df[, "species"] <- binomial_name
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
  threat_data <- species_data$threats$code
  # if species is LC this will return NULL - this will then populate the threat matrix with 0
  # for this species
  
  # restrict to second-level threat classification (after Stewart et al 2025)
  threat_data <- stringr::str_extract(threat_data, "[^_]*_[^_]*") # keep first two levels and remove second underscore and everything after
  # keep only unique second-level threat classifications
  threat_data <- unique(threat_data)
  
  # initialise df for data (following IUCN threat classification scheme)
  # https://www.iucnredlist.org/resources/threat-classification-scheme
  threats_df <- matrix(NA, nrow = 1, ncol = length(second_ord_threats) + 3)
  colnames(threats_df) <- c("species", "iucn_cat", "assessment_year", second_ord_threats)
  threats_df[, "species"] <- binomial_name
  threats_df[, "iucn_cat"] <- species_data$red_list_category$code
  threats_df[, "assessment_year"] <- species_data$year_published
  
  threats_df[, second_ord_threats] <- unlist(lapply(colnames(threats_df), function(col){
    if(!(col %in% c("species", "iucn_cat", "assessment_year")) & col %in% threat_data){
      threats_df[, col] <- 1
    } else if(!(col %in% c("species", "iucn_cat", "assessment_year")) & !(col %in% threat_data)){
      threats_df[, col] <- 0
    } 
    
  }))
  
  # 2 second delay makes API work better - recommended by IUCN
  Sys.sleep(2)
  
  # save CSV if required
  if(!is.null(csv_writepath)){
    
    handle_csv(iucn_species = iucn_species, binomial_name = binomial_name, csv_writepath = csv_writepath, threats_df = threats_df)
    
  }
  
  
  return(threats_df)
  
}

debug(get_threat_data)

# get threat data for all species - with a 2023 cutoff
cutoff_year <- 2023
filename <- paste("iucn_threat_matrix", cutoff_year, "cutoff_year.csv", sep = "_")
writepath <- here::here(
  "3_output_data", filename
)
allspec_threat_data <- lapply(iucn_species[608:1000], get_threat_data, latest = FALSE, year_cutoff = cutoff_year, csv_writepath = writepath)
threat_matrix <- do.call(rbind, allspec_threat_data)

# save as CSV
filename <- paste("iucn_threat_matrix", cutoff_year, "cutoff_year.csv", sep = "_")
write.csv(
  threat_matrix,
  here::here(
    "3_output_data", filename
  )
)

# get threat data for all species - latest data
allspec_threat_data <- lapply(iucn_species, get_threat_data, latest = TRUE)
threat_matrix <- do.call(rbind, allspec_threat_data)

# save as CSV
filename <- paste0("iucn_threat_matrix_latest_", Sys.Date(), ".csv")
write.csv(
  threat_matrix,
  here::here(
    "3_output_data", filename
  )
)




# Match taxonomies ----






rl_species <- rl_sp$assessments

species_list=unique(rl_species$taxon_scientific_name)
find_length=function(sp){
  return(length(stringr::str_split(sp," ")[[1]]))
}
split=sapply(species_list, find_length)


get_species=function(sp){
  return(paste(str_split(sp," ")[[1]][1],str_split(sp," ")[[1]][2], sep=" "))
}

probs=species_list[split!=2]
prob_res=sapply(probs, get_species)

sp_list_res=c(species_list_whole=species_list[split==2],
              prob_res[!prob_res %in% species_list])

threats_skeleton=rredlist::rl_threats(species_list[1], key = api)$result[1,]
threats_skeleton[1,]="NA"


# try extracting threat data for a single species
sp_data <- rredlist::rl_species_latest(genus = "Daubentonia", species = "madagascariensis", key = api)
threat_data <- sp_data$threats$code
# restrict to second-level threat classification (after Stewart et al 2025)
threat_data <- stringr::str_extract(threat_data, "[^_]*_[^_]*") # keep first two levels and remove second underscore and everything after
# keep only unique second-level threat classifications
threat_data <- unique(threat_data)

# now extract species data for all avian species
