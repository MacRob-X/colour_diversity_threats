# dev for Chapter 2 - colour loss mechanisms/threats to colour diversity

# 12/12/2025 -----
# Obtain threat info from IUCN redlist
# Bastardised/updated code from Stewart et al 2025

# Functions ----

# create or append CSV of threat data
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

# function to deal with species that have no IUCN data (e.g. Cacicus leucoramphus)
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
# These are probably mostly invalid taxa, e.g. Calliphlox_evelynae was moved to Nesophlox_evelynae 
# and subsequently split into two taxa, N. evelynae and N. lyrura). We will keep these species in 
# but put in the notes that they are now reassigned to other taxa
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

# get specific threat data for all species
# note that this gets the latest data - this may be different to the 2023-1 IUCN Red List version
# I used in Chapter 1 for some species
# There seem to be a lot of species that were reassessed in 2024
# Probably the easiest way is to use the latest IUCN Red List data and rerun Chapter 1 global 
# diversity loss analysis with the newest IUCN version - this should show that the results are 
# near-identical and it's fine to use the newest version of the threats as well
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
  
  # restrict to second-level threat classification (after Stewart et al 2025)
  #threat_data <- stringr::str_extract(threat_data, "[^_]*_[^_]*") # keep first two levels and remove second underscore and everything after
  # keep only unique second-level threat classifications
  #threat_data <- unique(threat_data)
  
  # initialise df for data (following IUCN threat classification scheme)
  # https://www.iucnredlist.org/resources/threat-classification-scheme
  assessment_year <- species_data$year_published
  iucn_cat <- species_data$red_list_category$code
  threats_df <- cbind(binomial_name, threat_data, assessment_year, iucn_cat, notes)
  # for second-order threats only
  # threats_df <- matrix(NA, nrow = 1, ncol = length(iucn_threat_types) + 3)
  # colnames(threats_df) <- c("species", "iucn_cat", "assessment_year", iucn_threat_types)
  # threats_df[, "species"] <- binomial_name
  # threats_df[, "iucn_cat"] <- species_data$red_list_category$code
  # threats_df[, "assessment_year"] <- species_data$year_published
  
  # threats_df[, iucn_threat_types] <- unlist(lapply(colnames(threats_df), function(col){
  #   if(!(col %in% c("species", "iucn_cat", "assessment_year")) & col %in% threat_data){
  #     threats_df[, col] <- 1
  #   } else if(!(col %in% c("species", "iucn_cat", "assessment_year")) & !(col %in% threat_data)){
  #     threats_df[, col] <- 0
  #   } 
  #   
  # }))
  
  # 2 second delay makes API work better - recommended by IUCN
  Sys.sleep(2)
  
  # save CSV if required
  if(!is.null(csv_writepath)){
    
    handle_csv(iucn_species = iucn_species, binomial_name = binomial_name, csv_writepath = csv_writepath, threats_df = threats_df)
    
  }
  
  
  return(threats_df)
  
}

# Workflow ----

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
  "03_output_data", "iucn_data_class_aves.RDS"
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
# second_ord_threats <- stringr::str_extract(iucn_threat_types$code, "[^_]*_[^_]*")
# second_ord_threats <- unique(second_ord_threats[!is.na(second_ord_threats)])


# debug(get_threat_data)

# set up threat skeleton - we will use this for species which have no threats
# use a species which has threats, but set them to NA
# let's use the assessment ID for the 2023 assessment of Aburria aburri
am_assess <- rredlist::rl_species(genus = "Aburria", species = "aburri", key = api)$assessments
a_a_assess_id <- am_assess[am_assess$year_published == "2023", "assessment_id"]
threats_skeleton <- rredlist::rl_assessment(id = a_a_assess_id, key = api)$threats[1, ]
threats_skeleton[1, ] = NA

# get threat data for all species - with a 2023 cutoff
# this will write to a new CSV if the specified filename and path doesn't exist yet, or append if
# it already exists
cutoff_year <- 2023
filename <- paste("iucn_threat_matrix", cutoff_year, "cutoff_year.csv", sep = "_")
writepath <- here::here(
  "03_output_data", filename
)
allspec_threat_data <- lapply(iucn_species[1501:length(iucn_species)], get_threat_data, latest = FALSE, year_cutoff = cutoff_year, csv_writepath = writepath, threats_skeleton = threats_skeleton)
threat_matrix <- do.call(rbind, allspec_threat_data)

# troubleshoot individual species
# Problem species: 606; 1163
# 606 (Apalis_karamojae): no data pre-2023
# 1163 (Cacicus leucoramphus): no data at all
# 1214 (Calliphlox_evelynae):
debug(get_threat_data)
get_threat_data(iucn_species[1214], latest = F, threats_skeleton = threats_skeleton, year_cutoff = cutoff_year)
try_res <- try(rredlist::rl_species_latest(genus = "Cacicus", species = "leucoramphus", key = api))
rredlist::rl_species(genus = "Cacicus", species = "leucoramphus", key = api)
# Cacicus leucoramphus doesn't seem to have any assessment data even though it's in the big IUCN
# species list - it looks like it's now conspecific with C. chrysonotus, but they're split in 
# the colour data so ideally I'd like to get the C. leucoramphus assessment if possible
# for now let's add a "try" to the species data retrieval

# get threat data for all species - latest data
filename <- paste("iucn_threat_matrix_latest_", Sys.Date(), ".csv")
writepath <- here::here(
  "03_output_data", filename
)
allspec_threat_data <- lapply(iucn_species, get_threat_data, latest = TRUE, csv_writepath = writepath, threats_skeleton = threats_skeleton)
threat_matrix <- do.call(rbind, allspec_threat_data)




# Match taxonomies ----


# Regress distance to centroid against threat type ----

library(dplyr)
library(ggplot2)

# load threat matrix
cutoff_year <- 2023
filename <- paste("iucn_threat_matrix", cutoff_year, "cutoff_year.csv", sep = "_")
threat_matrix <- read.csv(
  here::here(
    "03_output_data", filename
  )
)

# load colour pattern space (created in Chapter 1 - patch-pipeline)
colspace_path <- "G:/My Drive/patch-pipeline/2_Patches/3_OutputData/Neognaths/2_PCA_ColourPattern_spaces/1_Raw_PCA/Neognaths.matchedsex.patches.250716.PCAcolspaces.rds"
colour_space <- readRDS(colspace_path)[["lab"]][["x"]]

# add second-order threat codes
threat_matrix$second_ord_code <- stringr::str_extract(threat_matrix$code, "[^_]*_[^_]*")

# assign second-order IUCN threat types to grouped 'driver of extinction' categories
# same system as Stewart et al 2025 Nat Ecol Evol - grouping is provided in Supplementary Dataset 1
# of that paper
# There are many threats that aren't assigned to one of these groups - this is because these
# threats were non-significant in predicting IUCN threat level in Stewart et al 2025

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

threat_matrix <- threat_matrix |> 
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

# Plot 
threat_matrix |>
  ggplot(aes(x = ex_driver)) + 
  geom_bar()


# calculate distance to centroid from colourspace
centr_dists <- dispRity::dispRity(colour_space, metric = dispRity::centroids)$disparity[[1]][[1]]
centr_dists <- data.frame(
  species = sapply(strsplit(rownames(colour_space), split = "-"), "[", 1),
  sex = sapply(strsplit(rownames(colour_space), split = "-"), "[", 2),
  centr_dists = centr_dists
  )

# add distance to centroid onto threat matrix
threat_centr <- threat_matrix |> 
  full_join(centr_dists, by = join_by("binomial_name" == "species"))

# remove duplicates based on second-order code
threat_centr_clean <- threat_centr |> 
  distinct(second_ord_code, binomial_name, sex, .keep_all = TRUE)

# boxplot of centroid distances by threat type
threat_centr_clean |> 
  filter(
    !is.na(sex)
  ) |> 
  ggplot(aes(x = ex_driver, y = centr_dists, fill = ex_driver)) + 
  geom_boxplot() + 
  facet_grid(~ sex)

# ANOVA to check if mean centroid distances of each extinction driver are different
an_m <- aov(centr_dists ~ ex_driver, data = threat_centr_clean[threat_centr_clean$sex == "M", ])
summary(an_m)
# post-hoc Tukey test
lm_m <- lm(centr_dists ~ ex_driver, data = threat_centr_clean[threat_centr_clean$sex == "M", ])
av_m <- aov(lm_m)
TukeyHSD(av_m)

log_mod_f <- lm(centr_dists ~ ex_driver - 1, data = threat_centr_clean[threat_centr_clean$sex == "F", ])
summary(log_mod_f)









# Junk ------------------------------------

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
