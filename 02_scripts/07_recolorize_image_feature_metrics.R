## Calculate image feature metrics using colour mapping (recolorize) and adjacency analysis (recolorize, pavo)
## Following https://doi.org/10.1093/evolut/qpaf242 
## https://cran.r-project.org/web/packages/recolorize/vignettes/Introduction.html
## https://github.com/christophercooney/avian-colour-loss/blob/main/code/R/01_Recolorize_image_feature_metrics_v9.R


# Clear environment
rm(list=ls())

# Load libraries ----
library(recolorize)

# Read in cone catch functions -----
# maps pixels in NEF image to cone stimulation values
source("G:/My Drive/patch-pipeline/3_SharedScripts/00_Cone_mapping_functions_v8.R") # pixel vals need to be in range 0-100

## CHOOSE PARAMETERS ##
# Remove extinct (EX) species?
remove_extinct <- TRUE

# Load data ----

# Load patch data
patches <- readRDS("G:/My Drive/patch-pipeline/2_Patches/1_InputData/patches.250716.fixed.rds")

# load list of species with male/female colour data and threat data, with DD spp removed (and
# EX if specified)
spp_list_filename <- if(remove_extinct){
  "species_list_remove_EX.RDS"
} else if(!remove_extinct){
  "species_list_keep_EX.RDS"
}
spp_list <- readRDS(
  here::here(
    "03_output_data", "06_pc_axes_threats", spp_list_filename
  )
)

# Analysis ----

# remove specimens of unknown sex
patches <- patches[patches$sex %in% c("M", "F"), ]

# aggregate region values for sexes within species
patches <- aggregate(cbind(vR, vG, vB, uR, uG, uB) ~ species + sex + region, patches, mean)

# restrict patch data to only species in species list (with male/female colour data and threat data)
patches <- patches[patches$species %in% spp_list, ]

# convert patch data to sRGB
xyz <- fun.camA.ideal.ciexyz.int2poly1.expanded(vR = patches$vR*100, vG = patches$vG*100, vB = patches$vB*100)
srgb <- convertColor(xyz, from = "XYZ", to = "sRGB", to.ref.white = "E") # 'Equal energy'
# 2026-10-14 - need to rerun this with the above Equal energy parameter and see if it fixes 
# discrepancy with patch-pipeline method of generating sRGB values (which goes XYZ->Lab->sRGB, 
# rather than XYZ->sRGB)

# add rownames for species/sex/region ID
rownames(srgb) <- paste(patches$species, patches$sex, patches$region, sep = "-")
colnames(srgb) <- c("r", "g", "b")

## Convert patch data into image array to facilitate input into recolorize functions


# function to find factors of focal number (to create image containing all patch observations)
find_factor_pairs <- function(target_number) {
  
  # Error handling: Ensure input is a positive integer
  if (target_number <= 0 || target_number %% 1 != 0) {
    stop("Please enter a positive integer.")
  }
  
  # 1. Create a sequence from 1 to the square root of the target
  # We only need to check up to sqrt because factors repeat after that.
  possible_factors <- 1:floor(sqrt(target_number))
  
  # 2. Find which numbers divide evenly (Modulus operator %% gives the remainder)
  # If remainder is 0, it is a factor.
  first_factors <- possible_factors[target_number %% possible_factors == 0]
  
  # 3. Calculate the corresponding second factors
  second_factors <- target_number / first_factors
  
  # 4. Combine into a neat data frame
  results <- data.frame(
    Factor_1 = first_factors,
    Factor_2 = second_factors
  )
  
  return(results)
}

fac_pairs <- find_factor_pairs(nrow(srgb))

x <- tail(fac_pairs$Factor_1, 1)
y <- tail(fac_pairs$Factor_2, 1)

r <- matrix(srgb[,1], nrow = x, ncol = y)
g <- matrix(srgb[,2], nrow = x, ncol = y)
b <- matrix(srgb[,3], nrow = x, ncol = y)

global_patch_img <- array(c(r, g, b), dim = c(x, y, 3))

# initial cluster of patch data using recolorise
init_cc <- recolorize(global_patch_img, method = "histogram", bins = 6, color_space = "Lab")

# refine clustering
refin_cc <- recluster(init_cc, cutoff = 30, plot_hclust = T, color_space = "Lab", color_space_fit = "Lab")

# save for later use
colour_class_filename <- if(remove_extinct){
  "colour_classes_remove_EX.RDS"
} else if(!remove_extinct){
  "colour_classes_keep_EX.RDS"
}
saveRDS(
  refin_cc,
  here::here(
    "03_output_data", "07_recolorize_image_feature_metrics", colour_class_filename
  )
)

# plot colour categories
plot_filename <- if(remove_extinct){
  "colour_classes_remove_EX.pdf"
} else if(!remove_extinct){
  "colour_classes_keep_EX.pdf"
}
pdf(
  here::here(
    "04_output_plots", "07_recolorize_image_feature_metrics", plot_filename
  ), 
  width = 12, height = 2)
par(mar=c(1,1,1,1))
cols <- rgb(refin_cc$centers[,1], refin_cc$centers[,2], refin_cc$centers[,3], maxColorValue = 1)
plot(1:length(cols), rep(1,length(cols)), axes = F, pch = 16, col = cols, cex = 6, xlab="", ylab="")
text(1:length(cols), rep(1,length(cols)), labels = 1:16)
text(1:length(cols), rep(1,length(cols)), labels = refin_cc$sizes, pos = 1, offset = 2)
dev.off()

# generate colour category names
# generate colour category names
col_names <- refin_cc$centers
colnames(col_names) <- c("r","g","b")
col_names <- data.frame(parameter.name=paste0("p_", 1:nrow(col_names)), col_names)
col_names$colour_name_rxm <- c("Blue",
                      "Charcoal grey",
                      "Teal",
                      "Dark brown",
                      "Light brown",
                      "Beige",
                      "Sand",
                      "Dark red",
                      "Brick red",
                      "Olive green",
                      "Ochre",
                      "Gold",
                      "Orange")
colour_filename <- if(remove_extinct){
  "colour_classes_remove_EX.csv"
} else if(!remove_extinct){
  "colour_classes_keep_EX.csv"
}
write.csv(
  col_names,
  here::here(
    "03_output_data", "07_recolorize_image_feature_metrics", colour_filename
  )
  , row.names = F
)

## Use adjacency and boundary strength analysis to derive image feature metrics ----

# First need to get patch representation of each species/sex combo as image

# get species, sex and region columns
srgb_form <- as.data.frame(srgb)
srgb_form$species <- sapply(
  strsplit(
    rownames(srgb_form),
    split = "-"),
  "[",
  1
)
srgb_form$sex <- sapply(
  strsplit(
    rownames(srgb_form),
    split = "-"),
  "[",
  2
)
srgb_form$region <- sapply(
  strsplit(
    rownames(srgb_form),
    split = "-"),
  "[",
  3
)

# convert region to factor with correct order
## WHERE TO PUT WINGS? NO REAL CORRECT ANSWER
regions_ord <- c(
  "thr", "bre", "bel", # ventral
  "cov", "fli", # wings
  "cro", "nap",  "man", "rum", "tai" # dorsal 
  )
srgb_form$region <- factor(srgb_form$region, levels = regions_ord)
# order by sex, then species, then region
srgb_form <- with(srgb_form, srgb_form[order(srgb_form$sex, srgb_form$species, srgb_form$region), ])

# split into individual dataframes for each species/sex combo
srgb_ind <- split(
  srgb_form,
  list(srgb_form$species, srgb_form$sex),
  drop = TRUE
)


# reshape each dataframe into image array
patch_imgs <- lapply(
  srgb_ind,
  function(df){
    # message(df$species[1], "-", df$sex[1])
    spp_sex_cols <- which(colnames(df) %in% c("species", "sex", "region"))
    mat <- as.matrix(df[, -spp_sex_cols])
    r <- matrix(mat[, "r"], nrow = 5, ncol = 2)
    g <- matrix(mat[, "g"], nrow = 5, ncol = 2)
    b <- matrix(mat[, "b"], nrow = 5, ncol = 2)
    img <- array(c(r, g, b), dim = c(5, 2, 3))
  }
)

# Impose global colour classes on each patch image

# load colour classes
col_classes <- readRDS(
  here::here(
    "03_output_data", "07_recolorize_image_feature_metrics", colour_class_filename
  )
)


# function to recolorize and run adjacncy analysis
run_rc <- function(patch_img, col_classes, write_to_file = FALSE, spp_sex = NULL, recol_outpath = NULL) {
  
  try ({
    # write sRGB version
    if(write_to_file){
      png::writePNG(
        patch_img,
        target = paste0(recol_outpath, "/", "original", "/", gsub(".", "-", spp_sex, fixed = T), ".png")
      )
    }
    img <- suppressMessages(
      imposeColors(patch_img, centers = col_classes$centers, color_space = "Lab", plotting = F)
      )
    
    if(write_to_file){
      recolorize_to_png(img, filename = paste0(recol_outpath, "/", "recolorized", "/", gsub(".", "-", spp_sex, fixed = T), ".png"))
    }
    
  })
  
  return(img)
}


recol_imgs <- lapply(
  1:length(patch_imgs),
  function(img_num){
    
    cat('\r', "Image ", img_num, "of", length(patch_imgs))
    
    spp_sex <- names(patch_imgs[img_num])
    
    image <- patch_imgs[[img_num]]
    
    run_rc(
      image,
      col_classes = col_classes,
      write_to_file = TRUE,
      spp_sex = spp_sex,
      recol_outpath = "C:/Users/bop23rxm/Documents/recolorize_patch_images"
    )
    
  }
)
