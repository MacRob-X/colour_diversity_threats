## Calculate image feature metrics using colour mapping (recolorize) and adjacency analysis (recolorize, pavo)
## Following https://doi.org/10.1093/evolut/qpaf242 
## https://cran.r-project.org/web/packages/recolorize/vignettes/Introduction.html
## https://github.com/christophercooney/avian-colour-loss/blob/main/code/R/01_Recolorize_image_feature_metrics_v9.R


# Clear environment
rm(list=ls())

# Load libraries ----
library(recolorize)
library(ggplot2)
library(scico)

# Read in cone catch functions -----
# maps pixels in NEF image to cone stimulation values
source("G:/My Drive/patch-pipeline/3_SharedScripts/00_Cone_mapping_functions_v8.R") # pixel vals need to be in range 0-100

## CHOOSE PARAMETERS ##
# Remove extinct (EX) species?
remove_extinct <- TRUE

# Load data ----

# Load patch sRGB data
orig_srgb <- readRDS("G:/My Drive/patch-pipeline/2_Patches/3_OutputData/Aves/1_RawColourspaces/Aves.matchedsex.patches.250716.prePCAcolspaces.rds")[["sRGB"]]

# Load patch Lab data (for luminance values)
lab <- readRDS("G:/My Drive/patch-pipeline/2_Patches/3_OutputData/Aves/1_RawColourspaces/Aves.matchedsex.patches.250716.prePCAcolspaces.rds")[["lab"]]

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

# load colour space
clade <- "Aves"
colspace_path <- paste0("G:/My Drive/patch-pipeline/2_Patches/3_OutputData/", clade, "/2_PCA_ColourPattern_spaces/1_Raw_PCA/", clade, ".matchedsex.patches.250716.PCAcolspaces.rds")
colour_space <- readRDS(colspace_path)[["lab"]][["x"]]
colour_space <- as.data.frame(colour_space)

# load UMAP space
umap_path <- paste0("G:/My Drive/patch-pipeline/2_Patches/3_OutputData/", clade, "/2_PCA_ColourPattern_spaces/2_UMAP/", clade, ".matchedsex.patches.nn.25.mindist.0.1.lab.UMAP.rds")
umap <- readRDS(umap_path)[["layout"]]
colnames(umap) <- c("UMAP1", "UMAP2")

# Analysis ----

## Prepare data ----

# get long format sRGB data
colnames(orig_srgb) <- gsub("sRGB.", "", colnames(orig_srgb))
orig_srgb$species <- sapply(
  strsplit(
    rownames(orig_srgb),
    split = "-"),
  "[",
  1
)
orig_srgb$sex <- sapply(
  strsplit(
    rownames(orig_srgb),
    split = "-"),
  "[",
  2
)

orig_srgb_long <- tidyr::pivot_longer(
  orig_srgb,
  cols = r.bel:b.thr,
  names_to = "channel.region",
  values_to = "value", 
)
colsave <- colnames(orig_srgb_long)
split <- strsplit(orig_srgb_long$channel.region, ".", fixed = T)
split <- do.call(rbind, split)
orig_srgb_long <- cbind(orig_srgb_long, split)
colnames(orig_srgb_long) <- c(colsave, "channel", "region")
orig_srgb_long <- orig_srgb_long[, which(!colnames(orig_srgb_long) %in% "channel.region")]
orig_srgb_long <- tidyr::pivot_wider(
  orig_srgb_long,
  names_from = "channel",
  values_from = "value"
)
orig_srgb_long <- as.data.frame(orig_srgb_long)
rm(colsave)

# restrict patch data to only species in species list (with male/female colour data and threat data)
orig_srgb_long <- orig_srgb_long[orig_srgb_long$species %in% spp_list, ]

rownames(orig_srgb_long) <- paste(
  orig_srgb_long$species, 
  orig_srgb_long$sex, 
  orig_srgb_long$region,
  sep = "-")

# get format for colour clustering
srgb <- orig_srgb_long[, c("r", "g", "b")]

## Identify broad colour classes ----

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
col_names <- data.frame(
  parameter_name = paste0("p_", 1:nrow(col_names)),
  class_number = 1:nrow(col_names),
  col_names
  )
# assign colour names (note: subjective)
col_names$colour_name_rxm <- c("Dark Cerulean",
                      "Teal-Grey",
                      "Charcoal Blue",
                      "Blue Slate",
                      "Moss Green",
                      "Dark Charcoal Brown",
                      "Dark Taupe",
                      "Light Taupe",
                      "Mid Beige",
                      "Deep Terracotta",
                      "Cinnamon Terracotta",
                      "Olive Brown",
                      "Golden Olive",
                      "Gold",
                      "Ochre"
                      )
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



## Impose colour classes on patch data ----

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
colour_class_filename <- if(remove_extinct){
  "colour_classes_remove_EX.RDS"
} else if(!remove_extinct){
  "colour_classes_keep_EX.RDS"
}
col_classes <- readRDS(
  here::here(
    "03_output_data", "07_recolorize_image_feature_metrics", colour_class_filename
  )
)


# function to recolorize and save as PNG
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


# run across each patch image
recol_imgs <- lapply(
  1:length(patch_imgs),
  function(img_num){
    
    cat('\r', "Image ", img_num, "of", length(patch_imgs))
    
    spp_sex <- names(patch_imgs[img_num])
    
    image <- patch_imgs[[img_num]]
    
    recolorized_imgs <- run_rc(
      image,
      col_classes = col_classes,
      write_to_file = TRUE,
      spp_sex = spp_sex,
      #recol_outpath = here::here("04_output_plots", "07_recolorize_image_feature_metrics", "recolorize_patch_images")
      recol_outpath = "C:/Users/bop23rxm/Documents/recolorize_patch_images"
    )
    
    return(recolorized_imgs)
  }
)

# add species/sex names
names(recol_imgs) <- names(patch_imgs)

# save for later
recol_filename <- if(remove_extinct){
  "recolorize_images_remove_EX.RDS"
} else if(!remove_extinct){
  "recolorize_images_keep_EX.RDS"
}
saveRDS(
  recol_imgs,
  here::here(
    "03_output_data", "07_recolorize_image_feature_metrics", recol_filename
  )
)


## Use adjacency and boundary strength analysis to derive image feature metrics ----

# reload recolorized patch images if needed
recol_filename <- if(remove_extinct){
  "recolorize_images_remove_EX.RDS"
} else if(!remove_extinct){
  "recolorize_images_keep_EX.RDS"
}
recol_imgs <- readRDS(
  here::here(
    "03_output_data", "07_recolorize_image_feature_metrics", recol_filename
  )
) 

# load colour class numbers and names
colour_filename <- if(remove_extinct){
  "colour_classes_remove_EX.csv"
} else if(!remove_extinct){
  "colour_classes_keep_EX.csv"
}
col_names <- read.csv(
  here::here(
    "03_output_data", "07_recolorize_image_feature_metrics", colour_filename
  )
)



## Compute metrics -----

### Mean, max, variance luminance ----
mean_lum <- rowMeans(lab[, grep("L.", colnames(lab))]) # subset to columns with luminance values
max_lum <- apply(
  lab[, grep("L.", colnames(lab))],
  1,
  max
)
var_lum <- apply(
  lab[, grep("L.", colnames(lab))],
  1,
  var
)
# combine into single dataframe
lum <- data.frame(
  mean_lum,
  max_lum,
  var_lum
)

### Contrast between dorsal and ventral mean luminance (countershading) ----
# for dorsal, use patch luminance for crown, nape, mantle, rump
mean_dorsal <- rowMeans(lab[, which(colnames(lab) %in% c("L.cro", "L.nap", "L.man", "L.rum"))])
# for ventral, use throat, breast, belly
mean_ventral <- rowMeans(lab[, which(colnames(lab) %in% c("L.thr", "L.bre", "L.man", "L.bel"))])
# get difference in dorsal and ventral luminance
# positive values indicate darker dorsal, lighter ventral (i.e. countershading)
# negative values indicate lighter dorsal, darker ventral (i.e. reverse countershading)
countershading <- data.frame(countershading = mean_ventral - mean_dorsal)

### Colourfulness  ----

# number of colour classes

# function to count colour classes in a single recolorized image
count_col_classes <- function(recolorized_img){
  
  n_cols <- length(
    unique(
      as.vector(recolorized_img$pixel_assignments)
    )
  )
  
  return(n_cols)
  
}


# apply across all images
num_colour_classes <- unlist(
  lapply(
    recol_imgs,
    count_col_classes
  )
)

# mean pairwise distance between patches in Lab space
lab_rn <- lab
lab_rn$row_names <- rownames(lab)
lab_long <- tidyr::pivot_longer(
  lab_rn,
  !row_names,
  names_to = "ch.patch",
  values_to = "value"
)
lab_long$ch <- sapply(
  strsplit(
    lab_long$ch.patch,
    split = ".",
    fixed = TRUE
    ),
  "[",
  1
)
lab_long$patch <- sapply(
  strsplit(
    lab_long$ch.patch,
    split = ".",
    fixed = TRUE
  ),
  "[",
  2
)
lab_long <- tidyr::pivot_wider(
  lab_long,
  id_cols = c("row_names", "patch"),
  names_from = ch,
  values_from = value
  
)
mean_var_pairwise_dists <- lapply(
  unique(lab_long$row_names),
  function(spp_sex){
    
    dat <- as.matrix(lab_long[lab_long$row_names == spp_sex, c("L", "a", "b")])
    
    pairwise_dists <- vegan::vegdist(dat, method = "euclidean")
    
    mean_pairwise_dist <- mean(pairwise_dists)
    var_pairwise_dist <- var(pairwise_dists)
    
    res <- data.frame(
      mean_pairwise_dist,
      var_pairwise_dist
    )
    return(res)
    
  }
)
names(mean_var_pairwise_dists) <- unique(lab_long$row_names)
mean_var_pairwise_dists <- do.call(rbind, mean_var_pairwise_dists)

colourfulness <- as.data.frame(num_colour_classes)
rownames(colourfulness) <- gsub(".", "-", rownames(colourfulness), fixed = TRUE)
colourfulness <- merge(
  colourfulness,
  mean_var_pairwise_dists,
  by = "row.names"
)
rownames(colourfulness) <- colourfulness$Row.names
colourfulness$Row.names <- NULL

### Proportion of each colour class ----

# function to calculate proportion of a single colour class in a single image
prop_col_class <- function(recolorized_img, col_class_number){
  
  # ensure colour class number is numeric
  col_class_number <- as.numeric(col_class_number)
  
  n_class <- length(
    which(
      as.vector(recolorized_img$pixel_assignments) == col_class_number
    )
  )
  
  prop_class <- n_class / 10 # total number of patches in each patch representation

  return(prop_class)
  
}

# apply across all images for a single colour class
# use class 10 (browny red)
class_10_prop <- unlist(
  lapply(
    recol_imgs,
    prop_col_class,
    col_class_number = 10
  )
)
hist(class_10_prop)
head(sort(class_10_prop[which(class_10_prop > 0)], decreasing = TRUE), 20)


# apply across all images for all colour classes
class_props_list <- lapply(
  col_names$class_number,
  function(class_num){
    props <- unlist(
      lapply(
        recol_imgs,
        prop_col_class,
        col_class_number = class_num
      )
    )
    return(props)
  }
)
names(class_props_list) <- col_names$parameter_name

# inspect the highest proportion species of a given colour class
col_class_param <- "p_1"
head(
  sort(
    class_props_list[[col_class_param]][which(class_props_list[[col_class_param]] > 0)], 
    decreasing = TRUE
  ), 
  10
)

# convert to dataframe
class_props <- as.data.frame(do.call(cbind, class_props_list))
rownames(class_props) <- gsub(".", "-", rownames(class_props), fixed = TRUE)


# Combine metrics together in one dataframe with PC and UMAP axes


# merge UMAP data to raw colour space
colsp_umap <- merge(
  colour_space,
  umap,
  by = "row.names"
)
rownames(colsp_umap) <- colsp_umap$Row.names
colsp_umap$Row.names <- NULL

# merge to mean and max luminance
colsp_metrics <- merge(
  colsp_umap,
  lum,
  by = "row.names"
)
rownames(colsp_metrics) <- colsp_metrics$Row.names
colsp_metrics$Row.names <- NULL

# merge to countershading
colsp_metrics <- merge(
  colsp_metrics,
  countershading,
  by = "row.names"
)
rownames(colsp_metrics) <- colsp_metrics$Row.names
colsp_metrics$Row.names <- NULL

# merge to colour class proportions
colsp_metrics <- merge(
  colsp_metrics,
  class_props,
  by = "row.names"
)
rownames(colsp_metrics) <- colsp_metrics$Row.names
colsp_metrics$Row.names <- NULL

# merge to colourfulness
colsp_metrics <- merge(
  colsp_metrics,
  colourfulness,
  by = "row.names"
)
rownames(colsp_metrics) <- colsp_metrics$Row.names
colsp_metrics$Row.names <- NULL

# Save as CSV
metrics_filename <- if(remove_extinct){
  paste("handcrafted_colour_metrics", "remove_EX.csv", sep = "_")
} else if(!remove_extinct){
  paste("handcrafted_colour_metrics", "keep_EX.csv", sep = "_")
}
write.csv(
  colsp_metrics,
  here::here(
    "03_output_data", "07_recolorize_image_feature_metrics",
    metrics_filename
  )
)


# Plot metrics on PC/UMAP axes ----
# easier to do parts of this with tidyverse
library(dplyr)


# Plot metric on colour space axes
# choose metric to plot
plot_metr <- "var_lum"
# choose axes to plot
ax_1 <- "UMAP1"
ax_2 <- "UMAP2"
# sort data by metric first (makes plot interpretation easier)
# for countershading, this can be positive or negative, so we sort by the absolute value
# this will plot the most extremely countershaded or reverse-countershaded species on top
# to plot the least countershaded (i.e. smallest absolute countershading values) on top,
# add rev() around order() in the below line
sorted_dat <- colsp_metrics[order(abs(colsp_metrics[[plot_metr]])), ]
metric_plot <- ggplot(sorted_dat, aes(x = get(ax_1), y = get(ax_2), fill = get(plot_metr))) + 
  geom_point(
    shape = 21,
    colour = "transparent"
  ) + 
  labs(x = ax_1, y = ax_2, fill = plot_metr) + 
  scale_fill_continuous(trans = "reverse") +
  theme_minimal() +
  coord_fixed(ratio = 1) +
  guides(fill = guide_colourbar(reverse = TRUE))

# adjust colour palette to suit metric (e.g. countershading needs a divergent palette)
if(plot_metr == "countershading"){
  max_val <- max(abs(colsp_metrics[[plot_metr]]))
  metric_plot <- metric_plot + 
    scale_fill_scico(
      palette = "vik", 
      limits = c(-max_val, max_val),
      trans = "reverse"
      )
} else if(plot_metr %in% c("mean_lum", "max_lum")){
  metric_plot <- metric_plot + 
    scale_fill_viridis_c(
      limits = c(0, 100)
    )
} else if(plot_metr == "colourfulness"){
  metric_plot <- metric_plot + 
    scale_fill_viridis_c()
}
metric_plot




# Plot all proportions of colour classes on axes, faceted by colour class
# Choose axes to plot
ax_1 <- "UMAP1"
ax_2 <- "UMAP2"

# create plot for each colour class (have to do it this way rather than faceting to preserve
# highlighting and ordering non-zero proportion values)
col_class_plots <- lapply(
  col_names[["parameter_name"]],
  function(col_class_param){
    
    # plot (set all zero-value points to be very transparent and below non-zero values)
    # first order by increasing colour class proportion, to make plot interpretation easier (
    # plots the higher values on top)
    sorted_dat <- colsp_metrics[order(colsp_metrics[[col_class_param]]), ]
    class_plot <- ggplot() + 
      geom_point(
        data = colsp_metrics[colsp_metrics[, col_class_param] == 0, ], 
        aes(x = get(ax_1), y = get(ax_2)), 
        shape = 21,
        alpha = 0.1,
        colour = "lightgrey",
        fill = "lightgrey"
      ) + 
      geom_point(
        data = sorted_dat[sorted_dat[[col_class_param]] > 0, ], 
        aes(x = get(ax_1), y = get(ax_2), fill = get(col_class_param)), 
        shape = 21,
      #  size = 2,
        colour = "transparent"
      ) + 
      labs(fill = col_class_param) + 
      scale_fill_continuous(trans = "reverse", limits = c(0, 1)) + 
      theme_minimal() +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      coord_fixed(ratio = 1) +
      guides(fill = guide_colourbar(reverse = TRUE)) # reverse legend direction
    
    return(class_plot)
    
  }
)

# plot all of them on one plot
combined_plot <- ggpubr::annotate_figure(
  ggpubr::ggarrange(
    plotlist = col_class_plots,
    labels = col_names[["parameter_name"]],
    font.label = list(size = 10),
    common.legend = TRUE, legend = "right"
  ),
  bottom = ggpubr::text_grob(ax_1),
  left = ggpubr::text_grob(ax_2, rot = 90)
)

# save plot
class_prop_plot_filename <- if(remove_extinct){
    paste("colour_class_proportions", ax_1, ax_2, "remove_EX.png", sep = "_")
  } else if(!remove_extinct){
    paste("colour_class_proportions", ax_1, ax_2, "keep_EX.png", sep = "_")
    }
ggsave(
  filename = here::here(
    "04_output_plots", "07_recolorize_image_feature_metrics",
    class_prop_plot_filename
  ),
  plot = combined_plot,
  bg = "white",
  width = 21, height = 21, units = "cm"
)


# same for colourfulness

colour_space <- merge(
  colsp_umap,
  colourfulness,
  by = "row.names"
)
ggplot(colour_space[colour_space$colourfulness_res > 2, ], aes(x = PC7, y = PC8, colour = colourfulness_res, fill = colourfulness_res)) + 
  geom_point(alpha = 0.5)
