# PC-specific threat abatement using MCMCglmm ----
# Robert MacDonald

# Clear environment
rm(list=ls())

# Load libraries ----
library(dplyr)
library(ggplot2)
library(MCMCglmm)


## EDITABLE CODE ##
# Use latest IUCN assessment data or use most recent assessment data pre- specified cutoff year?
latest <- TRUE
# If not using latest assessment data, specify a cutoff year. Set to NULL if using latest.
cutoff_year <- NULL
# Clade to focus on ("Aves", "Neognaths", "Neoaves", "Passeriformes")
clade <- "Aves"
# Colour space to use (e.g. "lab" [CIELAB], "srgb" [sRGB], "xyz" [TCSxyz], "jndxyzlum" [JND with xyz and luminance])
space <- "lab"
# Work with male or female colour data?
sex <- "M"

# Load data ----

# load colour pattern space (created in Chapter 1 - patch-pipeline)
colspace_path <- paste0("G:/My Drive/patch-pipeline/2_Patches/3_OutputData/", clade, "/2_PCA_ColourPattern_spaces/1_Raw_PCA/", clade, ".matchedsex.patches.250716.PCAcolspaces.rds")
colour_space <- readRDS(colspace_path)[[space]][["x"]]

# Load raw species extinction matrices (1=extant, 0=extinct) for 100% abatement scope
extinctions = read.csv(
  here::here("01_input_data", "all_pamat_all_threats_1000_17122024.csv"), 
  row.names = 1
) # 100% scope

# Load phylogeny (distribution of Hackett backbone BirdTree trees)
phy <- ape::read.tree(
  here::here(
    "01_input_data", "First10_AllBirdsHackett1.tre"
  )
)


# Analysis ----

# get species and sex with colourspace values
colour_space_sppsex <- data.frame(species = sapply(strsplit(rownames(colour_space), split = "-"), "[", 1),
                                  sex = sapply(strsplit(rownames(colour_space), split = "-"), "[", 2), 
                                  colour_space)

# trim extinctions data to only species inluded in colour dataset
ext <- extinctions[, colnames(extinctions) %in% colour_space_sppsex$species]

# trim colour space data to only species included in extinctions dataset
colour_space_sppsex <- colour_space_sppsex[colour_space_sppsex$species %in% colnames(ext), ]

# remove the 'full' (extant) row which contains total species count
ext = ext[rownames(ext) != "full", ]

# get the number of species
n_species <- ncol(ext)

# Extract the threat code (e.g., 'hab', 'exp', 'none') from the rownames
codes = stringr::str_split_i(as.character(rownames(ext)), "_", 1)
threats <- unique(codes)[unique(codes) != "none"]

# Get column names (species), excluding the scenario name (X) and the new Code column
species = colnames(ext)

# Group by the threat code and sum the species status across all 1000 iterations for each species
# The resulting matrix 'totsum' shows, for each species (columns), the number of times it remained EXTANT 
# in the 1000 simulations of each abatement scenario (rows)
survival_numbers <- lapply(
  c("none", threats),
  function(code){
    row_nums <- which(codes == code)
    dat <- ext[row_nums, ]
    dat <- as.data.frame(t(colSums(dat)))
    dat$code <- code
    dat <- dat[, c(length(dat), 1:(length(dat)-1))]
    return(dat)
  }
)
survival_numbers <- do.call(rbind, survival_numbers)


# Calculate Extinctions Averted per Species, per Threat

# We calculate the number of times a species remained extant under the 'none' (no abatement) scenario
# then we subtract this number from the number of times a species remained extant under each
# of the abatement scenarios
# This gives us the number of times a species was extant in each scenario ON TOP OF the number
# of times it survives under no abatement
# Effectively this is P(survival under abatement) - P(survival under no abatement)
# which we could call the extra survival probability added by the abatement
# A negative number means the species is MORE likely to go extinct under the abatement scenario compared
# to the baseline scenario, a positive number means species are LESS likely to go extinct under 
# specific abatement than in the baseline scenario 
# Calculation: Scenario_Count - Baseline_Count
# The baseline 'none' row is repeated 6 times - one for each threat
averted_exts <- survival_numbers
averted_exts[averted_exts$code != "none", -1] <- survival_numbers[survival_numbers$code != "none", -1] -
  rbind(survival_numbers[survival_numbers$code == "none", -1],
        survival_numbers[survival_numbers$code == "none", -1],
        survival_numbers[survival_numbers$code == "none", -1],
        survival_numbers[survival_numbers$code == "none", -1],
        survival_numbers[survival_numbers$code == "none", -1],
        survival_numbers[survival_numbers$code == "none", -1])

# Transpose the data frame for merging (to have species as rows and threat codes as columns)
averted_exts_t = t(averted_exts) |> 
  as.data.frame() |> 
  tibble::rownames_to_column("species")
# Assign the first row (the codes) as the column names
colnames(averted_exts_t) = c("species", averted_exts_t[1, -1]) 
# Remove the first row (the old column names) and the 'all' and 'res' codes (only keeping the 6 drivers and 'none')
averted_exts_t = averted_exts_t[-1, ]
# Get the threat column names
threat_cols <- colnames(averted_exts_t)[colnames(averted_exts_t) %in% threats]
# Get the INVERSE of the number of times species extinction was averted
num_extincts <- data.frame(
  species = averted_exts_t$species,
  1000 - apply(averted_exts_t[, threat_cols], 2, as.numeric)
  )
# Conceptually, this is essentially centring the aversion data around 1000
# So if the number for a particular species in a specific abatement scenario is lower than 1000,
# that means it's more likely to be SAVED from extinction by the abatement, but if the number is
# higher than 1000 it's more likely to go extinct under the abatement (compared to the baseline 
# scenario)
# So later when we run the MCMCglmm, a positive estimate of each scenario means that threat (e.g. exp) 
# THREATENS the positive values of the PC, while a negative estimate means that threat threatens the
# negative values of that PC
# So, for example, if a species has a positive association of exp (hunting) with PC1, this means
# that hunting threatens light birds more than dark birds, because luminance is positively
# loaded on PC1

# Merge the PCA colourspace scores (PCs 1-7 only) with the extinction numbers data
dat <- merge(colour_space_sppsex[, c("species", "sex", paste0("PC", 1:3))], num_extincts, by = "species")

# Run MCMCglmm ----

# set seed for reproducibility
set.seed(42)

#G priors for each response variable
prior = list(G = list(G1 = list(V = diag(3), nu = 2)))

# trim the tree to match the dataset
tree <- phy[[1]]
phy <- lapply(phy, drop.tip, tip = setdiff(tree$tip.label, dat$species))


# Set up a dummy run

i <- 1 # this is arbitrary

# pull a tree from the distribution (number i - starts at one)
tree <- phy[[i]] 

# test if tree is ultrametric (all tips equidistant from root)
is.ultrametric(tree)

# check if tree is rooted
is.rooted(tree)

# invert covariance matrix for use by MCMCglmm 
ainv <- inverseA(tree, scale = T)$Ainv 

# Dummy run

# Model Formula:
# cbind(PC1,PC2,PC3) ~ 
#   trait          : Intercept that varies across the 3 response variables (PC1, PC2, PC3).
#   trait:hab + ...: Allows the slope of each threat effect to vary for each response variable (PC1, PC2, PC3).
#   -1             : Suppresses the global intercept.
mod <- MCMCglmm(
  cbind(PC1, PC2, PC3) ~ trait + trait:hab + trait:exp + trait:cli + trait:inv + trait:dis - 1,
  data = dat, 
  family = rep("gaussian", 3), 
  prior = prior, 
  random = ~ us(trait):species, 
  rcov= ~ us(trait):units,
  ginverse = list(species = ainv),
  nitt=1050, 
  thin=1, 
  burnin=50
  )

# check convergence visually
plot(mod)

# inspect results
summary(mod)

#set up a structure that we'll populate with the real model
Final.mod <- mod 
# using 10 here allows us to select 10 samples from each tree run - we can up this if need be
Final.mod$VCV[((i-1)*10+1):(i*10), ] <- mod$VCV[1:10, ] # [VCV is posterior distrib of covariance matrices]
Final.mod$Sol[((i-1)*10+1):(i*10), ] <- mod$Sol[1:10, ] # [Sol is posterior distrib of MME solutions (???) - includes fixed effects]
Final.mod$Liab[((i-1)*10+1):(i*10), ] <- mod$Liab[1:10, ] # [Liab is posterior distrib of latent variables]

nsamp.l <- nrow(mod$VCV)
units_cols <- colnames(mod$VCV)[grepl(".units", colnames(mod$VCV))]
species_cols <- colnames(mod$VCV)[grepl(".species", colnames(mod$VCV))]
start1.l = list(R = mod$VCV[nsamp.l, units_cols], G = list(G1 = mod$VCV[nsamp.l, species_cols])) 

# Create model file to save
save(
  Final.mod,
  file = here::here(
    "03_output_data", "multiresponse_mod.Rdata"
  )
)


# Set up final model - loop over 10 trees
for(i in 1:10){
  
  # pull the tree from the distribution
  tree <- phy[[i]]
  
  # invert covariance matrix for use by MCMCglmm 
  ainv <- inverseA(tree, scale = T)$Ainv 
  
  # run model
  mod <- MCMCglmm(
    cbind(PC1, PC2, PC3) ~ trait + trait:hab + trait:exp + trait:cli + trait:inv + trait:dis - 1,
    data = dat, 
    family = rep("gaussian", 3), 
    prior = prior, 
    random = ~ us(trait):species, 
    rcov= ~ us(trait):units,
    ginverse = list(species = ainv),
    nitt=1050, # need 10 samples per tree, to fit with the structure we set up
    thin=100, 
    burnin=50
  )
  
  print(i) #print which tree you're on (for your sanity as the loop runs)
  
  # populate final model
  Final.mod$VCV[((i-1)*10+1):(i*10), ]<-mod$VCV[1:10,] # [VCV is posterior distrib of covariance matrices]
  Final.mod$Sol[((i-1)*10+1):(i*10), ]<-mod$Sol[1:10,] # [Sol is posterior distrib of MME solutions (???) - includes fixed effects]
  Final.mod$Liab[((i-1)*10+1):(i*10), ]<-mod$Liab[1:10,] # [Liab is posterior distrib of latent variables]
  
  # set up start of next model
  nsamp.l <- nrow(mod$VCV)
  units_cols <- colnames(mod$VCV)[grepl(".units", colnames(mod$VCV))]
  species_cols <- colnames(mod$VCV)[grepl(".species", colnames(mod$VCV))]
  start1.l = list(R = mod$VCV[nsamp.l, units_cols], G = list(G1 = mod$VCV[nsamp.l, species_cols])) 
  
  # save (for sanity)
  if(i > 99){
    save(
      Final.mod,
      file = here::here(
        "03_output_data", "multiresponse_mod.Rdata"
      )
    )
  }
  
}

plot(Final.mod)
summary(Final.mod)


# Select only significant threats and rerun
# Only exp and cli significant
i <- 1 # this is arbitrary

# pull a tree from the distribution (number i - starts at one)
tree <- phy[[i]] 

# invert covariance matrix for use by MCMCglmm 
ainv <- inverseA(tree, scale = T)$Ainv 

# Dummy run

# Model Formula:
# cbind(PC1,PC2,PC3) ~ 
#   trait          : Intercept that varies across the 3 response variables (PC1, PC2, PC3).
#   trait:hab + ...: Allows the slope of each threat effect to vary for each response variable (PC1, PC2, PC3).
#   -1             : Suppresses the global intercept.
mod.2 <- MCMCglmm(
  cbind(PC1, PC2, PC3) ~ trait + trait:exp + trait:cli - 1,
  data = dat, 
  family = rep("gaussian", 3), 
  prior = prior, 
  random = ~ us(trait):species, 
  rcov= ~ us(trait):units,
  ginverse = list(species = ainv),
  nitt=1050, 
  thin=1, 
  burnin=50
)

# check convergence visually
plot(mod.2)

# inspect results
summary(mod.2)

#set up a structure that we'll populate with the real model
Final.mod.2 <- mod.2 
# using 10 here allows us to select 10 samples from each tree run - we can up this if need be
Final.mod.2$VCV[((i-1)*10+1):(i*10), ] <- mod.2$VCV[1:10, ] # [VCV is posterior distrib of covariance matrices]
Final.mod.2$Sol[((i-1)*10+1):(i*10), ] <- mod.2$Sol[1:10, ] # [Sol is posterior distrib of MME solutions (???) - includes fixed effects]
Final.mod.2$Liab[((i-1)*10+1):(i*10), ] <- mod.2$Liab[1:10, ] # [Liab is posterior distrib of latent variables]

nsamp.l <- nrow(mod.2$VCV)
units_cols <- colnames(mod.2$VCV)[grepl(".units", colnames(mod.2$VCV))]
species_cols <- colnames(mod.2$VCV)[grepl(".species", colnames(mod.2$VCV))]
start1.l = list(R = mod.2$VCV[nsamp.l, units_cols], G = list(G1 = mod.2$VCV[nsamp.l, species_cols])) 

# Create model file to save
save(
  Final.mod.2,
  file = here::here(
    "03_output_data", "multiresponse_mod.2.Rdata"
  )
)

# let's time it
start_time <- Sys.time()

# Set up final model - loop over 10 trees
for(i in 1:10){
  
  # pull the tree from the distribution
  tree <- phy[[i]]
  
  # invert covariance matrix for use by MCMCglmm 
  ainv <- inverseA(tree, scale = T)$Ainv 
  
  # run model
  mod.2 <- MCMCglmm(
    cbind(PC1, PC2, PC3) ~ trait + trait:exp + trait:cli - 1,
    data = dat, 
    family = rep("gaussian", 3), 
    prior = prior, 
    random = ~ us(trait):species, 
    rcov= ~ us(trait):units,
    ginverse = list(species = ainv),
    nitt=1050, # need 10 samples per tree, to fit with the structure we set up
    thin=100, 
    burnin=50
  )
  
  print(i) #print which tree you're on (for your sanity as the loop runs)
  
  # populate final model
  Final.mod.2$VCV[((i-1)*10+1):(i*10), ]<-mod.2$VCV[1:10,] # [VCV is posterior distrib of covariance matrices]
  Final.mod.2$Sol[((i-1)*10+1):(i*10), ]<-mod.2$Sol[1:10,] # [Sol is posterior distrib of MME solutions (???) - includes fixed effects]
  Final.mod.2$Liab[((i-1)*10+1):(i*10), ]<-mod.2$Liab[1:10,] # [Liab is posterior distrib of latent variables]
  
  # set up start of next model
  nsamp.l <- nrow(mod.2$VCV)
  units_cols <- colnames(mod.2$VCV)[grepl(".units", colnames(mod.2$VCV))]
  species_cols <- colnames(mod.2$VCV)[grepl(".species", colnames(mod.2$VCV))]
  start1.l = list(R = mod.2$VCV[nsamp.l, units_cols], G = list(G1 = mod.2$VCV[nsamp.l, species_cols])) 
  
  # save (for sanity)
  if(i > 99){
    save(
      Final.mod.2,
      file = here::here(
        "03_output_data", "multiresponse_mod.2.Rdata"
      )
    )
  }
  
}

# end time
end_time <- Sys.time()

save(
  Final.mod.2,
  file = here::here(
    "03_output_data", "multiresponse_mod.2.Rdata"
  )
)

plot(Final.mod.2)
summary(Final.mod.2)


# Re-load model (all predictors)
load(
  here::here(
    "03_output_data", "multiresponse_mod.Rdata"
  )
)

plot(Final.mod)
summary(Final.mod)




#### Plotting ----

# Based on code from Stewart et al 2025 Nat Ecol Evol - script "figure_4a_multiresponse_posterior_plots.R"

# Function to get posterior values of interest in a format for plotting

gen_posts = function(run, n_pcs) {
  
  # Extract posterior samples for the fixed effects
  # run$Sol contains the posterior samples for the fixed effects (Solutions)
  posts = as.data.frame(run$Sol[])
  
  # set intercept names
  intercept_cols <- paste0("traitPC", 1:n_pcs)
  
  # convert o long format to give a single posterior sample for a specific term
  # we also remove the intercepts (traitPCX) as we don't care about these
  long = posts %>%
    dplyr::select(- intercept_cols) |> 
    pivot_longer(cols = colnames(posts[(n_pcs + 1):ncol(posts)]))
  
  # Get the PC and threat names from the results
  long$PC = as.factor(gsub("trait", "", as.factor(str_split_i(long$name, ":", 1))))
  long$threat = as.factor(str_split_i(long$name, ":", 2))
  
  # Order the threat scenarios for plotting purposes
  order = c("hab", "exp", "cli","inv","dis")
  long$threat = factor(long$threat, levels = c(order), ordered = TRUE)
  
  # Kerry's original function reverses the order of PC1 as it makes more sense to have high body size
  # haveing high PC1 values. We don't need to do this
  # long[long$PC=="PC1", "value"]=long[long$PC=="PC1", "value"]*(-1)
  
  # Determine significance - Kerry uses a pMCMC < 0.1 significance level, not sure why not a 0.05
  sig_var = rownames(summary(run)$solutions)[summary(run)$solutions[, "pMCMC"] < 0.1]
  # Create binary flag - 1 if significant, 0 if not
  long$sig=ifelse(long$name %in% sig_var, 1, 0)
  
  # get summary groupings to determine if mean effect is positive or not
  sum = long |> 
    group_by(threat, PC) |> 
    summarise(
      pos = as.numeric(mean(value) > 0),
      sig = unique(sig)
      ) |> 
    as.data.frame()
  
  # Deine logic for colouring
  # -1 = significant and negative
  # 0 = non-significant (will be grey)
  # 1 = significant and positive
  sum$col=ifelse(sum$sig==1 & sum$pos==0, -1, 0)
  sum$col=ifelse(sum$sig==1 & sum$pos==1, 1, sum$col)
  
  # merge summary stats back into main long df
  long = merge(long, sum, by = c("PC", "threat", "sig"))
  return(long)
}

# Load model
load(
  here::here(
    "03_output_data", "multiresponse_mod.Rdata"
  )
)

# get posterior estimates in format for plotting
posterior_formatted <- gen_posts(Final.mod, n_pcs = 3)

# Create labels for plot facets
threat_labs <- c(
  hab = "Habitat", inv = "Invasive", cli = "Climate", exp = "Hunting", dis = "Disturbance"
)

# define custom palette - make non-significant results grey
cols =   c(
    "not_sig" = "grey90", 
    "inv" = brewer.pal(n=6, "Oranges")[3], 
    "cli" = brewer.pal(n=6, "Greens")[3], 
    "exp" = brewer.pal(n=6, "Reds")[3], 
    "hab" = brewer.pal(n=6, "Purples")[3],
    "dis" = brewer.pal(n=6, "Blues")[3]
    )


# create fill variable - this determines which colour each density plot gets
# first get threat names
posterior_formatted$other_col = str_split_i(posterior_formatted$name, ":", 2)
# if effect is non-significant, make value "not_sig"
posterior_formatted[posterior_formatted$sig == 0, "other_col"] = "not_sig"
# order the colour factor to match the 'cols' vector defined above
posterior_formatted$other_col = factor(posterior_formatted$other_col, ordered = TRUE, 
                         levels = c("not_sig", "inv", "cli", "exp", "hab"))

# Plot
p1 = ggplot(posterior_formatted, aes(x = value)) +
  # Apply the custom color palette defined above
  scale_fill_manual(values = cols) +
  
  # Create the Grid: Rows = Threats, Columns = pPCs
  # switch="y" moves the Y-axis labels to the left/right for better readability
  facet_grid(threat ~ PC, labeller = labeller(threat = threat_labs), switch = "y") +
  
  # Minimal theme with heavy customization
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),       # Remove Y-axis numbers (density height is relative)
    axis.text.x = element_text(size = 14), # Increase X-axis text size
    axis.title.y = element_blank(),      # Remove Y-axis title
    legend.title = element_blank(),      # No legend title
    legend.position = "none",            # Hide legend (colors are intuitive via row labels)
    panel.grid.major = element_blank(),  # Remove grid lines
    panel.grid.minor = element_blank(),
   # strip.text.y.left = element_text(size=16, angle=0, hjust=0.9), # Adjust Row Labels
  # strip.placement.x = "outside",
  #  strip.text.x = element_text(size=18, hjust=0.4, vjust=6),      # Adjust Column Labels
   # text = element_text(family="Leelawadee UI Semilight", size=18), # Specific font
    plot.margin = unit(c(2, 0.2, 0.2, 0.2), "cm"), # add thick margin at top for images of PC loadings
    panel.background = element_rect(fill="white", colour="NA"),
    plot.background = element_rect(fill="white", colour="NA")
  ) +
  
  # DRAWING CUSTOM AXES (Manual aesthetic choices)
  # These segments draw a "box" or axis lines manually.
  # WARNING: The yend=1785 is hardcoded! If your data density is lower/higher, 
  # these lines might look wrong or disappear.
  geom_segment(x=0, xend=0, y=0, yend=110, lwd=0.3) + # Vertical line at 0
#  geom_segment(x=-0.00123, xend=-0.00123, y=0, yend=100, col="black", lwd=0.3) + # Left border
#  geom_segment(x=0.001745, xend=0.001745, y=0, yend=100, col="black", lwd=0.3) + # Right border
#  geom_segment(x=-0.00123, xend=0.001745, y=1785, yend=100, col="black", lwd=0.3) + # Top border
  
  # PLOTTING THE DISTRIBUTIONS
  # 1. The Outline (Black line)
  geom_density(aes(alpha = 1), fill = NA, colour = "black") +
  # 2. The Fill (Colored based on significance/threat)
  geom_density(aes(fill = other_col), 
               colour = "NA", 
               alpha = 0.8) +
  
  labs(x = "Posterior values")
  
  # Zoom settings
  # clip="off" allows elements (like text) to exist outside the plot area
#  coord_cartesian(xlim = c(-0.0011, 0.0015), ylim = c(50, 1500), clip = "off") +
  
  # Manual X-axis breaks
 # scale_x_continuous(breaks = c(-0.001, 0, 0.001))

p1

