# Calculate disproportionality of colour pattern diversity loss under abatement 
# of each Extinction Driver relative to species richness loss
# Based on code in extended_figure6_intercepts_plot.R from Stewart et al 2025 Nat Ecol Evol
# Robert MacDonald
# 2026-11-10

# Clear environment
rm(list=ls())


# Load libraries ----
library(dplyr)
library(ggplot2)
library(lme4)
library(RColorBrewer)

## EDITABLE CODE ##
# Use latest IUCN assessment data or use most recent assessment data pre- specified cutoff year?
latest <- TRUE
# If not using latest assessment data, specify a cutoff year. Set to NULL if using latest.
cutoff_year <- NULL
# Clade to focus on ("Aves", "Neognaths", "Neoaves", "Passeriformes")
clade <- "Aves"
# Colour space to use (e.g. "lab" [CIELAB], "srgb" [sRGB], "xyz" [TCSxyz], "jndxyzlum" [JND with xyz and luminance])
space <- "lab"
# remove extinction drivers with negligible contribution to colour pattern diversity loss?
remove_negligible <- FALSE

# Load data ----

# Load colour pattern diversity and species richness loss in each of 8000 simulations
# 1000 sims per scenario
div_filename <- "colpat_div_loss_sims_ext_drivers.csv"
div <- read.csv(
  here::here(
    "03_output_data", "04_threat_abatements", div_filename
  )
)

# Prepare data ----

# add iteration number (i.e. the Monte Carlo simulation/MCMC iteration number - 1 to 1000)
div <- div %>% 
  mutate(
    iteration = stringr::str_split_i(X, "_", 2)
  )

# discard unneeded columns and rename others
div <- div %>% 
  select(
    X, code, iteration, cd_loss_avoided_abs, sr_loss_avoided_abs
  ) %>% 
  rename(
    scenario = X,
    cd_loss_av = cd_loss_avoided_abs,
    sr_loss_av = sr_loss_avoided_abs
  )

# Remove 'all' and 'none' scenario - we're not interested in abatement of all threats, 
# or of no threats, only the individual drivers
# (following Stewart et al 2025)
div <- div %>% 
  filter(
    !(code %in% c("all", "none"))
  )

# remove negligible extinction drivers
if(remove_negligible){
  # Remove Climate Change, Invasive Species and Pollution abatement scenarios - calculation of Cohen's D
  # for driver-specific abatement scenarios shows that these have negligible contribution to colour-pattern
  # diversity loss
  div <- div %>% 
    filter(
      !(code %in% c("cli", "pol"))
    )
}


# inspect summary
div %>% 
  group_by(code) %>% 
  summarise(
    mean_cd_loss_av = mean(cd_loss_av),
    sd_cd_loss_av = sd(cd_loss_av),
    mean_sr_loss_av = mean(sr_loss_av),
    sd_sr_loss_av = sd(sr_loss_av)
  )

# Data analysis ----

# Run full model, with all individual drivers, species loss, iteration, and interactions between each 
# driver and species loss as predictors and colour pattern diversity loss as response
lme_full <- lmer(cd_loss_av ~ code * sr_loss_av + (1 | iteration), data = div, na.action = "na.fail")

# identify best model from all combinations of explanatory variables/interactions (using MuMIn::dredge)
MuMIn::dredge(lme_full)
# best model doesn't include interactions

# Run best model
lme <- lmer(cd_loss_av ~ code + sr_loss_av + (1 | iteration), data = div)
# note that the baseline value for code is taken as the first alphabetically
# inspect
summary(lme)

# Pairwise posthoc cmoparison between predictors
em <- emmeans::emmeans(lme, specs = ~ code, pbkrtest.limit = 10000)
# inspect
summary(em)

# Plotting ----

# extract intercepts for each extinction driver
intercepts <- summary(lme)[["coefficients"]][, "Estimate"]
intercepts[2:length(intercepts)-1] <- intercepts[1] + intercepts[2:length(intercepts)-1] # to compare to baseline intercept
names(intercepts) <- gsub("code", "", names(intercepts))
names(intercepts)[1] <- sort(unique(div$code))[[1]] # get first code alphabetically for baseline
intercepts <- intercepts %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("code")
colnames(intercepts) <- c("code", "value")

# define order of explanatory variables
if(remove_negligible){
  # for restricted model (without negligible drivers)
  var_order <- c("hab", "exp", "inv", "dis", "sr_loss_av")
  var_labels <- c("Habitat", "Hunting", "Invasive", "Disturbance", "sr_loss_av")
} else if(!remove_negligible){
  # for full model
  var_order <- c("hab", "exp", "cli", "inv", "dis", "pol", "sr_loss_av")
  var_labels <- c("Habitat", "Hunting", "Climate", "Invasive", "Disturbance", "Pollution", "sr_loss_av")
}

# convert codes to ordered factor
intercepts$code <- factor(
  intercepts$code,
  ordered = TRUE,
  levels = var_order,
  labels = var_labels
)

div$code <- factor(
  div$code,
  ordered = TRUE,
  levels = var_order,
  labels = var_labels
)

# define colour palette
col_pal <- c(
  brewer.pal(n = 7,"Purples")[3], 
  brewer.pal(n = 7,"Reds")[3], 
  brewer.pal(n = 7, "Blues")[3],
  brewer.pal(n = 7,"Oranges")[3],
  brewer.pal(n = 7,"Greens")[3], 
  "#f09fb1"
  )

## Plot
points_plot <- div %>%
  ggplot(aes(x = sr_loss_av, y = cd_loss_av, col = code)) +
  geom_hline(yintercept = 0, col = "grey50", alpha = 0.2)+
  geom_point(alpha = 0.2)+
  geom_abline(
    slope = intercepts[intercepts$code == "sr_loss_av", "value"],
    intercept = intercepts[intercepts$code == "Habitat", "value"],
    colour = col_pal[1], lwd = 1
    ) +
  geom_abline(
    slope = intercepts[intercepts$code == "sr_loss_av", "value"], 
    intercept = intercepts[intercepts$code == "Hunting", "value"], 
    colour = col_pal[2], lwd = 1
    ) +
  geom_abline(
    slope = intercepts[intercepts$code == "sr_loss_av", "value"],
    intercept = intercepts[intercepts$code == "Invasive", "value"],
    colour = col_pal[4], lwd = 1
  ) +
  geom_abline(
    slope = intercepts[intercepts$code == "sr_loss_av", "value"], 
    intercept = intercepts[intercepts$code == "Disturbance", "value"], 
    colour = col_pal[5], lwd = 1
    ) +
  theme_minimal() +
  labs(x = "Number of species extinctions avoided", y = "Colour-pattern diversity loss avoided", col = "Threat") +
  scale_colour_manual(values = col_pal, name = "") +
  guides(colour = guide_legend(override.aes  =  list(alpha  =  1, size = 2))) +
  theme(text = element_text(family = "Leelawadee UI Semilight", size = 15),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_line(colour = "grey90"),
        axis.line = element_line(colour = "grey20", linewidth = 0.5),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        legend.position = c(0.86,0.18),
        legend.background = element_rect(fill = "white", colour = "NA"),
        legend.title = element_blank(),
        plot.background  =  element_rect(fill = "white", colour = "NA"))

if(!remove_negligible){
  points_plot <- points_plot + 
    geom_abline(
      slope = intercepts[intercepts$code == "sr_loss_av", "value"],
      intercept = intercepts[intercepts$code == "Climate", "value"],
      colour = col_pal[3], lwd = 1
      ) +
    geom_abline(
      slope = intercepts[intercepts$code == "sr_loss_av", "value"],
      intercept = intercepts[intercepts$code == "Pollution", "value"],
      colour = col_pal[6], lwd = 1
      )
}
points_plot

# Plot values of intercepts for each extinction driver
intercept_plot <- intercepts %>% 
  subset(code != "sr_loss_av") %>% 
  ggplot(aes(x = code, y = value, fill = code)) +
  geom_col(linewidth = 0.25, colour = "black") +
  scale_fill_manual(values = col_pal) +
  theme_minimal() +
  theme(text = element_text(family = "Leelawadee UI Semilight", size = 15),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(colour = "grey20", linewidth = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(angle = 25, hjust = 0.75, vjust = 1),
        plot.background  =  element_rect(fill = "white", colour = "NA"),
        legend.position = "none",
        strip.text = element_text(size = 14),
        panel.spacing = grid::unit(0, "lines")) +
  labs(y = "Intercept",x = "Driver of extinction" ) +
  geom_hline(yintercept = 0)
intercept_plot

# Put plots together into figure
fig <- ggpubr::ggarrange(points_plot, intercept_plot, widths = c(4, 3.5), labels = c("a)", "b)"))

fig_filename <- if(remove_negligible){
  "ext_drivers_intercepts_remove_neglig_drivers.svg"
} else if(!remove_negligible){
  "ext_drivers_intercepts_all_drivers.svg"
}
ggsave(
  filename = here::here(
    "04_output_plots", "04_threat_abatements",
    fig_filename
    ),
  plot = fig,
  device = "svg",
  width = 9, height = 5.21
)
