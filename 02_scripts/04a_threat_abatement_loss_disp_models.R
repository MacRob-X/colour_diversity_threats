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
library(MuMIn)
library(lme4)
library(RColorBrewer)

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
    X, code, iteration, cd_loss_avoided_pc, sr_loss_avoided_abs
  ) %>% 
  rename(
    scenario = X,
    cd_loss_av = cd_loss_avoided_pc,
    sr_loss_av = sr_loss_avoided_abs
  )

# Remove 'all' and 'none' scenario - we're not interested in abatement of all threats, 
# or of no threats, only the individual drivers
# (following Stewart et al 2025)
div <- div %>% 
  filter(
    !(code %in% c("all", "none"))
  )

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
dredge(lme_full)
# best model doesn't include interactions

# Run best model
lme <- lmer(cd_loss_av ~ code + sr_loss_av + (1 | iteration), data = div)
# note that the baseline value for code is taken as 'cli' automatically (alphabetically)
# inspect
summary(lme)

# Pairwise posthoc cmoparison between predictors
em <- emmeans::emmeans(lme, specs = ~ code, pbkrtest.limit = 10000)
# inspect
summary(em)

# Plotting ----

# extract intercepts for each extinction driver
intercepts <- summary(lme)[["coefficients"]][, "Estimate"]
intercepts[2:5] = intercepts[1] + intercepts[2:5] # to compare to baseline (cli) intercept
names(intercepts) <- gsub("code", "", names(intercepts))
names(intercepts)[1] <- "cli"
intercepts <- intercepts %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("code")
colnames(intercepts) <- c("code", "value")

# define order of explanatory variables
var_order <- c("hab", "exp", "cli", "inv", "dis", "pol", "sr_loss_av")
intercepts$code <- factor(
  intercepts$code,
  ordered = TRUE,
  levels = var_order,
  labels = c("Habitat", "Hunting", "Climate", "Invasive", "Disturbance", "Pollution", "sr_loss_av")
)

div$code <- factor(
  div$code,
  ordered = TRUE,
  levels = var_order,
  labels = c("Habitat", "Hunting", "Climate", "Invasive", "Disturbance", "Pollution", "sr_loss_av")
)

# define colour palette
palette2 <- c(rev(c("#f09fb1", brewer.pal(n = 7, "Blues")[3], brewer.pal(n = 7,"Oranges")[3],
                 brewer.pal(n = 7,"Greens")[3], brewer.pal(n = 7,"Reds")[3], brewer.pal(n = 7,"Purples")[3])))


## Pl
points_plot <- div %>%
  ggplot(aes(x = sr_loss_av, y = cd_loss_av, col = code)) +
  geom_hline(yintercept = 0, col = "grey50", alpha = 0.2)+
  geom_point(alpha = 0.2)+
  geom_abline(
    slope = intercepts[intercepts$code == "sr_loss_av", "value"], 
    intercept = intercepts[intercepts$code == "Habitat", "value"], 
    colour = palette2[1], lwd = 1
    ) +
  geom_abline(
    slope = intercepts[intercepts$code == "sr_loss_av", "value"], 
    intercept = intercepts[intercepts$code == "Hunting", "value"], 
    colour = palette2[2], lwd = 1
    ) +
  geom_abline(
    slope = intercepts[intercepts$code == "sr_loss_av", "value"], 
    intercept = intercepts[intercepts$code == "Climate", "value"], 
    colour = palette2[3], lwd = 1
    )+
  geom_abline(
    slope = intercepts[intercepts$code == "sr_loss_av", "value"], 
    intercept = intercepts[intercepts$code == "Invasive", "value"], 
    colour = palette2[4], lwd = 1
    )+
  geom_abline(
    slope = intercepts[intercepts$code == "sr_loss_av", "value"], 
    intercept = intercepts[intercepts$code == "Disturbance", "value"], 
    colour = palette2[5], lwd = 1
    )+
  geom_abline(
    slope = intercepts[intercepts$code == "sr_loss_av", "value"], 
    intercept = intercepts[intercepts$code == "Pollution", "value"], 
    colour = palette2[6], lwd = 1
    ) +
  theme_minimal() +
  labs(x = "Number of species extinctions avoided", y = "Colour-pattern diversity loss avoided (%)", col = "Threat") +
  scale_colour_manual(values = palette2, name = "") +
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
points_plot

# Plot values of intercepts for each extinction driver
intercept_plot <- intercepts %>% 
  subset(code != "sr_loss_av") %>% 
  ggplot(aes(x = code, y = value, fill = code)) +
  geom_col(linewidth = 0.25, colour = "black") +
  scale_fill_manual(values = palette2) +
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
fig <- ggpubr::ggarrange(points_plot, intercept_plot, widths = c(4, 3.5))

ggsave(
  filename = here::here(
    "04_output_plots", "04_threat_abatements",
    "ext_drivers_intercepts.svg"),
  plot = fig,
  device = "svg",
  width = 9, height = 5.21
)
