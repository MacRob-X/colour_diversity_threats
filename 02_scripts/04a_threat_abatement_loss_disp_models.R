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
    X, code, iteration, cd_pc_loss, sr_loss_abs
  ) %>% 
  rename(
    scenario = X,
    cd_loss = cd_pc_loss,
    sr_loss = sr_loss_abs
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
    mean_cd_loss = mean(cd_loss),
    sd_cd_loss = sd(cd_loss),
    mean_sr_loss = mean(sr_loss),
    sd_sr_loss = sd(sr_loss)
  )

# Data analysis ----

# Run full model, with all individual drivers, species loss, iteration, and interactions between each 
# driver and species loss as predictors and colour pattern diversity loss as response
lme_full <- lmer(cd_loss ~ code * sr_loss + (1 | iteration), data = div, na.action = "na.fail")

# identify best model from all combinations of explanatory variables/interactions (using MuMIn::dredge)
dredge(lme_full)
# best model doesn't include interactions

# Run best model
lme <- lmer(cd_loss ~ code + sr_loss + (1 | iteration), data = div)
# inspect
summary(lme)

# Pairwise posthoc cmoparison between predictors
em <- emmeans::emmeans(lme, specs = ~ code, pbkrtest.limit = 10000)
# inspect
summary(em)

# Plotting ----

# extract intercepts for each extinction driver
intercepts <- summary(lme)[["coefficients"]][, "Estimate"]
intercepts[2:5]=intercepts[1]+intercepts[2:5] # to compare to overall intercept?
names(intercepts) <- gsub("code", "", names(intercepts))
intercepts <- intercepts %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("code")
colnames(intercepts) <- c("code", "value")

# define order of explanatory variables
var_order <- c("hab", "exp", "cli", "inv", "dis", "pol", "sr_loss")
intercepts$code <- factor(
  intercepts$code,
  ordered = TRUE,
  levels = var_order,
  labels = c("Habitat", "Hunting", "Climate", "Invasive", "Disturbance", "Pollution", "sr_loss")
)

div$code <- factor(
  div$code,
  ordered = TRUE,
  levels = var_order,
  labels = c("Habitat", "Hunting", "Climate", "Invasive", "Disturbance", "Pollution", "sr_loss")
)

# define colour palette
palette2 <- c(rev(c("#f09fb1", brewer.pal(n = 7, "Blues")[3], brewer.pal(n = 7,"Oranges")[3],
                 brewer.pal(n = 7,"Greens")[3], brewer.pal(n = 7,"Reds")[3], brewer.pal(n = 7,"Purples")[3])))


## I DON'T UNDERSTAND THIS PLOT - need to go over it tomorrow
## Some confusion between SR loss AVOIDED and SR loss
plot_1 <- div %>%
  ggplot(aes(x=sr_loss, y=cd_loss, col=code))+
  geom_hline(yintercept=0, col="grey50", alpha=0.2)+
  geom_point(alpha=0.2)+
  geom_abline(slope=intercepts[intercepts$code=="sr_loss","value"], intercept=intercepts[intercepts$code=="Habitat","value"], colour=palette2[1], lwd=1)+
  geom_abline(slope=intercepts[intercepts$code=="sr_loss","value"], intercept=intercepts[intercepts$code=="Hunting","value"], colour=palette2[2], lwd=1)+
  geom_abline(slope=intercepts[intercepts$code=="sr_loss","value"], intercept=intercepts[intercepts$code=="Climate","value"], colour=palette2[3], lwd=1)+
  geom_abline(slope=intercepts[intercepts$code=="sr_loss","value"], intercept=intercepts[intercepts$code=="Invasive","value"], colour=palette2[4], lwd=1)+
  geom_abline(slope=intercepts[intercepts$code=="sr_loss","value"], intercept=intercepts[intercepts$code=="Disturbance","value"], colour=palette2[5], lwd=1)+
  geom_abline(slope=intercepts[intercepts$code=="sr_loss","value"], intercept=intercepts[intercepts$code=="Pollution","value"], colour=palette2[6], lwd=1)+
  theme_minimal()+
  labs(x="Number of species extinctions avoided", y="Functional richness loss avoided (%)", col="Threat")+
  scale_colour_manual(values=palette2, name="")+
  guides(colour=guide_legend(override.aes = list(alpha = 1, size=2)))+
  theme(text=element_text(family="Leelawadee UI Semilight", size=15),
        panel.grid.major.y=element_line(colour="grey90"),
        panel.grid.minor.y=element_line(colour="grey90"),
        axis.line=element_line(colour="grey20", linewidth=0.5),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13),
        axis.title=element_text(size=16),
        legend.position=c(0.86,0.18),
        legend.background=element_rect(fill="white", colour="NA"),
        legend.title=element_blank(),
        plot.background = element_rect(fill="white", colour="NA"))
plot_1
