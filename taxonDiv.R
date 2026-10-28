# =============================================================================
# TAXONOMIC COMMUNITY ATTRIBUTES — DUNG BEETLES
# =============================================================================

# --- 0. Load packages --------------------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(glmmTMB)
library(emmeans)
library(multcomp)
library(multcompView)
library(ggplot2)
library(DHARMa)

# --- 1. Data import and preparation -----------------------------------------
dataDB <- read.csv("dataCASBufferDB.csv", header = TRUE, sep = ";")
dataDB$vegcover <- factor(dataDB$vegcover,
                          levels = c("SF", "CAS", "PA"),
                          labels = c("Secondary forest",
                                     "Cacao agroforestry system",
                                     "Pasture"))

# 1a. Aggregate at site level (sum across traps)
general_data <- dataDB %>%
  group_by(locality, vegcover, site) %>%
  summarise(
    richness        = n_distinct(species),
    abund_total     = sum(abund, na.rm = TRUE),
    biom_total      = sum(biom,  na.rm = TRUE),
    length_weighted = sum(length * biom, na.rm = TRUE) / sum(biom, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(site = as.integer(site))

# 1b. Full sampling design (fill zeros for combinations without captures)
design <- tibble::tribble(
  ~locality, ~vegcover,
  "EH", "Cacao agroforestry system",
  "ET", "Pasture",
  "LP", "Secondary forest",
  "LP", "Cacao agroforestry system",
  "LV", "Secondary forest",
  "LV", "Pasture",
  "LV", "Cacao agroforestry system",
  "SE", "Secondary forest",
  "SE", "Cacao agroforestry system",
  "SR", "Secondary forest",
  "SR", "Pasture",
  "SR", "Cacao agroforestry system",
  "VN", "Pasture"
)

general_data_complete <- design %>%
  crossing(site = 1:4) %>%
  left_join(general_data, by = c("locality", "vegcover", "site")) %>%
  mutate(
    across(c(richness, abund_total, biom_total), ~ replace_na(., 0)),
    length_weighted = replace_na(length_weighted, NA)
  )

# 1c. Subset for body length (only sites with captures)
data_length <- general_data_complete %>%
  filter(!is.na(length_weighted) & length_weighted > 0)

# 1d. Ensure factor order (least → most disturbed)
veg_order <- c("Secondary forest", "Cacao agroforestry system", "Pasture")
general_data_complete <- general_data_complete %>%
  mutate(vegcover = factor(vegcover, levels = veg_order))
data_length <- data_length %>%
  mutate(vegcover = factor(vegcover, levels = veg_order))

# --- 2. Modelling ------------------------------------------------------------
# 2a. Species richness — Poisson
model_rich <- glmmTMB(richness ~ vegcover + (1 | locality),
                      data = general_data_complete, family = poisson)
sim_rich <- simulateResiduals(model_rich, n = 250)
plot(sim_rich)
testUniformity(sim_rich)
testOutliers(sim_rich)
summary(model_rich)
drop1(model_rich, test = "Chisq")

# 2b. Total abundance — Negative binomial type 1
model_abund <- glmmTMB(abund_total ~ vegcover + (1 | locality),
                       data = general_data_complete, family = nbinom1)
sim_abund <- simulateResiduals(model_abund, n = 250)
plot(sim_abund)
testUniformity(sim_abund)
testOutliers(sim_abund)
summary(model_abund)
drop1(model_abund, test = "Chisq")

# 2c. Total biomass — Tweedie
model_biom <- glmmTMB(biom_total ~ vegcover + (1 | locality),
                      data = general_data_complete,
                      family = tweedie(link = "log"))
sim_biom <- simulateResiduals(model_biom, n = 250)
plot(sim_biom)
testUniformity(sim_biom)
testOutliers(sim_biom)
summary(model_biom)
drop1(model_biom, test = "Chisq")

# 2d. Body length — t-Student on log scale (best fit among 6 distributions)
model_length <- glmmTMB(log(length_weighted) ~ vegcover + (1 | locality),
                        data = data_length, family = t_family(link = "identity"))
sim_length <- simulateResiduals(model_length, n = 250)
plot(sim_length)
testUniformity(sim_length)
testOutliers(sim_length)
summary(model_length)
drop1(model_length, test = "Chisq")

library(emmeans)

emm_length <- emmeans(model_length, ~ vegcover)
contrastes_length <- pairs(emm_length, adjust = "tukey", reverse = TRUE)
contrastes_length

# --- 3. Post-hoc and compact letter display ----------------------------------
cld_rich  <- cld(emmeans(model_rich,  ~ vegcover, type = "response"),
                 Letters = letters, adjust = "bonferroni")
cld_abund <- cld(emmeans(model_abund, ~ vegcover, type = "response"),
                 Letters = letters, adjust = "bonferroni")
cld_biom  <- cld(emmeans(model_biom,  ~ vegcover, type = "response"),
                 Letters = letters, adjust = "bonferroni")
cld_length <- cld(emmeans(model_length, ~ vegcover, type = "response"),
                  Letters = letters, adjust = "bonferroni")

# --- 4. Figures --------------------------------------------------------------
# Colour palette (original style)
cols <- c("Secondary forest"          = "forestgreen",
          "Cacao agroforestry system" = "darkorange4",
          "Pasture"                   = "springgreen")

upper_whisker <- function(x) quantile(x, 0.75, na.rm = TRUE) + 0.5 * IQR(x, na.rm = TRUE)

theme_tax <- theme_classic() +
  theme(
    axis.title.x     = element_text(margin = margin(t = 10),
                                    family = "serif", size = 20, face = "bold"),
    axis.title.y     = element_text(margin = margin(r = 10),
                                    family = "serif", size = 20, face = "bold"),
    axis.text        = element_text(family = "serif", size = 18),
    axis.line        = element_line(colour = "black"),
    panel.grid.major = element_line(color = "gray95", linetype = "dashed"),
    legend.position  = "none"
  )

# Helper to create a panel with letters
make_panel <- function(data, var, cld_obj, ylab) {
  pos <- data %>%
    group_by(vegcover) %>%
    summarise(ymax = upper_whisker(.data[[var]]), .groups = "drop")
  cld_df <- as.data.frame(cld_obj) %>%
    left_join(pos, by = "vegcover") %>%
    mutate(.group = trimws(.group))
  ggplot(data, aes(x = vegcover, y = .data[[var]], fill = vegcover)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21) +
    geom_text(data = cld_df,
              aes(x = vegcover, y = ymax * 1.08, label = .group),
              inherit.aes = FALSE, family = "serif", size = 8) +
    scale_fill_manual(values = cols) +
    labs(x = NULL, y = ylab) +
    theme_tax
}

# Individual plots
p_rich  <- make_panel(general_data_complete, "richness",    cld_rich,  "Species richness") +
  labs(x = NULL)+
  theme(axis.text.x = element_blank())
p_abund <- make_panel(general_data_complete, "abund_total", cld_abund, "Total abundance") +
  labs(x = NULL)+
  theme(axis.text.x = element_blank())
p_biom  <- make_panel(general_data_complete, "biom_total",  cld_biom,  "Total biomass (mg)") +
  labs(x = NULL)
p_len <- make_panel(data_length, "length_weighted", cld_length, "Body length (mm)") +
  labs(x = NULL) +
  scale_y_continuous(limits = c(14, NA))


library(patchwork)
fig_tax <- (p_rich + p_abund + p_biom + p_len) +
  plot_layout(ncol = 2) + # Organiza los 4 gráficos en 2 columnas (2x2)
  plot_annotation(
    caption = "Vegetation cover",
    theme   = theme(plot.caption = element_text(hjust  = 0.5,
                                                family = "serif",
                                                size   = 20, face = "bold"))
  )

ggsave("figures/Fig_tax.svg", fig_tax, width = 17, height = 12, dpi = 600)
