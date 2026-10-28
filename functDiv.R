# =============================================================================
# FUNCTIONAL DIVERSITY — DUNG BEETLES
# =============================================================================

# --- 0. Load packages --------------------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)  
library(FD)     
library(glmmTMB)
library(emmeans)
library(multcomp)
library(multcomp)
library(ggplot2)
library(DHARMa)

# --- 1. Data import ---------------------------------------------------------
dataDB <- read.csv("dataCASBufferDB.csv", header = TRUE, sep = ";")
dataDB$vegcover <- factor(dataDB$vegcover,
                          levels = c("SF", "CAS", "PA"),
                          labels = c("Secondary forest",
                                     "Cacao agroforestry system",
                                     "Pasture"))

# --- 2. Community matrix ----------------------------------------------------
comm_site <- dataDB %>%
  group_by(locality, vegcover, site, species) %>%
  summarise(abund = sum(abund), .groups = "drop") %>%
  mutate(siteID = paste(locality, vegcover, site, sep = "_")) %>%
  pivot_wider(names_from = species, values_from = abund, values_fill = 0)

metadata <- comm_site %>% distinct(siteID, locality, vegcover)
comm_mat <- as.data.frame(comm_site[, setdiff(colnames(comm_site),
                                              c("locality", "vegcover", "site", "siteID"))])
rownames(comm_mat) <- comm_site$siteID

# --- 3. Trait matrix with allometric shape -----------------------------------
# 3a. Species‑level means
sp_traits <- dataDB %>%
  group_by(species) %>%
  summarise(length_m = mean(length, na.rm = TRUE),
            mass_m   = mean(biom,   na.rm = TRUE),
            .groups  = "drop")

# 3b. Allometric model and shape residual
allo_fit <- lm(log(mass_m) ~ log(length_m), data = sp_traits)
sp_traits$shape_res <- residuals(allo_fit)   # body shape

# 3c. Build trait matrix: length, shape, guild
traits <- sp_traits %>%
  dplyr::select(species, length = length_m, shape = shape_res) %>%
  left_join(dataDB %>% group_by(species) %>% 
              summarise(guild = first(strategy), .groups = "drop"),
            by = "species") %>%
  column_to_rownames("species")

traits$length <- c(scale(traits$length))   # standardize
traits$shape  <- c(scale(traits$shape))
traits$guild  <- as.factor(traits$guild)

# 3d. Keep only shared species
common_sp <- intersect(colnames(comm_mat), rownames(traits))
comm_mat2 <- as.matrix(comm_mat[, common_sp])
traits2   <- traits[common_sp, ]

# --- 4. Functional diversity indices -----------------------------------------
fd <- dbFD(traits2, comm_mat2, calc.FRic = FALSE, calc.FDiv = FALSE,
           calc.CWM = TRUE, calc.FGR = FALSE)   # FDis and RaoQ are computed automatically

fd_results <- data.frame(
  siteID = rownames(comm_mat2),
  FDis  = fd$FDis,
  RaoQ  = fd$RaoQ
) %>% left_join(metadata, by = "siteID")

# Correlation between FDis and RaoQ (for supplementary justification)
cor_fdis_raoq <- cor(fd_results$FDis, fd_results$RaoQ, method = "pearson")
print(paste("Pearson r between FDis and RaoQ:", round(cor_fdis_raoq, 4)))

# --- 5. Community‑weighted mean (CWM) of raw traits --------------------------
traits_raw <- sp_traits %>%
  dplyr::select(species, length_raw = length_m, shape_raw = shape_res) %>%
  left_join(dataDB %>% group_by(species) %>% 
              summarise(guild = first(strategy), .groups = "drop"),
            by = "species") %>%
  column_to_rownames("species")
traits_raw2 <- traits_raw[common_sp, ]

comm_rel <- comm_mat2 / rowSums(comm_mat2)
cwm_raw <- as.data.frame(as.matrix(comm_rel) %*% 
                           as.matrix(traits_raw2[, c("length_raw", "shape_raw")])) %>%
  rownames_to_column("siteID") %>%
  left_join(metadata, by = "siteID")

# --- 6. GLMMs for RaoQ and CWMs ----------------------------------------------
# 6a. RaoQ — Gaussian
model_RaoQ <- glmmTMB(RaoQ ~ vegcover + (1 | locality),
                      data = fd_results, family = gaussian)
sim_RaoQ <- simulateResiduals(model_RaoQ, n = 250)
plot(sim_RaoQ)
testUniformity(sim_RaoQ)
testOutliers(sim_RaoQ)
summary(model_RaoQ)
drop1(model_RaoQ, test = "Chisq")

# 6b. CWM length — log‑normal
model_cwm_len <- glmmTMB(log(length_raw) ~ vegcover + (1 | locality),
                         data = cwm_raw, family = gaussian)
sim_cwm_len <- simulateResiduals(model_cwm_len, n = 250)
plot(sim_cwm_len)
testUniformity(sim_cwm_len)
testOutliers(sim_cwm_len)
summary(model_cwm_len)
drop1(model_cwm_len, test = "Chisq")

# 6c. CWM shape — Gaussian identity
model_cwm_shape <- glmmTMB(shape_raw ~ vegcover + (1 | locality),
                           data = cwm_raw, family = gaussian)
sim_cwm_shape <- simulateResiduals(model_cwm_shape, n = 250)
plot(sim_cwm_shape)
testUniformity(sim_cwm_shape)
testOutliers(sim_cwm_shape)
summary(model_cwm_shape)
drop1(model_cwm_shape, test = "Chisq")

# --- 7. Post‑hoc letters -----------------------------------------------------
cld_RaoQ  <- cld(emmeans(model_RaoQ,  ~ vegcover, type = "response"),
                 Letters = letters, adjust = "bonferroni")
cld_len   <- cld(emmeans(model_cwm_len, ~ vegcover, type = "response"),
                 Letters = letters, adjust = "bonferroni")
cld_shape <- cld(emmeans(model_cwm_shape, ~ vegcover, type = "response"),
                 Letters = letters, adjust = "bonferroni")

# --- 8. Guild proportions analysis ------------------------------------------
# 8a. Compute relative abundance of each guild per site
# (rollers excluded because virtually absent)
guild_abund <- dataDB %>%
  group_by(locality, vegcover, site, strategy) %>%
  summarise(total = sum(abund), .groups = "drop") %>%
  filter(strategy != "Roller") %>%
  group_by(locality, vegcover, site) %>%
  mutate(prop = total / sum(total)) %>%
  ungroup()

# Focus on dwellers (the most contrasting guild)
tunneler_data <- guild_abund %>%
  filter(strategy == "Tunneler") %>%
  mutate(siteID = paste(locality, vegcover, site, sep = "_")) %>%
  # keep only the columns we need – avoid duplicate vegcover from metadata
  dplyr::select(siteID, vegcover, locality, prop)

# 8b. Model proportion of tunnelers (ordbeta handles 0s and 1s)
model_tunneler <- glmmTMB(prop ~ vegcover + (1 | locality),
                         data = tunneler_data, family = ordbeta(link = "logit"))
sim_tunneler <- simulateResiduals(model_tunneler, n = 250)
plot(sim_tunneler)
testUniformity(sim_tunneler)
testOutliers(sim_tunneler)
summary(model_tunneler)
drop1(model_tunneler, test = "Chisq")

medias_shape <- emmeans(model_tunneler, ~ vegcover)
contrastes_shape <- pairs(medias_shape, adjust = "tukey")
print(contrastes_shape)

cld_tunneler <- cld(emmeans(model_tunneler, ~ vegcover, type = "response"),
                   Letters = letters, adjust = "bonferroni")

# --- 9. Figures --------------------------------------------------------------
cols <- c("Secondary forest"          = "forestgreen",
          "Cacao agroforestry system" = "darkorange4",
          "Pasture"                   = "springgreen")

upper_whisker <- function(x) quantile(x, 0.75, na.rm = TRUE) + 0.5 * IQR(x, na.rm = TRUE)

theme_fun <- theme_classic() +
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

make_fun_panel <- function(data, var, cld_obj, ylab) {
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
    theme_fun
}

p_tunneler <- make_fun_panel(tunneler_data, "prop", cld_tunneler,
                            "Proportion of tunnelers")+
  labs(x = "Vegetation cover")+
  theme(axis.title.x = element_blank())

p_RaoQ  <- make_fun_panel(fd_results, "RaoQ", cld_RaoQ,
                          "Rao's quadratic entropy (RaoQ)") +
  labs(x = "Vegetation cover")+
  theme(axis.title.x = element_blank())

p_len   <- make_fun_panel(cwm_raw, "length_raw", cld_len,
                          "CWM body length (mm)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

p_shape <- make_fun_panel(cwm_raw, "shape_raw", cld_shape,
                          "CWM body shape (residual)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

# Display all figures

library(patchwork)
fig_fun <- (p_len + p_shape + p_RaoQ + p_tunneler) +
  plot_layout(ncol = 2) + # Organiza los 4 gráficos en 2 columnas (2x2)
  plot_annotation(
    caption = "Vegetation cover",
    theme   = theme(plot.caption = element_text(hjust  = 0.5,
                                                family = "serif",
                                                size   = 20, face = "bold"))
  )

ggsave("figures/Fig_fun.svg", fig_fun, width = 17, height = 12, dpi = 600)

# =============================================================================
# END OF SCRIPT
# =============================================================================