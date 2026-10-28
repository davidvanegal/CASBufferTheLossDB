# =============================================================================
# Rank-abundance distribution — Dung beetles by vegetation cover
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(BiodiversityR)   
library(indicspecies)    
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(ggpubr)
library(viridis)
library(openxlsx)

# --- Data --------------------------------------------------------------------
dataDB <- read.csv("dataCASBufferDB.csv", header = TRUE, sep = ";")
dataDB$vegcover <- factor(dataDB$vegcover,
                          levels = c("SF", "CAS", "PA"),
                          labels = c("Secondary forest",
                                     "Cacao agroforestry system",
                                     "Pasture"))

# =============================================================================
# 1. BUILD MATRICES FROM dataDB
# =============================================================================

# --- 1a. Site-level matrix (locality × vegcover × site) ----------------------
# Aggregate trap-level records to site level
comm_site <- dataDB %>%
  group_by(locality, vegcover, site, species) %>%
  summarise(abund = sum(abund), .groups = "drop") %>%
  mutate(siteID = paste(locality, vegcover, site, sep = "_")) %>%
  pivot_wider(names_from  = species,
              values_from = abund,
              values_fill = 0)

# Site metadata
metadata <- comm_site %>%
  distinct(siteID, locality, vegcover)

# --- 1b. Cover-level matrix (vegcover × species) -----------------------------
# Aggregate all sites within each vegetation cover type
sp_mat_long <- dataDB %>%
  group_by(vegcover, species) %>%
  summarise(abund = sum(abund), .groups = "drop")

sp_mat <- sp_mat_long %>%
  pivot_wider(names_from  = species,
              values_from = abund,
              values_fill = 0) %>%
  tibble::column_to_rownames("vegcover")

# Ensure cover order is preserved
sp_mat <- sp_mat[c("Secondary forest", "Cacao agroforestry system", "Pasture"), ]

# Environment table with grouping factor (required by BiodiversityR)
env <- data.frame(
  VegetationCover = factor(
    rownames(sp_mat),
    levels = c("Secondary forest", "Cacao agroforestry system", "Pasture")
  ),
  row.names = rownames(sp_mat)
)

# IndVal group codes (1 = SF, 2 = Cacao, 3 = Pasture)
indval_mat  <- sp_mat
cover_group <- as.integer(env$VegetationCover)

# =============================================================================
# 2. RANK-ABUNDANCE CURVES — BiodiversityR::rankabundance()
# =============================================================================

ra_sf <- rankabundance(sp_mat, y = env,
                       factor = "VegetationCover",
                       level  = "Secondary forest") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Species") %>%
  mutate(VegetationCover = "Secondary forest",
         number          = as.integer(rank))

ra_ca <- rankabundance(sp_mat, y = env,
                       factor = "VegetationCover",
                       level  = "Cacao agroforestry system") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Species") %>%
  mutate(VegetationCover = "Cacao agroforestry system",
         number          = as.integer(rank))

ra_pa <- rankabundance(sp_mat, y = env,
                       factor = "VegetationCover",
                       level  = "Pasture") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Species") %>%
  mutate(VegetationCover = "Pasture",
         number          = as.integer(rank))

# Combine and keep species with proportion > 0
outtotal <- bind_rows(ra_sf, ra_ca, ra_pa) %>%
  filter(proportion > 0) %>%
  mutate(
    VegetationCover = factor(
      VegetationCover,
      levels = c("Secondary forest", "Cacao agroforestry system", "Pasture")
    )
  )

# =============================================================================
# 3. RANK-ABUNDANCE PLOT
# =============================================================================

cols <- c("Secondary forest"          = "forestgreen",
          "Cacao agroforestry system" = "darkorange4",
          "Pasture"                   = "springgreen")

plot_ra <- ggplot(outtotal,
                  aes(x = number, y = proportion,
                      color = VegetationCover, shape = VegetationCover)) +
  geom_point(size = 4) +
  geom_line(linewidth = 1, linetype = 3) +
  scale_color_manual(values = cols) +
  scale_x_continuous("Species rank") +
  scale_y_continuous("Proportional abundance") +
  facet_grid(cols = vars(VegetationCover)) +
  labs(color = "Vegetation cover", shape = "Vegetation cover") +
  theme(
    axis.title.x      = element_text(margin = margin(t = 10),
                                     family = "serif", size = 18, face = "bold"),
    axis.title.y      = element_text(margin = margin(r = 10),
                                     family = "serif", size = 18, face = "bold"),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.border      = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "gray95", linetype = "dashed"),
    axis.text         = element_text(family = "serif", size = 14, face = "bold"),
    axis.line         = element_line(colour = "black"),
    legend.position   = "top",
    legend.text       = element_text(family = "serif", size = 14),
    legend.title      = element_text(family = "serif", size = 18, face = "bold"),
    strip.background  = element_blank(),
    strip.text.x      = element_blank()
  )

# =============================================================================
# 4. THEORETICAL MODEL FITTING — radfit(), AIC-based selection
# Models compared: Null, Preemption, Log-Normal, Zipf, Mandelbrot
# =============================================================================

rad_sf <- radfit(sp_mat["Secondary forest", ])
rad_ca <- radfit(sp_mat["Cacao agroforestry system", ])
rad_pa <- radfit(sp_mat["Pasture", ])

# Best model per cover (lowest AIC)
best_sf <- names(which.min(AIC(rad_sf)))
best_ca <- names(which.min(AIC(rad_ca)))
best_pa <- names(which.min(AIC(rad_pa)))

# AIC comparison table
make_aic_tbl <- function(rad_obj, cover_name) {
  as.data.frame(AIC(rad_obj)) %>%
    tibble::rownames_to_column("model") %>%
    mutate(VegetationCover = cover_name, .before = 1) %>%
    rename(AIC = 3) %>%
    mutate(best = AIC == min(AIC))
}

aic_tbl <- bind_rows(
  make_aic_tbl(rad_sf, "Secondary forest"),
  make_aic_tbl(rad_ca, "Cacao agroforestry system"),
  make_aic_tbl(rad_pa, "Pasture")
)

# --- 4b. Extract predictions from the best model per cover -------------------

get_model_preds <- function(rad_obj, cover_name) {
  best          <- names(which.min(AIC(rad_obj)))
  best_model    <- rad_obj$models[[best]]
  fitted_vals   <- fitted(best_model)
  n             <- sum(rad_obj$y)
  fitted_sorted <- sort(fitted_vals, decreasing = TRUE)

  data.frame(
    number          = seq_along(fitted_sorted),
    prop_predicted  = as.numeric(fitted_sorted) / n * 100,
    VegetationCover = cover_name,
    model           = best
  )
}

preds_df <- bind_rows(
  get_model_preds(rad_sf, "Secondary forest"),
  get_model_preds(rad_ca, "Cacao agroforestry system"),
  get_model_preds(rad_pa, "Pasture")
) %>%
  mutate(
    VegetationCover = factor(
      VegetationCover,
      levels = c("Secondary forest", "Cacao agroforestry system", "Pasture")
    )
  )

# --- 4c. Combined plot: observed points + fitted curve -----------------------

# Label top 5 species by proportional abundance per cover
labels_df <- outtotal %>%
  group_by(VegetationCover) %>%
  slice_max(proportion, n = 5) %>%
  ungroup()

plot_ra_models <- ggplot() +
  geom_point(data = outtotal,
             aes(x = number, y = proportion,
                 color = VegetationCover, shape = VegetationCover),
             size = 4) +
  geom_line(data = outtotal,
            aes(x = number, y = proportion, color = VegetationCover),
            linewidth = 0.8, linetype = "dotted") +
  geom_line(data = preds_df,
            aes(x = number, y = prop_predicted),
            color = "black", linewidth = 1.2, linetype = "solid") +
  geom_text(data = labels_df,
            aes(x = number, y = proportion,
                label   = Species,
                color   = "black"),
            size        = 4,
            family      = "serif",
            fontface    = "bold.italic",
            show.legend = FALSE,
            hjust       = -0.15,
            vjust       = 0.5) +
  scale_color_manual(values = cols) +
  scale_x_continuous("Species rank") +
  scale_y_log10("Proportional abundance (log scale)",
                labels = scales::label_number(accuracy = 0.01)) +
  coord_cartesian(clip = "off") +
  facet_grid(cols = vars(VegetationCover)) +
  labs(color = "Vegetation cover", shape = "Vegetation cover") +
  theme(
    axis.title.x      = element_text(margin = margin(t = 10),
                                     family = "serif", size = 18, face = "bold"),
    axis.title.y      = element_text(margin = margin(r = 10),
                                     family = "serif", size = 18, face = "bold"),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.border      = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "gray95", linetype = "dashed"),
    axis.text         = element_text(family = "serif", size = 14, face = "bold"),
    axis.line         = element_line(colour = "black"),
    legend.position   = "top",
    legend.text       = element_text(family = "serif", size = 14),
    legend.title      = element_text(family = "serif", size = 18, face = "bold"),
    strip.background  = element_blank(),
    strip.text.x      = element_blank()
  )

ggsave("RA/RangeAbundance_withModels.svg",
       plot_ra_models, width = 15, height = 8, dpi = 600)

# =============================================================================
# 5. PIELOU'S EVENNESS (J') PER VEGETATION COVER
# =============================================================================

jevenness_vals <- diversityresult(sp_mat, index = "Jevenness", method = "each site")
shannon_vals   <- diversityresult(sp_mat, index = "Shannon",   method = "each site")
richness_vals  <- diversityresult(sp_mat, index = "richness",  method = "each site")

pielou_df <- data.frame(
  VegetationCover = factor(
    rownames(sp_mat),
    levels = c("Secondary forest", "Cacao agroforestry system", "Pasture")
  ),
  richness = richness_vals[, 1],
  shannon  = round(shannon_vals[, 1],   3),
  pielou_J = round(jevenness_vals[, 1], 3)
)

cat("\n--- Pielou's evenness (J') per vegetation cover ---\n")
print(pielou_df)

# Pielou bar plot
plot_pielou <- ggplot(pielou_df,
                      aes(x = VegetationCover, y = pielou_J,
                          fill = VegetationCover)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_text(aes(label = pielou_J),
            vjust = -0.5, family = "serif", size = 5) +
  scale_fill_viridis_d(option = "D") +
  scale_y_continuous("Pielou's evenness (J')",
                     limits = c(0, 1.1),
                     breaks = seq(0, 1, 0.25)) +
  labs(x = NULL) +
  theme(
    axis.title.y     = element_text(margin = margin(r = 10),
                                    family = "serif", size = 16, face = "bold"),
    axis.text        = element_text(family = "serif", size = 13, face = "bold"),
    axis.text.x      = element_text(angle = 15, hjust = 1),
    axis.line        = element_line(colour = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.border     = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95", linetype = "dashed"),
    legend.position  = "none"
  )

# =============================================================================
# 6. INDICATOR SPECIES ANALYSIS — indicspecies::multipatt()
# =============================================================================

cat("\n--- Indicator Species Analysis (IndVal, 999 permutations) ---\n")

set.seed(123)
indval_result <- multipatt(
  indval_mat,
  cluster = cover_group,
  func    = "IndVal.g",         # corrected for unequal group sizes
  control = how(nperm = 999)
)

cat("\nIndicator species summary (all results):\n")
summary(indval_result, alpha = 1)   # alpha = 1 shows all species

# Extract results for export and plotting
# NOTE: with n = 3 sampling units, permutation p-values have limited power
# (minimum resolution = 1/3). Exclusivity (s columns) is used as criterion.
indval_tbl <- indval_result$sign %>%
  tibble::rownames_to_column("Species") %>%
  rename(IndVal = stat) %>%
  mutate(
    IndVal  = round(IndVal, 3),
    n_groups = s.1 + s.2 + s.3,
    VegetationCover = case_when(
      s.1 == 1 & s.2 == 0 & s.3 == 0 ~ "Secondary forest",
      s.1 == 0 & s.2 == 1 & s.3 == 0 ~ "Cacao agroforestry system",
      s.1 == 0 & s.2 == 0 & s.3 == 1 ~ "Pasture",
      TRUE                             ~ "Generalist"
    ),
    exclusive = n_groups == 1
  ) %>%
  arrange(VegetationCover, desc(IndVal))

cat("\n--- Exclusive indicator species (present in only one cover type) ---\n")
print(indval_tbl %>% filter(exclusive))

# =============================================================================
# 7. EXPORT
# =============================================================================

ggsave("RA/RangeAbundance.svg", plot_ra, width = 13, height = 7, dpi = 600)

# Excel workbook — one sheet per analysis
wb <- createWorkbook()

addWorksheet(wb, "RadfitAIC")
writeDataTable(wb, "RadfitAIC", aic_tbl,    tableStyle = "TableStyleMedium9")
setColWidths(wb,  "RadfitAIC", cols = 1:ncol(aic_tbl),    widths = "auto")

addWorksheet(wb, "Pielou")
writeDataTable(wb, "Pielou",    pielou_df,  tableStyle = "TableStyleMedium9")
setColWidths(wb,  "Pielou",    cols = 1:ncol(pielou_df),  widths = "auto")

addWorksheet(wb, "IndVal")
writeDataTable(wb, "IndVal",    indval_tbl, tableStyle = "TableStyleMedium9")
setColWidths(wb,  "IndVal",    cols = 1:ncol(indval_tbl), widths = "auto")

saveWorkbook(wb, "RA/RankAbundance_Results.xlsx", overwrite = TRUE)
