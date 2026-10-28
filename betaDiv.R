# =============================================================================
# BETA DIVERSITY & TAXONOMIC–FUNCTIONAL DECOUPLING — DUNG BEETLES
# =============================================================================

# --- 0. Load packages --------------------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library(BAT)
library(glmmTMB)
library(emmeans)
library(multcomp)     
library(DHARMa)
library(ggplot2)
library(ggforce)
library(pairwiseAdonis)

# --- 1. Data import ---------------------------------------------------------
dataDB <- read.csv("dataCASBufferDB.csv", header = TRUE, sep = ";")
dataDB$vegcover <- factor(dataDB$vegcover,
                          levels = c("SF", "CAS", "PA"),
                          labels = c("Secondary forest",
                                     "Cacao agroforestry system",
                                     "Pasture"))
dataDB <- rename(dataDB, guild = strategy)

# Base theme
theme_base <- theme_classic() +
  theme(
    axis.title.x      = element_text(margin = margin(t = 10),
                                     family = "serif", size = 20, face = "bold"),
    axis.title.y      = element_text(margin = margin(r = 10),
                                     family = "serif", size = 20, face = "bold"),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.border      = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "gray95", linetype = "dashed"),
    axis.text         = element_text(family = "serif", size = 16),
    axis.line         = element_line(colour = "black"),
    legend.position   = "none",
    legend.text       = element_text(family = "serif", size = 16),
    legend.title      = element_text(family = "serif", size = 20, face = "bold"),
    strip.background  = element_blank()
  )

# Strategy theme: adds facet strip styling
theme_strategy <- theme_base +
  theme(strip.background = element_rect(fill = "grey90", color = NA),
        strip.text       = element_text(size = 10))

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

# --- 3. Trait matrix (length + shape + guild) --------------------------------
sp_traits <- dataDB %>%
  group_by(species) %>%
  summarise(length_m = mean(length, na.rm = TRUE),
            mass_m   = mean(biom,   na.rm = TRUE), .groups = "drop")
allo_fit <- lm(log(mass_m) ~ log(length_m), data = sp_traits)
sp_traits$shape_res <- residuals(allo_fit)

traits <- sp_traits %>%
  dplyr::select(species, length = length_m, shape = shape_res) %>%
  left_join(dataDB %>% group_by(species) %>%
              summarise(guild = first(guild), .groups = "drop"), by = "species") %>%
  column_to_rownames("species")
traits$length <- c(scale(traits$length))
traits$shape  <- c(scale(traits$shape))
traits$guild  <- as.factor(traits$guild)

common_sp <- intersect(colnames(comm_mat), rownames(traits))
comm_mat2 <- as.matrix(comm_mat[, common_sp])
traits2   <- traits[common_sp, ]

# --- 4. BAT: dissimilarity matrices + components (substitution & nestedness) --
tax_beta <- beta(comm_mat2, abund = TRUE, func = "sorensen", comp = TRUE)
tax_dist_total <- as.dist(tax_beta$Btotal)
tax_sub <- as.dist(tax_beta$Brepl)
tax_nest <- as.dist(tax_beta$Brich)

fun_beta <- beta(comm_mat2, tree = traits2, abund = TRUE, func = "sorensen", comp = TRUE)
fun_dist_total <- as.dist(fun_beta$Btotal)
fun_sub <- as.dist(fun_beta$Brepl)
fun_nest <- as.dist(fun_beta$Brich)

# --- 5. NMDS, PERMANOVA, PERMDISP (taxonomic) --------------------------------
nmds_tax <- metaMDS(tax_dist_total, k = 2, trymax = 100)
permanova_tax <- adonis2(tax_dist_total ~ vegcover, data = metadata, permutations = 9999)
disp_tax <- betadisper(tax_dist_total, metadata$vegcover)
print(permanova_tax)
print(anova(disp_tax))
print(permutest(disp_tax))

nmds_scores_tax <- as.data.frame(scores(nmds_tax))
nmds_scores_tax$siteID <- rownames(nmds_scores_tax)
nmds_scores_tax <- left_join(nmds_scores_tax, metadata, by = "siteID")

cols <- c("Secondary forest" = "forestgreen",
          "Cacao agroforestry system" = "darkorange4",
          "Pasture" = "springgreen")

p_nmds_tax <- ggplot(nmds_scores_tax,
                     aes(x = NMDS1, y = NMDS2, color = vegcover, fill = vegcover)) +
  geom_point(aes(shape = vegcover), size = 3, stroke = 1) +
  geom_mark_ellipse(alpha = 0.15, expand = 0) +
  labs(fill = "Vegetation cover", color = "Vegetation cover",
       shape = "Vegetation cover") + 
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_shape_manual(values = c(21, 22, 24)) +
  theme_strategy +
  theme(legend.position = "top")

# --- 6. NMDS, PERMANOVA, PERMDISP (functional) -------------------------------
nmds_fun <- metaMDS(fun_dist_total, k = 2, trymax = 100)
permanova_fun <- adonis2(fun_dist_total ~ vegcover, data = metadata, permutations = 9999)
disp_fun <- betadisper(fun_dist_total, metadata$vegcover)
print(permanova_fun)
print(anova(disp_fun))
print(permutest(disp_fun))

nmds_scores_fun <- as.data.frame(scores(nmds_fun))
nmds_scores_fun$siteID <- rownames(nmds_scores_fun)
nmds_scores_fun <- left_join(nmds_scores_fun, metadata, by = "siteID")

p_nmds_fun <- ggplot(nmds_scores_fun,
                     aes(x = NMDS1, y = NMDS2, color = vegcover, fill = vegcover)) +
  geom_point(aes(shape = vegcover), size = 3, stroke = 1) +
  geom_mark_ellipse(alpha = 0.15, expand = 0) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_shape_manual(values = c(21, 22, 24)) +
  theme_strategy +
  theme(legend.position = "top")

# --- 7. Pairwise PERMANOVA (taxonomic and functional) ------------------------
cat("\n=== Pairwise PERMANOVA: Taxonomic ===\n")
pair_tax <- pairwise.adonis2(tax_dist_total ~ vegcover, data = metadata,
                             permutations = 9999, p.adjust = "bonferroni")
print(pair_tax)

cat("\n=== Pairwise PERMANOVA: Functional ===\n")
pair_fun <- pairwise.adonis2(fun_dist_total ~ vegcover, data = metadata,
                             permutations = 9999, p.adjust = "bonferroni")
print(pair_fun)

# --- 8. Build dyadic dataset -------------------------------------------------
tax_mat_total <- as.matrix(tax_dist_total)
tax_mat_sub   <- as.matrix(tax_sub)
tax_mat_nest  <- as.matrix(tax_nest)
fun_mat_total <- as.matrix(fun_dist_total)
fun_mat_sub   <- as.matrix(fun_sub)
fun_mat_nest  <- as.matrix(fun_nest)
sites <- rownames(comm_mat2)

dyad <- expand.grid(Site1 = sites, Site2 = sites, stringsAsFactors = FALSE) %>%
  filter(Site1 < Site2) %>%
  rowwise() %>%
  mutate(
    Cover1   = metadata$vegcover[match(Site1, metadata$siteID)],
    Cover2   = metadata$vegcover[match(Site2, metadata$siteID)],
    TaxDist  = tax_mat_total[Site1, Site2],
    TaxSub   = tax_mat_sub[Site1, Site2],
    TaxNest  = tax_mat_nest[Site1, Site2],
    FunDist  = fun_mat_total[Site1, Site2],
    FunSub   = fun_mat_sub[Site1, Site2],
    FunNest  = fun_mat_nest[Site1, Site2]
  ) %>%
  ungroup()

# --- 9. Categorise pair types ------------------------------------------------
dyad <- dyad %>%
  mutate(
    PairType = case_when(
      Cover1 == "Secondary forest" & Cover2 == "Secondary forest" ~ "SF-SF",
      Cover1 == "Cacao agroforestry system" & Cover2 == "Cacao agroforestry system" ~ "CAS-CAS",
      Cover1 == "Pasture" & Cover2 == "Pasture" ~ "PA-PA",
      (Cover1 == "Secondary forest" & Cover2 == "Cacao agroforestry system") |
        (Cover1 == "Cacao agroforestry system" & Cover2 == "Secondary forest") ~ "SF-CAS",
      (Cover1 == "Secondary forest" & Cover2 == "Pasture") |
        (Cover1 == "Pasture" & Cover2 == "Secondary forest") ~ "SF-PA",
      (Cover1 == "Cacao agroforestry system" & Cover2 == "Pasture") |
        (Cover1 == "Pasture" & Cover2 == "Cacao agroforestry system") ~ "CAS-PA",
      TRUE ~ NA_character_
    ),
    PairType = factor(PairType, levels = c("SF-SF", "CAS-CAS", "PA-PA",
                                           "SF-CAS", "SF-PA", "CAS-PA"))
  )

# --- 10. Transform dissimilarities for beta regression -----------------------
N_pairs <- nrow(dyad)
dyad <- dyad %>%
  mutate(
    across(c(TaxDist, TaxSub, TaxNest, FunDist, FunSub, FunNest),
           ~ (.x * (N_pairs - 1) + 0.5) / N_pairs,
           .names = "{.col}_t")
  )

# --- 11. Beta GLMM for taxonomic total dissimilarity -------------------------
model_tax <- glmmTMB(TaxDist_t ~ PairType + (1 | Site1) + (1 | Site2),
                     data = dyad, family = beta_family(link = "logit"))
sim_tax <- simulateResiduals(model_tax, n = 250)
plot(sim_tax)
testUniformity(sim_tax)
testOutliers(sim_tax)
summary(model_tax)
drop1(model_tax, test = "Chisq")

# --- 12. Beta GLMM for functional total dissimilarity ------------------------
model_fun <- glmmTMB(FunDist_t ~ PairType + (1 | Site1) + (1 | Site2),
                     data = dyad, family = beta_family(link = "logit"))
sim_fun <- simulateResiduals(model_fun, n = 250)
plot(sim_fun)
testUniformity(sim_fun)
testOutliers(sim_fun)
summary(model_fun)
drop1(model_fun, test = "Chisq")

# --- 13. Contrasts of interest (total dissimilarities) -----------------------
emm_tax <- emmeans(model_tax, ~ PairType, type = "response")
emm_fun <- emmeans(model_fun, ~ PairType, type = "response")

cat("\n--- Taxonomic contrasts ---\n")
contrast(emm_tax, method = list(
  "SF-CAS vs SF-PA"  = c(0, 0, 0, 1, -1, 0),
  "SF-CAS vs CAS-PA" = c(0, 0, 0, 1, 0, -1),
  "SF-PA vs CAS-PA"  = c(0, 0, 0, 0, 1, -1)
), adjust = "bonferroni")

cat("\n--- Functional contrasts ---\n")
contrast(emm_fun, method = list(
  "SF-CAS vs SF-PA"  = c(0, 0, 0, 1, -1, 0),
  "SF-CAS vs CAS-PA" = c(0, 0, 0, 1, 0, -1),
  "SF-PA vs CAS-PA"  = c(0, 0, 0, 0, 1, -1)
), adjust = "bonferroni")

# --- 14. Post‑hoc letters for between‑cover pairs (total dissimilarity) -------
pairs_plot <- c("SF-CAS", "SF-PA", "CAS-PA")
emm_tax_sub <- emm_tax[which(levels(dyad$PairType) %in% pairs_plot)]
emm_fun_sub <- emm_fun[which(levels(dyad$PairType) %in% pairs_plot)]
cld_tax_sub <- cld(emm_tax_sub, Letters = letters, adjust = "bonferroni")
cld_fun_sub <- cld(emm_fun_sub, Letters = letters, adjust = "bonferroni")
print(cld_tax_sub)
print(cld_fun_sub)

# --- 15. Boxplots of total dissimilarities -----------------------------------
dyad_sub <- dyad %>% filter(PairType %in% pairs_plot)

add_letters <- function(data, value_col, cld_obj) {
  pos <- data %>%
    group_by(PairType) %>%
    summarise(
      ymax = quantile(.data[[value_col]], 0.75, na.rm = TRUE) +
        0.5 * IQR(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
  cld_df <- as.data.frame(cld_obj) %>%
    rename(PairType = PairType) %>%
    left_join(pos, by = "PairType") %>%
    mutate(.group = trimws(.group))
  cld_df
}

p_dyad_tax <- ggplot(dyad_sub, aes(x = PairType, y = TaxDist, fill = PairType)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_text(data = add_letters(dyad_sub, "TaxDist", cld_tax_sub),
            aes(x = PairType, y = ymax * 1.05, label = .group),
            inherit.aes = FALSE, family = "serif", size = 6) +
  scale_fill_manual(values = c("SF-CAS" = "darkorange4",
                               "SF-PA"  = "forestgreen",
                               "CAS-PA" = "springgreen")) +
  labs(x = "Site pair type", 
       y = "Taxonomic dissimilarity\nBray–Curtis") +
  theme_strategy +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 10, b = 10, l = 15))

p_dyad_fun <- ggplot(dyad_sub, aes(x = PairType, y = FunDist, fill = PairType)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_text(data = add_letters(dyad_sub, "FunDist", cld_fun_sub),
            aes(x = PairType, y = ymax * 1.05, label = .group),
            inherit.aes = FALSE, family = "serif", size = 6) +
  scale_fill_manual(values = c("SF-CAS" = "darkorange4",
                               "SF-PA"  = "forestgreen",
                               "CAS-PA" = "springgreen")) +
  labs(x = "Site pair type", 
       y = "Functional dissimilarity\nSørensen") + 
  theme_strategy +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 10, b = 10, l = 15))

# --- 16. Proportion of nestedness: between-pair comparisons per dimension -----
dyad <- dyad %>%
  mutate(PropTaxNest = TaxNest / (TaxSub + TaxNest),
         PropFunNest = FunNest / (FunSub + FunNest))

# Re‑create dyad_sub with the new columns
dyad_sub <- dyad %>% filter(PairType %in% pairs_plot)

# Ordinal beta models
model_tax_nest_ord <- glmmTMB(PropTaxNest ~ PairType + (1 | Site1) + (1 | Site2),
                              data = dyad, family = ordbeta(link = "logit"))
sim_tax_ord <- simulateResiduals(model_tax_nest_ord, n = 250)
plot(sim_tax_ord)
testUniformity(sim_tax_ord)
testOutliers(sim_tax_ord)

model_fun_nest_ord <- glmmTMB(PropFunNest ~ PairType + (1 | Site1) + (1 | Site2),
                              data = dyad, family = ordbeta(link = "logit"))
sim_fun_ord <- simulateResiduals(model_fun_nest_ord, n = 250)
plot(sim_fun_ord)
testUniformity(sim_fun_ord)
testOutliers(sim_fun_ord)

emm_prop_tax <- emmeans(model_tax_nest_ord, ~ PairType, type = "response")
emm_prop_fun <- emmeans(model_fun_nest_ord, ~ PairType, type = "response")

cat("\n--- Taxonomic nestedness contrasts ---\n")
contrast(emm_prop_tax, method = list(
  "SF-CAS vs SF-PA"  = c(0, 0, 0, 1, -1, 0),
  "SF-CAS vs CAS-PA" = c(0, 0, 0, 1, 0, -1),
  "SF-PA vs CAS-PA" = c(0, 0, 0, 0, 1, -1)
), adjust = "bonferroni")

cat("\n--- Functional nestedness contrasts ---\n")
contrast(emm_prop_fun, method = list(
  "SF-CAS vs SF-PA"  = c(0, 0, 0, 1, -1, 0),
  "SF-CAS vs CAS-PA" = c(0, 0, 0, 1, 0, -1),
  "SF-PA vs CAS-PA" = c(0, 0, 0, 0, 1, -1)
), adjust = "bonferroni")

# Post‑hoc letters for the between‑cover pairs (for the faceted plot)
emm_prop_tax_sub <- emm_prop_tax[which(levels(dyad$PairType) %in% pairs_plot)]
emm_prop_fun_sub <- emm_prop_fun[which(levels(dyad$PairType) %in% pairs_plot)]
cld_prop_tax <- cld(emm_prop_tax_sub, Letters = letters, adjust = "bonferroni")
cld_prop_fun <- cld(emm_prop_fun_sub, Letters = letters, adjust = "bonferroni")

# Prepare the data frame for the faceted plot
prop_plot_data <- dyad_sub %>%
  pivot_longer(c(PropTaxNest, PropFunNest),
               names_to = "Dimension", values_to = "PropNestedness") %>%
  mutate(Dimension = recode(Dimension,
                            PropTaxNest = "Taxonomic",
                            PropFunNest = "Functional"))

# Combine CLD results with dimension info
cld_prop_combined <- bind_rows(
  as.data.frame(cld_prop_tax) %>% mutate(Dimension = "Taxonomic"),
  as.data.frame(cld_prop_fun) %>% mutate(Dimension = "Functional")
) %>% rename(PairType = PairType) %>%
  left_join(prop_plot_data %>%
              group_by(Dimension, PairType) %>%
              summarise(ymax = quantile(PropNestedness, 0.75, na.rm = TRUE) +
                          0.5 * IQR(PropNestedness, na.rm = TRUE),
                        .groups = "drop"),
            by = c("Dimension", "PairType")) %>%
  mutate(.group = trimws(.group))

# Faceted boxplot – emphasises between‑pair comparisons within each dimension
p_prop_nest <- ggplot(prop_plot_data,
                      aes(x = PairType, y = PropNestedness, fill = PairType)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_text(data = cld_prop_combined,
            aes(x = PairType, y = ymax * 1.05, label = .group),
            inherit.aes = FALSE, family = "serif", size = 5) +
  scale_fill_manual(values = c("SF-CAS" = "lightblue",
                               "SF-PA"  = "cornflowerblue",
                               "CAS-PA" = "darkblue")) +
  facet_wrap(~ Dimension, ncol = 1, scales = "free_x") +
  labs(x = "Site pair type", y = "Proportion of nestedness") +
  theme_strategy +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        strip.text = element_blank(), 
        legend.position = "none")

# --- 17. Mantel -----------------------------------------------
cat("\n=== Mantel test ===\n")
mantel(tax_dist_total, fun_dist_total, method = "pearson", permutations = 9999)

# --- 18. Display all plots ---------------------------------------------------
p_nmds_tax <- p_nmds_tax + 
  theme(axis.title.x    = element_blank())

p_nmds_fun <- p_nmds_fun + 
  theme(legend.position = "none")

library(patchwork)
fig_tax2 <- (p_nmds_tax + p_nmds_fun) +
  plot_layout(ncol = 1)+
  plot_annotation(
    theme   = theme(plot.caption = element_text(hjust  = 0.5,
                                                family = "serif",
                                                size   = 20, face = "bold")))

ggsave("figures/Fig_tax2.svg", fig_tax2, width = 13, height = 12, dpi = 600)

p_dyad_tax <- p_dyad_tax +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

fig_tax3 <- (p_dyad_tax + p_dyad_fun) +
  plot_layout(ncol = 1)

ggsave("figures/Fig_tax3.svg", fig_tax3, width = 10, height = 8, dpi = 600)

print(p_prop_nest)

ggsave("figures/Fig_nest.svg", p_prop_nest, width = 10, height = 8, dpi = 600)
# =============================================================================
# END OF SCRIPT
# =============================================================================
