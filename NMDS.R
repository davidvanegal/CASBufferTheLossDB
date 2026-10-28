# =============================================================================
# NMDS analysis — Dung beetles
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggforce)
library(indicspecies)
library(pairwiseAdonis)

# --- Data --------------------------------------------------------------------
dataDB <- read.csv("dataCASBufferDB.csv", header = TRUE, sep = ";")
dataDB$vegcover <- factor(dataDB$vegcover,
                          levels = c("SF", "CAS", "PA"),
                          labels = c("Secondary forest",
                                     "Cacao agroforestry system",
                                     "Pasture"))

# =============================================================================
# 1. COMMUNITY MATRIX — aggregated at site level
# One row per locality × vegcover × site (traps are pseudoreplicates)
# =============================================================================

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

# Species-only matrix (numeric)
comm_mat <- as.data.frame(comm_site[, setdiff(colnames(comm_site),
                                              c("locality", "vegcover",
                                                "site", "siteID"))])
rownames(comm_mat) <- comm_site$siteID
comm_mat <- as.matrix(comm_mat)
mode(comm_mat) <- "numeric"

# Metadata vectors for vegan
meta_vec <- metadata$vegcover
meta_df  <- as.data.frame(metadata)

# =============================================================================
# 2. TRANSFORMATIONS AND DISTANCES
# =============================================================================

# Hellinger transformation + Bray-Curtis (recommended for abundance data)
comm_hell <- decostand(comm_mat, method = "hellinger")
dist_bray <- vegdist(comm_hell, method = "bray")

# Presence-absence + Jaccard (composition only, no abundance weighting)
comm_pa   <- decostand(comm_mat, method = "pa")
dist_jacc <- vegdist(comm_pa, method = "jaccard")

# =============================================================================
# 3. NMDS
# =============================================================================

set.seed(42)
nmds <- metaMDS(comm_hell, distance = "bray", k = 2,
                trymax = 200, autotransform = FALSE)

cat("NMDS stress (k = 2):", round(nmds$stress, 4), "\n")
stressplot(nmds)

# Variance explained per axis (proportion of total dispersion)
scores_nmds <- scores(nmds, display = "sites")
var_axes    <- apply(scores_nmds, 2, var)
prop_axes   <- var_axes / sum(var_axes)
print(prop_axes)

# Site scores with metadata
scores_sites        <- as.data.frame(scores(nmds, display = "sites"))
scores_sites$siteID <- rownames(scores_sites)
scores_sites        <- left_join(scores_sites, metadata, by = "siteID")

# Group centroids
centroids <- scores_sites %>%
  group_by(vegcover) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = "drop")

# Spider segments: site to group centroid
segments <- left_join(scores_sites, centroids, by = "vegcover",
                      suffix = c(".pt", ".cent"))

# =============================================================================
# 4. FIGURE
# =============================================================================

shapes_nmds <- c("Secondary forest"          = 21,
                 "Cacao agroforestry system" = 22,
                 "Pasture"                   = 24)

cols <- c("Secondary forest"          = "forestgreen",
          "Cacao agroforestry system" = "darkorange4",
          "Pasture"                   = "springgreen")

theme_nmds <- theme_classic() +
  theme(
    axis.title      = element_text(size = 18, family = "serif", face = "bold"),
    axis.text       = element_text(size = 16, family = "serif"),
    legend.text     = element_text(size = 16, family = "serif"),
    legend.title    = element_text(size = 18, family = "serif", face = "bold"),
    legend.position = "top"
  )

p_nmds <- ggplot(scores_sites,
                 aes(x = NMDS1, y = NMDS2,
                     fill = vegcover, color = vegcover)) +
  # Spider lines: site to centroid
  geom_segment(data = segments,
               aes(x = NMDS1.pt, y = NMDS2.pt,
                   xend = NMDS1.cent, yend = NMDS2.cent),
               linetype = "dashed", linewidth = 0.4, alpha = 0.6) +
  # Confidence ellipses
  ggforce::geom_mark_ellipse(expand = 0, alpha = 0.15, show.legend = FALSE) +
  # Individual sites
  geom_point(aes(shape = vegcover), size = 3) +
  # Centroids (larger symbols)
  geom_point(data = centroids,
             aes(x = NMDS1, y = NMDS2, shape = vegcover),
             size = 5) +
  scale_fill_manual(values  = cols) +
  scale_color_manual(values = cols) +
  scale_shape_manual(values = shapes_nmds) +
  scale_x_continuous("NMDS1", expand = expansion(mult = c(0.05, 0.12))) +
  labs(x     = "NMDS1",
       y     = "NMDS2",
       color = "Vegetation cover",
       shape = "Vegetation cover",
       fill  = "Vegetation cover") +
  annotate("text",
           x = Inf, y = -Inf,
           hjust = 1.1, vjust = -1.5,
           label = paste0("Stress = ", round(nmds$stress, 3)),
           size = 6, family = "serif") +
  theme_nmds

print(p_nmds)

ggsave("NMDS/NMDS_DungBeetles.png", p_nmds,
       width = 12, height = 7, dpi = 600, bg = "white")
ggsave("NMDS/NMDS_DungBeetles.svg", p_nmds,
       width = 12, height = 7, dpi = 600, bg = "white")

# =============================================================================
# 5. BETADISPER — homogeneity of multivariate dispersion
# =============================================================================

bd <- betadisper(dist_bray, meta_vec)

print(anova(bd))
print(permutest(bd, permutations = 999))
print(TukeyHSD(bd))

# Dispersion boxplot
boxplot(bd$distances ~ meta_vec,
        xlab = "Vegetation cover",
        ylab = "Distance to centroid",
        col  = cols[levels(meta_vec)])

# =============================================================================
# 6. PERMANOVA — global test (Bray-Curtis and Jaccard)
# =============================================================================

adonis_bray <- adonis2(dist_bray ~ vegcover, data = meta_df, permutations = 9999)
print(adonis_bray)

adonis_jacc <- adonis2(dist_jacc ~ vegcover, data = meta_df, permutations = 9999)
print(adonis_jacc)

# =============================================================================
# 7. PAIRWISE PERMANOVA (Bray-Curtis and Jaccard)
# =============================================================================

pair_bray <- pairwise.adonis2(dist_bray ~ vegcover, data = meta_df,
                              permutations = 9999)
print(pair_bray)

pair_jacc <- pairwise.adonis2(dist_jacc ~ vegcover, data = meta_df,
                              permutations = 9999)
print(pair_jacc)

# =============================================================================
# 8. SIMPER — species contributions to dissimilarity
# =============================================================================

sim_all <- simper(comm_mat, meta_vec)
print(summary(sim_all))

# =============================================================================
# 9. INDICATOR SPECIES ANALYSIS
# =============================================================================

ind <- multipatt(comm_mat, meta_vec, func = "r.g",
                 control = how(nperm = 999))
print(summary(ind))
