# =============================================================================
# INDICATOR SPECIES & SIMPER — DUNG BEETLES
# =============================================================================

# --- 0. Load packages --------------------------------------------------------
library(dplyr)
library(tidyr)
library(vegan)
library(indicspecies)
library(ggplot2)

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
comm_mat <- as.matrix(comm_site[, setdiff(colnames(comm_site),
                                          c("locality", "vegcover", "site", "siteID"))])
rownames(comm_mat) <- comm_site$siteID
meta_vec <- metadata$vegcover

# --- 3. SIMPER ---------------------------------------------------------------
sim_all <- simper(comm_mat, meta_vec, permutations = 999)
# Average contribution across pairwise comparisons
contrib_list <- lapply(sim_all, function(x) data.frame(species = x$species, average = x$average))
avg_contrib <- bind_rows(contrib_list) %>%
  group_by(species) %>%
  summarise(mean_contrib = mean(average), .groups = "drop") %>%
  arrange(desc(mean_contrib))

# Plot top 10
top_sim <- head(avg_contrib, 10)
ggplot(top_sim, aes(x = reorder(species, mean_contrib), y = mean_contrib)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  coord_flip() +
  labs(x = NULL, y = "Average contribution to dissimilarity")+
  theme_classic(base_family = "serif") +
  theme(axis.text = element_text(size = 14))

# --- 4. Indicator species (IndVal.g) -----------------------------------------

# Run IndVal
ind <- multipatt(comm_mat, meta_vec, func = "r.g",
                 control = how(nperm = 999))

# Identify the cover columns
cover_cols <- grep("^s\\.", colnames(ind$sign), value = TRUE)

# Create a clean data frame of significant indicators
ind_df <- data.frame(
  species = rownames(ind$sign),
  stat    = ind$sign$stat,
  p.value = ind$sign$p.value,
  cover   = cover_cols[apply(ind$sign[, cover_cols], 1, which.max)],
  stringsAsFactors = FALSE
) %>%
  filter(p.value <= 0.05) %>%
  mutate(cover = gsub("^s\\.", "", cover),
         cover = factor(cover, levels = c("Secondary forest",
                                          "Cacao agroforestry system",
                                          "Pasture")))

# Print to verify
print(ind_df)

# --- Add body length information for functional interpretation ---------------
sp_length <- dataDB %>%
  group_by(species) %>%
  summarise(body_length_mm = mean(length, na.rm = TRUE), .groups = "drop")

ind_df <- left_join(ind_df, sp_length, by = "species")

# --- Plot --------------------------------------------------------------------
cols_ind <- c("Secondary forest"          = "forestgreen",
              "Cacao agroforestry system" = "darkorange4",
              "Pasture"                   = "springgreen")

p_ind <- ggplot(ind_df, aes(x = reorder(species, stat), y = stat, fill = cover)) +
  geom_col(alpha = 0.85) +
  geom_text(aes(label = paste0(round(body_length_mm, 1), " mm")),
            hjust = -0.1, size = 4, family = "serif") +
  scale_fill_manual(values = cols_ind, name = "Indicator for") +
  coord_flip() +
  labs(x = NULL,
       y = "Indicator value (r.g)",
       title = "Significant indicator species (IndVal.g)")+
  theme_strategy +
  theme(
    axis.text       = element_text(size = 14),
    axis.title      = element_text(size = 16, face = "bold"),
    legend.position = "top"
  )

print(p_ind)
ggsave("figures/indspecies.svg", p_ind, width = 11, height = 8, dpi = 600)