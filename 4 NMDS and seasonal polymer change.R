# polymer composition comparison

library(tidyr)
library(ggplot2)
library(dplyr)
library(data.table)
library(readr)
library(vegan)
library(patchwork)



###_____Prepare data

setwd("C:/Users/rachelk/OneDrive - Desert Research Institute/Desktop/Mekong_data_share")

#load in dry data
dry_data_335_cor <- read_csv("Match_category_335_final_blank_corrected_dry_season.csv")
dry_data_100_cor <- read_csv("Match_category_100_final_blank_corrected_dry_season.csv")

#add season column
dry_data_100_cor <- dry_data_100_cor %>%
  mutate(season = "dry") %>%
  mutate(filter = "100um")

#add season column
dry_data_335_cor <- dry_data_335_cor %>%
  mutate(season = "dry") %>%
  mutate(filter = "335um")

#load in wet data
wet_data_100_cor <- read_csv("Match_category_100_final_blank_corrected_wet_season.csv")
wet_data_335_cor <- read_csv("Match_category_335_final_blank_corrected_wet_season.csv")

#add season column 
wet_data_100_cor <- wet_data_100_cor%>%
  mutate(season = "wet")%>%
  mutate(filter = "100um")

#add season column
wet_data_335_cor <- wet_data_335_cor %>%
  mutate(season = "wet") %>%
  mutate(filter = "335um")

# Combine the results
data_all_cor<- bind_rows(dry_data_100_cor, dry_data_335_cor, wet_data_100_cor, wet_data_335_cor)

# plastic counts per L for all sites
Plastic_counts_total <- data_all_cor%>%
  group_by(site, season) %>%
  summarize(total_count = sum(corrected_N), na.rm = TRUE)

# plastic counts per L for all sites
Plastic_counts_filter_total <- data_all_cor%>%
  group_by(site, season, filter) %>%
  summarize(total_count = sum(corrected_N), na.rm = TRUE)

percent_data <- data_all_cor %>%
  group_by(season, site) %>%
  summarise(total_count = sum(N_per_L), .groups = "drop") %>%
  left_join(
    data_all_cor %>%
      group_by(season, site, match_category_corrected) %>%
      summarise(total_polymer = sum(N_per_L), .groups = "drop"),
    by = c("season", "site")
  ) %>%
  mutate(percent = total_polymer / total_count)


####_____NMDS (non-metric multidimensional scaling)
#NMDS is a rank-based, distance-based method that preserves relative differences using Bray-Curtis dissimilarity
#(which accounts for both abundance and presence/absence of polymer types). NMDS also ignores double zeros so its a 
#good fit for data where polymers might not be present at all sites. Sites plotted closer together have more similar polymer profiles. 

#Make a site-by-polymer matrix
# Summarize counts per polymer type per site and season
nmds_data <- data_all_cor %>%
  group_by(site, season, match_category_corrected) %>%
  summarise(total_count = sum(N_per_L), .groups = "drop") %>%
  pivot_wider(names_from = match_category_corrected, values_from = total_count, values_fill = 0)

##__compute Bray-Curtis Dissimilarity
# Remove non-numeric columns
nmds_matrix <- nmds_data %>%
 dplyr::select(-site, -season) %>%
  as.matrix()

# Compute Bray-Curtis distance
bray_dist <- vegdist(nmds_matrix, method = "bray")


##___run NMDS
# Run NMDS (k=2 for a 2D plot)
nmds_result <- metaMDS(bray_dist, k = 4, trymax = 100)

# Check NMDS stress (should be < 0.2 for good fit)
nmds_result$stress

##____visualize the results
# Extract site and season info for plotting
nmds_sites <- nmds_data %>% 
  dplyr::select(site, season)

# Convert NMDS points to a tibble
nmds_points <- as_tibble(nmds_result$points) %>%
  bind_cols(nmds_sites)

# NMDS plot with ggplot
NMDS <- ggplot(nmds_points, aes(x = MDS1, y = MDS2, fill = season)) +
  geom_point(shape = 21, size = 4, color = "black", stroke = 1) +
  stat_ellipse(aes(color = season), show.legend = FALSE) +  # Hides ellipse legend
  theme_minimal() +
  theme(legend.position = "bottom", legend.margin = margin(t = 10, b = 10)) +
  labs(title = "", x = "NMDS1", y = "NMDS2", fill = "Season") +
  scale_fill_manual(values = c("tan", "forestgreen")) +
  scale_color_manual(values = c("tan", "forestgreen"))

NMDS
#ggsave("NMDS_seasonal_polymer_var_blank_corrected.jpg", width = 4, height = 4, dpi = 400)

##___ Test for statistical differences
# PERMANOVA - test is composition differs by Season/site (Pr (>F) is the P value). 
adonis2(bray_dist ~ season + site, data = nmds_data, permutations = 999)
# the model explains about 60% of the variability, but not enough to be significant. 
#so the seasonal difference in plastic composition for each site is not significant

#test if polymer composition varies significantly by season 
dispersion_test <- betadisper(bray_dist, nmds_data$season)
permutest(dispersion_test)
# The permutation dispersion test has a highly significant P value. This 
# means that the variability in polymer type is higher in the dry season

#Box plot visualizing the dispersion
boxplot(dispersion_test$distances ~ nmds_data$season, main = "Dispersion by Season",
        ylab = "Distance to Centroid", 
        xlab = "",
        col = c("tan", "forestgreen"))

# Similarity percentage analysis shows what polymer types contribute most to the 
#differences in microplastic composition between wet and dry seasons
simper_result <- simper(nmds_matrix, nmds_data$season)
summary(simper_result)
# PP, PET, and PE are significantly higher in the wet season compared to the dry season 
# and they account 


#####_____Bar plots showing main polymers associated with seasonal variation

# Filter data for significant polymer types from SIMPER
significant_polymers <- c("Polyterephthalates", "Polypropylene")

plot_data <- data_all_cor %>%
  group_by(site, season, match_category_corrected) %>%
  summarise(N_per_L = sum(N_per_L), .groups = "drop") 

filtered_data <- plot_data %>%
  filter(match_category_corrected %in% significant_polymers) %>%
  mutate(season = factor(season))  # Ensure season is categorical

# Create the boxplot with properly spread points
boxplot <- ggplot(filtered_data, aes(x = season, y = N_per_L, fill = match_category_corrected)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8)) +  
  geom_jitter(aes(fill = match_category_corrected), 
              shape = 21, color = "black", stroke = 0.5,
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), 
              alpha = 0.6, size = 2) +    
  scale_fill_manual(values = c("Polypropylene" = "#E69F00", 
                               "Polyterephthalates" = "#56B4E9" 
                               #"Polyethylenes" = "#009E73"
                               )) +
  labs(title = "",
       y = expression( "n L" ^ "-1"), 
       x = "Season", fill = "Polymer\ntype", color = "Polymer\ntype") +
  #ylim(0, 1) +  # Set y-axis limits
  theme_minimal() +
  theme(legend.position = "bottom") + 
  guides(fill = guide_legend(nrow = 2), color = guide_legend(nrow = 2))

boxplot

#ggsave("Plastic_types_driving_seasonal_var_blank_corrected.jpg", width = 7.5, height = 5, dpi = 400)

NMDS <- NMDS + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
boxplot <- boxplot + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

combined_plot <- (NMDS | boxplot) + 
  plot_layout(widths = c(1, 1)) + 
  plot_annotation(tag_levels = 'A', tag_prefix = "")

combined_plot
#ggsave("NMDS_and_polymer_change_corrected.jpg", width = 8, height = 5, dpi = 400)

