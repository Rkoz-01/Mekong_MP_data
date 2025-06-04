# polymer size distribution

library(tidyr)
library(ggplot2)
library(dplyr)
library(data.table)
library(readr)
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

head(data_all_cor)

# Aggregate data by width and length categories
heatmap_data <- data_all_cor %>%
  group_by(size_W, size_L) %>%
  summarise(total_N = sum(N_per_L, na.rm = TRUE))%>%
  mutate(size_W = ifelse(is.na(size_W), "<50", size_W))



#organize categories
heatmap_data <- heatmap_data %>%
  mutate(size_W = factor(size_W, levels = c("<50", "50-100", "100-200", "200-300", "300-400","400-500", "500-600", "600-700")),
         size_L = factor(size_L, levels = c("<50", "50-100", "100-200", "200-300", "300-400", "400-500", "500-600", "600-700", "700-800", "800-900", "900-1000", ">1000")))




# Create the heatmap for all plastics together
ggplot(heatmap_data, aes(x = size_L, y = size_W, fill = total_N)) +
  geom_tile() +
  labs(x = "Length Category (µm)", y = "Width Category (µm)", title = "Microplastic Size Distribution") +
  theme_minimal() +
  scale_fill_viridis_c(option = "plasma", name = "frequency")+  # "magma", "plasma", "inferno", "cividis"
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("Size variation all.jpg", width = 8, height = 5, dpi = 400)

#create heatmap for each season

# Calculate percentages within each season
heatmap_season <- data_all_cor %>%
  group_by(size_W, size_L, season) %>%
  summarise(total_N = sum(N_per_L, na.rm = TRUE))%>%
  mutate(size_W = ifelse(is.na(size_W), "<50", size_W))

heatmap_season <- heatmap_season %>%
  group_by(season) %>%
  mutate(percent_N_per_L = (total_N / sum(total_N)) * 100) %>%
  ungroup()

# Find the global min/max of percent_N_per_L for a consistent scale
global_min <- min(heatmap_season$percent_N_per_L, na.rm = TRUE)
global_max <- max(heatmap_season$percent_N_per_L, na.rm = TRUE)

# Define a common fill scale
common_fill_scale <- scale_fill_viridis_c(
  option = "plasma", 
  name = "Percentage of Total (%)",
  limits = c(global_min, global_max)  # Ensure both plots use the same range
)

heatmap_season <- heatmap_season %>%
  mutate(size_W = factor(size_W, levels = c("<50", "50-100", "100-200", "200-300", "300-400","400-500", "500-600", "600-700")),
         size_L = factor(size_L, levels = c("<50", "50-100", "100-200", "200-300", "300-400", "400-500", "500-600", "600-700", "700-800", "800-900", "900-1000", ">1000")))


# Base plot function
create_plot <- function(data, title_label,remove_y_axis = FALSE) {
  ggplot(data, aes(x = size_L, y = size_W, fill = percent_N_per_L)) +
    geom_tile() +
    common_fill_scale +  # Apply the common scale
    labs(x = "Length Category (µm)", y = expression("Width Category (" * mu * "m)")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = if (remove_y_axis) element_blank() else element_text(),
      axis.title.y = if (remove_y_axis) element_blank() else element_text(),
      plot.title = element_text(size = 16, face = "bold", hjust = 0),
      panel.spacing = unit(1, "cm"),
      legend.position = "bottom"
    ) +
    ggtitle(title_label)
}

# Create plots for each season
plot_dry <- create_plot(filter(heatmap_season, season == "dry"), "C", remove_y_axis = FALSE)
plot_wet <- create_plot(filter(heatmap_season, season == "wet"), "D", remove_y_axis = TRUE)

# Combine plots with a shared, standardized legend
final_plot <- (plot_dry + plot_wet) + 
  plot_layout(ncol = 2, guides = "collect") &  
  theme(legend.position = "none")  

# Display final plot
final_plot

ggsave("size_vatiation_no_legend.jpg", width = 8, height = 4, dpi = 400)



