#Microplastics Color comparison

library(dplyr)
library(ggplot2)
library(readr)
library(forcats)
library(stringr)

setwd("C:/Users/rachelk/OneDrive - Desert Research Institute/Desktop/Mekong_data_share")

#load in MP color data  (these datasets are only plastics, no natural materials)
dry_MP_100_colors <- read_csv("100_dry_season_colors.csv")
dry_MP_335_colors <- read_csv("335_dry_season_colors.csv")
wet_MP_100_colors <- read_csv("100_wet_season_colors.csv")
wet_MP_335_colors <- read_csv("335_wet_season_colors.csv")


# Combine datasets and add object name
combined_data <- bind_rows(
  dry_MP_100_colors %>% mutate(object_name = "100um, Dry Season"),
  dry_MP_335_colors %>% mutate(object_name = "335um, Dry Season"),
  wet_MP_100_colors %>% mutate(object_name = "100um, Wet Season"),
  wet_MP_335_colors %>% mutate(object_name = "335um, Wet Season"),
)

# Group, summarize, and calculate percentages
summary_data <- combined_data %>%
  group_by(object_name, plastic_or_not, color) %>%
  summarise(total_N = sum(N, na.rm = TRUE), .groups = "drop") %>%
  group_by(object_name, plastic_or_not) %>%
  mutate(percentage = total_N / sum(total_N) * 100)

# Concatenate object_name and plastic_or_not
summary_data <- summary_data %>%
  filter(!is.na(plastic_or_not)) %>%
  filter(!is.na(color))%>%
  mutate(plastic_or_not = if_else(plastic_or_not == "not plastic", 
                                  "natural material", 
                                  "plastic"),
         object_name_plastic = str_c(object_name, plastic_or_not, sep = " - "))

summary_data <- summary_data %>%
  mutate(object_name = factor(object_name_plastic, 
                              levels = c("100um, Dry Season - natural material",
                                         "100um, Dry Season - plastic",
                                         "100um, Wet Season - natural material",
                                         "100um, Wet Season - plastic",
                                         "335um, Dry Season - natural material",
                                         "335um, Dry Season - plastic",
                                         "335um, Wet Season - natural material",  
                                         "335um, Wet Season - plastic")))

natural_and_plastic_counts<- summary_data %>%
  group_by(plastic_or_not)%>%
  summarise(total_N = sum(total_N, na.rm = TRUE), .groups = "drop")

# Define custom color mapping
color_mapping <- c(
  "black" = "black",
  "brown" = "chocolate4",
  "purple" = "purple",      
  "red" = "red",
  "translucent" = "burlywood1", # 'translucent' can use a light gray
  "white" = "white",
  "yellow" = "yellow",
  "blue" = "deepskyblue",
  "green" = "chartreuse4",
  "gold" = "gold",
  "silver" = "gray",
  "multi" = "deeppink"
)

# Reorder object_name based on the total percentage
summary_data <- summary_data %>%
  mutate(color = factor(color, levels = c("purple","silver", "gold", "yellow", "red", "multi", "green", "blue","white", "black","translucent","brown")))

# Create a horizontal bar plot with the custom color scale
ggplot(summary_data, aes(x = percentage, y = object_name_plastic, fill = color)) +
  geom_bar(stat = "identity") +
  labs(title = "",
       x = "Color counts as a percentage of total particles",
       y = "") +
  theme_minimal() +
  #coord_flip() +
  scale_fill_manual(values = color_mapping)

ggsave("Color_comparison_plastics_only.jpg", width = 8, height = 4.5, dpi = 400)



# fisher's exact test to compare the color distribution between plastic and natural materials for each size and season. 
# this test is good for data sets where there are lot of cells with small counts.   
fisher.test(table(summary_data$color, summary_data$object_name_plastic))


### color plot for just plastics and naturals

natural_and_plastic_counts_with_colors<- summary_data %>%
  group_by(plastic_or_not, color)%>%
  summarise(total_N = sum(total_N, na.rm = TRUE), .groups = "drop")%>%
  group_by(plastic_or_not) %>%
  mutate(percentage = total_N / sum(total_N) * 100)

# Create a horizontal bar plot with the custom color scale
ggplot(natural_and_plastic_counts_with_colors, aes(x = percentage, y = plastic_or_not, fill = color)) +
  geom_bar(stat = "identity") +
  labs(title = "",
       x = "Color counts as a percentage of total particles",
       y = "") +
  theme_minimal() +
  #coord_flip() +
  scale_fill_manual(values = color_mapping)+
  theme(
    legend.position = "none"
  )

ggsave("Color_comparison_plastics_only_small.jpg", width = 5, height = 2, dpi = 400)


fisher.test(table(natural_and_plastic_counts_with_colors$color, natural_and_plastic_counts_with_colors$plastic_or_not))
