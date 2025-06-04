#compare dry and wet season

library(tidyr)
library(ggplot2)
library(dplyr)
library(data.table)
library(NatParksPalettes)
library(readr)
library(lme4) #for glm models
library(MASS)
library(patchwork)
library(circlize)#chord diagrams
library(grid)
library(ggbreak)


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

data_all_cor<- data_all_cor%>%
  mutate(match_category_corrected = recode(match_category_corrected, 
                                           "Polyamides (polylactams)" = "Polyamides",
                                           "Polysiloxanes" = "Other plastic",
                                           "Polyacrylamides" = "Other plastic"))

#Dry season only 
data_all_cor_filtered_dry <- data_all_cor %>%
  filter(season == "dry") #%>%
  #filter(!match_category_corrected %in% c("Rubber", "Silicon Rubber"))

#wet season only 
data_all_cor_filtered_wet <- data_all_cor %>%
  filter(season == "wet") #%>%
  #filter(!match_category_corrected %in% c("Rubber", "Silicon Rubber"))


###_____________ wet season Average duplicate polymer counts
# Step 1: Group by site and match_category_corrected, then calculate the mean
averaged_SS1 <- data_all_cor_filtered_wet %>%
  filter(site %in% c("SS1", "SS1D")) %>%
  group_by(site, match_category_corrected) %>%
  summarise(
    N_per_L = mean(N_per_L, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(site = gsub("D$", "", site))  # Remove the 'D' from the site name

# Step 2: Do the same for M3 and M3D
averaged_M3 <- data_all_cor_filtered_wet %>%
  filter(site %in% c("M3", "M3D")) %>%
  group_by(site, match_category_corrected) %>%
  summarise(
    N_per_L = mean(N_per_L, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(site = gsub("D$", "", site))  # Remove the 'D' from the site name

# Step 3: Remove SS1, SS1D, M3, and M3D from the original dataset
data_cleaned <- data_all_cor_filtered_wet %>%
  filter(!site %in% c("SS1", "SS1D", "M3", "M3D"))

# Step 4: Combine the cleaned dataset with the averaged rows
averaged_wet_data <- bind_rows(data_cleaned, averaged_SS1, averaged_M3)


# plastic counts per L for all sites
Plastic_counts_all <- data_all_cor%>%
  group_by(site, season) %>%
  summarize(count_per_L = sum(N_per_L), na.rm = TRUE)

#average concentrations
average_by_season <- Plastic_counts_all %>%
  #filter(!(site == "SS1")) %>%
  group_by(season)%>%
  summarize(avg_NperL = mean(count_per_L), na.rm = TRUE)

# Define color-blind friendly colors for each category
cb_palette <- c(
  "Polypropylene" = "#E69F00",       # Orange
  "Other plastic" = "#56B4E9",       # Sky Blue
  "Polyethylenes" = "#009E73",       # Green
  "Polyterephthalates" = "#F0E442",  # Yellow
  "Polyamides" = "#0072B2", # Blue
  "Polyvinyl chlorides" = "#CC79A7", # Reddish Purple
  "Polystyrenes" = "#B3B3B3",        # light Gray
  "Polycarbonates" = "#882255",       # Dark Red
  "Rubber" = "#3B3B3B",    # dark gray for Rubber
  "Silicon Rubber" = "#8B4513"  # Brown for Silicon Rubber      
)

# Define the order of sites
site_order <- c("SK1", "SS1","SS1D", "M1", "M2","M3","M3D", "M4", "M5", "M6", "M7", "M8", "M9", "B1", "B2")

unique(data_all_cor$match_category_corrected)


####______dry season plots, percent and count combined


###____Make percent plot
# Add a percent column for each process

dry_season_polymer_percent <- data_all_cor_filtered_dry %>%
  group_by(site) %>%  # Group by site to get the total for each site
  mutate(total_N_per_L = sum(N_per_L)) %>%  # Calculate total N_per_L per site
  ungroup() %>%
  group_by(site, match_category_corrected) %>%
  summarise(percent = (sum(N_per_L) / unique(total_N_per_L)) * 100, .groups = "drop")  # Compute percentage


# Create bar graph for summary by material type (standard version)
dry_percent <- ggplot(data = dry_season_polymer_percent, aes(x = factor(site, levels = site_order), y = percent, fill = match_category_corrected)) +
  geom_bar(stat = "identity", width = 0.8) +
  #geom_text(aes(label = round(total_values, 1)), position = position_dodge(width = 0.8), vjust = -0.5, size = 3) + # Add rounded count values
  theme_minimal(base_size = 12) +
  labs(title = "",
       x = "Site",
       y = "Percent",
       fill = "Polymer Type"
  ) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = cb_palette) +
  geom_vline(xintercept = 9.5, linetype = "dashed", color = "black") #+ #adjust location of Phnom Penh line (between M7 and M8). 
  #annotate("text", x = 9.5, y = Inf, label = "Phnom Penh", vjust = 2, hjust = 1.15, color = "black")


###____ Make count (confidence interval) plot

dry_season_polymer_total_cor<- data_all_cor_filtered_dry %>%
  group_by(site) %>%  # Group by site to get the total for each site
  mutate(total_N_per_site = sum(corrected_N))  # Calculate total N_per_L per site

#load in dry data
dry_data__95CI_cor <- read_csv("95_CI_dry_season_cor.csv")

#add season column
dry_data__95CI_cor <- dry_data__95CI_cor %>%
  mutate(season = "dry")

dry_CI <- ggplot(dry_data__95CI_cor, aes(x = factor(Site, levels = site_order), y = PL)) +
  geom_point(position = position_dodge(width = 0.5), size = 1.5) +
  geom_errorbar(aes(ymin = `95_CIL_PL`, ymax = `95_CIU_PL`), 
                position = position_dodge(width = 0.5), width = 0.2) +
  theme_minimal(base_size = 12) +
  labs(x = "Site", y = expression("n L" ^ "-1")) +
  theme(
    legend.position = "bottom",
    axis.text.y.right = element_blank(),  
    axis.ticks.y.right = element_blank(), 
    axis.title.y.right = element_blank(),
    panel.spacing = unit(0.5, "cm")  # Adjust panel spacing to balance layout
  ) +
  geom_vline(xintercept = 9.5, linetype = "dashed", color = "black") +
  
  scale_y_break(c(8, 23), scales = "fixed", space = 0.2) +  # Increase space for a wider break
  scale_y_continuous(breaks = c(2, 4, 6, 8, 22, 24, 26, 28), expand = c(0, 0))   
  
  #annotate("text", x = 12, y = 26, label = "Phnom Penh", vjust = -0.5, hjust = 1.15, color = "black") +
  
  #coord_cartesian(clip = "off")  # Ensures text annotation is not clipped

dry_CI


# Combine the two plots vertically
combined_plot <- dry_CI /dry_percent + 
  plot_layout(widths = c(1, 1)) + 
  plot_annotation(tag_levels = 'A') +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) #+  # Ensure no additional margin

combined_plot

#ggsave("Plastic_concentrations_and_type_dry_season__corrected.jpg", width = 11, height = 7, dpi = 400)



###_________wet season polymers

###____Make percent plot
# Add a percent column for each process

wet_season_polymer_percent <- data_all_cor_filtered_wet %>%
  group_by(site) %>%  # Group by site to get the total for each site
  mutate(total_N_per_L = sum(N_per_L)) %>%  # Calculate total N_per_L per site
  ungroup() %>%
  group_by(site, match_category_corrected) %>%
  summarise(percent = (sum(N_per_L) / unique(total_N_per_L)) * 100, .groups = "drop")  # Compute percentage


# Create bar graph for summary by material type (standard version)
wet_percent <- ggplot(data = wet_season_polymer_percent, aes(x = factor(site, levels = site_order), y = percent, fill = match_category_corrected)) +
  geom_bar(stat = "identity", width = 0.8) +
  #geom_text(aes(label = round(total_values, 1)), position = position_dodge(width = 0.8), vjust = -0.5, size = 3) + # Add rounded count values
  theme_minimal(base_size = 12) +
  labs(title = "",
       x = "Site",
       y = "Percent",
       fill = "Polymer Type"
  ) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = cb_palette) +
  geom_vline(xintercept = 10.5, linetype = "dashed", color = "black") #+ #adjust location of Phnom Penh line (between M7 and M8). 
#annotate("text", x = 9.5, y = Inf, label = "Phnom Penh", vjust = 2, hjust = 1.15, color = "black")

wet_percent

###____ Make count (confidence interval) plot

wet_season_polymer_total_cor<- data_all_cor_filtered_wet %>%
  group_by(site) %>%  # Group by site to get the total for each site
  mutate(total_N_per_site = sum(corrected_N))  # Calculate total N_per_L per site

#load in wet data
wet_data__95CI_cor <- read_csv("95_CI_wet_season_cor.csv")

#add season column
wet_data__95CI_cor <- wet_data__95CI_cor %>%
  mutate(season = "wet")

wet_CI <- ggplot(wet_data__95CI_cor, aes(x = factor(Site, levels = site_order), y = PL)) +
  geom_point(position = position_dodge(width = 0.5), size = 1.5) +
  geom_errorbar(aes(ymin = `95_CIL_PL`, ymax = `95_CIU_PL`), 
                position = position_dodge(width = 0.5), # Ensure same dodge as points
                width = 0.2) + # Adjust bar width as needed
  theme_minimal(base_size = 12) +
  labs(
    title = "",
    x = "Site",
    y = expression("n L" ^ "-1")) +
  theme(
    legend.position = "bottom"
  ) +
  geom_vline(xintercept = 10.5, linetype = "dashed", color = "black") 
  #annotate("text", x = 13, y = Inf, label = "Phnom Penh", vjust = 2, hjust = 1.15, color = "black")



# Combine the two plots vertically
combined_plot <- wet_CI /wet_percent + 
  plot_layout(widths = c(1, 12)) + 
  plot_annotation(tag_levels = 'A') +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) #+  # Ensure no additional margin

combined_plot

ggsave("Plastic_concentrations_and_type_wet_season_corrected.jpg", width = 9, height = 7, dpi = 400)

###combine wet and dry season

wet_percent <- wet_percent + theme(legend.position = "none")

# Ensure dry_percent has the legend positioned at the bottom
dry_percent <- dry_percent + theme(legend.position = "none")

# Combine the plots
combined_plot_all <- free(dry_CI | wet_CI) / (dry_percent | wet_percent) + 
  plot_layout(widths = c(1, 2)) + 
  plot_annotation(tag_levels = 'A') +
  theme(
    plot.tag.position = c(0.02, 1),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))  # Ensure no extra margin


combined_plot_all

#ggsave("Plastic_concentrations_and_type_all_corrected.jpg", width = 12, height = 9, dpi = 400)


####_______dry sesason chord diagram

# Select relevant columns
dry_chord <- data_all_cor_filtered_dry %>%
  dplyr::select(match_category_corrected, morphology, N_per_L) %>%
  group_by(match_category_corrected, morphology) %>%
  summarise(N_per_L = sum(N_per_L, na.rm = TRUE), .groups = "drop")  # Summing if duplicates exist

# Convert to wide format for matrix creation
chord_data_dry <- dry_chord %>%
  tidyr::pivot_wider(names_from = morphology, values_from = N_per_L, values_fill = 0)

# Convert to matrix (excluding the first column which is category names)
matrix_data_dry <- as.matrix(chord_data_dry[, -1])
rownames(matrix_data_dry) <- chord_data_dry$match_category_corrected

Grid.col = c( "Polypropylene" = "#E69F00",      
              "Other plastic" = "#56B4E9",     
              "Polyethylenes" = "#009E73",      
              "Polyterephthalates" = "#F0E442", 
              "Polystyrenes" = "#B3B3B3",        
              "Polycarbonates" = "#882255",
              "Rubber" = "#3B3B3B",
              "fiber/filament" = "violet", 
              "fragment" = "navy")


#remove film category and minor plastics (PVC, Polyamides)
matrix_data_dry <- matrix_data_dry[, !(colnames(matrix_data_dry) %in% c("film", "sphere/spheroid"))]
matrix_data_dry_1 <- matrix_data_dry[!(rownames(matrix_data_dry) %in% c("Polyvinyl chlorides","Polyamides", "Polyamides (polylactams)", "Polyacrylamides", "Silicon Rubber")), ]
matrix_data_dry_2 <- matrix_data_dry[!(rownames(matrix_data_dry) %in% c("Polyvinyl chlorides", "Polyamides (polylactams)", "Other plastic")), ]


circos.clear()

# Adjust label orientation to be perpendicular to the circle
circos.par(track.margin = c(0.04, 0.04), gap.degree = 2) #gap degree changes the white space between each sector

# Set margins to zero for better fitting (if needed)
par(mar = c(4, 4, 4, 4), xpd = NA)

pdf("chord_diagram_dry_cor.pdf", height = 4.5, width = 4.5)
# Create the chord diagram 
chordDiagram(matrix_data_dry_1, 
             transparency = 0.5, 
             grid.col = Grid.col,
             annotationTrack = "grid")

# Adjust text direction to be perpendicular
circos.track(track.index = 1, panel.fun = function(x, y) {
  name <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  
  # Custom labels for plastics
  custom_labels <- c(
    "Polypropylene" = "PP",
    "Polyethylenes" = "PE",
    "Polyterephthalates" = "PET",
    "Polystyrenes" = "PS",
    "Polycarbonates" = "PC"
  )
  # Check if the label is in the custom list and replace it
  new_label <- ifelse(name %in% names(custom_labels), custom_labels[name], name)
  
  #text direction (dd) and adjusmtents (aa)
  theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
  #dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
  aa = c(.5, 0.5)
  if(theta < 90 || theta > 270)  aa = c(0, 0.5)
  
  #plot plastic labels
  circos.text(x=mean(xlim), y=1.7, labels=new_label, facing = "bending.inside", niceFacing = TRUE,  cex=0.6,  adj = aa)

},bg.border = NA)

# Close the PDF device
dev.off()

####_______wet season chord diagram

# Select relevant columns
wet_chord <- data_all_cor_filtered_wet %>%
  dplyr::select(match_category_corrected, morphology, N_per_L) %>%
  group_by(match_category_corrected, morphology) %>%
  summarise(N_per_L = sum(N_per_L, na.rm = TRUE), .groups = "drop")  # Summing if duplicates exist

# Convert to wide format for matrix creation
chord_data_wet <- wet_chord %>%
  tidyr::pivot_wider(names_from = morphology, values_from = N_per_L, values_fill = 0)

# Convert to matrix (excluding the first column which is category names)
matrix_data_wet <- as.matrix(chord_data_wet[, -1])
rownames(matrix_data_wet) <- chord_data_wet$match_category_corrected

Grid.col = c( "Polypropylene" = "#E69F00",      
              "Other plastic" = "#56B4E9",     
              "Polyethylenes" = "#009E73",      
              "Polyterephthalates" = "#F0E442", 
              "Polystyrenes" = "#B3B3B3",  
              "Polycarbonates" = "#882255",
              "Rubber" = "#3B3B3B",
              "fiber/filament" = "violet", 
              "fragment" = "navy",
              "film" = "darkslategray")


#remove film category and minor plastics (PVC, Polyamides)
matrix_data_wet <- matrix_data_wet[, !(colnames(matrix_data_wet) %in% "sphere/spheroid")]
matrix_data_wet_1 <- matrix_data_wet[!(rownames(matrix_data_wet) %in% c("Polysiloxanes", "Polyamides (polylactams)","Polyamides", "Polyvinyl chlorides", "Polycarbonates", "Silicon Rubber")), ]
#matrix_data_dry_2 <- matrix_data_dry[!(rownames(matrix_data_dry) %in% c("Polysiloxanes", "Polyamides (polylactams)", "Other plastic")), ]


circos.clear()

# Adjust label orientation to be perpendicular to the circle
circos.par(track.margin = c(0.04, 0.04), gap.degree = 2) #gap degree changes the white space between each sector

# Set margins to zero for better fitting (if needed)
par(mar = c(4, 4, 4, 4), xpd = NA)

pdf("chord_diagram_wet_cor.pdf", height = 4.5, width = 4.5)
# Create the chord diagram 
chordDiagram(matrix_data_wet_1, 
             transparency = 0.5, 
             grid.col = Grid.col,
             annotationTrack = "grid")

# Adjust text direction to be perpendicular
circos.track(track.index = 1, panel.fun = function(x, y) {
  name <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  
  # Custom labels for plastics
  custom_labels <- c(
    "Polypropylene" = "PP",
    "Polyethylenes" = "PE",
    "Polyterephthalates" = "PET",
    "Polystyrenes" = "PS",
    "Polyvinyl chlorides" = "PVC"
  )
  # Check if the label is in the custom list and replace it
  new_label <- ifelse(name %in% names(custom_labels), custom_labels[name], name)
  
  #text direction (dd) and adjusmtents (aa)
  theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
  #dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
  aa = c(.5, 0.5)
  if(theta < 90 || theta > 270)  aa = c(0, 0.5)
  
  #plot plastic labels
  circos.text(x=mean(xlim), y=1.7, labels=new_label, facing = "bending.inside", niceFacing = TRUE,  cex=0.6,  adj = aa)
  
},bg.border = NA)

# Close the PDF device
dev.off()

 
# Combine the two chord plots horixontally
library(magick)
library(pdftools)


# Read images
dry_cor_chord <- image_read("chord_diagram_dry_cor.jpg")
wet_cor_chord <- image_read("chord_diagram_wet_cor.jpg")

# Add letter annotations before combining
dry_cor_chord <- image_annotate(dry_cor_chord, "A", size = 50, gravity = "northwest", location = "+50+50")
wet_cor_chord <- image_annotate(wet_cor_chord, "B", size = 50, gravity = "northwest", location = "+50+50")  # Adjust as needed

# Combine images side by side
combined <- image_append(c(dry_cor_chord, wet_cor_chord), stack = FALSE)

# Display and save the result
combined

image_write(combined, "chord_dry+wet_cor.jpg")

