#compare dry and wet season

library(tidyr)
library(ggplot2)
library(dplyr)
library(data.table)
library(NatParksPalettes)
library(readr)
library(lme4) #for glm models
library(MASS)


#update.packages(ask = FALSE)

setwd("C:/Users/rachelk/OneDrive - Desert Research Institute/Desktop/Mekong_data_share")

#load in dry data
dry_data_100_cor <- read_csv("Match_category_100_final_blank_corrected_dry_season.csv")
dry_data_335_cor <- read_csv("Match_category_335_final_blank_corrected_dry_season.csv")

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
data_all_cor<- bind_rows(dry_data_335_cor, wet_data_335_cor, dry_data_100_cor, wet_data_100_cor)


# plastic counts per L for all sites
Plastic_counts_all <- data_all_cor%>%
  group_by(site, season) %>%
  summarize(count_per_L = sum(N_per_L), na.rm = TRUE)

# plastic counts per L for all sites
Plastic_counts_all_m3 <- data_all_cor %>%
  group_by(site, season) %>%
  summarize(count_per_m3 = round(sum(N_per_L, na.rm = TRUE) * 1000),
            .groups = "drop")

#Set duplicate names to the same as the primary sample and combine duplicates (average of counts per L)
Plastic_counts_all_clean <- Plastic_counts_all_m3 %>%
  mutate(site = gsub("D", "", site)) %>%
  group_by(site, season) %>%
  summarise(across(everything(), ~round(mean(.x, na.rm = TRUE))), .groups = "drop") # Take the mean of all numeric columns

#add river column 
Plastic_counts_all_clean <- Plastic_counts_all_clean %>%
  mutate(river = case_when(
    grepl("B", site) ~ "Bassac",
    grepl("M", site) ~ "Mekong",
    grepl("SK", site) ~ "Sekong",
    grepl("SS", site) ~ "Sesan",
    TRUE ~ "Unknown"  # In case no match is found
  ))

#wet/dry season comparison
seasonal_variation <- Plastic_counts_all_clean %>%
  dplyr::select(site, season, count_per_m3) %>%  # Ensure only necessary columns are used
  pivot_wider(names_from = season, values_from = count_per_m3) %>%
  mutate(seasonal_variation = wet - dry) %>%  # Compute seasonal variation
  filter(!site %in% c("SS1", "M7")) #M7 doesn't have a matching wet season measurement and SS1 is an outlier. 

#add river column 
seasonal_variation <- seasonal_variation %>%
  mutate(river = case_when(
    grepl("B", site) ~ "Bassac",
    grepl("M", site) ~ "Mekong",
    grepl("SK", site) ~ "Sekong",
    grepl("SS", site) ~ "Sesan",
    TRUE ~ "Unknown"  # In case no match is found
  ))


### waste management stats
pop_and_waste <- read_csv("river_buff_2km_dist_pop_clean.csv")

WC_by_area <- read_csv("river_buff_2km_dist_wc_by_area.csv")

#Group by river corridor length and sum pop without access to Waste collection
Total_pop <- pop_and_waste %>%
  group_by(site,length) %>%
  summarize(total_pop = sum(pop, na.rm = TRUE))

Total_area <- WC_by_area %>%
  group_by(site,length) %>%
  summarize(sqr_km = sum(sqr_km, na.rm = TRUE))

pop_density <- left_join(Total_area, Total_pop, by = c("site", "length"))
pop_density <- pop_density %>%
  mutate(pop_density = total_pop/sqr_km)

#Group by river corridor length and sum pop without access to Waste collection
WC_by_area <- WC_by_area %>%
  group_by(site,length,WC_acces) %>%
  summarize(sqr_km = sum(sqr_km, na.rm = TRUE))

# create waste collection percent weighted by area. 
WC_by_area_weighted <- WC_by_area %>%
  group_by(site, length) %>%
  summarize(weighted_WC = sum(WC_acces * sqr_km) / sum(sqr_km), .groups = 'drop')

#make a column with the number of people who do not have access to waste collection
pop_and_waste <- pop_and_waste %>%
  mutate(pop_without_WC = pop * (1 - WC_access_percent))

#Group by river corridor length and sum pop without access to Waste collection
pop_no_WC<- pop_and_waste %>%
  group_by(site,length) %>%
  summarize(total_pop_without_WC = sum(pop_without_WC, na.rm = TRUE))

#Group by river corridor length and sum pop
pop_and_plastic<- pop_and_waste %>%
  group_by(site,length) %>%
  summarize(total_pop = sum(pop, na.rm = TRUE))

#join pop and weighted waste collection tables
pop_summary <- left_join(pop_no_WC, WC_by_area_weighted, by = c("site", "length"))
pop_all_summary <-left_join(pop_and_plastic, WC_by_area_weighted, by = c("site", "length"))


# Create separate tables for each unique value in 'length'
pop_summary_2 <- pop_summary %>% filter(length == 2)
pop_summary_5 <- pop_summary %>% filter(length == 5)
pop_summary_10 <- pop_summary %>% filter(length == 10)
pop_summary_20 <- pop_summary %>% filter(length == 20)

pop_all_summary_2 <- pop_all_summary %>% filter(length == 2)
pop_all_summary_5 <- pop_all_summary %>% filter(length == 5)
pop_all_summary_10 <- pop_all_summary %>% filter(length == 10)
pop_all_summary_20 <- pop_all_summary %>% filter(length == 20)

#create separate plastic counts for wet and dry season
Plastic_counts_all_wet <- Plastic_counts_all_m3 %>% filter(season == "wet")

#Set duplicate names to the same as the primary sample
Plastic_counts_all_wet <- Plastic_counts_all_wet %>%
  mutate(site = gsub("D", "", site))

Plastic_counts_all_dry <- Plastic_counts_all_m3 %>% filter(season == "dry")

# Create combined dry season tables
plastic_pop_2_dry <- left_join(Plastic_counts_all_dry, pop_summary_2, by = "site")
plastic_pop_5_dry <- left_join(Plastic_counts_all_dry, pop_summary_5, by = "site")
plastic_pop_10_dry <- left_join(Plastic_counts_all_dry, pop_summary_10, by = "site")
plastic_pop_20_dry <- left_join(Plastic_counts_all_dry, pop_summary_20, by = "site")

plastic_pop_all_2_dry <- left_join(Plastic_counts_all_dry, pop_all_summary_2, by = "site")
plastic_pop_all_5_dry <- left_join(Plastic_counts_all_dry, pop_all_summary_5, by = "site")
plastic_pop_all_10_dry <- left_join(Plastic_counts_all_dry, pop_all_summary_10, by = "site")
plastic_pop_all_20_dry <- left_join(Plastic_counts_all_dry, pop_all_summary_20, by = "site")


# Create combined wet season tables
plastic_pop_2_wet <- left_join(Plastic_counts_all_wet, pop_summary_2, by = "site")
plastic_pop_5_wet <- left_join(Plastic_counts_all_wet, pop_summary_5, by = "site")
plastic_pop_10_wet <- left_join(Plastic_counts_all_wet, pop_summary_10, by = "site")
plastic_pop_20_wet <- left_join(Plastic_counts_all_wet, pop_summary_20, by = "site")

plastic_pop_all_2_wet <- left_join(Plastic_counts_all_wet, pop_all_summary_2, by = "site")
plastic_pop_all_5_wet <- left_join(Plastic_counts_all_wet, pop_all_summary_5, by = "site")
plastic_pop_all_10_wet <- left_join(Plastic_counts_all_wet, pop_all_summary_10, by = "site")
plastic_pop_all_20_wet <- left_join(Plastic_counts_all_wet, pop_summary_20, by = "site")



#demonstration of poisson distribution for plastic data
ggplot(Plastic_counts_all_m3, aes(x=count_per_m3))+
  geom_histogram(binwidth = 1000)


# Fit a Poisson GLM
glm_model <- glm(count_per_m3 ~ weighted_WC, 
                 data = plastic_pop_5_dry, 
                 family = poisson)

summary(glm_model)


# Compute overdispersion ratio
overdispersion_ratio <- sum(residuals(glm_model, type = "pearson")^2) / df.residual(glm_model)
overdispersion_ratio 

#data is over dispersed so switching to a negative binomial model 

# Negative Binomial model 
library(MASS)
glm_nb <- glm.nb(count_per_m3 ~ weighted_WC, data = plastic_pop_20_dry)
summary(glm_nb)

par(mfrow = c(2, 2))
plot(glm_nb, which = 1)  # Residuals vs. Fitted
plot(glm_nb, which = 2)  # Q-Q Plot
plot(glm_nb, which = 3)  # Scale-Location
plot(glm_nb, which = 4)  # Cook’s Distance

par(mfrow = c(1,1))

#### models for weighted waste collection availability. Change data source to view various plots
plot(plastic_pop_20_wet$weighted_WC, plastic_pop_20_wet$count_per_m3)
plot(plastic_pop_20_dry$weighted_WC, plastic_pop_20_dry$count_per_m3)


ggplot(data = plastic_pop_20_dry, aes(x= weighted_WC, y = count_per_m3))+
  geom_jitter(width = 0.05, height = 0.05)+
  geom_smooth(method = 'glm.nb')+
  labs(title = "count_per_M3 ~ weighted_WC, data:plastic_pop_20_dry", x = "Weighted WC", y = "Count per L") 
  
# Run negative binomial regression models, change data source each model
nb_model <- glm.nb(count_per_m3 ~ weighted_WC, 
                   data = plastic_pop_20_wet)

# View model summary
summary(nb_model)


### WC and pop in 4km buffer
WC_and_pop_4km <- read_csv("River_buffer_4km_pop.csv")

WC_and_pop_4km_clean <- WC_and_pop_4km %>%
  group_by(site,length,WC_acces, SUM) %>%
  summarize(sqr_km = sum(sqr_km, na.rm = TRUE))

#add average pop density by square km and weighted average of waste collection 
WC_and_pop_summary_4km<- WC_and_pop_4km_clean %>%
  group_by(site, length) %>%
  summarize(
    avg_pop_per_km_sqr = mean(SUM / sqr_km, na.rm = TRUE),
    weighted_WC_access = sum(WC_acces * sqr_km, na.rm = TRUE) / sum(sqr_km, na.rm = TRUE)
  )

# Create separate tables for each unique value in 'length'
WC_and_pop_4km_2 <- WC_and_pop_summary_4km %>% filter(length == 2)
WC_and_pop_4km_5 <- WC_and_pop_summary_4km %>% filter(length == 5)
WC_and_pop_4km_10 <- WC_and_pop_summary_4km %>% filter(length == 10)
WC_and_pop_4km_20 <- WC_and_pop_summary_4km %>% filter(length == 20)
WC_and_pop_4km_40 <- WC_and_pop_summary_4km %>% filter(length == 40)

# Create combined dry season tables
plastic_WC_pop_2_dry <- left_join(Plastic_counts_all_dry, WC_and_pop_4km_2, by = "site")
plastic_WC_pop_5_dry <- left_join(Plastic_counts_all_dry, WC_and_pop_4km_5, by = "site")
plastic_WC_pop_10_dry <- left_join(Plastic_counts_all_dry, WC_and_pop_4km_10, by = "site")
plastic_WC_pop_20_dry <- left_join(Plastic_counts_all_dry, WC_and_pop_4km_20, by = "site")
plastic_WC_pop_40_dry <- left_join(Plastic_counts_all_dry, WC_and_pop_4km_40, by = "site")

# Create combined wet season tables
plastic_WC_pop_2_wet <- left_join(Plastic_counts_all_wet, WC_and_pop_4km_2, by = "site")
plastic_WC_pop_5_wet <- left_join(Plastic_counts_all_wet, WC_and_pop_4km_5, by = "site")
plastic_WC_pop_10_wet <- left_join(Plastic_counts_all_wet, WC_and_pop_4km_10, by = "site")
plastic_WC_pop_20_wet <- left_join(Plastic_counts_all_wet, WC_and_pop_4km_20, by = "site")
plastic_WC_pop_40_wet <- left_join(Plastic_counts_all_wet, WC_and_pop_4km_40, by = "site")


# Run negative binomial regression models, change data source each model
nb_model <- glm.nb(count_per_m3 ~ weighted_WC_access, 
                     data = plastic_WC_pop_40_dry)
              

# View model summary
summary(nb_model)

ggplot(data = plastic_WC_pop_10_dry, aes(x= weighted_WC_access, y = count_per_m3))+
  geom_jitter(width = 0.05, height = 0.05)+
  geom_smooth(method = 'glm.nb')+
  labs(title = "count_per_L ~ weighted_WC, data:plastic_pop_10_dry, 4km", x = "Weighted WC", y = "Count per L") 




#____seasonal variation
### waste management stats
pop_and_waste <- read_csv("river_buff_2km_dist_pop_clean.csv")

WC_by_area <- read_csv("river_buff_2km_dist_wc_by_area.csv")

#Create population density table
Total_pop <- pop_and_waste %>%
  group_by(site,length) %>%
  summarize(total_pop = sum(pop, na.rm = TRUE))

Total_area <- WC_by_area %>%
  group_by(site,length) %>%
  summarize(sqr_km = sum(sqr_km, na.rm = TRUE))

pop_density <- left_join(Total_area, Total_pop, by = c("site", "length"))
pop_density <- pop_density %>%
  mutate(pop_density = total_pop/sqr_km)

# Create separate tables for each unique value in 'length'
pop_density_2 <- pop_density %>% filter(length == 2)
pop_density_5 <- pop_density %>% filter(length == 5)
pop_density_10 <- pop_density %>% filter(length == 10)
pop_density_20 <- pop_density %>% filter(length == 20)

# Join pop density and plastic seasonality tables
seasonality_pop_2 <- left_join(seasonal_variation, pop_density_2, by = "site")
seasonality_pop_5 <- left_join(seasonal_variation, pop_density_5, by = "site")
seasonality_pop_10 <- left_join(seasonal_variation, pop_density_10, by = "site")
seasonality_pop_20 <- left_join(seasonal_variation, pop_density_20, by = "site")



ggplot(data = seasonality_pop_2, aes(x= pop_density, y = dry))+
  geom_jitter(width = 0.05, height = 0.05)+
  geom_smooth(method = 'glm')+
  ggtitle("Plastic counts vs population density, 2km corridor length")+
  theme_minimal()

ggplot(data = seasonality_pop_5, aes(x= pop_density, y = seasonal_variation, color = river))+
  geom_jitter(width = 0.05, height = 0.05)+
  geom_smooth(method = 'glm')+
  ggtitle("Seasonal variation vs population density")+
  theme_minimal()

# Run glm regression models, plastic ~ total pop, change data source to test each dataset
glm_model <- glm(seasonal_variation ~ pop_density, 
                   data = seasonality_pop_20
                 %>% filter(river == "Mekong"))

# View model summary
summary(glm_model)




#__seasonal variation test
Plastic_counts_335 <- data_all_cor %>%
  filter(filter != "100um") %>%
  group_by(season, site) %>%
  summarise(N_per_L = sum(N_per_L, na.rm = TRUE))


Plastic_counts_wide <- Plastic_counts_all_clean |> 
  pivot_wider(names_from = season, values_from = count_per_m3)

Plastic_counts_Mekong <- Plastic_counts_wide %>% filter(river == "Mekong")
Plastic_counts_tributaries <- Plastic_counts_wide %>% filter(river != "Mekong")



# Run Wilcoxon signed-rank test to test for variation by season
wilcox.test(Plastic_counts_wide$wet, Plastic_counts_wide$dry, paired = TRUE)

#compare variation between sites above and below urban areas
pairwise.wilcox.test(
  Plastic_counts_all_clean %>%
    filter(site %in% c("M5", "M8")) %>%
    pull(count_per_m3),
  
  Plastic_counts_all_clean %>%
    filter(site %in% c("M5", "M8")) %>%
    pull(site),
  
  p.adjust.method = "fdr"
)

head(Plastic_counts_all_clean)

shapiro.test(Plastic_counts_all_clean$count_per_m3) # confirms that data is non-parametric

#check for significant difference between seasons
model <-glm.nb(count_per_m3 ~ season, Plastic_counts_all_clean)
summary(model)
