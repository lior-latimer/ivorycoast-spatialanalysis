#load packages
# Load in safe order to avoid S4 class conflicts
library(sp)       # required for automap
library(gstat)    # dependency for automap
library(automap)  # variogram + kriging
library(sf)
library(terra)
library(mapview)
library(Hmisc)
library(geoR)
library(ggplot2)
library(dplyr)
library(RiskMap)
library(shiny)
library(tmap)
library(gt)
library(tidyverse)
library(spdep)


#read files
opendef <- read.csv('adm2_opendef.csv')
rdt <- read.csv('cluster_rdt.csv')
civ_dhs <- st_read('civ_DHS.shp')

#Descriptive statistics
#Gather urban/rural data
urban <- rdt %>% 
  group_by(adm1name) %>% 
  summarise(urban_prop = mean(urban)) %>% 
  mutate(adm1name = str_to_lower(adm1name))

#Aggregate open defecation data by region
region_data <- civ_dhs %>% 
  group_by(DHS_REGION) %>% 
  summarise(n_clusters = sum(N_CLUSTERS),
            n_households = sum(N_HOUSEHOL),
            n_opendef = sum(N_OPENDEF),
            avg_prop_od = mean(prop_od))

#join urban data to region data 
region_data <- full_join(region_data, urban, by = c('DHS_REGION' = "adm1name")) %>% 
  st_drop_geometry() %>% 
  select(DHS_REGION, n_clusters, n_households, n_opendef, avg_prop_od, urban_prop) %>% 
  rename(urban_prop = urban_prop)

#Make table one 
region_data %>% 
  mutate(DHS_REGION = str_to_title(DHS_REGION)) %>%   # capitalise region names
  gt(rowname_col = "DHS_REGION") %>% 
  tab_header(title = "Table 1 - Côte d'Ivoire DHS 2021: Regional Summary") %>% 
  fmt_number(
    columns = c(n_clusters, avg_prop_od, urban_prop),
    decimals = 3
  ) %>% 
  fmt_integer(
    columns = c(n_clusters, n_households, n_opendef)
  ) %>% 
  cols_label(
    n_clusters = "Clusters (n)",
    n_households  = "Households (n)",
    n_opendef     = "Open Defecation (n)",
    avg_prop_od   = "Avg Open Defecation Proportion",
    urban_prop    = "Proportion of Urban Clusters"
  ) %>% 
  grand_summary_rows(
    columns = c(n_clusters, n_households, n_opendef),
    fns = list(
      Median ~ round(median(.), 0)
    )
  ) %>%
  grand_summary_rows(
    columns = c(avg_prop_od, urban_prop),
    fns = list(
      Median ~ round(median(.), 3)
    ))

#Histogram of open defacation rate
civ_dhs %>%
  ggplot(aes(prop_od)) +
  geom_histogram(bins = 10, color = "white", fill = "seagreen") +
  labs(
    title = "Figure 1: Distribution of Open Defecation Rates",
    subtitle = "Côte d'Ivoire DHS clusters",
    x = "Proportion Open Defecation",
    y = "Number of districts"
  ) +
  scale_x_continuous(labels = scales::percent_format()) +  # show as % instead of decimals
  theme_minimal() +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 11),
    plot.caption  = element_text(color = "grey60", size = 9),
    axis.title    = element_text(size = 11),
    panel.grid.minor = element_blank()  # remove minor gridlines
  )

#Choropleth Map 
tm_shape(civ_dhs) +
  tm_polygons(
    fill = "prop_od",
    fill.scale = tm_scale_intervals(
      breaks = c(0, 0.10, 0.25, 0.50, 1),
      labels = c("0–.1", ".1-.25", ".25–.5", ".5–1"),
      values = "brewer.yl_or_rd"),
    fill.legend = tm_legend(title = "Rate of Open Defecation")) +
  tm_borders() +
  tm_compass(position = c("right", "top")) +
  tm_scalebar() + 
  tm_title("Figure 2: Rate of Open Defecation\nCôte d'Ivoire", size = 0.9)

#Define which regions are neighbours to one another based on queen contiguity
neighbours_dist <- poly2nb(civ_dhs)  

#Visualization of neighbors 
district_centroids <- civ_dhs %>%
  st_centroid() 
plot(civ_dhs$geometry, lwd = 0.2)
plot(neighbours_dist, district_centroids$geometry, 
     col = 'red', cex = 0.1, lwd = 1, add = TRUE)

#Update the neighbour pattern to a row standardised weights matrix.
weights_queen <- nb2listw(neighbours_dist, style = "W", zero.policy = T) 

# complete a Global Moran's I test
district_moran <- moran.test(civ_dhs$prop_od, listw = weights_queen, zero.policy = T)
moran_slope <- round(lm(lag_scaled_od ~ scaled_od, data = civ_dhs)$coefficients[2], 3)


# scale the OD parameter
civ_dhs$scaled_od <- scale(civ_dhs$prop_od)
# create a lag vector from the neighbour list and the scaled OD values
civ_dhs$lag_scaled_od <- lag.listw(weights_queen, 
                                   civ_dhs$scaled_od, zero.policy = T)

# plot the output
ggplot(data = civ_dhs, aes(x = scaled_od, y = lag_scaled_od)) +
  geom_smooth(method = "lm", se = FALSE, colour = "#e05c5c", linetype = "dashed", linewidth = 0.8) +
  geom_point(size = 2, alpha = 0.7, colour = "grey30") +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey20") +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey20") +
  annotate("text", 
           x = 2, y = -0.6,  # adjust position as needed
           label = paste0("Moran's I (slope) = ", moran_slope),
           colour = "red4", size = 3.5) +
  labs(
    title = "Figure 3 - Moran's I Scatterplot: Open Defecation",
    subtitle = "Spatial lag vs. standardised OD rate",
    x = "Raw Open Defecation Rate (scaled)",
    y = "Spatial Lag (scaled OD)"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(colour = "grey40", size = 10),
    axis.title       = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90")
  )

# Calculate the local Moran statistic for each region using queen contiguity
local_moran <- localmoran(civ_dhs$prop_od, 
                          weights_queen, # our weights object
                          zero.policy = TRUE, 
                          alternative = "two.sided")

# Replace the column names to make it clearer
colnames(local_moran) <- c("local_I", "expected_I", 
                           "variance_I", "z_statistic", "p_value")


# Join Local Moran's I information back on to the spatial data for plotting
local_moran <- cbind(civ_dhs, local_moran) %>%
  select(scaled_od, lag_scaled_od, p_value, geometry)

# Reset the bbox as it can become incompatible after data manipulation.
local_moran <- st_as_sf(local_moran)
local_moran <- st_transform(local_moran, 4326)

# plot the data to look for hot and cold spot clusters based on p_values
tm_shape(st_transform(local_moran, crs = 4326)) +
  tm_polygons(
    fill = "p_value",
    fill.scale = tm_scale_intervals(
      style = "fixed",
      breaks = c(1e-20, 0.01, 0.05, 1),
      labels = c("p < 0.01", "p < 0.05", "p ≥ 0.05"),
      values = "-brewer.yl_or_rd"
    ),
    fill.legend = tm_legend(
      title = "Significance",
      title.fontface = "bold",
      frame = FALSE
    ),
    lwd = 0.2,
    col = "grey30"
  ) +
  tm_compass(position = c(0.02, 0.18)) +
  tm_scalebar() +
  tm_title("Figure 4 - Local Moran's I: Spatial Clusters\nCôte d'Ivoire", size = 0.9)

# Local Moran Map at .05 Significance
# Join our local moran outputs to the main data and select columns for analysis
local_map <- civ_dhs %>%
  cbind(., local_moran) %>% # 
  select(scaled_od, lag_scaled_od, p_value, geometry)

# Run the map_maker function  
local_map <- map_maker(local_map)

# Plot the number of each cluster type
table(local_map$cluster)

tm_shape(st_transform(local_map, crs = 4326)) +
  tm_polygons(
    fill = "cluster",
    col = "grey30",
    lwd = 0.2,
    fill.scale = tm_scale_categorical(
      values = c(
        "high-high" = "#d7191c",
        "low-low"   = "#2c7bb6",
        "low-high"  = "#fdae61",
        "high-low"  = "#abd9e9",
        "p > 0.05"  = "#d3d3d3"
      )
    ),
    fill.legend = tm_legend(
      title = "Cluster Type",
      title.fontface = "bold",
      frame = FALSE
    )
  ) +
  tm_compass(position = c(0.85, 0.98), size = 1.5) +
  tm_scalebar() +
  tm_title("Figure 5 - Local Moran's I: Cluster Types (p < 0.05)\nCôte d'Ivoire", size = 0.9) 

# Local Moran Map at .01 Significance 
# Run the map_maker function  
local_map2 <- map_maker2(local_map)

# Plot the number of each cluster type
table(local_map2$cluster)

# Plot the data to look for hot and cold spot clusters based on p_values
tm_shape(st_transform(local_map2, crs = 4326)) +
  tm_polygons(
    fill = "cluster",
    col = "grey30",
    lwd = 0.2,
    fill.scale = tm_scale_categorical(
      values = c(
        "high-high" = "#d7191c",
        "low-low"   = "#2c7bb6",
        "low-high"  = "#fdae61",
        "high-low"  = "#abd9e9",
        "p > 0.01"  = "#d3d3d3"
      )
    ),
    fill.legend = tm_legend(
      title = "Cluster Type",
      title.fontface = "bold",
      frame = FALSE
    )
  ) +
  tm_compass(position = c(0.85, 0.98), size = 1.5) +
  tm_scalebar() +
  tm_title("Figure 6 - Local Moran's I: Cluster Types (p < 0.01)\nCôte d'Ivoire", size = 0.9) 






