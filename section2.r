#Get positive rdt rate
rdt <- rdt %>% 
  mutate(pos_rate = n_rdt_pos / n_child)

#Obtain summary statistics 
summary(rdt$pos_rate)

rdt_summary <- rdt %>% 
  group_by(adm1name) %>% 
  summarise(n_clusters = n_distinct(dhs_cluster),
            n_households = sum(n_hholds),
            n_children = sum(n_child),
            prop_urban = mean(urban),
            mean_pos_rate = mean(pos_rate)
  )

rdt_summary %>%
  mutate(adm1name = str_to_title(adm1name)) %>%
  gt(rowname_col = "adm1name") %>%
  tab_header(title = "Table 2 - Côte d'Ivoire DHS 2021: Regional Malaria RDT Summary") %>%
  fmt_number(
    columns = c(prop_urban, mean_pos_rate),
    decimals = 3
  ) %>%
  fmt_integer(
    columns = c(n_clusters, n_households, n_children)
  ) %>%
  cols_label(
    n_clusters    = "Clusters (n)",
    n_households  = "Households (n)",
    n_children    = "Children Tested (n)",
    prop_urban    = "Proportion of Urban Clusters",
    mean_pos_rate = "Mean RDT Positivity Rate"
  ) %>%
  grand_summary_rows(
    columns = c(n_clusters, n_households, n_children),
    fns = list(
      Median ~ round(median(.), 0)
    )
  ) %>%
  grand_summary_rows(
    columns = c(prop_urban, mean_pos_rate),
    fns = list(
      Median ~ round(median(.), 3)
    )
  )

#histogram of positive rate
rdt %>%
  ggplot(aes(pos_rate)) +
  geom_histogram(bins = 15, color = "white", fill = "red4") +
  theme_bw() +
  labs(
    title = "Figure 7: Distribution of Malaria RDT Positivity Rates",
    subtitle = "DHS Survey Clusters, Côte d'Ivoire 2021",
    x = "RDT Positivity Rate",
    y = "Number of Clusters"
  ) +
  scale_x_continuous(labels = scales::percent_format()) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    plot.caption = element_text(size = 8, color = "grey40"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

#Boxplot Rural vs Urban
rdt %>%
  mutate(urban = factor(urban, levels = c(0, 1), 
                        labels = c("Rural", "Urban"))) %>%
  ggplot(aes(x = urban, y = pos_rate, fill = urban)) +
  geom_boxplot(color = "grey30", alpha = 0.8, outlier.shape = 21,
               outlier.fill = "white", outlier.color = "grey30") +
  scale_fill_manual(values = c("Rural" = "red4", "Urban" = "steelblue4")) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw() +
  labs(
    title = "Figure 8: Malaria RDT Positivity Rate by Settlement Type",
    subtitle = "DHS Survey Clusters, Côte d'Ivoire 2021",
    x = NULL,
    y = "RDT Positivity Rate"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    plot.caption = element_text(size = 8, color = "grey40"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

#Make Choropleth map 
#Get region borders 
adm_region <- civ_dhs %>%
  group_by(DHS_REGION) %>%
  summarise(geometry = st_union(geometry))

#change dataframe to shapefile 
rdt_sf <- rdt %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  mutate(pos_rate = n_rdt_pos / n_child)

#create categorical variable
rdt_sf <- rdt_sf %>%
  mutate(urban = factor(urban, levels = c(0, 1),
                        labels = c("Rural", "Urban")),
         pos_rate_cat = cut(pos_rate,
                            breaks = c(0, 0.1, 0.50, 1),
                            labels = c("Low (0-10%)", 
                                       "Moderate (10-50%)", "High (>50%)"),
                            include.lowest = TRUE
         ))

#join dataframes 
rdt_summary <- rdt_summary %>% 
  mutate(adm1name = str_to_lower(adm1name))
adm_region_pos <- adm_region %>%
  left_join(rdt_summary, by = c("DHS_REGION" = "adm1name"))

#make map 
tm_shape(adm_region_pos) +
  tm_polygons(
    fill = "mean_pos_rate",
    fill.scale = tm_scale_continuous(values = "YlOrRd"),
    fill.legend = tm_legend(title = "Regional Mean\nRDT Positivity Rate"),
    col = "grey40",
    lwd = 1.2
  ) +
  tm_shape(rdt_sf) +
  tm_dots(
    fill = "pos_rate_cat",
    fill.scale = tm_scale_categorical(
      values = c("Low (0-10%)" = "dodgerblue",
                 "Moderate (10-50%)" = "yellow",
                 "High (>50%)" = "white")
    ),
    fill.legend = tm_legend(title = "Cluster RDT\nPositivity Rate"),
    size = 0.15,
    fill_alpha = 0.9
  ) +
  tm_title("Figure 9: Malaria RDT Positivity Rate by Region and Cluster") 

# Semivariogram prework steps

# Transform the prevalence to the empirical logit scale
rdt_sf$logitp <- log((rdt_sf$n_rdt_pos + 0.5)/(rdt_sf$n_child - rdt_sf$n_rdt_pos + 0.5))

hist(rdt_sf$logitp, xlab = " ", main = "Figure 10: Empirical logit",
     border = "white", col = "seagreen3")

ggplot(rdt_sf, aes(x = logitp)) +
  geom_histogram(bins = 8, fill = "seagreen3", color = "white") +
  theme_bw() +
  labs(
    title = "Figure 10: Empirical Logit of RDT Positivity Rate",
    subtitle = "DHS Survey Clusters, Côte d'Ivoire 2021",
    x = "Empirical Logit",
    y = "Number of Clusters"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    plot.caption = element_text(size = 8, color = "grey40"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )


#Empirical semivariogram
empirical <- ggvario(coords = st_coordinates(rdt_sf), data = rdt_sf$logitp) +
  labs(title = "Figure 11: Empirical Semivariogram of Logit RDT Positivity Rate")

print(update(plot(empirical), 
             main = "Figure 12: Experimental Variogram and Fitted Model",
             xlab = "Distance",
             ylab = "Semi-variance"))


#fit the theoretical semivariogram
rdt_spatial <- as(rdt_sf, "Spatial")  # sf -> SpatialPointsDataFrame

variofit <- automap::autofitVariogram(
  formula = logitp ~ 1,
  input_data = rdt_spatial,
  model = c("Sph", "Exp", "Gau")
)
plot(variofit, 
     main = "Figure 12: Experimental Variogram and Fitted Model\nLogit RDT Positivity Rate, Côte d'Ivoire 2021")

print(update(plot(variofit), 
             main = "Figure 12: Experimental Variogram and Fitted Model",
             xlab = "Distance",
             ylab = "Semi-variance"))


#get range
practicalRange(cov.model = "matern", phi = 19)

#Kriging
# Load the shapefile and plot it
map <- st_crop(civ_dhs, rdt_sf)

# reproject to 3857
rdt_sf_proj <- st_transform(rdt_sf, crs = 3857)
map_proj <- st_transform(map, crs = 3857)

# regenerate pred_locations in 3857
pred_locations <- st_make_grid(map_proj, cellsize = 5000, what = "centers")

# verify
plot(map_proj$geometry, col = "white")
plot(pred_locations, cex = 0.01, col = "red", pch = 19, add = TRUE)

# run kriging
orkriging <- autoKrige(
  formula = logitp ~ 1,
  input_data = as(rdt_sf_proj, "Spatial"),
  new_data = as(pred_locations, "Spatial"),
  model = c("Sph", "Exp", "Gau"),
  verbose = FALSE
)

plot(orkriging)

# First extract predictions on the logit scale 
pred <- (orkriging$krige_output$var1.pred)

# Convert them to the prevalence scale
pred_prevalence <- plogis(pred)

# compare predictions to data
summary(pred_prevalence)
summary(rdt_sf$pos_rate)

pred_stdev <- orkriging$krige_output$var1.stdev
summary(pred_stdev)

# Use standard deviation to calculate 95% CI for prevalence 
lower <- plogis(pred - 1.96 * pred_stdev)
summary(lower)
upper <- plogis(pred + 1.96 * pred_stdev)
summary(upper)

# approximate standard deviation of prevalnece 
sd_prev <- (upper - lower) / (2 * 1.96)

# confidence interval width
ci_width <- upper - lower   # CI width is a useful measure of error as it is on the same scale as the prediction (prevalence) and makes intuitive sense

#  Create a dataframe with coordinates, prevalence, st dev values, LCI, UCI, and CI width
pred_df <- data.frame(
  st_coordinates(pred_locations),
  prevalence = pred_prevalence,
  sd_prev = sd_prev,
  lower = lower,
  upper = upper,
  ci_width = ci_width
)
summary(pred_df)
# Convert dataframe to raster
pred_rast <- terra::rast(pred_df, type = "xyz", crs = st_crs(rdt_sf_proj)$wkt)
pred_rast <- crop(pred_rast, map_proj, mask = TRUE)


# force correct source CRS then reproject
crs(pred_rast) <- "EPSG:3857"
pred_rast_4326 <- terra::project(pred_rast, "EPSG:4326")
summary(pred_rast_4326)

# verify
ext(pred_rast_4326)  # should show degrees now

# now crop
pred_rast_4326 <- crop(pred_rast_4326, map, mask = TRUE)
plot(pred_rast_4326[["prevalence"]], 
     col = viridis::viridis(100),
     main = "Interpolated Malaria Prevalence Surface")
plot(st_geometry(map), add = TRUE, border = "white", lwd = 0.5)

plot(pred_rast_4326[["ci_width"]], 
     col = rev(RColorBrewer::brewer.pal(9, "RdYlGn")),
     main = "Prediction Uncertainty (95% CI Width)",
     legend = TRUE)
plot(st_geometry(civ_dhs), add = TRUE, border = "white", lwd = 0.5)

summary(pred_rast_4326$prevalence)

#Part i: Define endemic threshold 
high_risk <- pred_rast_4326[["prevalence"]] > 0.539
plot(high_risk, main = "High Risk Areas (Prevalence > 54%)")
plot(st_geometry(civ_dhs), add = TRUE, border = "white", lwd = 0.5)

#Part ii: 
pop <- terra::rast("civ_pop_2021_CN_1km_R2025A_UA_v1.tif")  
pop <- terra::project(pop, crs(pred_rast_4326))

plot(pop)
summary(pop)

# resample population to match pred_rast_4326 resolution
pop_resampled <- terra::resample(pop, pred_rast_4326[["prevalence"]], 
                                 method = "bilinear")

# recreate high risk layer
high_risk <- pred_rast_4326[["prevalence"]] > 0.54

# mask population to high risk areas only
pop_high_risk <- terra::mask(pop_resampled, high_risk, maskvalue = 0)

# national totals
total_pop <- global(pop_resampled, sum, na.rm = TRUE)
total_pop_risk <- global(pop_high_risk, sum, na.rm = TRUE)

cat("Total population:", round(total_pop$sum), "\n")
cat("Population at high risk:", round(total_pop_risk$sum), "\n")
cat("Proportion at high risk:", 
    round(total_pop_risk$sum / total_pop$sum, 3), "\n")

# 
# create summary dataframe
national_summary <- data.frame(
  level = "Côte d'Ivoire (National)",
  total_pop = round(total_pop$sum),
  pop_at_risk = round(total_pop_risk$sum),
  prop_at_risk = total_pop_risk$sum / total_pop$sum
)

adm1_vect <- terra::vect(adm_region)

# extract total population by region
pop_total_adm1 <- terra::extract(pop_resampled, adm1_vect,
                                 fun = sum, na.rm = TRUE)

# extract high risk population by region
pop_risk_adm1 <- terra::extract(pop_high_risk, adm1_vect,
                                fun = sum, na.rm = TRUE)

# combine into summary dataframe
adm1_summary <- data.frame(
  level = adm_region$DHS_REGION,
  total_pop = round(pop_total_adm1[, 2]),
  pop_at_risk = round(pop_risk_adm1[, 2])
) %>%
  mutate(
    level = str_to_title(level),
    pop_at_risk = replace_na(pop_at_risk, 0),
    prop_at_risk = pop_at_risk / total_pop
  )

# Final GT Table
adm1_summary %>%
  arrange(desc(prop_at_risk)) %>%
  bind_rows(national_summary) %>%
  gt() %>%
  tab_header(
    title = "Table 3: Population at High Malaria Risk (Prevalence > 54%)",
    subtitle = "Côte d'Ivoire, 2021"
  ) %>%
  fmt_integer(columns = c(total_pop, pop_at_risk)) %>%
  fmt_percent(columns = prop_at_risk, decimals = 1) %>%
  cols_label(
    level = "Region",
    total_pop = "Total Population",
    pop_at_risk = "Population at High Risk (n)",
    prop_at_risk = "Population at High Risk (%)"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = level == "Côte d'Ivoire (National)"
    )
  )

adm1_map <- adm_region %>%
  mutate(DHS_REGION = str_to_title(DHS_REGION)) %>%
  left_join(adm1_summary, by = c("DHS_REGION" = "level"))

p1 <- tm_shape(pred_rast_4326[["prevalence"]]) +
  tm_raster(
    col.scale = tm_scale_continuous(values = "viridis"),
    col.legend = tm_legend(title = "Predicted\nPrevalence")
  ) +
  tm_shape(adm1_map) +
  tm_borders(col = "white", lwd = 0.8) +
  tm_title("Figure 13a: Malaria Prevalence Surface")

p2 <- tm_shape(pred_rast_4326[["ci_width"]]) +
  tm_raster(
    col.scale = tm_scale_continuous(values = "brewer.yl_or_rd"),
    col.legend = tm_legend(title = "95% CI Width")
  ) +
  tm_shape(adm1_map) +
  tm_borders(col = "white", lwd = 0.8) +
  tm_title("Figure 13b: Prediction Uncertainty")

tmap_arrange(p1, p2, ncol = 2)


# create high risk layer
high_risk <- pred_rast_4326[["prevalence"]] > 0.54

# convert population raster to points
pop_points <- as.data.frame(pop_resampled, xy = TRUE) %>%
  rename(population = 3) %>%
  filter(!is.na(population), population > 0) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)

# create log population raster
pop_log <- log1p(pop_resampled)  # log1p handles zeros safely

tm_shape(pop_log) +
  tm_raster(
    col.scale = tm_scale_continuous(values = "greys"),
    col.legend = tm_legend(title = "Population\nper cell (log)")
  ) +
  tm_shape(high_risk) +
  tm_raster(
    col.scale = tm_scale_continuous(values = c("white", "red")),
    col.legend = tm_legend(show = FALSE),
    col_alpha = 0.4
  ) +
  tm_shape(adm1_map) +
  tm_borders(col = "grey20", lwd = 1) +
  tm_shape(st_union(adm1_map)) +
  tm_borders(col = "black", lwd = 1.5) +
  tm_title("Figure 14: Population at High Malaria Risk (Prevalence > 54%)") +
  tm_add_legend(
    type = "polygons",
    fill = "red",
    labels = "High Risk (>54%)"
  )

#prediction uncertainty with clusters
tm_shape(pred_rast_4326[["ci_width"]]) +
  tm_raster(
    col.scale = tm_scale_continuous(values = "YlOrRd"),
    col.legend = tm_legend(title = "95% CI Width")
  ) +
  tm_shape(rdt_sf) +
  tm_dots(
    fill = "black",
    size = 0.1,
    fill_alpha = 0.6
  ) +
  tm_shape(adm1_map) +
  tm_borders(col = "white", lwd = 0.8) +
  tm_title("Figure 15: Prediction Uncertainty and Existing Cluster Locations") +
  tm_add_legend(
    type = "polygons",
    fill = "black",
    labels = "Existing DHS Cluster"
  )

#high risk with clusters
tm_shape(high_risk) +
  tm_raster(
    col.scale = tm_scale_continuous(values = c("lightyellow", "red4")),
    col.legend = tm_legend(show = FALSE)
  ) +
  tm_shape(rdt_sf) +
  tm_dots(
    fill = "black",
    size = 0.1,
    fill_alpha = 0.6
  ) +
  tm_shape(adm1_map) +
  tm_borders(col = "grey40", lwd = 0.8) +
  tm_title("Figure 16: High Risk Areas and Existing Cluster Locations") +
  tm_add_legend(
    type = "polygons",
    fill = "red4",
    labels = "High Risk (>54%)"
  ) +
  tm_add_legend(
    type = "polygons",
    fill = "black",
    labels = "Existing DHS Cluster"
  )








