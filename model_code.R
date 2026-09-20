
# Importing packages ------------------------------------------------------

pacman::p_load(
  fastverse,
  data.table,
  dplyr,
  ggplot2,
  stringr,
  scales,
  patchwork,
  sf,
  RColorBrewer,
  lubridate,
  INLA,
  inlabru,
  spdep,
  purrr,
  recipes
)
source("R/functions.R")

# Importing data ----------------------------------------------------------

adm1 <- st_read("shapefiles/hdx_adm1.shp", quiet = T)
adm2 <- st_read("shapefiles/hdx_adm2.shp", quiet = T)
adm3 <- st_read("shapefiles/hdx_adm3.shp", quiet = T)
covs <- fread("data/covs.csv")
all_df <- fread("data/all_df.csv")

# Small arrangement of data ----------------------------------------------

#' Source of endemic counties is from WHO:
#' https://www.who.int/publications/i/item/who-wer10004-21-32 and from the
#' Kenya VL strategic plan see: 
#' http://guidelines.health.go.ke:8000/media/KSPC-OF-LEISHMANIASIS-STRATEGY-2021-2025.pdf
#' on page 5.
#' 
dm2_recipe <- all_df |>
  mutate(
    IsEndemic = fifelse(
      county %chin% 
        c(
          "Baringo",
          "Garissa",
          "Isiolo",
          "Kajiado",
          "Kitui",
          "Marsabit",
          "Mandera",
          "Tharaka-Nithi",
          "Turkana",
          "Wajir",
          "West Pokot"
        ),
      "Yes",
      "No"
    ) |> 
      as.factor(),
    IsEndemic = relevel(IsEndemic, ref = "No")
  ) |>
  recipe(county + subcounty + ward + N + pop + prev ~ .) |>
  step_center(all_double_predictors()) |> 
  step_scale(all_double_predictors()) |> 
  step_dummy(all_nominal_predictors()) |> 
  prep()

dm <- juice(dm2_recipe) |> 
  mutate(b0 = 1, wardID = 1:n(), present = fifelse(N == 0, 0, 1))

# Creating neighborhood structure ----------------------------------------

h_temp <- poly2nb(adm3)

# Plotting the two areas with no neighbours
mi_plot <- adm3 |> 
  filter(subcounty == 'Mbita') |> 
  ggplot() +
  geom_sf(aes(fill = ward)) +
  geom_sf_text(aes(label = ward)) + 
  theme_void()

faza_plot <- adm3 |> 
  filter(subcounty == 'Lamu East') |> 
  ggplot() +
  geom_sf(aes(fill = ward)) +
  geom_sf_text(aes(label = ward)) + 
  theme_void()

# Temporarily for assigining
temp_adm3 <- adm3 |> 
  mutate(id = 1:n())
temp_adm3[temp_adm3$ward %in% c("Rusinga Island"), ]
temp_adm3[temp_adm3$ward %in% c("Basuba"), ]

# Manually assigning neighbours to those missing
h_temp[[64]] <- as.integer(c(65))
h_temp[[65]] <- as.integer(c(h_temp[[65]], 64)) 

h_temp[[146]] <- as.integer(c(148))
h_temp[[148]] <- as.integer(c(h_temp[[148]], 146)) 

# Finding the disjoint subgraphs ------------------------------------------

#' Connecting the 3 disjoint connected graphs. To see why we are doing this
#' see `kenya_disjoint_subgraphs.png` in the plots folder

# How many disjoint connected graphs do we have after solving for lack of neighbours?
comp <- n.comp.nb(h_temp)
comp$nc # 3 dijsoint connected graphs

# Getting coordinates and mainland Wards
coords <- st_coordinates(st_centroid(adm3))
main_comp_wards <- which(comp$comp.id == 1)

# We start with component 2 disjoint connected graph with wards:
temp_adm3 <- adm3 |> 
  mutate(component = comp$comp.id, id = 1:n()) 

comp2_wards <- temp_adm3[temp_adm3$component == 2, ]

# We pick any ward from the component 
island_comp2_ward <- comp2_wards[1,]$id  
dists_comp2 <- sqrt((coords[island_comp2_ward, 1] - coords[main_comp_wards, 1])^2 + 
                      (coords[island_comp2_ward, 2] - coords[main_comp_wards, 2])^2)
best_mainland_comp2 <- main_comp_wards[which.min(dists_comp2)]

h_temp[[island_comp2_ward]] <- as.integer(sort(unique(c(h_temp[[island_comp2_ward]], best_mainland_comp2))))
h_temp[[best_mainland_comp2]] <- as.integer(sort(unique(c(h_temp[[best_mainland_comp2]], island_comp2_ward))))

# Component 3 disjoint connected graph
comp3_wards <- temp_adm3[temp_adm3$component == 3, ]

# Pick any ward from component 3
island_comp3_ward <- comp3_wards[1, ]$id  

# Find nearest mainland ward
dists_comp3 <- sqrt((coords[island_comp3_ward, 1] - coords[main_comp_wards, 1])^2 + 
                      (coords[island_comp3_ward, 2] - coords[main_comp_wards, 2])^2)
best_mainland_comp3 <- main_comp_wards[which.min(dists_comp3)]

# Create bidirectional connection
h_temp[[island_comp3_ward]] <- as.integer(sort(unique(c(h_temp[[island_comp3_ward]], best_mainland_comp3))))
h_temp[[best_mainland_comp3]] <- as.integer(sort(unique(c(h_temp[[best_mainland_comp3]], island_comp3_ward))))


# Final check
comp_final <- n.comp.nb(h_temp)
comp_final$nc

# Converting to an INLA object
nb2INLA("H.adj", h_temp)
H <- inla.read.graph(filename = "H.adj")

# Plotting the neighborhood -----------------------------------------------

centroids <- st_centroid(st_geometry(adm3))
centroid_coords <- st_coordinates(centroids)

# Build edges 
edges <- do.call(rbind, lapply(1:length(h_temp), function(i) {
  if (length(h_temp[[i]]) == 0) return(NULL)
  data.frame(
    x = centroid_coords[i, 1],
    y = centroid_coords[i, 2],
    xend = centroid_coords[h_temp[[i]], 1],
    yend = centroid_coords[h_temp[[i]], 2],
    row.names = NULL  
  )
}))

first_order_nb <- ggplot() +
  geom_sf(data = adm3, fill = "white", color = "gray50") +
  geom_sf(data = adm3, aes(fill = ward, color = ward), show.legend = F) +
  geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend),
               color = "black", linewidth = 0.3) +
  geom_point(data = as.data.frame(centroid_coords),
             aes(X, Y), shape = 21, size = 1) +
  theme_void()
ggsave(
  plot = first_order_nb,
  filename = "plots/first_order_nb.png",
  dpi = 5e2,
  width = 2 * 5,
  height = 2 * 5,
  units = "in",
  bg = 'white'
)
# Testing for local spatial correlation -----------------------------------

# The non spatial formula
formula_non_spatial <- N ~ 0 + b0  + prec + temp + soil_moisture  + 
  clay_content + time_to_treatment + prop_itn_access + prop_livestock_ownership  + 
  prop_lowest_wealth_quantile + number_per_room + prop_wall_material 

glm_mod <- glm(
  formula = formula_non_spatial,
  data = dm,
  family = "poisson",
  offset = log(dm$pop)
)
summary(glm_mod)
# Extracting the residuals
resid_zip <- residuals(glm_mod,  type = "pearson")

# List of neighbours
lw <- nb2listw(h_temp, style = "W", zero.policy = TRUE)

# Global Moran's I
moran.test(resid_zip, lw, zero.policy = TRUE)
global_moran <- data.table(
  `Moran I statistic` = moran.test(resid_zip, lw, zero.policy = TRUE)$estimate[1],
  expectation = moran.test(resid_zip, lw, zero.policy = TRUE)$estimate[2],
  variance = moran.test(resid_zip, lw, zero.policy = TRUE)$estimate[3]
)

local_Is = localmoran(resid_zip, lw, zero.policy = TRUE)

local_moran <- data.table(
  county = dm$county,
  subcounty = dm$subcounty,
  ward = dm$ward,
  local_I = local_Is[, 1],
  p_value = local_Is[, 5],
  z_score = local_Is[, 4],
  cat_local_I = attr(local_Is,"quadr")[, 1]
) |>
  merge(data.table(adm3)) |>
  st_as_sf()

# Plotting for LISA
palette_clusters <- c(
  "High-High" = "#d7191c",  # Red
  "High-Low"  = "#fdae61",  # Orange
  "Low-High"  = "#abd9e9",  # Light Blue
  "Low-Low"   = "#2c7bb6",  # Blue
  "Not Significant" = "#ffffff"  # White
)

plot_local_I <- local_moran |>
  ggplot() +
  geom_sf(aes(fill = cat_local_I), col = "grey80", lwd = 0.001) +
  scale_fill_manual(values = palette_clusters) +
  theme_void() +
  theme(
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black")
  ) +
  labs(fill = "")
ggsave(
  plot = plot_local_I,
  filename = "plots/ward_lisa.png",
  dpi = 5e2,
  width = 1 * 5,
  height = 1 * 5,
  units = "in",
  bg = NULL
)

# Zero inflated models ----------------------------------------------------

# The formula and the prior
# The prior
prior <- list(prec = list(prior = "pc.prec", param = c(0.5 / 0.31, 0.01)),
              phi = list(prior = "pc", param = c(0.5, 2 / 3)))

formula <- N ~ 0 + b0  + prec + temp + soil_moisture +  
  clay_content + time_to_treatment + prop_itn_access + prop_livestock_ownership  + 
  prop_lowest_wealth_quantile + number_per_room + prop_wall_material +
  f(wardID, model = "bym2", graph = H, hyper = prior, scale.model = T) 

# The hurdle model
zap <- inla(
  formula = formula,
  family = "zeroinflatedpoisson0",
  E = dm$pop,
  data = dm,
  control.fixed = list(mean = 0, prec = 1),
  control.compute = list(config = TRUE, return.marginals.predictor = TRUE, cpo = T, waic = T, dic = T),
  num.threads	= '7:1',
  verbose = F
)
zap <- inla.rerun(zap)

# Zero Inflated Poisson Model
zip <- inla(
  formula = formula,
  family = "zeroinflatedpoisson1",
  E = dm$pop,
  data = dm,
  control.fixed = list(mean = 0, prec = 1),
  control.compute = list(config = TRUE, return.marginals.predictor = TRUE, cpo = T, waic = T, dic = T),
  num.threads	= '7:1',
  verbose = F
)
zip <- inla.rerun(zip)

# The joint model ---------------------------------------------------------

comps <- ~
  intercept_count(1) +
  precipitation_count(prec, model = "linear") +
  temperature_count(temp, model = "linear" ) +
  soil_moisture_count(soil_moisture, model = "linear" ) +
  clay_content_count(clay_content, model = "linear" ) +
  time_to_treatment_count(time_to_treatment, model = "linear" ) +
  prop_itn_access_count(prop_itn_access, model = "linear" ) +
  prop_livestock_ownership_count(prop_livestock_ownership, model = "linear" ) +
  prop_lowest_wealth_quantile_count(prop_lowest_wealth_quantile, model = "linear" ) +
  number_per_room_count(number_per_room, model = "linear" ) +
  prop_wall_material_count(prop_wall_material, model = "linear" ) +
  bym2_count(wardID, model = "bym2", graph = H, scale.model = T, hyper = prior) +
  intercept_occurence(1) +
  precipitation_occurence(prec, model = "linear") +
  temperature_occurence(temp, model = "linear" ) +
  soil_moisture_occurence(soil_moisture, model = "linear" ) +
  clay_content_occurence(clay_content, model = "linear" ) +
  time_to_treatment_occurence(time_to_treatment, model = "linear" ) +
  prop_itn_access_occurence(prop_itn_access, model = "linear" ) +
  prop_livestock_ownership_occurence(prop_livestock_ownership, model = "linear" ) +
  prop_lowest_wealth_quantile_occurence(prop_lowest_wealth_quantile, model = "linear" ) +
  number_per_room_occurence(number_per_room, model = "linear" ) +
  prop_wall_material_occurence(prop_wall_material, model = "linear") +
  bym2_occurence(wardID, model = "bym2", graph = H, scale.model = T, hyper = prior) 

truncated_poisson_obs <- bru_obs(
  family = "nzpoisson",
  data = dm[dm$N > 0, ],
  formula = N ~ 0 +
    intercept_count +
    precipitation_count +
    temperature_count +
    soil_moisture_count +
    clay_content_count +
    time_to_treatment_count +
    prop_itn_access_count +
    prop_livestock_ownership_count +
    prop_lowest_wealth_quantile_count +
    number_per_room_count +
    prop_wall_material_count +
    bym2_count,
  E = pop
)

occurence_obs <- bru_obs(
  family = "binomial",
  data = dm,
  formula = present ~ 0 +
    intercept_occurence +
    precipitation_occurence +
    temperature_occurence +
    soil_moisture_occurence +
    clay_content_occurence +
    time_to_treatment_occurence +
    prop_itn_access_occurence +
    prop_livestock_ownership_occurence +
    prop_lowest_wealth_quantile_occurence +
    number_per_room_occurence +
    prop_wall_material_occurence +
    bym2_occurence
) 

joint_mod <- bru(
  comps,
  truncated_poisson_obs,
  occurence_obs,
  options = list(
    control.fixed = list(mean = 0, prec = 1),
    control.compute = list(cpo = T, waic = T, dic = T)
  )
)

# Best model --------------------------------------------------------------

best_model <- list(zap, zip, joint_mod) |>
  setNames(c(
    "Zero-Adjusted Poisson (ZAP) model",
    "Zero-Inflated Poisson (ZIP) model",
    "Joint model"
  )) |>
  imap(\(x, y) {
    data.table(
      model = y,
      WAIC = x$waic$waic,
      DIC = x$dic$dic,
      logarithimic_score = -mean(log(x$cpo$cpo), na.rm = T)
    )
  }) |>
  rbindlist() |>
  _[, map_if(.SD, is.double, \(z) round(z, 2))]
fwrite(best_model, "data/model_comparison_metrics.csv", row.names = F)

# The joint model performed better than all the other 2 models
mod <- joint_mod
saveRDS(joint_mod, "data/joint_mod.rds")

# Predicting and generating samples ---------------------------------------

# Predicting
pred.mod <- predict(
  mod,
  dm,
  ~ {
    occurence_prob <-
      plogis(
        intercept_occurence +
          precipitation_occurence +
          temperature_occurence +
          soil_moisture_occurence +
          clay_content_occurence +
          time_to_treatment_occurence +
          prop_itn_access_occurence +
          prop_livestock_ownership_occurence +
          prop_lowest_wealth_quantile_occurence +
          number_per_room_occurence +
          prop_wall_material_occurence +
          bym2_occurence
      )
    lambda <- exp(
      intercept_count +
        precipitation_count +
        temperature_count +
        soil_moisture_count +
        clay_content_count +
        time_to_treatment_count +
        prop_itn_access_count +
        prop_livestock_ownership_count +
        prop_lowest_wealth_quantile_count +
        number_per_room_count +
        prop_wall_material_count +
        bym2_count
    )
    expect_param <- occurence_prob * lambda * pop
    expect <- expect_param / (1 - exp(-lambda * pop))
    variance <- expect * (1 - exp(-lambda * pop) * expect)
    list(
      occurence = occurence_prob,
      lambda = lambda,
      expect = expect,
      variance = variance,
      obs_prob = (1 - occurence_prob) * (N == 0) +
        (N > 0) * occurence_prob * dpois(N, expect_param) /
        (1 - dpois(0, expect_param))
    )
  },
  n.samples = 1000
)


# Generating posterior samples
posterior.samples <- generate(
  mod,
  dm,
  ~ {
    occurence_prob <-
      plogis(
        intercept_occurence +
          precipitation_occurence +
          temperature_occurence +
          soil_moisture_occurence +
          clay_content_occurence +
          time_to_treatment_occurence +
          prop_itn_access_occurence +
          prop_livestock_ownership_occurence +
          prop_lowest_wealth_quantile_occurence +
          number_per_room_occurence +
          prop_wall_material_occurence +
          bym2_occurence
      )
    lambda <- exp(
      intercept_count +
        precipitation_count +
        temperature_count +
        soil_moisture_count +
        clay_content_count +
        time_to_treatment_count +
        prop_itn_access_count +
        prop_livestock_ownership_count +
        prop_lowest_wealth_quantile_count +
        number_per_room_count +
        prop_wall_material_count +
        bym2_count
    )
    expect_param = occurence_prob * lambda * pop
    expect = expect_param / (1 - exp(-lambda * pop))
    variance = expect * (1 - exp(-lambda * pop) * expect)
    
    tibble(
      occurence = occurence_prob,
      lambda = lambda,
      expect = expect,
      variance = variance,
      obs_prob = (1 - occurence_prob) * (N == 0) +
        (N > 0) * occurence_prob * dpois(N, expect_param) /
        (1 - dpois(0, expect_param))
    )
    },
  n.samples = 1000,
  num.threads = "4:1"
)

posterior.df <- posterior.samples |> 
  imap(\(x, y){
    matrix(round(x$expect))
  }) |> 
  reduce(cbind)

# The forest plot --------------------------------------------------------

coeff_dfs <- tidy.inla(mod, exp = T) |> 
  filter(!terms %in% c("intercept_count", "intercept_occurence")) |> 
  mutate(model = fifelse(str_detect(terms, "count"), "Number of cases", "Occurence of VL")) |> 
  mutate(
    model = factor(model, levels = c("Occurence of VL", "Number of cases")),
    terms = case_when(
      terms == "ndvi_count" ~ "NDVI",
      terms == "precipitation_count" ~ "Precipitation",
      terms == "temperature_count" ~ "Temperature",
      terms == "soil_moisture_count" ~ "Soil moisture",
      terms == "clay_content_count" ~ "Soil clay content (0-20cm deep)",
      terms == "time_to_treatment_count" ~ "Travel time to treatment centers",
      terms == "prop_floor_material_count" ~ "Proportion of population living in households with earth/dung floor material",
      terms == "prop_itn_access_count" ~ "Proportion of population with ITN access",
      terms == "prop_livestock_ownership_count" ~ "Proportion of population living in households that own livestock",
      terms == "number_per_room_count" ~ "Number of people per sleeping room",
      terms == "prop_lowest_wealth_quantile_count" ~ "Proportion of population in the lowest wealth quantile",
      terms == "prop_wall_material_count" ~ "Proportion of the population living in households with bamboo or mud wall materials",
      terms == "ndvi_occurence" ~ "NDVI",
      terms == "precipitation_occurence" ~ "Precipitation",
      terms == "temperature_occurence" ~ "Temperature",
      terms == "soil_moisture_occurence" ~ "Soil moisture",
      terms == "clay_content_occurence" ~ "Soil clay content (0-20cm deep)",
      terms == "time_to_treatment_occurence" ~ "Travel time to treatment centers",
      terms == "prop_floor_material_occurence" ~ "Proportion of population living in households with earth/dung floor material",
      terms == "prop_itn_access_occurence" ~ "Proportion of population with ITN access",
      terms == "prop_livestock_ownership_occurence" ~ "Proportion of population living in households that own livestock",
      terms == "number_per_room_occurence" ~ "Number of people per sleeping room",
      terms == "prop_lowest_wealth_quantile_occurence" ~ "Proportion of population in the lowest wealth quantile",
      terms == "prop_wall_material_occurence" ~ "Proportion of the population living in households with bamboo or mud wall materials"
    ) |>
      str_wrap(30)
  )
fwrite(coeff_dfs, "data/forest_plot_estimates.csv", row.names = F)

# Number of cases
forest_plot_count <- coeff_dfs |>
  filter(model == "Number of cases") |> 
  ggplot(aes(x = mean, y = reorder(terms, mean))) +
  geom_point(color = "#0072B2", size = 3) +
  geom_errorbar(
    aes(xmin = `0.025quant`, xmax = `0.975quant`),
    width = .2,
    linewidth = 1,
    color = "#0072B2"
  ) +
  geom_vline(xintercept = 1, color = "#009E73", linetype = 2, linewidth = 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(color = "black", size = 17, hjust = .5, face = "bold"),
    axis.title = element_text(color = "black", size = 15),
    axis.text = element_text(color = "black", size = 15)
  ) +
  labs(x = "Incidence Risk Ratio", y = "Variable")

ggsave(
  plot = forest_plot_count,
  filename = "plots/forest_plot_count.png",
  dpi = 5e2,
  width = 3 * 5,
  height = 2 * 5,
  units = "in",
  bg = 'white'
)

# Forest plot for occurence
forest_plot_occurence <- coeff_dfs |> 
  filter(model == "Occurence of VL") |> 
  ggplot(aes(x = mean, y = reorder(terms, mean))) +
  geom_point(color = "#0072B2", size = 3) +
  geom_errorbar(
    aes(xmin = `0.025quant`, xmax = `0.975quant`),
    width = .2,
    linewidth = 1,
    color = "#0072B2"
  ) +
  geom_vline(xintercept = 1, color = "#009E73", linetype = 2, linewidth = 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(color = "black", size = 17, hjust = .5, face = "bold"),
    axis.title = element_text(color = "black", size = 15),
    axis.text = element_text(color = "black", size = 15)
  ) +
  labs(x = "Incidence Risk Ratio", y = "Variable")

ggsave(
  plot = forest_plot_occurence,
  filename = "plots/forest_plot_occurence.png",
  dpi = 5e2,
  width = 3 * 5,
  height = 2 * 5,
  units = "in",
  bg = 'white'
)

# Combined forest plot
combined_forest_plot <- coeff_dfs |>
  mutate(terms = factor(
    terms,
    levels = coeff_dfs |>
      filter(model == unique(model)[2]) |>
      arrange(mean) |>
      pull(terms)
  )) |>
  ggplot(aes(x = mean, y = terms)) +
  geom_point(color = "#0072B2", size = 3) +
  geom_errorbar(
    aes(xmin = `0.025quant`, xmax = `0.975quant`),
    width = 0.2,
    linewidth = 1,
    color = "#0072B2"
  ) +
  facet_wrap(~model) +
  geom_vline(xintercept = 1, color = "#009E73", linetype = 2, linewidth = 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(color = "black", size = 17, hjust = 0.5, face = "bold"),
    axis.title = element_text(color = "black", size = 15),
    axis.text = element_text(color = "black", size = 15),
    strip.text = element_text(color = "black", size = 17)
  ) +
  labs(x = "Effect size (ratio)", y = "Variable")
ggsave(
  plot = combined_forest_plot,
  filename = "plots/combined_forest_plot.png",
  dpi = 5e2,
  width = 4 * 5,
  height = 2 * 5,
  units = "in",
  bg = NULL
)

# Probability of exceedance ----------------------------------------------------

ward_exceedance1 <- t(apply(posterior.df, 1, \(x) c(
  round(mean(x)), quantile(x, probs = c(0.025, 0.5, 0.975)), length(which(x >= 1500)) / 1000
)))

ward_exceedance2 <- t(apply(posterior.df, 1, \(x) c(
  round(mean(x)), quantile(x, probs = c(0.025, 0.5, 0.975)), length(which(x >= mean(ward_exceedance1[,1]))) / 1000
)))

colnames(ward_exceedance1) <- c(
  "mean_cases",
  "0.025quant",
  "0.5quant",
  "0.975quant",
  "thresh_1500"
)

colnames(ward_exceedance2) <- c(
  "mean_cases",
  "0.025quant",
  "0.5quant",
  "0.975quant",
  "thresh_kenya_average"
) 
ward_exceedance2 <- data.frame(ward_exceedance2)

ward_exceedance <- cbind(
  data.table(dm)[, .(county, subcounty, ward, N, pop)], 
  ward_exceedance1
  ) |>
  as.data.table() |> 
  _[, let(
    mean_cases = fifelse(mean_cases == 0, NA, mean_cases),
    thresh_1500 = fifelse(thresh_1500 == 0, NA, thresh_1500),
    thresh_kenya_average = fifelse(ward_exceedance2$thresh_kenya_average == 0, NA, ward_exceedance2$thresh_kenya_averag)
  )]

# Risk at ward level ------------------------------------------------------

dm_out <- data.table(copy(pred.mod$expect)) |> 
  _[, let(
    risk = round(mean/pop * 1e5, 1),
    prev = mean/pop,
    occurence = (pred.mod$occurence$mean),
    ci_lower = round(q0.025/pop * 1e5, 1),
    ci_upper = round(q0.975/pop * 1e5, 1)
  )] |> 
  _[, let(uncertainty = ci_upper - ci_lower)] |> 
  _[, let(
    risk_cat = case_when(
      risk == 0 ~ "0",
      risk > 0 & risk < 5 ~ "1-5",
      risk >= 5 & risk < 17 ~ "5-17",
      risk >= 17 & risk < 51 ~ "17-51",
      risk >= 51 & risk < 93 ~ "51-93",
      risk >= 93 & risk < 191 ~ "93-191",
      risk >= 191 & risk < 258 ~ "191-258",   # 97.5-98th percentile
      risk >= 258 & risk < 788 ~ "258-788",   # 98-99th percentile
      risk >= 788 & risk < 957 ~ "788-957",   # 99-99.5th percentile
      risk >= 957 ~ ">=957",                    # 99.5-100th percentile
      TRUE ~ "Not categorized"
    ),
    
    occurence_cat = case_when(
      occurence < 0.01 ~ "<0.01",
      occurence >= 0.01 & occurence < 0.05 ~ "0.01-0.05",
      occurence >= 0.05 & occurence < 0.15 ~ "0.05-0.15",
      occurence >= 0.15 & occurence < 0.30 ~ "0.15-0.30",
      occurence >= 0.30 & occurence < 0.50 ~ "0.30-0.50",
      occurence >= 0.50 & occurence < 0.70 ~ "0.50-0.70",
      occurence >= 0.70 & occurence < 0.90 ~ "0.70-0.90",
      occurence >= 0.90 ~ ">=0.90",
      TRUE ~ NA_character_
    )
  )] |> 
  _[, let(
    risk_cat = factor(
      risk_cat, 
      levels = c("0", "1-5", "5-17", "17-51", "51-93", "93-191", "191-258", "258-788", "788-957", ">=957", "Not categorized"),
      ordered = TRUE
    ),
    occurence_cat = factor(
      occurence_cat,
      levels = c("<0.01", "0.01-0.05", "0.05-0.15", "0.15-0.30", "0.30-0.50", "0.50-0.70", "0.70-0.90", ">=0.90"),    
      ordered = TRUE
    ),
    risk = fifelse(risk == 0, NA_real_, risk),
    uncertainty = fifelse(uncertainty == 0, NA_real_, uncertainty)
  )] |>
  merge(adm3) |> 
  st_as_sf()
fwrite(data.table(dm_out)[, !"geometry"], "data/modelling_results.csv", row.names = F)

#' This are the quantiles we used to create the categories above
# quantile(dm_out$risk, na.rm = TRUE,
#          probs = c(0, 0.5, 0.7, 0.9, 0.95, 0.975, 0.98, 0.99, 0.995, 1))
#
# For the probability of occurence
# quantile(dm_out$occurence, 
#          probs = c(0, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 1),
#          na.rm = TRUE)

# Continous risk plot
pred_ward_incidence_cont_plot <- dm_out |> 
  ggplot() +
  geom_sf(aes(fill = risk), lwd = 0.1, color = "grey80") +
  geom_sf(data =  adm1, color = "black", fill = NA) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, na.value = "white") +
  theme_void() +
  theme(
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black", size = 8, hjust = .5)
    ) +
  labs(fill = "Predicted VL incidence\n(per 100,000 person years)") +
  guides(fill = guide_colorbar( barwidth = 1, barheight = 8))

ggsave(
  plot = pred_ward_incidence_cont_plot,
  filename = "plots/pred_ward_incidence_cont_plot.png",
  dpi = 5e2,
  width = 1 * 5,
  height = 1 * 5,
  units = "in",
  bg = "white"
)

# Categorical risk plot
pred_ward_incidence_cat_plot <- dm_out |> 
  ggplot() +
  geom_sf(aes(fill = risk_cat), lwd = 0.1, color = "grey80") +
  geom_sf(data =  adm1, color = "black", fill = NA) +
  scale_fill_manual(
    values = c("white", brewer.pal(9, "YlOrRd")),
    labels = c(
      "0", "1-5", "5-17", "17-51", "51-93", "93-191", "191-258", "258-788", "788-957", expression(paste("" >= 957))
    )
  ) +
  theme_void() +
  theme(
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black", size = 8, hjust = .5)
    ) +
  labs(fill = "Predicted incidence\n(per 100,000 person years)") 

ggsave(
  plot = pred_ward_incidence_cat_plot,
  filename = "plots/pred_ward_incidence_cat_plot.png",
  dpi = 5e2,
  width = 1 * 5,
  height = 1 * 5,
  units = "in",
  bg = "white"
)

# Continous predicted probability of occurence
pred_ward_occurence_cont_plot <- dm_out |> 
  ggplot() +
  geom_sf(aes(fill = occurence), lwd = 0.1, color = "grey80") +
  geom_sf(data =  adm1, color = "black", fill = NA) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, na.value = "white") +
  theme_void() +
  theme(
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black", size = 8, hjust = .5)
  ) +
  labs(fill = "Probability of VL occurence") +
  guides(fill = guide_colorbar( barwidth = 1, barheight = 8))

ggsave(
  plot = pred_ward_occurence_cont_plot,
  filename = "plots/pred_ward_occurence_cont_plot.png",
  dpi = 5e2,
  width = 1 * 5,
  height = 1 * 5,
  units = "in",
  bg = "white"
)

# Continous predicted probability of occurence
pred_ward_occurence_cat_plot <- dm_out |> 
  ggplot() +
  geom_sf(aes(fill = occurence_cat), lwd = 0.1, color = "grey80") +
  geom_sf(data =  adm1, color = "black", fill = NA) +
  scale_fill_manual(
    values = c("white", brewer.pal(7, "YlOrRd")),
    labels = c(
      "<0.01", "0.01-0.05", "0.05-0.15", "0.15-0.30", "0.30-0.50", "0.50-0.70", "0.70-0.90", expression(paste("" >= 0.90))
    )
  ) +
  theme_void() +
  theme(
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black", size = 8, hjust = .5)
  ) +
  labs(fill = "Probability of VL occurence") 

ggsave(
  plot = pred_ward_occurence_cat_plot,
  filename = "plots/pred_ward_occurence_cat_plot.png",
  dpi = 5e2,
  width = 1 * 5,
  height = 1 * 5,
  units = "in",
  bg = "white"
)

# Uncertainty
ward_uncertainty_plot <- dm_out |> 
  ggplot() +
  geom_sf(aes(fill = uncertainty), color = "grey80", lwd = 0.1) +
  geom_sf(data = adm1, fill = NA, color = "black") +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, na.value = "white") +
  theme_void() +
  theme(
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black", size = 8, hjust = .5)
    ) +
  labs(fill = "Uncertainty") +
  guides(fill = guide_colorbar( barwidth = 1, barheight = 8))  

ggsave(
  plot = ward_uncertainty_plot,
  filename = "plots/ward_uncertainty_plot.png",
  dpi = 5e2,
  width = 1 * 5,
  height = 1 * 5,
  units = "in",
  bg = "white"
)

# Combining the categorical plots for occurence and incidence
cat_incidence_occurence <- wrap_plots(
  pred_ward_occurence_cat_plot,
  pred_ward_incidence_cat_plot,
  ncol = 2
) +
  plot_annotation(
    title = "",
    tag_levels = list(c("(a)", "(b)"))
  ) &
  theme(
    legend.key.size = unit(1.5, "lines"),
    plot.title = element_text(size = 16, hjust = 0.5),
    plot.tag = element_text(size = 16, color = "black"),
    plot.tag.position = c("topleft"),
    plot.tag.location = "panel",
    # Add transparent background elements
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )
ggsave(
  plot = cat_incidence_occurence,
  filename = "plots/cat_incidence_occurence.png",
  dpi = 5e2,
  width = 2 * 5,
  height = 1 * 5,
  units = "in",
  bg = "white"
)

# Some summaries
data.table(dm_out)[, .N, by = risk_cat] |> 
  _[order(N)] |> 
  _[, let(total = sum(N))] |> 
  _[, let(prop = round(N/total * 100, 1))] |> 
  _[, !"total"]

data.table(dm_out)[, .N, by = occurence_cat] |> 
  _[order(N)] |> 
  _[, let(total = sum(N))] |> 
  _[, let(prop = round(N/total * 100, 1))] |> 
  _[, !"total"]

# Probability of exceedance -----------------------------------------------

# 1. Predicted number of cases
pred_ward_cases_plot <- merge(ward_exceedance, adm3) |> 
  st_as_sf() |> 
  ggplot() +
  geom_sf(aes(fill = mean_cases), color = "grey80", lwd = 0.1) +
  geom_sf(data = adm1, fill = NA, color = "black") +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, na.value = "white") +
  theme_void() +
  theme(
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black", hjust = .5, size = 9)
  ) +
  labs(fill = 'Predicted number\nof VL cases') +
  guides(fill = guide_colorbar( barwidth = 1, barheight = 8))
ggsave(
  plot = pred_ward_cases_plot,
  filename = "plots/pred_ward_cases_plot.png",
  dpi = 5e2,
  width = 1 * 5,
  height = 1 * 5,
  units = "in",
  bg = "white"
)

 # 2. Probability of exceeding 1500 cases
exceed_1500_cases_plot <- merge(ward_exceedance, adm3) |> 
  st_as_sf() |> 
  ggplot() +
  geom_sf(aes(fill = thresh_1500), color = "grey80", lwd = 0.1) +
  geom_sf(data = adm1, fill = NA, color = "black") +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, na.value = "white") +
  theme_void() +
  theme(
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black", hjust = .5, size = 9)
  ) +
  labs(fill = 'Probability of\nexceeding 1500 cases') +
  guides(fill = guide_colorbar( barwidth = 1, barheight = 8))
ggsave(
  plot = exceed_1500_cases_plot,
  filename = "plots/exceed_1500_cases_plot.png",
  dpi = 5e2,
  width = 1 * 5,
  height = 1 * 5,
  units = "in",
  bg = "white"
)

 # 2. Probability of exceeding kenya average cases
exceed_kenya_average_cases_plot <- merge(ward_exceedance, adm3) |> 
  st_as_sf() |> 
  ggplot() +
  geom_sf(aes(fill = thresh_kenya_average), color = "grey80", lwd = 0.1) +
  geom_sf(data = adm1, fill = NA, color = "black") +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, na.value = "white") +
  theme_void() +
  theme(
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black", hjust = .5, size = 9)
  ) +
  labs(fill = "Probability of\nexceeding country's average cases") +
  guides(fill = guide_colorbar( barwidth = 1, barheight = 8))
ggsave(
  plot = exceed_kenya_average_cases_plot,
  filename = "plots/exceed_kenya_average_cases_plot.png",
  dpi = 5e2,
  width = 1 * 5,
  height = 1 * 5,
  units = "in",
  bg = "white"
)

