# Load necessary libraries
library(spdep)
library(sf)
library(leaflet)
library(tidyverse)
library(dplyr)
library(RColorBrewer)

# Read and preprocess species list
species_list <- read.csv('./data/input/DVW_prior.csv', sep="\t") %>%
  mutate(across(where(is.character), ~utf8::utf8_encode(.x) %>% trimws()))

# Read and transform spatial data
vla_1km <- st_read("./data/spatial/grids/vla_1km.geojson")
DVW_sectoren <- st_read(dsn="./data/spatial/DVW_SECTOREN_BASIS.shp") %>%
  st_transform(st_crs(vla_1km)) %>%
  st_make_valid()

DVW_1km <- st_intersection(vla_1km, DVW_sectoren)%>%
  as.data.frame() %>% 
  select(CELLCODE, DSTRCT, TYPE)

# Aggregate features by 'DSTRCT' and remove empty features
DVW_sectoren_agg <- DVW_sectoren %>%
  filter(!st_is_empty(.)) %>%
  group_by(DSTRCT) %>%
  summarize(geometry = st_union(geometry), .groups = 'drop')

# Find neighbors for each district
neighbors <- poly2nb(DVW_sectoren_agg)
polygon_names <- DVW_sectoren_agg$DSTRCT

# Map neighbors to district names
neighbor_list <- map(1:length(neighbors), ~list(district = polygon_names[.x], neighbors = polygon_names[neighbors[[.x]]]))

# Function to get neighbors of a given district
get_neighbors <- function(district) {
  district_index <- match(district, polygon_names)
  if (!is.na(district_index)) {
    neighbor_list[[district_index]]$neighbors
  } else {
    character(0)
  }
}

# Prepare IAS_Cube_DSTRCT data frame (assuming it's already loaded)
IAS_Cube <- read_csv(file = "./data/input/be_DVW_species_cube.csv") %>%
  filter(year >= 2010) %>% 
  filter(speciesKey %in% species_list$speciesKey) %>% 
  filter(eea_cell_code %in% vla_1km$CELLCODE) %>%
  filter(!is.na(speciesKey))

IAS_Cube_DSTRCT <- IAS_Cube %>%
left_join(DVW_1km, by = c("eea_cell_code" = "CELLCODE")) 

# Function to check species presence in district and neighbors
check_species_presence <- function(species, district) {
  species_in_district <- any(IAS_Cube_DSTRCT$speciesKey == species & IAS_Cube_DSTRCT$DSTRCT == district, na.rm = TRUE)
  if (!species_in_district) {
    neighbors <- get_neighbors(district)
    neighbors_with_species <- neighbors[neighbors %in% IAS_Cube_DSTRCT$DSTRCT[IAS_Cube_DSTRCT$speciesKey == species]]
    if (length(neighbors_with_species) == 0) NA else neighbors_with_species
  } else {
    NA
  }
}

# Apply function over each speciesKey and district combination
results <- expand.grid(speciesKey = unique(IAS_Cube_DSTRCT$speciesKey), DSTRCT = unique(IAS_Cube_DSTRCT$DSTRCT)) %>%
  mutate(neighbors_with_species = map2(speciesKey, DSTRCT, ~check_species_presence(.x, .y))) %>%
  filter(!is.na(neighbors_with_species) & !is.na(neighbors_with_species))

# Display the results
print(results)
results$neighbors <- 1

# Calculate absent/present species:
# Assuming IAS_Cube_DSTRCT and complete_species_list are already loaded

# Create all possible combinations of speciesKey and DSTRCT
all_combinations <- expand.grid(speciesKey = na.omit(species_list$speciesKey), 
                                DSTRCT = unique(IAS_Cube_DSTRCT$DSTRCT))

# Function to check presence
check_presence <- function(species, district) {
  if (any(IAS_Cube_DSTRCT$speciesKey == species & IAS_Cube_DSTRCT$DSTRCT == district, na.rm = TRUE)) {
    return(0)  # Present
  } else {
    return(1)  # Absent
  }
}

# Apply the function to each row
all_combinations$absent <- mapply(check_presence, all_combinations$speciesKey, all_combinations$DSTRCT)

# Ensure that no NA values in key columns before merging
all_combinations <- na.omit(all_combinations)

# Merge all_combinations with results
final_df <- merge(all_combinations, results, by = c("speciesKey", "DSTRCT"), all = TRUE)

# [Function for check_species_absence_in_neighbors and its application goes here]

# Merge with species_list to add verbatim_name
final_df <- merge(final_df, species_list[, c("speciesKey", "species")], by = "speciesKey", all.x = TRUE)

# Replace NA values where appropriate
final_df[is.na(final_df)] <- 0

# Reorder columns to make verbatim_name the first column
final_df <- final_df[c("species", setdiff(names(final_df), "species"))]

# View the updated dataframe
print(final_df)

# Functie om te controleren of een soort aanwezig is in een district maar niet in de naburige districten
check_species_absence_in_neighbors <- function(species, district) {
  species_in_district <- any(IAS_Cube_DSTRCT$speciesKey == species & IAS_Cube_DSTRCT$DSTRCT == district, na.rm = TRUE)
  if (species_in_district) {
    neighbors <- get_neighbors(district)
    absent_in_all_neighbors <- all(!neighbors %in% IAS_Cube_DSTRCT$DSTRCT[IAS_Cube_DSTRCT$speciesKey == species])
    if (absent_in_all_neighbors) TRUE else FALSE
  } else {
    FALSE
  }
}

# Pas deze functie toe over elk speciesKey en DSTRCT combinatie
results_with_absence <- expand.grid(speciesKey = unique(IAS_Cube_DSTRCT$speciesKey), DSTRCT = unique(IAS_Cube_DSTRCT$DSTRCT)) %>%
  mutate(absent_in_neighbors = mapply(check_species_absence_in_neighbors, speciesKey, DSTRCT))

# Voeg deze resultaten toe aan de final_df
final_df <- merge(final_df, results_with_absence, by = c("speciesKey", "DSTRCT"), all = TRUE)

final_df <- merge(final_df, species_list[, c("speciesKey", "species")], by = "speciesKey", all.x = TRUE)

# Herschik de kolommen indien nodig
final_df <- final_df[c("species", setdiff(names(final_df), "species"))]

# Toon de bijgewerkte dataframe
print(final_df)