library(sf)
library(leaflet)
library(tidyverse)
library(dplyr)
library(vecsets)
library(nngeo)
library(rgbif)

gc()

#####DATA PREPARATION: CHECKLIST & SPATIAL DATA ################################
#Get datasetKey from the online checklist Flemish Waterways
datasetKey <- "23e95da2-6095-4778-b893-9af18a310cb6"
DVW_list <- name_usage(datasetKey = datasetKey)

#Make a list of taxa from DVW-checklist from the dataframe
species_list <- DVW_list$data %>% 
  filter(origin == "SOURCE")


#get spatial data
vla_1km <- st_read("./data/spatial/grids/vla_1km.geojson")
Vlaanderen_grenzen <- st_read("./data/spatial/flanders_wgs84.geojson")
Provincies_grenzen <- st_read("./data/spatial/Provincies.geojson")
Waterlopen <- st_read("./data/spatial/Waterlopen.geojson")
DVW_sectoren <- st_read(dsn="./data/spatial/DVW_SECTOREN_BASIS.shp")
DVW_sectoren <- DVW_sectoren %>%
  st_transform(st_crs(vla_1km)) %>%
  st_make_valid()

#Construct polygons of the 17 separate districts
DVW_sectoren_extended <- DVW_sectoren %>%
  group_by(DSTRCT) %>%
  summarise(geometry = st_union(geometry)) %>%
  st_remove_holes() 

#Select all km² squares in Core Areas of DVW
KERN_1km <- vla_1km %>%
  filter(CELLCODE %in% (
    st_intersection(
      vla_1km, 
      DVW_sectoren %>% filter(TYPE == "KERN")
    )$CELLCODE
  ))

#upload DVW species cube
IAS_Cube <- read_csv(file = "./data/input/be_DVW_species_cube.csv") %>%
    filter(year >= 2000,
          speciesKey %in% species_list$nubKey,
          eea_cell_code %in% vla_1km$CELLCODE)

##### TEST DATA COMPLETENESS ###################################################
# Summarize IAS_Cube by speciesKey 
IAS_Cube_summarise <- IAS_Cube %>%
  group_by(speciesKey) %>%
  summarise(n_tot = sum(n, na.rm = TRUE))

# Perform the join operation once
joined_data <- species_list %>%
  left_join(IAS_Cube_summarise, by = c("nubKey" = "speciesKey"))

# Filter species not in cube
not_in_cube <- joined_data %>%
  filter(is.na(n_tot))
write_csv(not_in_cube, "./data/intermediate/taxonkeys_not_in_cube.csv")

# Filter species in cube
in_cube <- joined_data %>%
  filter(!is.na(n_tot))

# Create lists for present and not present species
List_species_present <- unique(in_cube$species)
List_species_not_present <- vsetdiff(unique(species_list$species), 
                                     List_species_present)

# Write the lists to output files
write(List_species_present, "./data/output/species_present.txt")
write(List_species_not_present, "./data/output/species_not_present.txt")

##### PLOT SPREADMAPS FOR PRESENT SPECIES ######################################

# write function to create species spreadmap
create_species_plot <- function(species_name) {
  print(species_name)
  
  GBIF_codes <- in_cube[in_cube$species == species_name,]$nubKey
  print(GBIF_codes)
  
  IAS_Cube_temp <- IAS_Cube %>%
    filter(speciesKey %in% GBIF_codes)
  print(IAS_Cube_temp)
  
  IAS_Cube_1km_spec <- merge(vla_1km,
                             IAS_Cube_temp,
                             by.x = 'CELLCODE',
                             by.y = 'eea_cell_code',
                             all.y = TRUE)
  
# Identify cells that are in KERN areas
  IAS_Cube_1km_spec$is_KERN <- IAS_Cube_1km_spec$CELLCODE %in% KERN_1km$CELLCODE
  
  temp_plot <- ggplot() +
    geom_sf(data = Vlaanderen_grenzen, fill = "#f0f0f0", size = 0.2) +
    geom_sf(data = Waterlopen, colour = "lightblue") +
    geom_sf(data = DVW_sectoren_extended, fill = NA, size = 0.05, colour = "green4") +
    # geom_sf(data = Provincies_grenzen, fill = NA, size = 0.2) +
    geom_sf(data = IAS_Cube_1km_spec, aes(color = factor(is_KERN))) +  # Color based on core area status
    scale_color_manual(
      values = c("orange", "red"), 
      labels = c("Vlaanderen", "DVW_Kern"),
      guide = guide_legend(override.aes = list(size = 1, shape = 16)))  +  # Set colors here
    theme_void() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    theme(plot.margin = grid::unit(c(0, 0, 0, 0), "mm")) +
    coord_sf()
  
  ggsave(temp_plot, file = paste0("./data/output/spread_maps/", gsub(" ", "_", species_name), ".png"), width = 15, height = 6.4, units = "cm", dpi = 200)
}

# Loop through all unique species
for (species_name in unique(in_cube$species)) {
  create_species_plot(species_name)
}