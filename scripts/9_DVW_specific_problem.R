library(sf)
library(leaflet)
library(tidyverse)
library(dplyr)
library(RColorBrewer)


rm(list=ls())

#Get species list
datasetKey <- "23e95da2-6095-4778-b893-9af18a310cb6"
DVW_list <- name_usage(datasetKey = datasetKey)

#Make a list of taxa from DVW-checklist from the dataframe
species_list <- DVW_list$data %>% 
  filter(origin == "SOURCE") %>%
  rename(sK = speciesKey) %>%
  rename(speciesKey = nubKey)

#r get shapes
vla_1km <- st_read("./data/spatial/grids/vla_1km.geojson")

DVW_sectoren <- st_read(dsn="./data/spatial/DVW_SECTOREN_BASIS.shp")
DVW_sectoren <-st_transform(DVW_sectoren,st_crs(vla_1km))
DVW_sectoren <- st_make_valid(DVW_sectoren)

DVW_1km <- st_intersection(vla_1km, DVW_sectoren)%>%
  as.data.frame() %>% 
  select(CELLCODE, DSTRCT, TYPE)

#r get cube data
IAS_Cube <- read_csv(file = "./data/input/be_DVW_species_cube.csv") %>%
  filter(year >= 2010) %>% 
  filter(speciesKey %in% species_list$speciesKey) %>% 
  filter(eea_cell_code %in% vla_1km$CELLCODE)

#numbers of 1m² in Flanders
Vlaanderen <- nrow(vla_1km)

# Aantal 1km² in de kerngebieden
n_kern <- DVW_1km %>%
  filter(TYPE == "KERN") %>%
  distinct(CELLCODE) %>%
  nrow()

kern_occupancy <- IAS_Cube %>%
  filter(eea_cell_code %in% DVW_1km$CELLCODE[DVW_1km$TYPE == "KERN"]) %>%
  distinct(speciesKey, eea_cell_code) %>%
  group_by(speciesKey) %>%
  summarise(num_hokken_kern = n(), .groups = 'drop')

# Berekenen van unieke bezette hokken in Vlaanderen per speciesKey
vlaanderen_occupancy <- IAS_Cube %>%
  distinct(speciesKey, eea_cell_code) %>%
  group_by(speciesKey) %>%
  summarise(num_hokken_vlaanderen = n(), .groups = 'drop')

# Berekenen van het totale aantal hokken in Vlaanderen en kerngebieden
total_hokken_vlaanderen <- length(unique(DVW_1km$CELLCODE))
total_hokken_kern <- length(unique(DVW_1km$CELLCODE[DVW_1km$TYPE == "KERN"]))

# Samenvoegen van de informatie in een nieuwe tabel
result <- vlaanderen_occupancy %>%
  left_join(kern_occupancy, by = "speciesKey") %>%
  mutate(
    total_hokken_vlaanderen = total_hokken_vlaanderen,
    total_hokken_kern = total_hokken_kern
  ) %>%
  replace_na(list(num_hokken_kern = 0)) # Zorg ervoor dat ontbrekende waarden worden vervangen door 0

# Het resultaat bekijken
print(result)



  

