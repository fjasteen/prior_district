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

#This data has more rows then KERN, since some eea-cells are within 2 districts!! 3900 iso 3647

#r get cube data
IAS_Cube <- read_csv(file = "./data/input/be_DVW_species_cube.csv") %>%
  filter(year >= 2010) %>% 
  filter(speciesKey %in% species_list$speciesKey) %>% 
  filter(eea_cell_code %in% vla_1km$CELLCODE)

IAS_Cube_summarise <- IAS_Cube %>% 
  group_by(speciesKey) %>% 
  summarise(n_tot = sum(n, na.rm = TRUE))


species_list_subset <- species_list %>%
  select(speciesKey, verbatim_name)

#Construct the dataframe with all cubes containing the speciesname, districtname year and number of occurrences.
IAS_Cube_DSTRCT <- IAS_Cube %>%
  left_join(DVW_1km, by = c("eea_cell_code" = "CELLCODE")) %>%
  left_join(species_list_subset, by = "speciesKey")


# Create the output directory if it doesn't exist
output_dir <- "./data/output/Districts/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

district_cells_count <- DVW_1km %>%
  group_by(DSTRCT) %>%
  summarise(cells_count = n_distinct(CELLCODE))

# Filter out NA values in DSTRCT and aggregate the data
filtered_data <- IAS_Cube_DSTRCT %>%
  filter(!is.na(DSTRCT)) %>%
  distinct(DSTRCT, verbatim_name, eea_cell_code) %>%
  group_by(DSTRCT, verbatim_name) %>%
  summarise(unique_squares = n()) %>%
  arrange(DSTRCT, desc(unique_squares))

filtered_data <- filtered_data %>%
  left_join(district_cells_count, by = "DSTRCT") %>%
  mutate(occupied_percentage = (unique_squares / cells_count) * 100)

# Identify the top 10 species based on unique square counts
top_species <- filtered_data %>%
  group_by(verbatim_name) %>%
  summarise(total_unique_squares = sum(unique_squares)) %>%
  arrange(desc(total_unique_squares)) %>%
  top_n(10, total_unique_squares) %>%
  pull(verbatim_name)

# Define distinct colors for the top 10 species
distinct_colors <- RColorBrewer::brewer.pal(min(10, length(top_species)), "Set3")
names(distinct_colors) <- top_species

# Define a color palette for the remaining species
remaining_species <- setdiff(unique(filtered_data$verbatim_name), top_species)
remaining_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(remaining_species))
names(remaining_colors) <- remaining_species

# Combine the two sets of colors
colors <- c(distinct_colors, remaining_colors)

# Get the list of unique districts
unique_districts <- unique(filtered_data$DSTRCT)

# Create and Save Plots Individually with Normalized Data
for (district in unique_districts) {
  district_data <- filter(filtered_data, DSTRCT == district)
  plot <- ggplot(district_data, aes(x = reorder(verbatim_name, -occupied_percentage), y = occupied_percentage, fill = verbatim_name)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors) +
    ylim(0, 100) +  # Set Y-axis limits
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(title = paste("District:", district),
         x = "Species", 
         y = "Percentage of Occupied 1 km² Cells")
  
  # Save the plot
  ggsave(filename = paste0(output_dir, district, ".png"), plot = plot, width = 10, height = 6)
}


# Calculate the number of unique 1km² cells within KERN of each district
district_cells_count_kern <- DVW_1km %>%
  filter(TYPE == "KERN") %>%
  group_by(DSTRCT) %>%
  summarise(cells_count = n_distinct(CELLCODE))

# Filter IAS_Cube_DSTRCT for entries where TYPE is KERN
filtered_data_kern <- IAS_Cube_DSTRCT %>%
  filter(TYPE == "KERN", !is.na(DSTRCT)) %>%
  distinct(DSTRCT, verbatim_name, eea_cell_code) %>%
  group_by(DSTRCT, verbatim_name) %>%
  summarise(unique_squares = n()) %>%
  arrange(DSTRCT, desc(unique_squares))

# Join with district_cells_count_kern and calculate occupied_percentage
filtered_data_kern <- filtered_data_kern %>%
  left_join(district_cells_count_kern, by = "DSTRCT") %>%
  mutate(occupied_percentage = (unique_squares / cells_count) * 100)

# Use the same color settings and top_species as before

# Get the list of unique districts for KERN
unique_districts_kern <- unique(filtered_data_kern$DSTRCT)

# Create and Save Plots Individually for KERN data
for (district in unique_districts_kern) {
  district_data <- filter(filtered_data_kern, DSTRCT == district)
  plot <- ggplot(district_data, aes(x = reorder(verbatim_name, -occupied_percentage), y = occupied_percentage, fill = verbatim_name)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors) +
    ylim(0, 100) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(title = paste("District:", district, " - KERN"),
         x = "Species", 
         y = "Percentage of Occupied 1 km² Cells")
  
  # Save the plot with "_KERN" in the filename
  ggsave(filename = paste0(output_dir, district, "_KERN.png"), plot = plot, width = 10, height = 6)
}
