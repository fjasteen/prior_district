library(sf)
library(leaflet)
library(tidyverse)
library(dplyr)
library(vecsets)
library(readr)
library(rgbif)


#Clear
rm(list=ls())
setwd("C:/Users/frederique_steen/Documents/GitHub/prior_district/")

datasetKey <- "23e95da2-6095-4778-b893-9af18a310cb6"
DVW_list <- name_usage(datasetKey = datasetKey)

#Make a list of taxa from DVW-checklist from the dataframe
species_list <- DVW_list$data %>% 
  filter(origin == "SOURCE") %>%
  rename(sK = speciesKey) %>%
  rename(speciesKey = nubKey)

#r get shapes
vla_1km <- st_read("./data/spatial/grids/vla_1km.geojson")
Vlaanderen_grenzen <- st_read("./data/spatial/flanders_wgs84.geojson")
Provincies_grenzen <- st_read("./data/spatial/Provincies.geojson")
Waterlopen <- st_read("./data/spatial/Waterlopen.geojson")
DVW_sectoren <- st_read(dsn="./data/spatial/DVW_SECTOREN_BASIS.shp")
DVW_sectoren <-st_transform(DVW_sectoren,st_crs(vla_1km))
DVW_sectoren <- st_make_valid(DVW_sectoren)

unique_DSTRCT <- unique(DVW_sectoren$DSTRCT)

#This scrips generates absent-present-sporadic for each KERNGEBIED for each District
#If removing   filter(TYPE == "KERN") you do it for the total of the district.

for (current_DSTRCT in unique_DSTRCT) {
  # Create directory for current DISTRICT
  dir_name <- paste0("./data/output/absent_sporadic_present/", current_DSTRCT)
  dir.create(dir_name, showWarnings = FALSE)

DSTRCT_EEA_ind <- DVW_sectoren %>% 
  filter(DSTRCT == current_DSTRCT) %>%
  filter(TYPE == "KERN") %>%
  st_union() %>%
  st_intersects(vla_1km)

DSTRCT_EEA <- vla_1km[unlist(DSTRCT_EEA_ind), ]

#Get Cube data
#Made a seperate Cube for the prior_list species + stored locally
IAS_Cube <- read_csv(file = "./data/input/be_DVW_species_cube.csv")

gc()
IAS_Cube_DSTRCT <- IAS_Cube %>%
  filter(year >= 2010) %>% 
  filter(speciesKey %in% species_list$speciesKey) %>% 
  filter(eea_cell_code %in% DSTRCT_EEA$CELLCODE)


#Test data completeness
IAS_Cube_summarise <- IAS_Cube_DSTRCT %>% 
  group_by(speciesKey) %>% 
  summarise(n_tot = sum(n, na.rm = TRUE))

not_in_cube <- species_list %>% 
  left_join(IAS_Cube_summarise, by = c("speciesKey")) %>% 
  filter(is.na(n_tot))

write_csv(not_in_cube, "./data/intermediate/taxonkeys_not_in_cube.csv")

in_cube <- species_list %>% 
  left_join(IAS_Cube_summarise, by = "speciesKey") %>% 
  filter(!is.na(n_tot))

not_in_cube$species %in% in_cube$species

List_species_present <- unique(in_cube$species)
List_species_not_present <- vsetdiff(unique(species_list$species) , List_species_present)
write(List_species_present, paste0(dir_name,"/species_present.txt"))
write(List_species_not_present, paste0(dir_name,"/species_not_present.txt"))


keep_all <- data.frame(year = c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2021)) ## Keep all species
keep_all$year <- as.factor(keep_all$year)                                  ## For compatibility with below

for (species_name in unique(in_cube$species)){
  print(species_name)
  GBIF_codes<- in_cube[which(in_cube$species==species_name),]$speciesKey
  print(GBIF_codes)
  IAS_Cube_temp <- IAS_Cube_DSTRCT %>%
    filter(speciesKey %in% GBIF_codes)
  print(IAS_Cube_temp)
  IAS_Cube_1km_spec <- merge(vla_1km,
                             IAS_Cube_temp,
                             by.x='CELLCODE',
                             by.y='eea_cell_code',
                             all.y=TRUE)
  IAS_Cube_sporadic <- IAS_Cube_1km_spec %>% 
    group_by(year)%>%
    summarize(n=n())
  
  IAS_Cube_sporadic$year <- as.factor(IAS_Cube_sporadic$year)
  
  data <- data.frame(year = c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021))
  
  data$year <- as.factor(data$year)
  
  x <- merge(data, IAS_Cube_sporadic)
  x <- replace(x, is.na(x), 0)
  
  if (nrow(x) == 0) {
    print("x is empty, skipping plotting.")
  } else {
    temp_plot <- ggplot(data=x, aes(x=year, y=n)) +
      geom_bar(stat="identity") +
      labs(y = expression ("Occupancy in"~km^2)) +
      scale_x_discrete(drop = FALSE)
    ggsave(temp_plot, file=paste0( dir_name, "/", species_name,".png"), dpi=200)
  }

  
  
  keep_all <- keep_all %>% 
    left_join(x[,1:2], by = "year") %>% 
    rename(!!species_name:=n)  ## Geef kolom naam van soort (m.o.: gebruik !! := want variabele)
  
}


#Summarize accross species.
showme <- as.data.frame(t(keep_all)[-1, , drop = FALSE])
colnames(showme) <- t(keep_all)[1,]

showme$status <- NA

# Loop through each row to assign a status
for (i in 1:nrow(showme)) {
  row_data <- showme[i, 1:(ncol(showme) - 1)]  # Exclude the "status" column
  
  # Check if all years are NA (Absent)
  if (all(is.na(row_data))) {
    showme$status[i] <- "Absent"
    next
  }
  
  # Check if the species wasn't seen in the last 5 years (2018-2022)
  last_five_years <- tail(row_data, 5)
  if (all(is.na(last_five_years))) {
    showme$status[i] <- "Sporadic"
    next
  }
  
  # Check for any NA after a non-NA value (Sporadic)
  non_na_indices <- which(!is.na(row_data))
  first_non_na <- min(non_na_indices)
  na_after_non_na = any(diff(non_na_indices) > 1)
  
  if (na_after_non_na || (max(non_na_indices) <= (ncol(showme) - 1 - 5))) {
    showme$status[i] <- "Sporadic"
    next
  }
  
  # If none of the conditions are met, the species is Present
  showme$status[i] <- "Present"
}

# Convert the "status" column to a factor for better data handling
showme$status <- as.factor(showme$status)

write.csv(showme, paste0(dir_name, "/absent-sporadic-present.csv"))

}

#Generate one table for each distrit:
parent_dir <- "./data/output/absent_sporadic_present"

# Get the list of directories (subsets for each district)
subdirs <- list.dirs(parent_dir, full.names = TRUE, recursive = FALSE)

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each directory
for (dir in subdirs) {
  # Define the path to the specific file
  file_path <- file.path(dir, "absent-sporadic-present.csv")
  
  # Read the table, assuming it's a CSV file
  table_data <- read_csv(file_path)
  
  # Assuming the 'status' column is always in the same position and the first column is the row names
  # Select only the first column and the 'status' column
  selected_data <- table_data %>% select(1, status)
  
  # Rename the 'status' column to the folder name, extracting from the directory path
  folder_name <- basename(dir)
  colnames(selected_data)[2] <- folder_name  # Assuming 'status' is the second column after row names
  
  # Combine the data by column, merging by the row names
  if (nrow(combined_data) == 0) {
    combined_data <- selected_data
  } else {
    combined_data <- full_join(combined_data, selected_data, by = colnames(selected_data)[1])
  }
}

write.csv(combined_data, "./data/output/combined_status_data.csv", row.names = TRUE)

#Reformat the combined table of all districts
# Load necessary libraries
library(gt)
library(dplyr)

# Assuming your data frame is named 'combined_data'
# Replace NA values with "Absent"
# Rename the first column to 'Species'
colnames(combined_data)[1] <- "Species"
combined_data[is.na(combined_data)] <- "Absent"

# Simplify the status if needed
combined_data <- combined_data %>%
  mutate(across(everything(), ~ case_when(
    . == "Present" ~ "P",
    . == "Sporadic" ~ "S",
    . == "Absent" ~ "A",
    TRUE ~ .
  )))

# Create a nicely formatted table using gt
gt_table <- combined_data %>%
  gt() %>%
  tab_header(
    title = "Species presence",
    subtitle = "across Different Districts"
  ) %>%
  data_color(
    # Exclude the Species column from coloring
    colors = scales::col_factor(
      palette = c("P" = "white", "S" = "orange", "A" = "darkred"),
      domain = c("P", "S", "A")
    )
  ) %>%
  cols_align(
    align = "center",
    columns = everything()
  ) 

# Print the table
print(gt_table)
