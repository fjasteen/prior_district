library(sf)
library(leaflet)
library(tidyverse)
library(dplyr)
library(vecsets)


#Clear
rm(list=ls())
setwd("C:/Users/frederique_steen/Documents/GitHub/prior_district/")


#get species list
#Get datasetKey from the online checklist Flemish Waterways
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

#Get Cube data
#Made a seperate Cube for the prior_list species + stored locally
IAS_Cube <- read_csv(file = "./data/input/be_DVW_species_cube.csv")

gc()
IAS_Cube <- IAS_Cube %>%
  filter(year >= 2010) %>% 
  filter(speciesKey %in% species_list$speciesKey) %>% 
  filter(eea_cell_code %in% vla_1km$CELLCODE)


#Test data completeness
IAS_Cube_summarise <- IAS_Cube %>% 
  group_by(speciesKey) %>% 
  summarise(n_tot = sum(n, na.rm = TRUE))

not_in_cube <- species_list %>% 
  left_join(IAS_Cube_summarise, by = c("speciesKey")) %>% 
  filter(is.na(n_tot))

write_csv(not_in_cube, "./data/intermediate/taxonkeys_not_in_cube.csv")

in_cube <- species_list %>% 
  left_join(IAS_Cube_summarise, by = "speciesKey") %>% 
  filter(!is.na(n_tot))

not_in_cube$Soort %in% in_cube$species

List_species_present <- unique(in_cube$Soort)
List_species_not_present <- vsetdiff(unique(species_list$Soort) , List_species_present)
write(List_species_present, "./data/output/species_present.txt")
write(List_species_not_present, "./data/output/species_not_present.txt")


keep_all <- data.frame(year = c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2021)) ## Keep all species
keep_all$year <- as.factor(keep_all$year)                                  ## For compatibility with below

for (species_name in unique(in_cube$species)){
  print(species_name)
  GBIF_codes<- in_cube[which(in_cube$species==species_name),]$speciesKey
  print(GBIF_codes)
  IAS_Cube_temp <- IAS_Cube %>%
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
  
  temp_plot<- ggplot(data=x, aes(x=year, y=n)) +
    geom_bar(stat="identity")+
    labs(y = expression ("Occupancy in"~km^2))+
    scale_x_discrete(drop = FALSE)
  
  ggsave(temp_plot, file=paste0("./data/output/absent-sporadic-present/plots/", species_name,".png"), dpi=200)
  
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

write.csv(showme, "./data/output/absent_sporadic_present/absent-sporadic-present.csv")
