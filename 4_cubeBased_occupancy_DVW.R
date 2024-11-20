library(sf)
library(leaflet)
library(tidyverse)
library(dplyr)
library(vecsets)

rm(list=ls())

#Get species list
species_list <- read.csv('./data/input/DVW_prior.csv', sep="\t")
species_list <- species_list %>%
  mutate_if(is.character, utf8::utf8_encode) %>%
  mutate_if(is.character, trimws)

#r get shapes
vla_1km <- st_read("./data/spatial/grids/vla_1km.geojson")

DVW_sectoren <- st_read(dsn="./data/spatial/DVW_SECTOREN_BASIS.shp")
DVW_sectoren <-st_transform(DVW_sectoren,st_crs(vla_1km))
DVW_sectoren <- st_make_valid(DVW_sectoren)

#Prepare "Kern"-regions
KERN <- DVW_sectoren %>%
  filter(TYPE == "KERN") %>%
  group_by(DSTRCT) %>%
  summarise(geometry = st_union(geometry),
            .groups = "drop") 

KERN <- st_intersection(vla_1km, KERN)%>%
  as.data.frame()

KERN_1km <- vla_1km %>% 
  filter(CELLCODE %in% KERN$CELLCODE)

KERN_1km <- left_join(KERN_1km, KERN[, c("CELLCODE", "DSTRCT")], by = "CELLCODE")

#This data has more rows then KERN, since some eea-cells are within 2 districts!! 3900 iso 3647

#r get cube data
IAS_Cube <- read_csv(file = "./data/input/be_DVW_species_cube.csv") %>%
  filter(year >= 2010) %>% 
  filter(speciesKey %in% species_list$speciesKey) %>% 
  filter(eea_cell_code %in% vla_1km$CELLCODE)

#r test data completeness
IAS_Cube_summarise <- IAS_Cube %>% 
  group_by(speciesKey) %>% 
  summarise(n_tot = sum(n, na.rm = TRUE))

not_in_cube <- species_list %>% 
  left_join(IAS_Cube_summarise, by = c("speciesKey" = "speciesKey")) %>% 
  filter(is.na(n_tot))

in_cube <- species_list %>% 
  left_join(IAS_Cube_summarise, by = c("speciesKey" = "speciesKey")) %>% 
  filter(!is.na(n_tot))

write.csv(not_in_cube, "./data/output/occupancy_DVW/Analyzed_species/species_present.txt")
write.csv(in_cube, "./data/output/occupancy_DVW/Analyzed_species/species_not_present.txt")

#r create dataframe for occupancy in kerngebieden}
df <- data.frame()

for (species_name in unique(in_cube$species)){
  
  print(species_name)
  GBIF_codes <- in_cube[which(in_cube$species==species_name),]$speciesKey
  print(GBIF_codes)
  IAS_Cube_temp <- IAS_Cube %>%
    filter(speciesKey %in% GBIF_codes)
  
  N_occ_spec <- as.numeric(n_distinct(IAS_Cube_temp$eea_cell_code)) # nodig, want zelfde cel voor verschillende jaren
  N_CORE <- as.numeric(n_distinct(KERN_1km$CELLCODE))
  N_occ_spec_in_CORE <- as.numeric(n_distinct(
    IAS_Cube_temp[which(IAS_Cube_temp$eea_cell_code %in% KERN_1km$CELLCODE),]$eea_cell_code))
  df <- rbind(df, data.frame(species_name,
                             N_occ_spec,
                             N_occ_spec_in_CORE,
                             N_occ_spec_in_CORE/N_occ_spec, #percentage van de bezette km hokken in Kern
                             round(100*(N_occ_spec_in_CORE/N_occ_spec), digits = 1),     # Idem, as percent (%), for report
                             N_occ_spec_in_CORE/N_CORE, #percentage van de totale CORE bezet door deze soort
                             round(100*(N_occ_spec_in_CORE/N_CORE), digits = 2)))       # Idem, as promille (‰), for report
}

colnames(df) <- c("species",
                  "N_occ",
                  "N_occ_in_CORE",
                  "share_in_CORE",           
                  "%_in_CORE",
                  "share_of_CORE",
                  "%_of_CORE")

write.csv(df, "./data/output/occupancy_DVW/occupancy_DVW.csv")

#make barcharts for occupancy in kerngebieden}
df <- df %>%
  arrange(N_occ)
temp_plot<-ggplot(df, aes(x = reorder(species, N_occ), y = N_occ))  +
  geom_bar(stat = "identity", fill = "magenta4") +
  coord_flip() +
  xlab("Species") +
  ylab("Aantal bezette km²-hokken") +
  ggtitle("Soorten volgens aantal bezette km²-hokken in Vlaanderen")

ggsave(temp_plot, file="./data/output/occupancy_DVW/occupancy_VL.png", width = 16, height = 16, units = "cm", dpi=200)

df <- df %>%
  arrange(N_occ_in_CORE)
temp_plot<-ggplot(df, aes(x = reorder(species, N_occ_in_CORE), y = N_occ_in_CORE))  +
  geom_bar(stat = "identity", fill = "magenta4") +
  coord_flip() +
  xlab("Species") +
  ylab("Aantal bezette km²-hokken") +
  ggtitle("Soorten volgens aantal bezette km²-hokken in kerngebieden")

ggsave(temp_plot, file="./data/output/occupancy_DVW/occupancy_CORE.png", width = 16, height = 16, units = "cm", dpi=200)


#Run analysis on DSTRCTS seperately
df_list <- list()
unique_districts <- unique(KERN_1km$DSTRCT)

for (i in unique_districts) {
  # Subset KERN_1km for the current district
  KERN_sub <- KERN_1km[KERN_1km$DSTRCT == i, ]
  
  # Initialize an empty data frame to store the results
  df <- data.frame()
  
  for (species_name in unique(in_cube$species)) {
    GBIF_codes <- in_cube[in_cube$species == species_name, ]$speciesKey
    IAS_Cube_temp <- IAS_Cube %>% filter(speciesKey %in% GBIF_codes)
    N_occ_spec <-
      as.numeric(n_distinct(IAS_Cube_temp$eea_cell_code))
    N_CORE <- as.numeric(n_distinct(KERN_sub$CELLCODE))
    N_occ_spec_in_CORE <-     as.numeric(n_distinct(IAS_Cube_temp[IAS_Cube_temp$eea_cell_code %in% KERN_sub$CELLCODE, ]$eea_cell_code))
    
    df <- rbind(
      df,
      data.frame(
        i,
        #District
        species_name,
        GBIF_codes,
        N_occ_spec,
        #aantal hokken over Vlaanderen
        N_occ_spec_in_CORE,
        #aantal hokken binnen dit district
        N_occ_spec_in_CORE / N_occ_spec,
        #aandeel km² hokken van Vlaanderen in dit district
        round(100 * (N_occ_spec_in_CORE / N_occ_spec), digits = 1),
        N_occ_spec_in_CORE / N_CORE,
        #aandeel van km² hokken binnen district
        round(100 * (N_occ_spec_in_CORE / N_CORE), digits = 2)
      )
    )
  }
  
  colnames(df) <-
    c(
      "District",
      "species",
      "speciesKey",
      "aantal km² in VL",
      "aantal km² in district",
      "aandeel km² in district",
      "% in district",
      "aandeel km² hokken van district",
      "% van district"
    )
  df_list[[i]] <- df
  
  # Write to CSV, named based on the district
  write.table(
    df,
    file = paste0(
      "./data/output/occupancy_DVW/Districts/occupancy_DVW_",
      i,
      ".csv"
    ),
    sep = ";",
    dec = ",",
    row.names = FALSE,
    quote = FALSE,
    fileEncoding = "UTF-8"
  )
}

occupancy_districts <- bind_rows(df_list)
write.table(
  occupancy_districts,
  file = paste0(
    "./data/output/occupancy_DVW/Districts/occupancy_District.csv"
  ),
  sep = ";",
  dec = ",",
  row.names = FALSE,
  quote = FALSE
)



#plots for every district
clean <- function(filename) {
  cleaned <- gsub("[[:space:]]", "", filename)  # Vervang spaties door underscores
  cleaned <- gsub("[(),]", "", cleaned)          # Verwijder komma's en haakjes
  return(cleaned)
}

unique_districts <- unique(occupancy_districts$District)
for (i in unique_districts) {
  sub_data <- occupancy_districts[occupancy_districts$District == i,]
  
  # Sorteer de soorten op basis van aantal km² in district
  sub_data <- sub_data[order(-sub_data$`aantal km² in district`), ]
  
  temp_plot <- ggplot(sub_data, aes(x = reorder(species, -`aantal km² in district`), y = `aantal km² in district`)) +
    geom_bar(stat = "identity", fill="magenta4") +
    coord_flip() +
    ggtitle(paste("Soorten in", i)) +
    xlab("Soort") +
    ylab("Aantal km² in district")
  ggsave(temp_plot, file=paste0("./data/output/occupancy_DVW/Per_District/",clean(i),".png"), width = 16, height = 16, units = "cm", dpi=200)
}

#plotten per soort
unique_species <- unique(occupancy_districts$species)

for (species_name in unique_species) {
  sub_data <- occupancy_districts[occupancy_districts$species == species_name,]
  sub_data <- sub_data %>%
    arrange(desc(`% in district`)) %>%
    mutate(District = factor(District, levels = unique(District)))
  
  # Barchart voor aandeel in district
  temp_plot <- ggplot(sub_data, aes(x = District, y = `% in district`)) +
    geom_bar(stat = "identity", fill="magenta4") +
    coord_flip() +
    ggtitle(paste("Aandeel van totaal km²-hokken", species_name, "in district")) +
    xlab("District") +
    ylab("Aandeel km² in district")
  
  ggsave(temp_plot, 
         file=paste0("../data/output/occupancy_DVW/Per_Soort/",species_name,".png"),
         width = 16, 
         height = 16, 
         units = "cm", dpi=200)
}


#r return ubiquitous species


# Extract unique combinations of species and districts where 'aantal km² in district' is greater than 0
unique_combinations <- occupancy_districts %>%
  filter(`aantal km² in district` > 0) %>%
  select(species, District)

# Calculate the total number of unique districts
total_districts <- length(unique(occupancy_districts$District))

# Calculate the number of districts each species is found in and then compute the percentage
result <- unique_combinations %>%
  group_by(species) %>%
  summarise(num_districts = n_distinct(District)) %>%
  mutate(`% of districts` = (num_districts / total_districts) * 100)

# Create a wider format dataframe indicating presence (1) or absence (0) of each species in each district
Full_result <- unique_combinations %>%
  mutate(present = 1) %>%
  spread(key = District, value = present, fill = 0) %>%
  left_join(result, by = "species") %>%
  arrange( by= `% of districts`)


# Print the result
print(result)




