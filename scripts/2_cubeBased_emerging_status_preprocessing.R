library(tidyverse) # To do data science
library(tidylog) # To provide feedback on dplyr functions
library(progress) # To add progress bars
library(here) # To find files
library(lubridate) # To work with dates
library(rgbif) # To get taxa from publisehd unified checklist
library(sf) # To perform spatial operations

#clear memory
rm(list=ls())

#NOTE: the species cube is (too) heavy to push to the git. Download locally

###### DATA PREPARATION ########################################################
# species list DVW 
#Get datasetKey from the online checklist Flemish Waterways
datasetKey <- "23e95da2-6095-4778-b893-9af18a310cb6"
DVW_list <- name_usage(datasetKey = datasetKey)

#Make a list of taxa from DVW-checklist from the dataframe
species_list <- DVW_list$data %>% 
  filter(origin == "SOURCE")

# eea_1km_vlaanderen
vla_1km <- st_read("./data/spatial/grids/vla_1km.geojson")

#read local DVW cube
df <- read_csv(file = "./data/input/be_DVW_species_cube.csv",
               col_types = cols(year = col_double(), eea_cell_code = col_character(),
                                speciesKey = col_double(), n = col_double(), 
                                min_coord_uncertainty = col_double()), na = "") %>%
  select(-min_coord_uncertainty) %>%
  rename(obs = n, taxonKey = speciesKey) %>%
  filter(taxonKey %in% species_list$nubKey,
         eea_cell_code %in% vla_1km$CELLCODE,
         year >= 2000)

###### MAKE BASELINE REFERENCE #################################################
# Read baseline occurrence data for standardisation
df_bl <- read_csv(
  file = "https://zenodo.org/records/3635510/files/be_classes_cube.csv?download=1",
  col_types = cols(year = col_double(), eea_cell_code = col_character(),
                   classKey = col_double(), n = col_double(), min_coord_uncertainty = col_double()), 
  na = ""
) %>%
  select(-min_coord_uncertainty) %>%
  rename(cobs = n) %>%
  filter(eea_cell_code %in% vla_1km$CELLCODE, year >= 2000)

# Read and process spatial
# EEA cells in DVW core area
DVW_sectoren <- st_read(dsn="./data/spatial/DVW_SECTOREN_BASIS.shp") %>%
  st_transform(st_crs(vla_1km)) %>%
  st_make_valid()

KERN_1km <- vla_1km %>%
  filter(CELLCODE %in% (
    st_intersection(
      vla_1km, 
      DVW_sectoren %>% filter(TYPE == "KERN")
    )$CELLCODE
  ))

#protected_areas
df_areas <- read_tsv(
  "https://raw.githubusercontent.com/trias-project/indicators/master/data/interim/intersect_EEA_ref_grid_protected_areas.tsv",
  na = "")

#Get a table of all 1km² cells and indicate whether they are in the core or extended
df_areas <- df_areas %>%
  mutate(DVW_CORE = ifelse(CELLCODE %in% KERN_1km$CELLCODE, TRUE, FALSE)) %>%
  select(CELLCODE, DVW_CORE) %>%
  filter(CELLCODE %in% vla_1km$CELLCODE)

# Add informative columns
taxon_key <- df %>%
  distinct(taxonKey) %>%
  pull()

pb <- progress_bar$new(total = length(taxon_key))

spec_names <- map_df(taxon_key,
                     function(k) {
                       pb$tick()
                       name_usage(key = k)$data
                     }) %>%
  select(taxonKey = key, canonicalName, scientificName, kingdomKey, classKey) %>%
  mutate(canonicalName = ifelse(is.na(canonicalName), scientificName, canonicalName))

spec_names <- spec_names %>%
  group_by(canonicalName) %>%
  add_tally() %>%
  ungroup() %>%
  mutate(canonicalName = if_else(n > 1, scientificName, canonicalName)) %>%
  select(-c(n, scientificName))

#Om enkel te vergelijken met Liliopsida, Magnoliopsida,Polypodiopsida
class_key <- spec_names %>%
  distinct(classKey) %>%
  filter(!is.na(classKey)) %>%
  pull()

pb <- progress_bar$new(total = length(class_key))

kingdom_class <- map_df(class_key,
                        function(x) {
                          pb$tick()
                          name_usage(key = x)$data
                        }) %>%
  select(classKey, class, kingdomKey, kingdom)

# add class & kingdom
spec_names <- spec_names %>%
  left_join(kingdom_class %>% distinct(.data$classKey, .data$class), by = "classKey") %>%
  left_join(kingdom_class %>% distinct(.data$kingdomKey, .data$kingdom), by = "kingdomKey")

spec_names %>% head(10)
spec_names %>% filter(is.na(classKey))

#extract coordinates from unique cellcodes
df_xy <- df %>%
  distinct(eea_cell_code) %>%
  bind_cols(tibble(x = unlist(str_extract_all(unique(df$eea_cell_code),
                                              pattern = "(?<=E)\\d+")),
                   y = unlist(str_extract_all(unique(df$eea_cell_code),
                                              pattern = "(?<=N)\\d+"))) %>%
              mutate_all(as.integer))
df_xy %>% head()



df_bl_xy <- df_bl %>%
  distinct(eea_cell_code) %>%
  bind_cols(tibble(x = unlist(str_extract_all(unique(df_bl$eea_cell_code),
                                              pattern = "(?<=E)\\d+")),
                   y = unlist(str_extract_all(unique(df_bl$eea_cell_code),
                                              pattern = "(?<=N)\\d+"))) %>%
              mutate_all(as.integer))
df_bl_xy %>% head()


# Add information about presence in DVW_CORE (join cube and KERN_1KM)
df <- df %>%
  left_join(KERN_1km, by = c("eea_cell_code" = "CELLCODE"))
df %>% head()

df_bl <-
  df_bl %>% left_join(KERN_1km, by = c("eea_cell_code" = "CELLCODE"))
df_bl %>% head()


# create time series
df_cc <- df %>%
  group_by(taxonKey) %>%
  summarize(begin_year = min(year)) %>%
  right_join(df %>% distinct(taxonKey, eea_cell_code), by = "taxonKey") %>%
  select(taxonKey, begin_year, eea_cell_code)
df_cc %>% head()

# create time series per eea_code and species
make_time_series <- function(eea_cell_code, taxonKey, begin_year, last_year ) {
  expand_grid(eea_cell_code = eea_cell_code,
              taxonKey = taxonKey,
              year = seq(from = begin_year, to = last_year))
}

# create timeseries slots
df_ts <- pmap_dfr(df_cc, .f = make_time_series, last_year = year(Sys.Date()))

## Add data
df_ts <- df_ts %>%
  left_join(df %>% select(taxonKey, year, eea_cell_code, obs), #add info
            by = c("taxonKey", "year", "eea_cell_code")) %>%
  left_join(df_areas %>% select(CELLCODE, DVW_CORE), #add membership to KERN
            by = c("eea_cell_code" = "CELLCODE")) %>%
  left_join(spec_names %>% select(taxonKey, classKey), #add classKey
            by = "taxonKey")

# Research effort correction
df_ts <- df_ts %>%
  left_join(df_bl %>% select(year, eea_cell_code, classKey, cobs), # add baseline data (at class level) 
            by = c("year", "eea_cell_code", "classKey")) %>%
  mutate(cobs = cobs - obs) %>% #diminished by obs of specific alien taxon
  replace_na(list(cobs = 0, obs = 0)) %>%
  mutate(core_cobs = if_else(cobs > 0, 1, 0),
         core_obs = if_else(obs > 0, 1, 0)) %>%
  select(taxonKey, year, eea_cell_code, obs, core_obs, cobs, core_cobs, classKey, DVW_CORE)


write_tsv(df_ts, file = "./data/intermediate/df_timeseries.tsv", na = "")
write_tsv(spec_names, file = "./data/intermediate/timeseries_taxonomic_info.tsv", na = "")
