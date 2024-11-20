
# 3: emerging
# 2: potentially emerging
# 1: unclear
# 0: not emerging

library(tidyverse) # To do data science
library(trias)

# clear memory
rm(list=ls())
setwd("C:/Users/frederique_steen/Documents/GitHub/prior_district/")


##### GET DATA #################################################################
#Get datasetKey from the online checklist Flemish Waterways
datasetKey <- "23e95da2-6095-4778-b893-9af18a310cb6"
DVW_list <- name_usage(datasetKey = datasetKey)

#Make a list of taxa from DVW-checklist from the dataframe
species_list <- DVW_list$data %>% 
  filter(origin == "SOURCE")

# read in preprocessing files
df_ts <- read_tsv("./data/intermediate/df_timeseries.tsv", na = "")

spec_names <- read_tsv(file = here::here("./data/intermediate/timeseries_taxonomic_info.tsv"), na = "") %>%
  select(taxonKey, canonicalName) %>%
  filter(taxonKey %in% df_ts$taxonKey)

# First & Last evaluation year
last_year <- lubridate::year(Sys.Date()) - 1
first_year <- last_year - 22

evaluation_years <- seq(first_year, last_year)
evaluation_years

#Take the subset of the data
df_ts_compact <- df_ts %>%
  group_by(taxonKey, year, classKey) %>%
  summarise(obs = sum(obs), cobs = sum(cobs), ncells = sum(core_obs), c_ncells = sum(core_cobs)) %>%
  ungroup() %>%
  left_join(spec_names, by = "taxonKey") %>% #add species
  filter(year <= last_year) %>%
  filter(year >= first_year)

df_ts_compact_DVW_CORE <- df_ts %>%
  filter(DVW_CORE == TRUE) %>%
  group_by(taxonKey, year, classKey) %>%
  summarise(obs = sum(obs), cobs = sum(cobs), ncells = sum(core_cobs),  c_ncells = sum(core_cobs)) %>%
  ungroup() %>%
  left_join(spec_names, by = "taxonKey")%>% #add species
  filter(year <= last_year) %>%
  filter(year >= first_year) 

appearing_taxa_to_remove <- df_ts_compact %>% #Remove all taxa that are not there in first_year
  group_by(taxonKey) %>%
  summarize(begin_year = min(year)) %>%
  filter(begin_year > min(evaluation_years)) %>%
  pull(taxonKey)

df_ts_compact <- df_ts_compact %>%
  filter(!taxonKey %in% appearing_taxa_to_remove)

df_ts_compact_DVW_CORE <- df_ts_compact_DVW_CORE %>%
  filter(!taxonKey %in% appearing_taxa_to_remove)

spec_names <- spec_names %>% #remove appearing taxa from spec_names
  filter(!taxonKey %in% appearing_taxa_to_remove)

# See https://trias-project.github.io/trias/reference/apply_decision_rules.html
# GAM output has priority (below), but these decision rules allow to output em_status for species with failed GAM output
# Only interested in occupancy. Can also be done for observations (... y_var = "obs" ...)

em_decision_rules_occupancy_FL <- map_dfr(evaluation_years,
                                          ~ apply_decision_rules(df = df_ts_compact, y_var = "ncells", eval_year = .))

em_decision_rules_occupancy_CORE <- map_dfr(evaluation_years,
                                            ~ apply_decision_rules(df = df_ts_compact_DVW_CORE, y_var = "ncells", eval_year = .))

em_decision_rules_occupancy_FL_2020 <- em_decision_rules_occupancy_FL %>% 
  filter(year == 2020) %>% 
  left_join(species_list, by= c('taxonKey'='nubKey'))

write_csv(em_decision_rules_occupancy_FL_2020, "./data/output/emerging_status/output_decision_rules_occupancy_FL_2020.csv")

em_decision_rules_occupancy_CORE_2020 <- em_decision_rules_occupancy_CORE %>% 
  filter(year == 2020) %>% 
  left_join(species_list, by= c('taxonKey'='nubKey'))

write_csv(em_decision_rules_occupancy_CORE_2020, "./data/output/emerging_status/output_decision_rules_occupancy_CORE_2020.csv")


# GAM for occupancy data in all of Flanders

dir_name <- "./data/output/emerging_status"
taxon_keys <- unique(df_ts_compact$taxonKey)
taxon_names <- unique(df_ts_compact$canonicalName)
gam_occupancy_BE <- map2(
  taxon_keys, taxon_names,
  function(t, n) {
    df_key <- df_ts_compact %>%
      filter(taxonKey == t) %>%
      filter(year <= max(evaluation_years))
    class_key <- unique(df_key[["classKey"]])
    if (!is.na(class_key)) {
      apply_gam(
        df = df_key,
        y_var = "ncells",
        eval_years = evaluation_years,
        type_indicator = "occupancy",
        taxon_key = t,
        name = n,
        baseline_var = "c_ncells",
        dir_name = dir_name,
        saveplot = FALSE,
        y_label = "occupancy (km2)"
      )
    } else {
      apply_gam(
        df = df_key,
        y_var = "ncells",
        eval_years = evaluation_years,
        type_indicator = "occupancy",
        taxon_key = t,
        name = n,
        dir_name = dir_name,
        saveplot = FALSE,
        y_label = "occupancy (km2)"
      )
    }
  }
)

names(gam_occupancy_BE) <- taxon_keys

for (tK in unique(df_ts_compact$taxonKey)) {
  tK <- as.character(tK)
  print(tK)
  if (!is.null(gam_occupancy_BE[[tK]]$plot)) {
    
    ##
    ggplot(data = gam_occupancy_BE[[tK]]$output, aes(x = year, y = ncells)) +
      theme_bw() +
      ggtitle("Vlaanderen") +
      xlab("jaar") + ylab("aantal km²-hokken") +
      theme(axis.text = element_text(size = 9), 
            axis.title = element_text(size = 9),
            legend.text = element_text(size = 9),
            legend.title= element_text(size=9),
            legend.position = "right",  # Position legend at the bottom
            legend.box = "horizontal",
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove margin around legend keys
            legend.box.margin = margin(t = -5, b = -5, r=-5, l=-5),
            legend.spacing.y= unit(0.2,"cm"),
            panel.background = element_rect(fill = "grey80"),  # Set background color to light grey
            panel.grid.major = element_line(color = "white", size= 0.2),  # Set major grid color to white
            panel.grid.minor = element_line(color = "white", size= 0.1) ) +
      geom_point(color = "black", size = 1) +                                             # Raw data
      geom_line(aes(x = .data$year, y = .data$fit), color = "black", size = 0.2) +       # The GAM line
      geom_point(aes(x = year, y = fit, color = as.factor(em_status)), shape = 4 , size = 1) +  # Points on the GAM line colored by em_status
      geom_ribbon(aes(ymax = ucl, ymin = lcl), fill = grey(0.5), alpha = 0.4) +           # Error on the GAM line
      scale_color_manual(
        values = c("darkgreen", "grey40", "red", "darkred"),  # Colors for 0, 1, 2, 3
        labels = c("niet opkomend", "onduidelijk", "potentieel opkomend", "opkomend"),  # Legend labels
        name = "") # Legend title (empty in this case)
    
    ggsave(filename = paste0("./data/output/emerging_status/",
                             gsub(" ","_",first(gam_occupancy_BE[[tK]]$output$canonicalName)),
                             "_occupancy_corrected_Flanders.jpeg"),
           width = 12, height = 7.2, units= "cm")
    ##
    
  }
  else {print('no plot')}
}


# GAM for occupancy data in CORE areas

dir_name <- "./data/output/emerging_status"
taxon_keys <- unique(df_ts_compact_DVW_CORE$taxonKey)
taxon_names <- unique(df_ts_compact_DVW_CORE$canonicalName)
gam_occupancy_CORE <- map2(
  taxon_keys, taxon_names,
  function(t, n) {
    df_key <- df_ts_compact_DVW_CORE %>%
      dplyr::filter(taxonKey == t)
    class_key <- unique(df_key[["classKey"]])
    if (!is.na(class_key)) {
      apply_gam(
        df = df_key,
        y_var = "ncells",
        eval_years = evaluation_years,
        type_indicator = "occupancy",
        taxon_key = t,
        name = n,
        baseline_var = "c_ncells",
        df_title = "DVW_Kern",
        dir_name = dir_name,
        saveplot = FALSE,
        y_label = "occupancy (km2)"
      )
    } else {
      apply_gam(
        df = df_key,
        y_var = "ncells",
        eval_years = evaluation_years,
        type_indicator = "occupancy",
        taxon_key = t,
        name = n,
        df_title = "DVW_Kern",
        dir_name = dir_name,
        saveplot = FALSE,
        y_label = "occupancy (km2)"
      )
    }
  }
)

names(gam_occupancy_CORE) <- taxon_keys

for (tK in unique(df_ts_compact_DVW_CORE$taxonKey)) {
  tK <- as.character(tK)
  print(tK)
  if (!is.null(gam_occupancy_CORE[[tK]]$plot)) {
    
    ## Zie hoger voor toelichtingen
    ggplot(data = gam_occupancy_CORE[[tK]]$output, aes(x = year, y = ncells)) +
      theme_bw() +
      ggtitle("Kerngebieden") +
      xlab("jaar") + ylab("aantal km²-hokken") +
      theme(axis.text = element_text(size = 9), 
            axis.title = element_text(size = 9),
            legend.text = element_text(size = 9),
            legend.title= element_text(size=9),
            legend.position = "right",  # Position legend at the bottom
            legend.box = "horizontal",
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove margin around legend keys
            legend.box.margin = margin(t = -5, b = -5, r=-5, l=-5),
            legend.spacing.y= unit(0.2,"cm"),
            panel.background = element_rect(fill = "grey80"),  # Set background color to light grey
            panel.grid.major = element_line(color = "white", size= 0.2),  # Set major grid color to white
            panel.grid.minor = element_line(color = "white", size= 0.1) ) +
      geom_point(color = "black", size = 1) +                                             # Raw data
      geom_line(aes(x = .data$year, y = .data$fit), color = "black", size = 0.2) +       # The GAM line
      geom_point(aes(x = year, y = fit, color = as.factor(em_status)), shape = 4 , size = 1) +  # Points on the GAM line colored by em_status
      geom_ribbon(aes(ymax = ucl, ymin = lcl), fill = grey(0.5), alpha = 0.4) +           # Error on the GAM line
      scale_color_manual(
        values = c("darkgreen", "grey40", "red", "darkred"),  # Colors for 0, 1, 2, 3
        labels = c("niet opkomend", "onduidelijk", "potentieel opkomend", "opkomend"),  # Legend labels
        name = "") # Legend title (empty in this case)
    
    ggsave(filename = paste0("./data/output/emerging_status/",
                             gsub(" ","_",first(gam_occupancy_CORE[[tK]]$output$canonicalName)),
                             "_occupancy_corrected_DVW_CORE.jpeg"),
           width = 12, height = 7.2, units = "cm")
    ##
  }
  else {print('no plot')}
}


# Save results GAM

method_em <- gam_occupancy_BE[[1]]$em_summary$method[1]

write_tsv(map_dfr(gam_occupancy_BE, function(x) {x$output}), na = "",
          path = paste0("./data/output/emerging_status/output_GAM_", method_em, "_occupancy_FL.tsv"))

method_em <- gam_occupancy_CORE[[1]]$em_summary$method[1]

write_tsv(map_dfr(gam_occupancy_CORE, function(x) {x$output}), na = "",
          path = paste0("./data/output/emerging_status/output_GAM_", method_em, "_occupancy_CORE.tsv"))

em_GAM_occ_FL <- map_dfr(gam_occupancy_BE, function(x) {x$em_summary})

em_GAM_occ_FL_2020 <- em_GAM_occ_FL %>% 
  filter(year == 2020) %>%
  left_join(species_list, by= c('taxonKey'='speciesKey'))

write_csv(em_GAM_occ_FL_2020, paste0("./data/output/emerging_status/output_GAM_",
                                     method_em, "_occupancy_FL_2020.csv"))

em_GAM_occ_CORE <- map_dfr(gam_occupancy_CORE, function(x) {x$em_summary})

em_GAM_occ_CORE_2020 <- em_GAM_occ_CORE %>% 
  filter(year == 2020) %>%
  left_join(species_list, by= c('taxonKey'='speciesKey'))

write_csv(em_GAM_occ_CORE_2020, paste0("./data/output/emerging_status/output_GAM_",
                                       method_em, "_occupancy_CORE_2020.csv"))

# Make master output

strip_decrules_FL <- em_decision_rules_occupancy_FL_2020 %>% 
  transmute(taxonKey, dec_rules_FL = em_status)

strip_decrules_CORE <- em_decision_rules_occupancy_CORE_2020 %>% 
  transmute(taxonKey, dec_rules_CORE = em_status)

strip_GAM_FL <- em_GAM_occ_FL_2020 %>% 
  transmute(taxonKey, GAM_FL = em_status)

strip_GAM_CORE <- em_GAM_occ_CORE_2020 %>% 
  transmute(taxonKey, GAM_CORE = em_status)

master <- species_list %>% 
  select("species", "scientificName", "speciesKey") %>% 
  full_join(strip_decrules_FL, by = c("speciesKey" = "taxonKey")) %>% 
  full_join(strip_decrules_CORE, by = c("speciesKey" = "taxonKey")) %>% 
  full_join(strip_GAM_FL, by = c("speciesKey" = "taxonKey")) %>%
  full_join(strip_GAM_CORE, by = c("speciesKey" = "taxonKey"))


# Consensus value: if GAM available, then GAM; if GAM not available, then dec_rules
master$consensus_FL <- master$GAM_FL
data.table::setDT(master)[is.na(consensus_FL), consensus_FL := dec_rules_FL] # Maybe other options possible

master$consensus_CORE <- master$GAM_CORE
data.table::setDT(master)[is.na(consensus_CORE), consensus_CORE := dec_rules_CORE] # Maybe other options possible

rep_string = c("3" = "opkomend", "2" = "potentieel opkomend", "1" = "onduidelijk", "0" = "niet opkomend")
master$FL <- str_replace_all(master$consensus_FL, rep_string)
master$CORE <- str_replace_all(master$consensus_CORE, rep_string)

rm(strip_decrules_FL, strip_decrules_CORE, strip_GAM_FL, strip_GAM_CORE, rep_string)
write_delim(master, "./data/output/emerging_status/output_em-status_all.csv", delim = ";")
