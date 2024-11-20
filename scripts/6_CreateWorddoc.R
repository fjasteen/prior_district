library(officer)
library(flextable)

# Read the species data
datasetKey <- "23e95da2-6095-4778-b893-9af18a310cb6"
DVW_list <- name_usage(datasetKey = datasetKey)

#Make a list of taxa from DVW-checklist from the dataframe
species_list <- DVW_list$data %>% 
  filter(origin == "SOURCE") %>%
  rename(sK = speciesKey) %>%
  rename(speciesKey = nubKey)



status <- read.csv("./data/output/absent_sporadic_present/absent-sporadic-present.csv")
status <- rename(status, canonicalName = X)
species_data <- left_join(DVW_prior, status, by = "canonicalName")
species_data <- species_data %>%
  mutate(status.y = replace_na(status.y, "Absent"))

# Loop through each species
for (species in species_data$canonicalName) {
  
  # Create a new Word document
  doc <- read_docx()
  
  # Add some text
  doc <- body_add_par(doc, paste0("Verspreiding van ", species), style = "heading 1")
  doc <- body_add_par(doc, "Dit is wat tekst.", style = "Normal")
  
  # Add a second title for presence
  doc <- body_add_par(doc, "Aanwezigheid tussen 2010-2023", style = "heading 2")
  
  status <- species_data %>%
    filter(canonicalName == .env$species) %>%
    pull(status.y)
  
  checkboxes <- ifelse(status == "Present", "\t\u2611 aanwezig \u2610 afwezig \u2610 sporadisch",
                       ifelse(status == "Sporadic", "\t\u2610 aanwezig \u2611 afwezig \u2610 sporadisch",
                              "\t\u2610 aanwezig \u2611 afwezig \u2610 sporadisch")
  )
  
  doc <- body_add_par(doc, checkboxes, style = "Normal")
  
  species <- gsub(" ", "_", species)
  
  # Add main image
  image_path <- paste0("./data/output/spread_maps/", species, ".png")
  if (file.exists(image_path)) {
    doc <- body_add_img(doc, src = image_path, width = 6, height = 2.72)
  }
  
  # Add a second title for additional figures
  doc <- body_add_par(doc, "Is de soort opkomend?", style = "heading 2")
  
  # Paths for additional images
  image_path1 <- paste0("./data/output/emerging_status/", species, "_occupancy_corrected_Flanders.jpeg")
  image_path2 <- paste0("./data/output/emerging_status/", species, "_occupancy_corrected_DVW_CORE.jpeg")
  
  if (file.exists(image_path1) && file.exists(image_path2)) {
    doc <- body_add_par(doc, "In Vlaanderen", style = "heading 3")
    doc <- body_add_img(doc, src = image_path1, width = 3, height = 1.36)
    doc <- body_add_par(doc, "In kerngebieden", style = "heading 3")
    doc <- body_add_img(doc, src = image_path2, width = 3, height = 1.36)
  } else if (file.exists(image_path1)) {
    doc <- body_add_par(doc, "In Vlaanderen", style = "heading 3")
        doc <- body_add_img(doc, src = image_path1, width = 3, height = 1.36)
  } else if (file.exists(image_path2)) {
    doc <- body_add_par(doc, "In kerngebieden", style = "heading 3")
    doc <- body_add_img(doc, src = image_path2, width = 3, height = 1.36)
  } else {
    message(paste("No additional images for ", species, " found. Skipping..."))
  }
  
  # Save the document
  doc_name <- paste0("./data/output/word/Document_for_", species, ".docx")
  print(doc, target = doc_name)
}
