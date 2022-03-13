
#----Setup----------------------------------------------------------------------
library(here)
library(tidyverse)
source(here("utils.R"))
load(here("data", "species_info.rda"))

# When was the last update?
load(here("data", "date.rda"))
pdats <- paste0(date, ":2025/01/01")

#----Create bioproject table for Brassicaceae-----------------------------------
load(here("data", "brassicaceae_projects.rda"))
species_df_brassicaceae <- species_info %>%
    filter(Family == "Brassicaceae")

retmax <- 30000
brassicaceae_old <- brassicaceae_projects
brassicaceae_new <- get_bp_table(species_df_brassicaceae$Species, pdats)

brassicaceae_projects <- rbind(
    brassicaceae_old, brassicaceae_new
)

#----Create bioproject table for Poaceae-----------------------------------
load(here("data", "poaceae_projects.rda"))
species_df_poaceae <- species_info %>%
    filter(Family == "Poaceae")

retmax <- 30000
poaceae_old <- poaceae_projects
poaceae_new <- get_bp_table(species_df_poaceae$Species, pdats)

poaceae_projects <- rbind(
    poaceae_old, poaceae_new
)

#----Create bioproject table for Fabaceae---------------------------------------
load(here("data", "fabaceae_projects.rda"))
species_df_fabaceae <- species_info %>%
    filter(Family == "Fabaceae")

retmax <- 30000
fabaceae_old <- fabaceae_projects
fabaceae_new <- get_bp_table(species_df_fabaceae$Species, pdats)

fabaceae_projects <- rbind(
    fabaceae_old, fabaceae_new
)

#----Create bioproject table for all other families-----------------------------
load(here("data", "otherfam_projects.rda"))
species_df_otherfam <- species_info %>%
    filter(!Family %in% c("Poaceae", "Fabaceae", "Brassicaceae"))

retmax <- 30000
otherfam_old <- otherfam_projects
otherfam_new <- get_bp_table(species_df_otherfam$Species, pdats)

otherfam_projects <- rbind(
    otherfam_old, otherfam_new
)

#----Save updated tables--------------------------------------------------------
brassicaceae_projects <- dplyr::distinct(brassicaceae_projects, .keep_all = TRUE)
fabaceae_projects <- dplyr::distinct(fabaceae_projects, .keep_all = TRUE)
poaceae_projects <- dplyr::distinct(poaceae_projects, .keep_all = TRUE)
otherfam_projects <- dplyr::distinct(otherfam_projects, .keep_all = TRUE)

save(brassicaceae_projects,
     file = here("data", "brassicaceae_projects.rda"),
     compress = "xz")

save(fabaceae_projects,
     file = here("data", "fabaceae_projects.rda"),
     compress = "xz")

save(poaceae_projects,
     file = here("data", "poaceae_projects.rda"),
     compress = "xz")

save(otherfam_projects,
     file = here("data", "otherfam_projects.rda"),
     compress = "xz")

# Finally, update date
date <- gsub("-", "/", Sys.Date())
save(date,
     file = here::here("data", "date.rda"),
     compress = "xz")





