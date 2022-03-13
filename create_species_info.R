
#----Set up---------------------------------------------------------------------
library(tidyverse)
library(rvest)
library(bears)
library(here)
source(here("utils.R"))

#----Create a data frame of included species and their info---------------------
# Only species with >50% of complete BUSCOs will be included
species_info <- readr::read_tsv(
    "https://raw.githubusercontent.com/almeidasilvaf/plantgenomes/master/data/nature_plants_data.txt",
    show_col_types = FALSE
) %>%
    mutate(Species = stringr::str_c(Genus, species, sep = " ")) %>%
    dplyr::rename(Complete_BUSCO = `BUSCO % complete (combined)`) %>%
    dplyr::select(Family, Species, Complete_BUSCO) %>%
    drop_na() %>%
    mutate(Complete_BUSCO = str_replace_all(
        Complete_BUSCO, 
        c("," = "\\.",
          "%" = ""
        ))
    ) %>%
    mutate(Complete_BUSCO = as.numeric(Complete_BUSCO)) %>%
    filter(Complete_BUSCO > 50) 

rnaseq_samples <- unlist(lapply(species_info$Species, count_available_samples))
species_info$RNAseq_samples <- rnaseq_samples
species_info <- species_info %>%
    filter(RNAseq_samples > 0) %>%
    arrange(Family, -RNAseq_samples)

readr::write_tsv(species_info, file = here("data", "species_info.tsv"))
