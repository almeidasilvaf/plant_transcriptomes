Original data - 1st version
================

# Data download

The data sets in `data/` are automatically updated every week. However,
an initial data set had to be created first, so new lines could be
appended after each week. The data sets were created with the following
code:

## species_info.tsv

``` r
library(tidyverse)
library(rvest)
library(bears)
library(here)

source(here("utils.R"))

#----Create a data frame of included species and their info---------------------
# Only species with >50% of complete BUSCOs will be included
species_df <- readr::read_tsv(
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

rnaseq_samples <- unlist(lapply(species_df$Species, count_available_samples))
species_df$RNAseq_samples <- rnaseq_samples
species_df <- species_df %>%
    filter(RNAseq_samples > 0) %>%
    arrange(Family, -RNAseq_samples)

retmax <- max(species_df$RNAseq_samples + 10)
```

## brassicaceae_projects.rda

``` r
species_df_brassicaceae <- species_df %>%
    filter(Family == "Brassicaceae")

#----Get BioProject tables------------------------------------------------------
pdats <- c(
  "2000:2017/12/31", # 2000-2017
  "2017/12/31:2018/07/31", # 2018.1
  "2018/07/31:2018/12/31", # 2018.2
  "2018/12/31:2019/07/31", # 2019.1
  "2019/07/31:2019/12/31", # 2019.2
  "2019/12/31:2020/07/31", # 2020.1
  "2020/07/31:2020/12/31", # 2020.2
  "2020/12/31:2021/07/31", # 2021.1
  "2021/07/31:2021/12/31", # 2021.2
  "2021/12/31:2022/07/31" # 2022
)


brassicaceae_2017 <- get_bp_table(species_df_brassicaceae$Species, pdats[1])
brassicaceae_2018_1 <- get_bp_table(species_df_brassicaceae$Species, pdats[2])
brassicaceae_2018_2 <- get_bp_table(species_df_brassicaceae$Species, pdats[3])
brassicaceae_2019_1 <- get_bp_table(species_df_brassicaceae$Species, pdats[4])
brassicaceae_2019_2 <- get_bp_table(species_df_brassicaceae$Species, pdats[5])
brassicaceae_2020_1 <- get_bp_table(species_df_brassicaceae$Species, pdats[6])
brassicaceae_2020_2 <- get_bp_table(species_df_brassicaceae$Species, pdats[7])
brassicaceae_2021_1 <- get_bp_table(species_df_brassicaceae$Species, pdats[8])
brassicaceae_2021_2 <- get_bp_table(species_df_brassicaceae$Species, pdats[9])
brassicaceae_2022 <- get_bp_table(species_df_brassicaceae$Species, pdats[10])
ls(pattern = "brassicaceae") # check if all objects exist

brassicaceae_projects <- rbind(
  brassicaceae_2017,
  brassicaceae_2018_1, brassicaceae_2018_2,
  brassicaceae_2019_1, brassicaceae_2019_2,
  brassicaceae_2020_1, brassicaceae_2020_2,
  brassicaceae_2021_1, brassicaceae_2021_2,
  brassicaceae_2022
)
```

## fabaceae_projects.rda

``` r
species_df_fabaceae <- species_df %>%
    filter(Family == "Fabaceae")

#----Get BioProject tables------------------------------------------------------
fabaceae_2017 <- get_bp_table(species_df_fabaceae$Species, pdats[1])
fabaceae_2018_1 <- get_bp_table(species_df_fabaceae$Species, pdats[2])
fabaceae_2018_2 <- get_bp_table(species_df_fabaceae$Species, pdats[3])
fabaceae_2019_1 <- get_bp_table(species_df_fabaceae$Species, pdats[4])
fabaceae_2019_2 <- get_bp_table(species_df_fabaceae$Species, pdats[5])
fabaceae_2020_1 <- get_bp_table(species_df_fabaceae$Species, pdats[6])
fabaceae_2020_2 <- get_bp_table(species_df_fabaceae$Species, pdats[7])
fabaceae_2021_1 <- get_bp_table(species_df_fabaceae$Species, pdats[8])
fabaceae_2021_2 <- get_bp_table(species_df_fabaceae$Species, pdats[9])
fabaceae_2022 <- get_bp_table(species_df_fabaceae$Species, pdats[10])
ls(pattern = "fabaceae") # check if all objects exist

fabaceae_projects <- rbind(
  fabaceae_2017,
  fabaceae_2018_1, fabaceae_2018_2,
  fabaceae_2019_1, fabaceae_2019_2,
  fabaceae_2020_1, fabaceae_2020_2,
  fabaceae_2021_1, fabaceae_2021_2,
  fabaceae_2022
)

save(fabaceae_projects,
     file = here("data", "fabaceae_projects.rda"),
     compress = "xz")
```

## poaceae_projects.rda

``` r
species_df_poaceae <- species_df %>%
    filter(Family == "Poaceae")

#----Get BioProject tables------------------------------------------------------
poaceae_2017 <- get_bp_table(species_df_poaceae$Species, pdats[1])
poaceae_2018_1 <- get_bp_table(species_df_poaceae$Species, pdats[2])
poaceae_2018_2 <- get_bp_table(species_df_poaceae$Species, pdats[3])
poaceae_2019_1 <- get_bp_table(species_df_poaceae$Species, pdats[4])
poaceae_2019_2 <- get_bp_table(species_df_poaceae$Species, pdats[5])
poaceae_2020_1 <- get_bp_table(species_df_poaceae$Species, pdats[6])
poaceae_2020_2 <- get_bp_table(species_df_poaceae$Species, pdats[7])
poaceae_2021_1 <- get_bp_table(species_df_poaceae$Species, pdats[8])
poaceae_2021_2 <- get_bp_table(species_df_poaceae$Species, pdats[9])
poaceae_2022 <- get_bp_table(species_df_poaceae$Species, pdats[10])
ls(pattern = "poaceae_") # check if all objects exist

poaceae_projects <- rbind(
  poaceae_2017,
  poaceae_2018_1, poaceae_2018_2,
  poaceae_2019_1, poaceae_2019_2,
  poaceae_2020_1, poaceae_2020_2,
  poaceae_2021_1, poaceae_2021_2,
  poaceae_2022
)

save(poaceae_projects,
     file = here("data", "poaceae_projects.rda"),
     compress = "xz")
```

## otherfam_projects.rda

``` r
species_df_other <- species_df %>%
    filter(!Family %in% c("Fabaceae", "Brassicaceae", "Poaceae"))

#----Get BioProject tables------------------------------------------------------
other_2017 <- get_bp_table(species_df_other$Species, pdats[1])
other_2018_1 <- get_bp_table(species_df_other$Species, pdats[2])
other_2018_2 <- get_bp_table(species_df_other$Species, pdats[3])
other_2019_1 <- get_bp_table(species_df_other$Species, pdats[4])
other_2019_2 <- get_bp_table(species_df_other$Species, pdats[5])
other_2020_1 <- get_bp_table(species_df_other$Species, pdats[6])
other_2020_2 <- get_bp_table(species_df_other$Species, pdats[7])
other_2021_1 <- get_bp_table(species_df_other$Species, pdats[8])
other_2021_2 <- get_bp_table(species_df_other$Species, pdats[9])
other_2022 <- get_bp_table(species_df_other$Species, pdats[10])
ls(pattern = "other_") # check if all objects exist

otherfam_projects <- rbind(
  other_2017,
  other_2018_1, other_2018_2,
  other_2019_1, other_2019_2,
  other_2020_1, other_2020_2,
  other_2021_1, other_2021_2,
  other_2022
)

save(otherfam_projects,
     file = here("data", "otherfam_projects.rda"),
     compress = "xz")
```
