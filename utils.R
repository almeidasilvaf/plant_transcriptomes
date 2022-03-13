
#' Count number of available RNA-seq on SRA for a given species
#'
#' @param species Character of species name.
#' 
#' @return Numeric scalar of number of samples
#' @importFrom rvest read_html html_nodes html_text
#' @examples 
#' species <- "Brassica napus"
#' count_available_samples(species)
count_available_samples <- function(species = NULL) {
    
    link <- paste0(
        "https://www.ncbi.nlm.nih.gov/sra/?term=%22", 
        gsub(" ", "+", species),
        "%22%5BORGN%5D+AND+%22RNA-seq%22%5BSTRA%5D"
    )
    page <- rvest::read_html(link)
    title <- rvest::html_nodes(page, ".title_and_pager") 
    number <- rvest::html_nodes(title, "h3.result_count")
    number <- rvest::html_text(number)
    number <- strsplit(number, " ")
    count <- vapply(number, function(x) { tail(x, 1) }, character(1))
    count <- as.numeric(count)
    if(length(count) == 0) { count <- 0 }
    return(count)
}


#' Create a BioProject-level summary table
#' 
#' @param metadata Data frame of sample metadata.
#' 
#' @return A data frame with the variables:
#' \itemize{
#'   \item BioProject
#'   \item N
#'   \item Tissue
#'   \item Study title
#'   \item Study abstract
#'   \item Pubmed
#' }
#' @importFrom dplyr add_count select rename distinct group_by filter 
#' summarise arrange
#' @importFrom stringr str_c
#' @noRd
create_project_table <- function(metadata = NULL) {
    table <- metadata %>%
        dplyr::filter(startsWith(BioProject, "PRJ")) %>%
        dplyr::add_count(BioProject) %>%
        dplyr::select(BioProject, n, Study_title, Study_abstract, Pubmed) %>%
        dplyr::rename(
            N = n, 
            `Study title` = Study_title,
            `Study abstract` = Study_abstract
        ) %>%
        dplyr::distinct()
    
    tissue_count <- metadata %>%
        dplyr::filter(startsWith(BioProject, "PRJ")) %>%
        group_by(BioProject, Tissue) %>%
        summarise(n = n()) %>%
        ungroup() %>%
        arrange(-n) %>%
        group_by(BioProject) %>%
        summarise(tissue_count = stringr::str_c(Tissue, ": ", n, 
                                                collapse = " | "))
    
    final_table <- dplyr::inner_join(table, tissue_count, 
                                     by = "BioProject") %>%
        dplyr::rename(Tissue = tissue_count) %>%
        dplyr::select(BioProject, N, Tissue, 
                      `Study title`, `Study abstract`, Pubmed)
    return(final_table)
}


get_bp_table <- function(species = NULL, pdat = NULL) {
    bp_table <- Reduce(rbind, lapply(species, function(x) {
        terms <- paste0(x, "[ORGN] AND RNA-seq[STRA] AND ", pdat, "[PDAT]")
        message("Working on species ", x)
        metadata <- bears::create_sample_info(terms, retmax = retmax)
        project_table <- NULL
        if(!is.null(metadata)) {
            project_table <- create_project_table(metadata)
            project_table$Species <- x
        }
        return(project_table)
    }))
    return(bp_table)
}
