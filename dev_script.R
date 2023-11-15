#######################################################################################
#### This script is used to test R code before adding it to the Rshiny application ####
#######################################################################################


#### Step 1 ####
#' Load in data, data downloaded from NCBI 
#'    - Use GEOquery (biocManager) package to load in series matrix data
#'    - Use Tidyverse to load in csv/tsv's
library(GEOquery)
library(tidyverse)

#GSE object, structures lots of information about GEO accession/data given to it
gse=getGEO(filename="data/GSE64810_series_matrix.txt")

#Sample information can be extracted from GSE object and then filtered
metadata <- gse@phenoData@data
md_filtered <- metadata %>%
    #select the relevant columns 
    dplyr::select(title, geo_accession, dplyr::starts_with('characteristics_ch1.')) %>%
    
    #rename column names from generic to specific
    rename(diagnosis         = characteristics_ch1.1, pmi               = characteristics_ch1.2,
           age_of_death      = characteristics_ch1.3, rin               = characteristics_ch1.4,
           total_reads       = characteristics_ch1.5, age_of_onset      = characteristics_ch1.6,
           duration          = characteristics_ch1.7, cag               = characteristics_ch1.8,
           vonsattel_grade   = characteristics_ch1.9, hv_striatal_score = characteristics_ch1.10,
           hv_cortical_score = characteristics_ch1.11) %>%
    
    #change column entries from "characteristic_name:characteristic_value" to "characteristic_value".
    #case_when to ignore the empty values in non Huntingtons patients rows for huntington's specific
    #     columns
    mutate(across(c(-1, -2), ~ case_when(. == "" ~ .,
                                         . != "" ~ str_split_i(., ":", -1)))) %>%
    
    #change columns that should be numeric to numeric, adds NA's where necessary (will give warning)
    mutate(across(c(-1:-3), as.numeric))
    



#### Questions for Adam ####
#' 13.5.3 Shiny Functionalities: Tab with content similar to that described in [Assignment 7] ?