#######################################################################################
#### This script is used to test R code before adding it to the Rshiny application ####
#######################################################################################


#### 15.3.1 Sample Information Exploration ####
#' Load in data, data downloaded from NCBI 
#'    - Use GEOquery (biocManager) package to load in series matrix data
#'    - Use Tidyverse to load in csv/tsv's
library(GEOquery)
library(tidyverse)

#GSE object, structures lots of information about GEO accession/data given to it
gse=getGEO(filename="data/GSE64810_series_matrix.txt")

#Sample information can be extracted from GSE object 
metadata <- gse@phenoData@data

#Filter the metadata to only relevant information and clean it up
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
    mutate(across(c(-1:-3), as.numeric)) %>%
    
    #get read of leading space in diagnosis column
    mutate(diagnosis = trimws(diagnosis)) %>%
    
    #change diagnosis column to a factor
    mutate_at('diagnosis', factor)
    
#relevel diagnosis column in md_filtered it so control is reference
md_filtered$diagnosis <- relevel(md_filtered$diagnosis, ref='Neurologically normal')
    

#### 13.5.1 OUTPUT TAB 1: SUMMARY TABLE ####
#Take metadata and generate summary table of it including column name | Type | Mean(sd) or distinct values
mean_sd_f <- function(df_col) {
    stringr::str_c(round(mean(df_col, na.rm=TRUE), 2), ' (+/-', round(sd(df_col, na.rm=TRUE), 2), ')')
}

md_summary <- (summarise_all(md_filtered, class))
md_summary <- as.data.frame(t(md_summary[,]))
colnames(md_summary)[1] = 'Type'
md_summary <- rownames_to_column(md_summary, var = "Column Name")

character_vector <- c('sample name', 'geo_accession', 'Neurologically Normal, Huntington\'s Disease')
numeric_vector <- md_filtered %>% 
    summarise(across(where(is.numeric), mean_sd_f))
summary_stats_vector <- c(character_vector, numeric_vector)

md_summary$`Mean(sd) or Description of Values` <- summary_stats_vector



#### Questions for Adam ####
