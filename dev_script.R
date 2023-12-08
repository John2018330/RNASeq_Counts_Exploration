#######################################################################################
#### This script is used to test R code before adding it to the Rshiny application ####
#######################################################################################

library(GEOquery)
library(tidyverse)

#### 13.5.1 Sample Information Exploration ####
#' Load in data, data downloaded from NCBI 
#'    - Use GEOquery (biocManager) package to load in series matrix data
#'    - Use Tidyverse to load in csv/tsv's

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


#### 13.5.1 OUTPUT TAB 3: DENSITY PLOTS
# Case when both HD and non-HD data exists
var <- 'total_reads'
md_dens_plot <- ggplot(md_filtered, aes(x=!!sym(var), fill=diagnosis)) + 
    geom_density(alpha=0.5) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.box.background = element_rect(color = 'black', linewidth = 0.8),
          legend.background = element_rect(fill=alpha('grey', 0.5))) +
    xlab(str_to_title(gsub('_', ' ', var))) +
    ylab('Density') + 
    theme(legend.position="bottom", legend.box = "horizontal") +
    guides(fill = guide_legend(title.position = 'top', title.hjust = 0.5)) + 
    scale_fill_discrete(name='Diagnosis') +
    scale_x_continuous(labels = scales::comma, expand = c(0,0)) +
    scale_y_continuous(expand=c(0,0))
    
# Case when only HD data exists
var <- 'hv_striatal_score'
md_dens_plot <- ggplot(md_filtered, aes(x=!!sym(var))) +
    geom_density(alpha=0.5, fill='#34BFC7') +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    xlab(paste(str_to_title(gsub('_', ' ', var)), 'of Huntington\'s patients/diagnosis')) +
    ylab('Density') +
    #scale_fill_discrete(name='Diagnosis') +
    scale_x_continuous(labels = scales::comma, expand = c(0,0)) +
    scale_y_continuous(expand=c(0,0))



#### 13.5.2 Counts Matrix Exploration INPUT
# input a csv, output a tibble with 2 new columns: percentile variance and num samples that are non zero for the gene
#dplyr::percent_rank()
load_normalized_counts <- function(counts_file) {
    #this needs to be separate because of the full_join in the below pipe
    norm_counts <- read_tsv(counts_file, show_col_types = FALSE)
    
    #pivoting longer and then grouping by the original row is WAAYY faster than doing this rowwise()...
    #use percent_rank() to generate percentile of a column
    norm_counts <- norm_counts %>%
        tidyr::pivot_longer(-GeneID) %>%
        group_by(GeneID) %>%
        summarize(var = var(value),
                  non_zero_count = sum(value != 0)) %>%
        full_join(y=norm_counts, by=join_by(GeneID)) %>%
        mutate(var_percentile = percent_rank(var), .after = var)
        
    return(norm_counts)
}

#tpm <- load_normalized_counts('data/GSE64810_norm_counts_TPM.tsv')
#### Questions for Adam ####
