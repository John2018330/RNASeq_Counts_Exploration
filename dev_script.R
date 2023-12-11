#######################################################################################
#### This script is used to test R code before adding it to the Rshiny application ####
#######################################################################################

library(GEOquery)
library(tidyverse)
library(gplots)
library(DESeq2)

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
    dplyr::rename(diagnosis         = characteristics_ch1.1, pmi               = characteristics_ch1.2,
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
        mutate(var_percentile = 100 * percent_rank(var), .after = var)
        
    return(norm_counts)
}

#tpm <- load_normalized_counts('data/GSE64810_norm_counts_TPM.tsv')


#### 13.5.2 Norm Counts Filter Summary
filter_norm_counts <- function(norm_counts, percentile, min_non_zeroes) {
    filtered <- norm_counts %>%
        dplyr::filter(var_percentile >= percentile & non_zero_count >= min_non_zeroes)
    return (filtered)
}

#filtered_counts <- filter_norm_counts(tpm, 95, 60)


#### 13.5.2 Norm Counts OUTPUT 2 Scatter plots w/ filter
# VARIANCE FILTER: Take normal counts tibble, filter for percentile variance & minimum samples w/ nonzero counts, and generate a scatter plot
filtered_norm_counts_variance_scatter <- function(norm_counts, percentile_filter) {
    # Make copy of original tibble w/ 0 variance filtered out
    norm_counts_wVariance <- norm_counts %>%
        filter(var != 0)
    
    # Take tibble and wrangle to add necessary columns
    variance_tibble <- norm_counts_wVariance %>%
        
        # Again perform get row-wise median count without using rowwise() because it's super slow
        tidyr::pivot_longer(cols = starts_with('GSM')) %>%
        group_by(GeneID) %>%
        summarize(median_count = median(value)) %>%
        full_join(y=norm_counts_wVariance, by=join_by(GeneID)) %>%
        
        # Ease for troubleshooting
        select(!starts_with('GSM')) %>% 
        
        # SUSPECT LINE!! Some genes have VERY HIGH variance but a 0 *MEDIAN* count; check out gene 28476
        filter(median_count != 0) %>%
        
        # Rank median (since it's also pretty skewed) and add column to color points by whether they pass the filtering threshold
        mutate(rank_median = min_rank(median_count),
               .after = median_count) %>%
        mutate(variance_filter_status = case_when(var_percentile >= percentile_filter ~ 'TRUE', 
                                                  .default = 'FALSE'),
               .after = var_percentile)
    
    # Color vector for plot
    color_map <- c('TRUE' = '#FFC107', 'FALSE' = '#004D40')
    
    # Make scatter plot with y axis log scale
    variance_scatter <- variance_tibble %>%
            ggplot(aes(x=rank_median, y=var, color=variance_filter_status)) +
            geom_point() + 
            theme_classic() +     
            xlab('Median (ranked)') +
            ylab('Variance on a Log10 Scale') + 
            ggtitle('Ranked median vs Variance on a Log10 Scale') + 
            
            # Change axis tick labels and scales
            scale_x_continuous(labels = scales::comma) +
            scale_y_continuous(trans = scales::log_trans(10), labels=scales::label_log(10)) + 
            
            # Legend and color and general text editing
            theme(plot.title = element_text(size=20, face='bold'), 
                  axis.title.x = element_text(size=14), axis.title.y = element_text(size=14),
                  legend.title = element_text(size=14), legend.text = element_text(size=14),
                  legend.box.background = element_rect(color = 'black', linewidth = 0.8),
                  legend.background = element_rect(fill=alpha('grey', 0.5)), 
                  legend.position="bottom", legend.box = "horizontal") +
            guides(color = guide_legend(title = paste('Percentile Variance > ', percentile_filter, '%', sep=''), 
                                        title.position = 'top', title.hjust = 0.5)) + 
            scale_color_manual(values=color_map)
    
    return (variance_scatter)
}

#filtered_norm_counts_variance_scatter(tpm, 99.9)
#test_scatter <- filtered_norm_counts_variance_scatter(tpm, 99.9)


# NONZERO FILTER
filtered_norm_counts_zero_scatter <- function(norm_counts, zero_filter) {
    # Used to inverse the non-zero count
    num_samples <- sum(startsWith(colnames(norm_counts), 'GSM'))
    
    # Add columns to tibble with required information for scatter plot
    zero_tibble <- norm_counts %>%
        
        # Again perform get row-wise median count without using rowwise() because it's super slow
        tidyr::pivot_longer(cols = starts_with('GSM')) %>%
        group_by(GeneID) %>%
        summarize(median_count = median(value)) %>%
        full_join(y=norm_counts, by=join_by(GeneID)) %>%
        
        # Plotting number of zeroes, not non-zeroes so inverse the count and add it as a column
        mutate(zero_count = num_samples - non_zero_count, .after = non_zero_count) %>%
        select(!starts_with('GSM')) %>%
    
        # Add status column for plotting if passes filter threshold and rank median count
        mutate(rank_median = dense_rank(median_count),
               .after = median_count) %>%
        mutate(zero_filter_status = case_when(non_zero_count >= zero_filter ~ 'TRUE', 
                                              .default = 'FALSE'),
           .after = zero_count)
    
    
    # Color vector for plot
    color_map <- c('TRUE' = '#FFC107', 'FALSE' = '#004D40')
    
    # Make scatter plot
    zero_scatter <- zero_tibble %>%
        ggplot(aes(x=rank_median, y=zero_count, color=zero_filter_status)) +
            geom_point() + 
            theme_classic() +     
            xlab('Median Count (ranked)') +
            ylab('Zero Count') +
            ggtitle('Ranked median vs. # of zero counts for a gene') +
            
            # Change axis tick labels and scales
            scale_x_continuous(labels = scales::comma) +

            # Legend and color editing
            theme(plot.title = element_text(size=20, face='bold'), 
                axis.title.x = element_text(size=14), axis.title.y = element_text(size=14),
                legend.title = element_text(size=14), legend.text = element_text(size=14),
                legend.box.background = element_rect(color = 'black', linewidth = 0.8),
                legend.background = element_rect(fill=alpha('grey', 0.5)),
                legend.position = "bottom", legend.box = "horizontal") +
            guides(color = guide_legend(title = paste('At least ', zero_filter, ' samples have non-zero counts', sep=''),
                                        title.position = 'top', title.hjust = 0.5)) +
            scale_color_manual(values = color_map)
    
    return (zero_scatter)
    
}

#filtered_norm_counts_zero_scatter(tpm, 60)
#test_zero_scatter <- filtered_norm_counts_zero_scatter(tpm, 60)


#### 13.5.2 Norm Counts OUTPUT 3 Heatmap
filtered_counts <- filter_norm_counts(tpm, 95, 60)

plot_heatmap <- function(filtered_counts) {
    
    filtered_matrix <- filtered_counts %>% 
        dplyr::select(dplyr::starts_with('GSM'))# %>%
        #mutate_all(log10) %>%
        #tidyr::drop_na() %>%
        #dplyr::filter_all(dplyr::all_vars(is.finite(.)))
    
    hm <- heatmap.2(t(filtered_matrix), trace='none', Colv = FALSE, 
                    dendrogram = 'row', scale='column', na.rm=TRUE)
    
    return (hm)
}

heatmap <- plot_heatmap(filtered_counts)



#### 13.5.2 PCA 
make_pca_plot <- function(filtered_metadata) {
    # Manual load complete dataset for PCA, remove zero variance genes
    intensity <- read.csv('data/GSE64810_norm_counts_TPM.tsv', sep='\t')
    intensity <- intensity[-c(1)]
    nonzero_rowvector <- rowVars(as.matrix(intensity[-1]))!=0 
    intensity <- intensity[nonzero_rowvector,]
    
    # Run PCA
    pca_results <- prcomp(scale(t(intensity)), center=FALSE, scale=FALSE)
    
    # Select metadata
    md_for_pca <- filtered_metadata %>%
        dplyr::select('geo_accession', 'diagnosis')
    
    # Make tibble with diagnosis
    pca_tibble <- as_tibble(pca_results$x, rownames=NA) %>%
        rownames_to_column('geo_accession') %>%
        left_join(y=md_for_pca, by=('geo_accession' = 'geo_accession')) %>%
        relocate(diagnosis, .after=geo_accession)
    
    # Get variance values for plot
    var <- pca_results$sdev**2
    prop_of_var <- (var/sum(var))
    prop_of_var <- signif(prop_of_var, digits=2) * 100
    
    # Color map for plot
    color_map <- c('Huntington\'s Disease' = '#FFC107', 'Neurologically normal' = '#004D40')
    
    # Construct plot
    plot <- ggplot(pca_tibble, aes(x=PC1,y=PC2, color = diagnosis)) +
        # Labels and themes
        geom_point(size=4) +
        theme_classic() +
        ggtitle('PCA of Count Data') +
        xlab(str_c('PC1: ', prop_of_var[1], '% variance')) +
        ylab(str_c('PC2: ', prop_of_var[2], '% variance')) +
        
        # Legend and font and colors
        theme(plot.title = element_text(size=20, face='bold'), 
              axis.title.x = element_text(size=14), axis.title.y = element_text(size=14),
              legend.title = element_text(size=14), legend.text = element_text(size=14),
              legend.box.background = element_rect(color = 'black', linewidth = 0.8),
              legend.background = element_rect(fill=alpha('grey', 0.5)),
              legend.position = "bottom", legend.box = "horizontal") +
        guides(color = guide_legend(title = 'Diagnosis',
                                    title.position = 'top', title.hjust = 0.5)) +
        scale_color_manual(values = color_map)
    
    return(plot)
    
}
make_pca_plot(md_filtered)
#pca_plot <- make_pca_plot(md_filtered)

