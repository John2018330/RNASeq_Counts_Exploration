################################
#'### RShiny Final Project #####
#'##############################
#'    BF591 Fall 2023 
#'    Johnathan Zhang
#'    Data from GEO accession: GSE64810
#'       - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810



#### Packages and Options ####
library(shiny)
library(bslib)
library(tidyverse)
library(colourpicker)
library(GEOquery)
library(DT)
library(gplots)
library(DESeq2)
library(fgsea)

# File upload limit
options(shiny.maxRequestSize=30*1024^2)

#### UI ####
ui <- fluidPage(
    titlePanel("BF591 Final Project"),
    h3('A Rshiny website to explore RNA-seq counts'),
    p(strong("Made by Johnathan Zhang")),
    
    tabsetPanel(
        
        #### 13.5.1 Sample Information Exploration ####
        tabPanel('Sample Metadata', 
                 sidebarLayout(
                     ### Upload metadata file
                     sidebarPanel(
                         fileInput('sample_info_csv', 
                                   paste0('Upload file: Sample information matrix in CSV format'))),
                     
                     ### Outputs
                     mainPanel(tabsetPanel(
                         ## Table containing a summary of metadata file
                         tabPanel('Summary', 'Summary of table that includes summary of 
                                  type and values in each column',
                                  DTOutput('metadata_summary')), 
                         
                         ## Full table of filtered metadata
                         tabPanel('Table', 'Tab with a data table displaying the sample information, 
                                  with sortable columns', 
                                  DTOutput('metadata')),
                         
                         ## Customizable density plots of variables in metadata
                         tabPanel('Plots', 
                                  sidebarLayout(
                                      sidebarPanel(
                                          radioButtons("metadata_plot_var", #ID
                                                       "Choose a variable to generate a density plot", 
                                                       c('pmi', 'age_of_death', 'rin', 'total_reads',
                                                         'age_of_onset', 'duration', 'cag', 
                                                         'vonsattel_grade', 'hv_striatal_score',
                                                         'hv_cortical_score'), 
                                                       selected='pmi'),
                                          actionButton('make_density_plot', label='Plot', 
                                                       style='width:100%')
                                          #style='width:100%; background: #5ECCAB'
                                      ),
                                      mainPanel(
                                          plotOutput(outputId = "metadata_density_plot", height='600px') #height='600px' 
                                      )))
                     ))
                 )),
        
        #### 13.5.2 Normalized Counts Matrix Exploration ####
        tabPanel('Counts Exploration', 
                 sidebarLayout(
                     ### Upload norm_counts_csv and filter sliders
                     sidebarPanel(
                         ## Upload the csv
                         fileInput('norm_counts_csv', 
                                   paste0('Normalized counts matrix, in CSV format, input sliders 1 and 2')),
                         
                         ## Slider for variance percentile (0-100%)
                         sliderInput("norm_counts_percentile_var", min = 0, max = 100, value = 99,  step = 0.1,
                                     label = "Select the percentile of genes to include based on variance across samples 
                                     (e.g., 90% means select genes with the top 10% variance values)"),
                         
                         ## Slider for number of non-zero counts
                         sliderInput('norm_counts_nonzeroes', min=0, max=69, value = 60, step = 1,
                                     label = 'Select how many samples must have a non-zero count value for the gene')),
                     
                     ### Output tabs: Summary, Plots, Heatmap, and PCA
                     mainPanel(tabsetPanel(
                         tabPanel('Summary', 
                                  'Tab with a table summarizing the effect of the filtering',
                                  DTOutput('filtered_norm_counts_table')),
                         
                         ## Scatter plot 
                         tabPanel('Scatter Plots', 
                                  'Tab with diagnostic scatter plots, where genes passing filters are darker, 
                                  and genes filtered out are lighter',
                                  fluidRow(splitLayout(
                                      cellWidths = c('50%', '50%'), 
                                      plotOutput(outputId = 'counts_scatter_variation', height = '500px'),
                                      plotOutput(outputId = 'counts_scatter_nonzeros', height = '500px')),
                                  )),
                         
                         ## Heatmap
                         tabPanel('Heatmap', 
                                  h3('Heatmap of filtered counts'),
                                  plotOutput(outputId = 'counts_heatmap', height='600px')),
                         
                         ## PCA
                         tabPanel('PCA', 
                                  sidebarLayout(
                                      sidebarPanel(
                                          radioButtons("pca_x", #ID
                                                       "Choose a Principal Component to plot on the x axis", 
                                                       c('PC1', 'PC2', 'PC3', 'PC4', 'PC5'), 
                                                       selected = 'PC1'),
                                          radioButtons('pca_y',
                                                       'Choose a Principal Component to plot on the y axis',
                                                       c('PC1', 'PC2', 'PC3', 'PC4', 'PC5'), 
                                                       selected = 'PC2'),
                                          actionButton('plot_pca', label='Plot', 
                                                       style='width:100%')),
                                      
                                       mainPanel('Tab with PCA plot of select principal components on entire dataset 
                                                 minus zero-variance genes. Please upload BOTH metadata/sample_info 
                                                 file (tab 1) and normalized counts',
                                                 plotOutput(outputId = 'counts_pca', height='600px'))
                            ))
                     ))
                 )),
        
        #### 13.5.3 Differential Expression ####
        tabPanel('Differential Expression', 
                 sidebarLayout(
                     ## Upload the DESeq2 results
                     sidebarPanel(
                         fileInput('deseq_csv', 
                                   paste0('Upload results of a differential expression analysis in CSV format'))
                     ),
                     
                     mainPanel(tabsetPanel(
                         ## Sortable table of DESeq results
                         tabPanel('Table', 
                                  'Tab with sortable table displaying differential expression results',
                                  DTOutput('deseq_table')),
                         
                         ## Customizable volcano plot of DE genes
                         tabPanel('Volcano Plot', 
                                  sidebarLayout(
                                      sidebarPanel(
                                          HTML(paste("A volcano plot can be generated with", "<b>", "\"log2 fold-change\"", 
                                                     "</b>", "on the x-axis and", "<b>", "\"p-adjusted\"", "</b>", 
                                                     "on the y-axis.", "<br/><br/>")),
                                          
                                          radioButtons("volcano_x", "Choose the column for the x-axis", 
                                                       c('baseMean', 'HD.mean', 'Control.mean', 'log2FoldChange', 
                                                         'lfcSE', 'stat', 'pvalue', 'padj'), 
                                                       selected='log2FoldChange'),
                                          radioButtons("volcano_y", "Choose the column for the y-axis", 
                                                       c('baseMean', 'HD.mean', 'Control.mean', 'log2FoldChange', 
                                                         'lfcSE', 'stat', 'pvalue', 'padj'),
                                                       selected='padj'),
                                          
                                          colourInput("volcano_base_color", "Base Point color", "#004D40"),
                                          colourInput("volcano_highlight_color", "Highlight Point color", "#FFCF56"),
                                          
                                          # tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar
                                          #                 {background: #5ECCAB}")),
                                          
                                          sliderInput("padj_threshold", min = -40, max = 0, 
                                                      label = "Select the magnitude of the p adjusted coloring:", 
                                                      value = -25, step = 1),
                            
                                          actionButton('make_volcano_plot', label='Plot', 
                                                       style='width:100%; background: #5ECCAB')
                                        ),
                                      
                                      mainPanel(plotOutput(outputId = "volcano", height='600px'),
                                      )
                                  ))
                        ))
                 )),
        
        #### 13.6.1 GSEA using fgsea
        tabPanel('GSEA',
                sidebarLayout(
                    
                    ### INPUT: Upload gmt file (and deseq file in previous tab)
                    sidebarPanel('Please upload gmt file here and DESeq results in previous tab',
                        fileInput('gmt_file', paste0('Upload gmt file here')),
                    ),
                    
                    ### OUTPUTS: NES Barplot, Data table, and scatterplot
                    mainPanel(tabsetPanel(
                        
                        ## Barplot of significant NES pathways
                        tabPanel('NES Barplot',
                                 sidebarLayout(
                                     sidebarPanel(
                                         sliderInput("fgsea_padj_threshold", min = -30, max = 0, 
                                                     label = "Select a threshold for padj for the barplot", 
                                                     value = -15, step = 1)
                                     ),
                                     
                                     mainPanel(
                                         plotOutput(outputId = 'fgsea_barplot', height='600px'))
                        )),
                        
                        ## Table of FGSEA results
                        tabPanel('FGSEA results',
                                 sidebarLayout(
                                     sidebarPanel(
                                         sliderInput("fgsea_padj_threshold_2", min = -30, max = 0, 
                                                     label = "Select a threshold for padj for the data table", 
                                                     value = -15, step = 1),
                                         
                                         radioButtons("nes_pathway", "Select positive, negative, or all NES pathways", 
                                                      c('Positive', 'Negative', 'Both'),
                                                      selected='Both'),
                                         
                                         downloadButton('download_fgsea', 'Download the data table to the right')
                                     ),
                                     
                                     mainPanel(
                                         DTOutput('fgsea_table')
                                     )
                        
                        )),
                        
                        ## Scatter plot of pathways and significance 
                        tabPanel('Pathway Scatter Plot',
                                 sidebarLayout(
                                     sidebarPanel(
                                         sliderInput("fgsea_padj_threshold_3", min = -30, max = 0, 
                                                     label = "Select a threshold for padj for the scatter plot", 
                                                     value = -15, step = 1)
                                     ),
                                     
                                     mainPanel(
                                         plotOutput(outputId = 'fgsea_scatter', height='600px')
                                     )
                        ))
                    ))
                ))
    )
)

#### Server side / Back end ####
server <- function(input, output, session) {

    #######################################################################
    #' 13.5.1 SAMPLE INFO EXPLORE INPUT FILE
    #' Function to take in file uploaded in Sample Information Exploration, 
    #' data wrangle it and return it as tibble
    
    load_sample_information <- reactive({
        # Don't run until file has been uploaded
        req(input$sample_info_csv)
        
        # GSE object, structures lots of information about GEO accession/data given to it
        #gse=getGEO(filename="data/GSE64810_series_matrix.txt")
        gse=getGEO(filename=input$sample_info_csv$datapath)
        
        # Sample information can be extracted from GSE object 
        metadata <- gse@phenoData@data
        
        # Filter the metadata to only relevant information and clean it up
        md_filtered <- metadata %>%
            # Select the relevant columns 
            dplyr::select(title, geo_accession, dplyr::starts_with('characteristics_ch1.')) %>%
            
            # Rename column names from generic to specific
            dplyr::rename(diagnosis         = characteristics_ch1.1, pmi               = characteristics_ch1.2,
                   age_of_death      = characteristics_ch1.3, rin               = characteristics_ch1.4,
                   total_reads       = characteristics_ch1.5, age_of_onset      = characteristics_ch1.6,
                   duration          = characteristics_ch1.7, cag               = characteristics_ch1.8,
                   vonsattel_grade   = characteristics_ch1.9, hv_striatal_score = characteristics_ch1.10,
                   hv_cortical_score = characteristics_ch1.11) %>%
            
            # Change column entries from "characteristic_name:characteristic_value" 
            # to "characteristic_value".
            # Case_when to ignore the empty values in non Huntingtons patients 
            # rows for huntington's specific columns
            dplyr::mutate(across(c(-1, -2), ~ case_when(. == "" ~ .,
                                                 . != "" ~ str_split_i(., ":", -1)))) %>%
            
            # Change columns that should be numeric to numeric, adds NA's where 
            # necessary (will give warning)
            dplyr::mutate(across(c(-1:-3), as.numeric)) %>%
        
            # Get read of leading space in diagnosis column
            dplyr::mutate(diagnosis = trimws(diagnosis)) %>%
            
            # Change diagnosis column to a factor
            dplyr::mutate_at('diagnosis', factor)
        
        # Relevel diagnosis column in md_filtered it so control is reference
        md_filtered$diagnosis <- relevel(md_filtered$diagnosis, ref='Neurologically normal')
        
        return (md_filtered)
    })
    
    #######################################################
    #' 13.5.1 SAMPLE INFO EXPLORE OUTPUT TAB 1 
    #' Render Data Table summary information about metadata
    
    # A helper function to find mean and SD of numeric columns and join them into a string
    mean_sd_f <- function(df_col) {
        stringr::str_c(round(mean(df_col, na.rm=TRUE), 2), ' (+/-', round(sd(df_col, na.rm=TRUE), 2), ')')
    }
    
    # Main function that takes in metadata tibble and summarizes it
    summarize_metadata <- function() {
        # Don't run until file has been uploaded
        #req(input$sample_info_csv)
        
        # Load in the metadata
        md_filtered <- load_sample_information()
        
        # Make first two columns of table: Col 1: Column Name), Col 2: Data type
        md_summary <- (summarise_all(md_filtered, class))                   # Class of column
        md_summary <- as.data.frame(t(md_summary[,]))                       # Transpose
        colnames(md_summary)[1] = 'Type'                                    # Rename 2nd column
        md_summary <- rownames_to_column(md_summary, var = "Column Name")   # Rename 1st column
        
        # Generate vectors for character/factor columns (hard coded) and then vector for numeric columns
        # using the helper function mean_sd_f. Combine the 2 and add to md_summary tibble
        character_vector <- c('sample name', 'geo_accession', 
                              'Neurologically Normal, Huntington\'s Disease')
        numeric_vector <- md_filtered %>% summarise(across(where(is.numeric), mean_sd_f))
        summary_stats_vector <- c(character_vector, numeric_vector)
        md_summary$`Mean(sd) or Description of Values` <- summary_stats_vector
        
        return (md_summary)
    }
    
    # Render the Data Table
    output$metadata_summary <- renderDT(
        summarize_metadata(), options = list(
            dom='ltp', # What to show in the data table (look up docs)
            pageLength = 20 # Default number of rows to show
        )
    )
    
    
    ##############################################
    #' 13.5.1 SAMPLE INFO EXPLORE OUTPUT TAB 2 
    #' Render the Data Table with sortable columns 
    
    output$metadata <- renderDT(
        load_sample_information(), options = list(
            pageLength = 50, # Default number of rows to show
            lengthMenu = c(5, 10, 20, 50, 100),  # Options for number of rows to show
            dom = 'ltp') # What to show in the data table (look up docs)
    )
    
    
    ###########################################################
    #' 13.5.1 SAMPLE INFO EXPLORE OUTPUT TAB 3
    #' Generate density plot of selected variable from metadata
    
    # Function to generate the plot, split into two options based on if variable is Huntington's specific 
    md_density_plot <- function(dataf, var) {
        # Variables that aren't HD specific
        non_hd_vars <- c('pmi', 'age_of_death', 'rin', 'total_reads')
        if (var %in% non_hd_vars) {
            md_filtered <- dataf
            
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
        
        # Other option is variable is HD specific (so no legend)
        } else {
            md_filtered <- dataf %>%
                drop_na(var)
            
            md_dens_plot <- ggplot(md_filtered, aes(x=!!sym(var))) +
                geom_density(alpha=0.5, fill='#34BFC7') +
                theme_classic() +
                theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank()) +
                xlab(paste(str_to_title(gsub('_', ' ', var)), 'of Huntington\'s patients/diagnosis')) +
                ylab('Density') +
                scale_x_continuous(labels = scales::comma, expand = c(0,0)) +
                scale_y_continuous(expand=c(0,0))
        }
        
        return (md_dens_plot)
        
    }
    
    # Render plot based on user input
    output$metadata_density_plot <- renderPlot ({
        # Require plot button to be pushed
        input$make_density_plot
        
        # Generate plot w/n isolate to require action button to be pushed
        isolate({
            md_density_plot(load_sample_information(), input$metadata_plot_var)
        })
    })
    
    
    #####################################################################################################
    #' 13.5.2 NORM COUNTS INPUT CSV
    #' Input a raw CSV and add columns that will be used for filtering/plotting and return it as a tibble
    #' - percentile variance: a percentile rank of how much variance a gene has 
    #' - non_zero_count:      a count of how many genes had non-zero counts
    
    load_normalized_counts <- reactive({
        # Don't run until file has been uploaded
        req(input$norm_counts_csv)
        counts_file <- input$norm_counts_csv$datapath
        
        # This needs to be separate because of the full_join in the below pipe
        norm_counts <- read_tsv(counts_file, show_col_types = FALSE)
        
        # Pivoting longer and then grouping by the original row is WAAYY faster than doing this rowwise()...
        # Use percent_rank() to generate percentile of a column
        norm_counts <- norm_counts %>%
            tidyr::pivot_longer(-GeneID) %>%
            group_by(GeneID) %>%
            summarize(var = var(value),
                      non_zero_count = sum(value != 0)) %>%
            full_join(y=norm_counts, by=join_by(GeneID)) %>%
            mutate(var_percentile = (100 * percent_rank(var)), .after = var)
        
        return(norm_counts)
    })
    
    
    #############################################################################################
    #' 13.5.2 NORM COUNTS OUTPUT 1 FILTER SUMMARIZATION
    #' Generate a table summarizing the effect of the filtering the normalized counts, including:
    #' number of samples, 
    #' total number of genes, 
    #' number and % of genes passing current filter, 
    #' number and % of genes not passing current filter
    
    #First make a function to filter raw normalized counts by variation percentile and number of non-zero count genes
    filter_norm_counts <- function(norm_counts, percentile, min_non_zeroes) {
        filtered <- norm_counts %>%
            dplyr::filter(var_percentile >= percentile & non_zero_count >= min_non_zeroes)
        return (filtered)
    }

    
    # Function to first filter the raw normalized counts by variation percentile and number of non-zero count genes,
    # then summarize the effects of filtering
    summarize_filtered_norm_counts <- function(norm_counts, percentile, min_non_zeroes) {
        # Generate filtered dataset
        filtered_norm_counts <- filter_norm_counts(norm_counts, percentile, min_non_zeroes)
        
        # Original dimensions of the dataset
        num_samples <- sum(startsWith(colnames(norm_counts), 'GSM'))
        num_genes   <- dim(norm_counts)[1]
        
        # Numbers based on filtered dataset
        num_filtered_genes         <- dim(filtered_norm_counts)[1]
        num_filtered_out_genes     <- num_genes - num_filtered_genes
        percent_filtered_genes     <- round((100 * num_filtered_genes / num_genes), 2)
        percent_filtered_out_genes <- 100 - percent_filtered_genes
        
        # Write filtered gene stats as 'Number (Percent%)'
        filtered_stat     <- paste(num_filtered_genes, '(', percent_filtered_genes, '%)', sep='')
        filtered_out_stat <- paste(num_filtered_out_genes, '(', percent_filtered_out_genes, '%)', sep='')
        
        # Generate columns of summary table
        criteria <- c('Number of samples', 'Total number of genes', 
                      'Number (%) of genes passing current filter', 'Number (%) of genes filtered out')
        summary  <- c(num_samples, num_genes, filtered_stat, filtered_out_stat)
        
        # Combine columns into a dataframe
        summary_table <- data.frame(criteria, summary)
        names(summary_table) <- c('Stat', 'Value')
        
        return (summary_table)
    }
    
    
    # Render the data table
    output$filtered_norm_counts_table <- renderDT(
            summarize_filtered_norm_counts(load_normalized_counts(), 
                                           input$norm_counts_percentile_var, 
                                           input$norm_counts_nonzeroes), 
            options = list(dom = 't')
    )
    
    #################################################
    #' 13.5.2 NORM COUNTS OUTPUT 2 SCATTER PLOTS
    #' Generate scatter plots of median count vs 
    #'  (1) variance and 
    #'  (2) number of zeroes
    #'  w/ filtering/coloring based on input sliders 
    
    #### Function to generate scatter plot of median count vs variance
    filtered_norm_counts_variance_scatter <- function(norm_counts, percentile_filter) {
        # Make copy of original tibble w/ 0 variance filtered out
        norm_counts_wVariance <- norm_counts %>%
            dplyr::filter(var != 0)
        
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
            dplyr::filter(median_count != 0) %>%
            
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
            xlab('Median count (ranked)') +
            ylab('Variance on a Log10 Scale') + 
            ggtitle('Ranked median count vs Variance on a Log10 Scale') + 
            
            # Change axis tick labels and scales
            scale_x_continuous(labels = scales::comma) +
            scale_y_continuous(trans = scales::log_trans(10), labels=scales::label_log(10)) + 
            
            # Legend and color and general text editing
            theme(plot.title = element_text(size=20, face='bold'), 
                  axis.title.x = element_text(size=14), axis.title.y = element_text(size=14),
                  legend.title = element_text(size=14), legend.text = element_text(size=14),
                  legend.box.background = element_rect(color = 'black', linewidth = 0.8),
                  legend.background = element_rect(fill=alpha('grey', 0.5))) +
            theme(legend.position="bottom", legend.box = "horizontal") +
            guides(color = guide_legend(title = paste('Percentile Variance > ', percentile_filter, '%', sep=''), 
                                        title.position = 'top', title.hjust = 0.5)) + 
            scale_color_manual(values=color_map)
        
        return (variance_scatter)
    }
    
    
    #### Function to generate scatter plot of median count vs zero counts for a gene
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
            ggtitle('Ranked median count vs. # of zero counts for a gene') +
            
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
    
    
    # Render filtered normalized counts variation scatter plot
    output$counts_scatter_variation <- renderPlot({
        
        # Require csv to be uploaded first
        input$norm_counts_csv
        
        # Generate plot to render
        filtered_norm_counts_variance_scatter(load_normalized_counts(),
                                              input$norm_counts_percentile_var)
    })
    
    # Render filtered normalized counts zero-count scatter plot
    output$counts_scatter_nonzeros <- renderPlot({
        
        # Require csv to be uploaded first
        input$norm_counts_csv
        
        # Generate plot
        filtered_norm_counts_zero_scatter(load_normalized_counts(),
                                          input$norm_counts_nonzeroes)
    })
    
    
    ###################################################################
    #' 13.5.2 NORM COUNTS OUTPUT 3 HEATMAP
    #' Tab with a clustered heatmap of counts remaining after filtering
    #' Consider enabling log-transforming counts for visualization
    
    # Function to generate heatmap from filtered counts data, 
    # but I feel heatmap's don't really make sense in this context
    plot_heatmap <- function(filtered_counts) {
        
        # Make counts tibble log transofmred and drop introduced NA's/Infinites
        filtered_matrix <- filtered_counts %>% 
            dplyr::select(dplyr::starts_with('GSM')) #%>%
            #mutate_all(log10) %>%
            #tidyr::drop_na() %>%
            #dplyr::filter_all(dplyr::all_vars(is.finite(.)))
        
        # Make heatmap plot
        hm <- heatmap.2(t(filtered_matrix), trace='none', Colv = FALSE, 
                        dendrogram = 'row', scale='column', na.rm=TRUE)
        
        return (hm)
    }
    
    # Render the heatmap
    output$counts_heatmap <- renderPlot ({
        
        # Require csv to be uploaded first
        input$norm_counts_csv
        
        filtered_counts <- (filter_norm_counts(load_normalized_counts(), 
                                               input$norm_counts_percentile_var, 
                                               input$norm_counts_nonzeroes))
        
        plot_heatmap(filtered_counts)
    })
    
    
    ###########################################
    #' 13.5.2 NORM COUNTS OUTPUT 4 PCA
    #' Make a PCA plot of the normalized counts 
    #' Incorporate % Variance explained
    
    # Function to make PCA scatter plot
    make_pca_plot <- function(filtered_metadata, pcx, pcy) {

        # Manual load complete dataset for PCA, remove zero variance genes
        intensity <- read.csv(input$norm_counts_csv$datapath, sep='\t')
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
        
        # PC number and variance explained as a variable
        pc_x_num <- str_sub(pcx, -1)
        pc_y_num <- str_sub(pcy, -1)
        pc_x_var <- prop_of_var[as.numeric(pc_x_num)]
        pc_y_var <- prop_of_var[as.numeric(pc_y_num)]
        
        # Construct plot
        plot <- ggplot(pca_tibble, aes(x=!!sym(pcx), y=!!sym(pcy), color = diagnosis)) +
            # Labels and themes
            geom_point(size=3) +
            theme_classic() +
            ggtitle('PCA of Count Data') +
            xlab(str_c(pcx, ': ', pc_x_var, '% variance')) +
            ylab(str_c(pcy, ': ', pc_y_var, '% variance')) +
            
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
    
    
    output$counts_pca <- renderPlot({
        # Don't run until both metadata and counts file has been uploaded
        req(input$norm_counts_csv)
        req(input$sample_info_csv)
        
        # Make plot update on pushing action button
        input$plot_pca
        
        # Run PCA plot and feed in metadata
        isolate ({
            make_pca_plot(load_sample_information(), input$pca_x, input$pca_y)
        })
    })
    
    
    #####################################################
    #' 13.5.3 Differential Expression INPUT DESEQ RESULTS
    #' Only need to turn results csv into table, 
    #' since deseq results already available
    
    load_deseq_results <- reactive({
        # Don't run until file has been uploaded
        req(input$deseq_csv)
        
        deseq_results <- readr::read_tsv(input$deseq_csv$datapath, show_col_types = FALSE) %>%
            dplyr::rename(Gene = 1)
        
        return (deseq_results)
    })
    
    
    ###########################################################
    #' 13.5.3 Differential Expression OUTPUT 1 table of results
    output$deseq_table <- renderDT(
        load_deseq_results(), options = list(
            pageLength = 50, # Default number of rows to show
            lengthMenu = c(5, 10, 20, 50, 100),  # Options for number of rows to show
            dom = 'ltpf') # What to show in the data table (look up docs)
    )
    
    
    ###############################################
    #' 13.5.3 DE OUTPUT 2 Customizable Volcano Plot
    
    # Function to make a volcano plot with many customizations
    volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
        volcano <- dataf %>%
            drop_na(padj) %>%
            
            # Column to color points by padj threshold
            mutate(volc_plot_status = case_when(
                padj <= 1*10^slider ~ 'True',
                padj > 1*10^slider ~ 'False')) %>%
            
            # Log transform for volcano effect
            mutate(padj = -log10(padj)) %>%
            
            # Plot
            ggplot(aes(x=!!sym(x_name), y=!!sym(y_name), color=volc_plot_status)) +
                geom_point() +
                theme_classic() +
                scale_colour_manual(values = c('True' = color2, 'False' = color1)) +
                
                # This is hard coded for padj unforunately
                guides(color=guide_legend(title=paste(c("padj <= 1*10^", slider), collapse=''))) +
                theme(legend.position="bottom", legend.box = "horizontal") 
                
        return(volcano)
    }
    
    
    output$volcano <- renderPlot({
        # Makes this renderPlot depend on the button w/ ID 'make_plot'
        input$make_volcano_plot      
        
        # Dont super understand what isolate does
        isolate({  
            volcano_plot(load_deseq_results(), input$volcano_x, input$volcano_y, input$padj_threshold,
                         input$volcano_base_color, input$volcano_highlight_color)
        })
    })
    
    
    ########################################################
    #' 13.6.1 FGSEA INPUTS
    #' Upload GMT file, choose appropriate ranking variable,
    #' and run FGSEA
    
    # Function to run fgsea (requires DESeq and gmt file)
    run_fgsea <- reactive({
        # Don't run until files has been uploaded
        req(input$deseq_csv)
        req(input$gmt_file)
        
        # Load in DESeq results again
        deseq_results <- load_deseq_results()
        
        # Make a ranked vector of genes ranked by their log 2 fold change values
        ranked_genes <- deseq_results$symbol
        ranked_variable <- deseq_results$log2FoldChange
        ranked_list <- stats::setNames(ranked_variable, ranked_genes)
       
        # Run fgsea using uploaded gmt file 
        hallmark_pathways_fgsea <- fgsea::gmtPathways(input$gmt_file$datapath)
        fgsea_results <- fgsea(hallmark_pathways_fgsea, ranked_list, minSize = 15, maxSize = 500) 
        fgsea_results$leadingEdge <- lapply(fgsea_results$leadingEdge, paste0, collapse=',')
        fgsea_results$leadingEdge <- as.character(fgsea_results$leadingEdge)
        
        write_csv(fgsea_results, file='data/fgsea_results.csv')
        
        return (fgsea_results)
    })
    
    # Function to make Hallmark pathways more readable and shorter
    fix_hallmark <- function(hallmark) {
        fixed <- str_sub(hallmark, 10)
        fixed <- str_replace_all(fixed, '_', ' ')
        return (fixed)
    }
    
    # Function to make horizontal barplot of significant pathways after fgsea (filtered by padj values)
    plot_fgsea_barplot <- function(fgsea_res, padj_threshold) {
        
        # Add a column to describe if pathway is upregulated or downregulated
        fgsea_res <- fgsea_res %>%
            mutate(NES_status = case_when(NES >= 0 ~ 'TRUE', 
                                          NES < 0 ~ 'FALSE'))
        
        # Color map for different +/- regulation of pathway
        color_map <- c('FALSE' = '#FFC107', 'TRUE' = '#004D40')
        
        # Make barplot after filtering for padj threshold
        fgsea_barplot <- fgsea_res %>%
            dplyr::filter(padj <= 1*10^(padj_threshold)) %>%
            dplyr::mutate(pathway = fix_hallmark(pathway)) %>%
            ggplot(aes(x = reorder(pathway, NES), y = NES, fill = NES_status)) + 
                geom_bar(stat = 'identity') + 
                theme_classic() +
                ggtitle('Barplot of fgsea NES') +
                #scale_x_discrete(labels = function(x){gsub("\\s", "\n", x)}) + # add new lines after each word to 
                                                                                # deal with long labels
                theme(axis.text = element_text(size = 12), axis.title.y = element_blank()) + 
                guides(fill = 'none') +
                scale_fill_manual(values = color_map) + 
                coord_flip() 
        
        return (fgsea_barplot)
    }
    
    # Render the horizontal barplot
    output$fgsea_barplot <- renderPlot({

        plot_fgsea_barplot(run_fgsea(), input$fgsea_padj_threshold)

    })
    
    
    ####################################################
    #' 13.6.1 FGSEA OUTPUT 2 DATA TABLE OF FGSEA RESULTS 
    #' WITH FILTERS & DOWNLOAD
    
    # Function to generate filtered fgsea results tibble based on inputs
    filter_fgsea_results <- function(fgsea_res, padj_threshold, sign) {
        
        # Generate tibble with filtered padj only, fix long text of pathways
        padj_filtered <- fgsea_res %>%
            dplyr::filter(padj <= 1*10^(padj_threshold)) %>%
            dplyr::mutate(pathway = fix_hallmark(pathway))
        
        # If/else to select for positive vs negative (or both) enriched pathways
        if (sign == 'Positive') {
            sign_filtered <- padj_filtered %>%
                dplyr::filter(NES >= 0)
            
        } else if (sign == 'Negative') {
            sign_filtered <- padj_filtered %>%
                dplyr::filter(NES <= 0)
            
        } else if (sign == 'Both') {
            return (padj_filtered)
        }
        
        return (sign_filtered)
    }

    
    # Render filtered fgsea results tibble
    output$fgsea_table <- renderDT(
        filter_fgsea_results(run_fgsea(), input$fgsea_padj_threshold_2, input$nes_pathway),
        options = list(
            pageLength = 10, # Default number of rows to show
            lengthMenu = c(5, 10, 20, 50)  # Options for number of rows to show
        )
    )
    
    # Download button for fgsea_results
    
    output$download_fgsea <- downloadHandler(
        filename = function() {'filtered_fgsea_results.csv'},
        content = function(file) {
            write.csv(
                filter_fgsea_results(run_fgsea(), input$fgsea_padj_threshold_2, input$nes_pathway), 
                file)
        }
    )
    
    
    ########
    #' 13.6.1 OUTPUT 3 SCATTER PLOT
    
    # Make plot
    fgsea_scatter <- function(fgsea_res, padj_threshold) {
        
        # Add column to see if padj values meet threshold and log transform
        filtered_fgsea <- fgsea_res %>%
            dplyr::mutate(padj_status = case_when(padj <= 1*10^(padj_threshold) ~ 'TRUE', 
                                                  .default = 'FALSE')) %>%
            dplyr::mutate(transformed_padj = -log10(padj))
        
        # Color vector for plot
        color_map <- c('TRUE' = '#FFC107', 'FALSE' = '#004D40')
        
        # Scatter plot
        scatter <- filtered_fgsea %>%
            ggplot(aes(x=NES, y=transformed_padj, color=padj_status)) + 
            geom_point(size=3) + 
            theme_classic() +
            ggtitle('NES vs. -log10 Adjusted P Value') +
            ylab('-log10(padj)') + 
            
            # Legend/text elements
            theme(plot.title = element_text(size=20, face='bold'), 
                  axis.title.x = element_text(size=14), axis.title.y = element_text(size=14),
                  legend.title = element_text(size=14), legend.text = element_text(size=14),
                  legend.box.background = element_rect(color = 'black', linewidth = 0.8),
                  legend.background = element_rect(fill=alpha('grey', 0.5)),
                  legend.position = "bottom", legend.box = "horizontal") +
            
            # Edit legend
            guides(color = guide_legend(title = paste('padj <= 1*10^', padj_threshold),
                                        title.position = 'top', title.hjust = 0.5)) +
            
            # Set colors
            scale_color_manual(values = color_map)
        
        return (scatter)
    }

    # Render plot
    output$fgsea_scatter <- renderPlot(
        
        # Generate scatter plot 
        fgsea_scatter(run_fgsea(), input$fgsea_padj_threshold_3)
        
    )
        
}

#Run the application
shinyApp(ui = ui, server = server)
