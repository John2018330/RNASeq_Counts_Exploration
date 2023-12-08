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

options(shiny.maxRequestSize=30*1024^2)


#### UI / Front end ####
ui <- fluidPage(
    titlePanel("BF591 Final Project"),
    h3('A Rshiny website to explore RNA-seq results'),
    p(strong("Made by Johnathan Zhang")),
    
    tabsetPanel(
        
        #### 13.5.1 Sample Information Exploration 
        tabPanel('Samples', 
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
        
        #### 13.5.2 Normalized Counts Matrix Exploration
        tabPanel('Counts Exploration', 
                 sidebarLayout(
                     ### Upload norm_counts_csv and filter sliders
                     sidebarPanel(
                         ## Upload the csv
                         fileInput('norm_counts_csv', 
                                   paste0('Normalized counts matrix, in CSV format, input sliders 1 and 2')),
                         
                         ## Slider for variance percentile (0-100%)
                         sliderInput("norm_counts_percentile_var", min = 0, max = 100, value = 90,  step = 0.5,
                                     label = "Select the percentile of genes to include based on variance across samples 
                                     (e.g., 90% means select genes with the top 10% variance values)"),
                         
                         ## Slider for number of non-zero counts
                         sliderInput('norm_counts_nonzeroes', min=0, max=69, value = 10, step = 1,
                                     label = 'Select how many samples must have a non-zero count value')),
                     
                     ### Output tabs: Summary, Plots, Heatmap, and PCA
                     mainPanel(tabsetPanel(
                         tabPanel('Summary', 'Tab with text or a table summarizing the effect of the
                                  filtering'),
                         tabPanel('Plots', 'Tab with diagnostic scatter plots, where genes passing 
                                  filters are darker, and genes filtered out are lighter'),
                         tabPanel('Heatmap', 'Tab with a clustered heatmap of counts remaining after
                                  filtering'),
                         tabPanel('PCA', 'Tab with a scatter plot of principal component analysis
                                  projections')
                     ))
                 )),
        
        tabPanel('Differential Expression', 
                 sidebarLayout(
                     sidebarPanel('Results of a differential expression analysis in CSV format.'),
                     mainPanel(tabsetPanel(
                         tabPanel('Table', 'Tab with sortable table displaying differential expression
                                  results'),
                         tabPanel('unknown', 'Tab with content similar to that described in [Assignment
                                  7]')
                     ))
                 )),
        tabPanel('CYO Adventure', 
                 'either GSEA or Correlation Network Analysis')
    )
    
)

#### Server side / Back end ####
server <- function(input, output, session) {

    ########################################
    #' 13.5.1 SAMPLE INFO EXPLORE INPUT FILE
    #' Function to take in file uploaded in Sample Information Exploration, data wrangle it
    #' and return it as tibble
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
            rename(diagnosis         = characteristics_ch1.1, pmi               = characteristics_ch1.2,
                   age_of_death      = characteristics_ch1.3, rin               = characteristics_ch1.4,
                   total_reads       = characteristics_ch1.5, age_of_onset      = characteristics_ch1.6,
                   duration          = characteristics_ch1.7, cag               = characteristics_ch1.8,
                   vonsattel_grade   = characteristics_ch1.9, hv_striatal_score = characteristics_ch1.10,
                   hv_cortical_score = characteristics_ch1.11) %>%
            
            # Change column entries from "characteristic_name:characteristic_value" 
            # to "characteristic_value".
            # Case_when to ignore the empty values in non Huntingtons patients 
            # rows for huntington's specific columns
            mutate(across(c(-1, -2), ~ case_when(. == "" ~ .,
                                                 . != "" ~ str_split_i(., ":", -1)))) %>%
            
            # Change columns that should be numeric to numeric, adds NA's where 
            # necessary (will give warning)
            mutate(across(c(-1:-3), as.numeric)) %>%
        
            # Get read of leading space in diagnosis column
            mutate(diagnosis = trimws(diagnosis)) %>%
            
            # Change diagnosis column to a factor
            mutate_at('diagnosis', factor)
        
        # Relevel diagnosis column in md_filtered it so control is reference
        md_filtered$diagnosis <- relevel(md_filtered$diagnosis, ref='Neurologically normal')
        
        return (md_filtered)
    })
    
    ##########################################
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
    
    
    ##########################################
    #' 13.5.1 SAMPLE INFO EXPLORE OUTPUT TAB 2 
    #' Render the Data Table with sortable columns 
    output$metadata <- renderDT(
        load_sample_information(), options = list(
            pageLength = 50, # Default number of rows to show
            lengthMenu = c(5, 10, 20, 50, 100),  # Options for number of rows to show
            dom = 'ltp') # What to show in the data table (look up docs)
    )
    
    
    ##########################################
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
    
    
    #######################################
    #' 13.5.2 NORM COUNTS EXPLORE INPUT CSV
    #' Input a raw CSV and add columns that will be used for filtering/plotting and return it as a tibble
    #' - percentile variance: a percentile rank of how much variance a gene has 
    #' - non_zero_count:      a count of how many genes had non-zero counts
    
    load_normalized_counts <- function(counts_file) {
        # Don't run until file has been uploaded
        req(input$norm_counts_csv)
        
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
            mutate(var_percentile = percent_rank(var), .after = var)
        
        return(norm_counts)
    }
    
}



#Run the application
shinyApp(ui = ui, server = server)
