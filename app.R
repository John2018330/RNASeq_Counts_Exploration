################################
#'### RShiny Final Project #####
#'##############################
#'    BF591 Fall 2023 
#'    Johnathan Zhang
#'    Data from GEO accession: GSE64810
#'       - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810



####Packages####
library(shiny)
library(bslib)
library(tidyverse)
library(colourpicker)
library(GEOquery)
library(DT)


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
                         tabPanel('Summary', 'Summary of table that includes summary of 
                                  type and values in each column',
                                  DTOutput('metadata_summary')), 
                         
                         tabPanel('Table', 'Tab with a data table displaying the sample information, 
                                  with sortable columns', 
                                  DTOutput('metadata')),
                         
                         tabPanel('Plots', 'histograms, density plots, or violin plots 
                                  of continuous variables')
                     ))
                 )),
        
        tabPanel('Counts Exploration', 
                 sidebarLayout(
                     sidebarPanel('Normalized counts matrix, in CSV format, input sliders 1 and 2'),
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

    ####################
    #' 13.5.1 INPUT FILE
    #' Function to take in file uploaded in Sample Information Exploration, data wrangle it
    #' and return it as tibble
    load_sample_information <- reactive({
        #Don't run until file has been uploaded
        req(input$sample_info_csv)
        
        #GSE object, structures lots of information about GEO accession/data given to it
        #gse=getGEO(filename="data/GSE64810_series_matrix.txt")
        gse=getGEO(filename=input$sample_info_csv$datapath)
        
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
            
            #change column entries from "characteristic_name:characteristic_value" 
            #to "characteristic_value".
            #case_when to ignore the empty values in non Huntingtons patients 
            #rows for huntington's specific columns
            mutate(across(c(-1, -2), ~ case_when(. == "" ~ .,
                                                 . != "" ~ str_split_i(., ":", -1)))) %>%
            
            #change columns that should be numeric to numeric, adds NA's where 
            #necessary (will give warning)
            mutate(across(c(-1:-3), as.numeric)) %>%
        
            #get read of leading space in diagnosis column
            mutate(diagnosis = trimws(diagnosis)) %>%
            
            #change diagnosis column to a factor
            mutate_at('diagnosis', factor)
        
        #relevel diagnosis column in md_filtered it so control is reference
        md_filtered$diagnosis <- relevel(md_filtered$diagnosis, ref='Neurologically normal')
        
        return (md_filtered)
    })
    
    ######################
    #' 13.5.1 OUTPUT TAB 1 
    #' Render Data Table summary information about metadata
    
    # A helper function to find mean and SD of numeric columns and join them into a string
    mean_sd_f <- function(df_col) {
        stringr::str_c(round(mean(df_col, na.rm=TRUE), 2), ' (+/-', round(sd(df_col, na.rm=TRUE), 2), ')')
    }
    
    # Main function that takes in metadata tibble and summarizes it
    summarize_metadata <- function() {
        #Don't run until file has been uploaded
        req(input$sample_info_csv)
        
        # Load in the metadata
        md_filtered <- load_sample_information()
        
        # Make first two columns of table: Col 1: Column Name), Col 2: Data type
        md_summary <- (summarise_all(md_filtered, class))                   #class of column
        md_summary <- as.data.frame(t(md_summary[,]))                       #transpose
        colnames(md_summary)[1] = 'Type'                                    #rename 2nd column
        md_summary <- rownames_to_column(md_summary, var = "Column Name")   #rename 1st column
        
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
            dom='ltp', #what to show in the data table (look up docs)
            pageLength = 20 #default number of rows to show
        )
    )
    
    ######################
    #' 13.5.1 OUTPUT TAB 2 
    #' Render the Data Table with sortable columns 
    output$metadata <- renderDT(
        load_sample_information(), options = list(
            pageLength = 50, #default number of rows to show
            lengthMenu = c(5, 10, 20, 50, 100),  #options for number of rows to show
            dom = 'ltp') #what to show in the data table (look up docs)
    )
    
    
    
    
}



#Run the application
shinyApp(ui = ui, server = server)
