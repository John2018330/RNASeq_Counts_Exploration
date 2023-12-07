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
                         tabPanel('Summary', 'summary of table that includes summary of 
                                  type and values in each column'), 
                         tabPanel('Table', 'data table displaying the sample information, 
                                  with sortable columns'),
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

    #' Function to take in csv uploaded in 13.5.1 Sample Information Exploration and return it
    #' as tibble
    load_sample_information <- reactive({
        #Don't run until file has been uploaded
        req(input$sample_info_csv)
        
        #GSE object, structures lots of information about GEO accession/data given to it
        #gse=getGEO(filename="data/GSE64810_series_matrix.txt")
        gse=getGEO(filename=input$sample_info_csv$datapath)
        
        #Sample information can be extracted from GSE object 
        metadata <- gse@phenoData@data
        
        return (metadata)
        
        #req(input$csv)
        #return(read.csv(file = input$csv$datapath))
    })
    
    
}



#Run the application
shinyApp(ui = ui, server = server)
