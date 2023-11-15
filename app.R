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
        tabPanel('Samples', 
                 sidebarLayout(
                     sidebarPanel('input: Sample information matrix in CSV format'),
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

}



#Run the application
shinyApp(ui = ui, server = server)
