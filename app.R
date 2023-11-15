#### RShiny Final Project ####
#    BF591 Fall 2023 
#    Johnathan Zhang
#    Downloaded data from: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE64810



#Packages
library(shiny)
library(bslib)
library(tidyverse)
library(colourpicker)
library(GEOquery)

#gse=getGEO(filename="data/GSE64810_series_matrix.txt")
#metadata <- gse@phenoData@data

ui <- fluidPage(

)


server <- function(input, output, session) {

}



#Run the application
shinyApp(ui = ui, server = server)
