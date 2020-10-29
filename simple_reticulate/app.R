library(shinyalert) # pop-up error messages
library(shinycssloaders) # spinning plot loading icon
library(icesTAF) # contains mkdir
library(tidyverse)
library(reticulate)
library(shiny)
PYTHON_DEPENDENCIES = c('numpy')

# set up python environment
virtualenv_dir = Sys.getenv('VIRTUALENV_NAME')
python_path = Sys.getenv('PYTHON_PATH')
reticulate::use_virtualenv(virtualenv_dir, required = T)

# source local scripts
source("dye_calling_support_script.R")
reticulate::source_python('analyze_data.py')

# app begins here
server <- function(input, output) {
### -------------- Set up reactive data containers and read in data -------------- ###  
    # set up the reactive values container 
    values <- reactiveValues()
    
    # read in sample formatted data
    values$df_all <- readRDS("sample_data_three_dyes.rds") 
    
### -------------- Condition data for neural net -------------- ###
    # triggers when values$df_all appears or changes
    observeEvent(values$df_all, {
      # write instructions for the plot
        values$main_plot <- facet_wrap_linetyped2(values$df_all, title = "Plot test!", facets_wide = 6 )
        
      # reformat data for the NN model
      values$df_wide <- values$df_all %>%
                            mutate(type_f = factor(type, levels = c("protein", "buffer"))) %>%
                            select(well, channel_f, type_f,  Temperature, value)  %>%
                            arrange(type_f, channel_f, Temperature) %>%
                            pivot_wider( names_from = c(type_f, channel_f, Temperature), values_from = value )
    })
    
### -------------- Call hits using the neural net -------------- ###
    # triggers when values$df_wide appears or changes
    observeEvent(values$df_wide, {
      values$python_output <- analyze_data(values$df_wide)
    })
    
### -------------- Render GUI elements -------------- ###
    output$plot <- renderPlot( values$main_plot ) 
    output$which_python <- renderText({  values$python_output })
    output$raw_table <- renderTable({ values$df_wide  %>% head() })
}

# Define UI for application that draws a histogram
ui <- fluidPage(useShinyalert(),
                          plotOutput("plot"), # display plot
                          tableOutput('raw_table') # print df_wide, the data passed to the NN
                          
)

shinyApp(ui, server) # run the app
