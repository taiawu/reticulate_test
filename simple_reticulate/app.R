library(shinyalert) # pop-up error messages
library(shinycssloaders) # spinning plot loading icon
library(icesTAF) # contains mkdir
library(tidyverse)
library(reticulate)
library(shiny)
PYTHON_DEPENDENCIES = c('numpy', 'torch', 'pandas')

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
    
    #values$data_name <- "Exp0618--OLF-processed.rds"
    
    # read in sample formatted data
    values$df_all <- readRDS("20201029_olf_reprocessed.rds") 
    
### -------------- Condition data for neural net -------------- ###
    # triggers when values$df_all appears or changes
    observeEvent(values$df_all, {
      # write instructions for the plot
        #values$main_plot <- facet_wrap_linetyped2(values$df_all, title = "Plot test!", facets_wide = 6 )
        
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
      values$python_output <- classify('dsf_net_dict.pt', values$df_wide)
                              
      
      values$df_outcomes <- values$python_output %>% 
                              set_names(c("reject", "sensitive", "hit", "well")) %>%
                              pivot_longer(-well, names_to = "outcome", values_to = "log_prob") %>%
                              left_join( ., values$df_all %>% select(well, dye) %>% distinct(), by = "well") %>%
                              filter(dye != "Empty")
      
      values$df_outcomes_summary <- values$df_outcomes %>%
                                    group_by(well) %>%
                                    filter(log_prob == max(log_prob))
      
      
      values$df_outcomes_tally <- values$df_outcomes %>%
                                    group_by(well) %>%
                                    filter(log_prob == max(log_prob)) %>%
                                    ungroup() %>%
                                    count(outcome)
      
      values$hits <- values$df_outcomes_summary %>% filter(outcome == "hit") %>% arrange(dye) %>% pull(dye)
      values$sensitive <- values$df_outcomes_summary %>% filter(outcome == "sensitive")  %>% arrange(dye) %>% pull(dye)
      
      print(values$python_output %>% str()) 
    })
    

    
### -------------- Render GUI elements -------------- ###
   # output$plot <- renderPlot( print(facet_wrap_linetyped2(values$df_all, title = "Plot test!", facets_wide = 6 ))) 
    output$which_python <- renderText({  values$python_output })
    #output$raw_table <- renderTable({ values$df_wide %>% head() })
    output$prob_table <- renderTable({ values$df_outcomes_summary %>% head() })
    output$prob_table_tally <- renderTable({ values$df_outcomes_tally })
    output$hit_list <- renderText( values$hits)
    output$sensitive_list <- renderText(values$sensitive)
}

# Define UI for application that draws a histogram
ui <- fluidPage(useShinyalert(),
                         # plotOutput("plot"), # display plot
                         # tableOutput('raw_table'), # print df_wide, the data passed to the NN
                tableOutput('prob_table'),
                tableOutput('prob_table_tally'),
                verbatimTextOutput('hit_list'),
                verbatimTextOutput('sensitive_list')
)

shinyApp(ui, server) # run the app
