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

facet_no_y_theme <- theme( # for the first hit-calling plot, the most austere
  text = element_text(size = 8),
  axis.title = element_blank(), # don't label the axes
  axis.text.x = element_text(), # don't label the numbers on the axes
  axis.ticks = element_blank(), # dont have ticks
  legend.position = "right", # put a legent on the right of the plot
  plot.title = element_text(lineheight=.8, face="bold", size = 12), # have a title
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  strip.background = element_blank(),
  panel.spacing.x = unit(0.1, "lines"),
  panel.spacing.y = unit(0.1, "lines"),
  aspect.ratio = (1/1.618)
)

facet_wrap_linetyped2 <- function(df_melt, title, facets_wide) {
  p <- ggplot(df_melt, aes(x = Temperature, # temperature on X
                           y = value, # RFU on y
                           color = channel_f, # colored by the state
                           linetype = type,
                           group = dye_conc_type_channel # group means series, as in, this defines the unique data sets
  )) +
    geom_line(size = 0.3, alpha = 0.8) + # change the line type depending on the dye concentration # linetype = df_melt$conc #
    facet_wrap(~dye, scales = "free", ncol = facets_wide) +
    labs(title = title, color = "Channel") +
    theme_bw() +
    scale_color_manual( values = c("Cy5.5" = "#67000d", "Cy5" = "#a50f15", "ROX" = "#ef3b2c", "TAMRA" = "#f16913", "JOE" = "#74c476", "FAM" = "#2171b5")) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    facet_no_y_theme
  
  p # return the plot
}

# app begins here
server <- function(input, output) {
### -------------- Set up reactive data containers and read in data -------------- ###  
    # set up the reactive values container 
    values <- reactiveValues()
    
    #values$data_name <- "Exp0618--OLF-processed.rds"
    
    # read in sample formatted data
   values$df_all <- readRDS("Exp0855_2020-09-27_df_all.rds") #readRDS("Exp0933_2020-10-22_df_all.rds") #readRDS("Exp0950_hisGST_2020-11-11_df_all.rds") #readRDS("Exp0891_2020-09-28_df_all_SP0182.rds") # readRDS("Exp0900_2020-10-02_SP0185-Nsp1410__df_all.rds") # readRDS("20201029_olf_reprocessed.rds") 
   # values$df_all <- readRDS("Exp0945_2020-11-11_df_all_SP197_Nucleocapsid.rds") # readRDS("20201029_olf_reprocessed.rds") 
    
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
      values$python_output <- classify('dsf_net.pt', values$df_wide)
      print(values$python_output %>% str()) 
                              
      
      values$df_outcomes <- values$python_output %>% 
                              set_names(c("reject", "sensitive", "hit", "well")) %>%
                              pivot_longer(-well, names_to = "outcome", values_to = "log_prob") %>%
                              left_join( ., values$df_all %>% select(well, dye) %>% distinct(), by = "well") %>%
                              filter(dye != "Empty")
      
      values$df_outcomes_summary <- values$df_outcomes %>%
                                    group_by(well) %>%
                                    filter(log_prob == max(log_prob))
      
      values$df_plot_outcomes <-  values$df_outcomes %>%
        arrange(outcome, log_prob) %>%
        mutate(dye_f = factor(dye, levels = dye %>% unique)) 
      
      
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
    
    
    output$prob_plot <- renderPlot({
 
      values$df_plot_outcomes %>%
                ggplot(aes(x = dye_f, y = log_prob, color = outcome)) +
          geom_point(alpha = 0.5, size = 3) +
          theme(axis.text.x = element_blank())
                                })
    
    output$hover_data <- renderTable({
      brushedPoints(values$df_plot_outcomes, input$plot_brush) %>% select(dye, outcome, log_prob) %>% pivot_wider(names_from = outcome, values_from = log_prob)
    })
    
    
    
    output$download_plot <- downloadHandler(
      filename = function() {
        
        paste0("rank_orderd_screen_plot.pdf")
      },
      content = function(file) {
        
        dye_order <- values$df_plot_outcomes %>%
          filter(outcome == "hit") %>%
          arrange(log_prob) %>%
          mutate(dye_f = factor(dye, levels = dye %>% unique())) %>%
          pull(dye_f) %>% unique()
        
        ### make and save the primary plot
        values$df_all_ord <- values$df_all %>%
                              mutate(dye = factor(dye, levels = rev(dye_order)))

        values$main_plot <- facet_wrap_linetyped2(values$df_all_ord, title = "Ordered raw dye screen results, ranked by NN", facets_wide = 6 )
        
        values$plot_ratio <- values$df_all_ord %>%
          select(dye) %>%
          distinct() %>%
          nrow()
        
        ggsave(file, values$main_plot, width = 10, height =  1.5*(1/1.618 * values$plot_ratio/6 + 1), limitsize = FALSE )
        
      }
    )
    
}

# Define UI for application that draws a histogram
ui <- fluidPage(useShinyalert(),
                         # plotOutput("plot"), # display plot
                         # tableOutput('raw_table'), # print df_wide, the data passed to the NN
                plotOutput("prob_plot", brush = "plot_brush"),
                downloadButton("download_plot", "Download ordered raw data plot"),
                tableOutput("hover_data"),
                tableOutput('prob_table'),
                tableOutput('prob_table_tally'),
                verbatimTextOutput('hit_list'),
                verbatimTextOutput('sensitive_list')
)

shinyApp(ui, server) # run the app
