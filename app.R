## MAF Explorer
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#increase maximum file upload size to 500 MB
options(shiny.maxRequestSize = 500*1024^2)

library(shiny)
library(tidyverse)
library(maftools)

# make relevant functions

# Define UI for application
ui <- fluidPage(
    
    # Application title
    titlePanel("MAF Explorer"),
    
    # Sidebar with file input 
    sidebarLayout(
        sidebarPanel(
            fileInput(inputId = "maf_file", buttonLabel = "Upload...",
                      label = "Load MAF file", 
                      multiple = F, accept = ".maf"),
            selectInput(inputId = "sample", label = "Select tumor sample", choices = character()),
            downloadButton("report", "Generate Batch report")),
        mainPanel(
            tabsetPanel(
                tabPanel("Batch",
                         br(),
                         plotOutput("ss_histplot"),
                         br(),
                         plotOutput("ds_histplot"),
                         br(),
                         plotOutput("duprate_scatter", dblclick = "duprate_scatter_dblclick",
                                    brush = brushOpts(id = "duprate_scatter_brush", resetOnNew = TRUE)
                         ),
                         br(),
                         plotOutput("nreads_barplot")
                ),
                tabPanel("Sample",
                         br(),
                         tableOutput("sample_table_summary"),
                         br(),
                         plotOutput("sample_ss_histplot"),
                         br(),
                         plotOutput("sample_ds_histplot")
                )
            )
        )
    )
)

# Define server logic 
server <- function(input, output) {
    
    input_data<-reactive({
        maf <- input$maf_file
        req(maf)
        
        #check file extension
        ext<-tools::file_ext(input$maf_file$datapath)
        validate(need(ext == "maf", "Please upload a .maf file"))
        # switch(ext, tsv = vroom::vroom(input$qc_files$datapath, delim = "\t"),
        #        validate("Invalid file; Please upload a .csv or .tsv file"))
    })

    #load sample names into selectInput
    observeEvent(input$maf_file, {
        freezeReactiveValue(input, "sample")
        maf_df<-read.delim(input$maf_file$datapath, stringsAsFactors = F, comment.char = "#", header = T)
        updateSelectInput(inputId = "sample", 
                          choices = unique(maf_df$sample))
    })
    # output$distPlot <- renderPlot({
    #     # generate bins based on input$bins from ui.R
    #     x    <- faithful[, 2]
    #     bins <- seq(min(x), max(x), length.out = input$bins + 1)
    # 
    #     # draw the histogram with the specified number of bins
    #     hist(x, breaks = bins, col = 'darkgray', border = 'white',
    #          xlab = 'Waiting time to next eruption (in mins)',
    #          main = 'Histogram of waiting times')
    # })
}

# Run the application 
shinyApp(ui = ui, server = server)
