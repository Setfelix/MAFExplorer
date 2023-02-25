## MAF Explorer
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#increase shiny maximum file upload size to 500 MB
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
                      multiple = F, accept = ".maf", placeholder = "Please upload a MAF file"),
            selectInput(inputId = "sample", label = "Select tumor sample", choices = character(), multiple = F),
            downloadButton("report", "Generate Batch report")),
        mainPanel(
            tabsetPanel(
                tabPanel("Batch",
                         br(),
                         plotOutput("total_maf_summary"),
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
server <- function(input, output, session) {
    
    input_data<-reactive({
        req(input$maf_file)
        
        #check file extension
        ext<-tools::file_ext(input$maf_file$name)
        #print(ext)
        validate(need(ext == "maf", "Please upload a .maf file"))
        # switch(ext, 
        #        maf = vroom::vroom(input$maf_file$datapath, comment = "#", delim = "\t"),
        #        validate("Invalid file; Please upload a .maf file")
        #        )
    })
    
    # list for reactive values
    rv <- reactiveValues()

    #load sample names into selectInput
    observeEvent(input$maf_file, {
        freezeReactiveValue(input, "sample")
        maf_df<-read.delim(input$maf_file$datapath, stringsAsFactors = F, comment.char = "#", header = T)
        rv$maf<-maf_df
        #rv$samples<-rv$maf$
        #maf_df$Tumor_Sample_Barcode
        updateSelectInput(inputId = "sample", 
                          choices = unique(rv$maf$Tumor_Sample_Barcode))
    })
    output$total_maf_summary <- renderPlot({
       # draw the histogram with the specified number of bins
        maf_obj<-read.maf(input_data())
        plotmafSummary(maf_obj)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
