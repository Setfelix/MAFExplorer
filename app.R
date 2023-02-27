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
    theme = bslib::bs_theme(bootswatch = "spacelab"),
    # Application title
    titlePanel("MAF Explorer"),
    
    # Sidebar with file input 
    sidebarLayout(
        sidebarPanel(
            fileInput(inputId = "maf_file", buttonLabel = "Upload...",
                      label = "Load MAF file", 
                      multiple = F, accept = ".maf", placeholder = "Please upload a MAF file"),
            selectInput(inputId = "sample", label = "Select tumor sample", choices = character(), multiple = F),
            selectInput(inputId = "gene", label = "Select gene", choices = character(), multiple = F),
            checkboxGroupInput("var_class", label = "Variant class", choices = character()),
            checkboxGroupInput("var_type", label = "Variant type", choices = character()),
            downloadButton("report", "Generate Batch report"), width = 3),
        mainPanel(
            tabsetPanel(
                tabPanel("Batch",
                         #br(),
                         fluidRow(
                             column(width = 6,
                                h4("MAF Summary"),
                                plotOutput("batch_maf_summary")
                             ),
                             column(width = 6,
                                 h4("Top 20 mutated genes"),
                                 plotOutput("oncoplot")
                             )
                         ),
                         fluidRow(
                            column(width = 6,
                                   h4("Sample summary"),
                                   dataTableOutput("sample_summ") 
                                   ),
                            column(width = 6,
                                   h4("Gene summary"),
                                   dataTableOutput("gene_summ")
                            )
                         )
                         
                ),
                tabPanel("Sample",
                         br(),
                         tableOutput("sample_table_summary"),
                         br(),
                         plotOutput("sample_ss_histplot"),
                         br(),
                         plotOutput("sample_ds_histplot")
                ),
                tabPanel("Gene")
            )
        , width = 9)
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
        
        rv$batch_name<-str_remove(input$maf_file$name, pattern = fixed(".maf"))
        #rv$samples<-rv$maf$
        #maf_df$Tumor_Sample_Barcode
        
        #get sample names
        updateSelectInput(inputId = "sample", choices = unique(rv$maf$Tumor_Sample_Barcode))
        #get gene names
        updateSelectInput(inputId = "gene", choices = unique(rv$maf$Hugo_Symbol))
        #get variant classes
        updateCheckboxGroupInput(inputId = "var_class", choices = unique(rv$maf$Variant_Classification))
        #get variant type
        updateCheckboxGroupInput(inputId = "var_type", choices = unique(rv$maf$Variant_Type))
        
        #maf summary plot
        rv$maf_obj<-read.maf(rv$maf)
        output$batch_maf_summary <- renderPlot({
            # draw the histogram with the specified number of bins
            plotmafSummary(rv$maf_obj, addStat = 'median')
        })
        #oncoplot
        output$oncoplot <- renderPlot({
            oncoplot(maf = rv$maf_obj)
        })
        #sample summary
        output$sample_summ<-renderDataTable(getSampleSummary(x=rv$maf_obj), 
                                            options = list(paging = TRUE, pageLength = 5, scrollX= TRUE))
        #gene summary
        output$gene_summ<-renderDataTable(getGeneSummary(x=rv$maf_obj), 
                                          options = list(paging = TRUE, pageLength = 5, scrollX=TRUE))
        
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
