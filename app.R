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
#library(DT)

# make relevant functions, objects


# Define UI for application
ui <- fluidPage(
    theme = bslib::bs_theme(bootswatch = "spacelab"),
    # Application title
    titlePanel("MAF Explorer"),
    
    # Sidebar with file input 
    sidebarLayout(
        sidebarPanel(
            fileInput(inputId = "maf_file",  buttonLabel = "Upload...",
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
                                   DT::dataTableOutput("sample_summ") 
                                   ),
                            column(width = 6,
                                   h4("Gene summary"),
                                   DT::dataTableOutput("gene_summ")
                            )
                         )
                ),
                tabPanel("Sample"
                         ),
                tabPanel("Gene",
                         fluidRow(
                             column(width = 12,
                                    h4("Amino acid changes on protein structure"),
                                    plotOutput("gene_lollipop"),
                                    br(),
                                    DT::dataTableOutput("gene_data")
                             )
                         ) 
                )
            )
        , width = 9)
    )
)

# Define server logic 
#shinyServer
server <- function(input, output, session) {
   
    
    input_data_chk<-reactive({
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
    
    # list for reactive data values
    rv <- reactiveValues()
    
    #input_data()
    
    observeEvent(input$maf_file, {
        
        
        rv$maf<-read.delim(input$maf_file$datapath, stringsAsFactors = F, comment.char = "#", header = T)
        # print(head(rv$maf))
        #create maf object
        rv$maf_obj<-read.maf(rv$maf)
        
        rv$batch_name<-str_remove(input$maf_file$name, pattern = fixed(".maf"))
        
        #get sample names
        updateSelectInput(inputId = "sample", choices = unique(rv$maf_obj@data$Tumor_Sample_Barcode))
        #get gene names
        updateSelectInput(inputId = "gene", choices = unique(rv$maf_obj@data$Hugo_Symbol))
        #get variant classes
        updateCheckboxGroupInput(inputId = "var_class", choices = unique(rv$maf_obj@data$Variant_Classification))
        #get variant type
        updateCheckboxGroupInput(inputId = "var_type", choices = unique(rv$maf_obj@data$Variant_Type))
        
        #maf summary plot
        output$batch_maf_summary <- renderPlot({
           plotmafSummary(rv$maf_obj, addStat = 'median')
        })
        #oncoplot
        output$oncoplot <- renderPlot({
            oncoplot(maf = rv$maf_obj)
        })
        #sample summary
        output$sample_summ<-DT::renderDataTable({
            DT::datatable(getSampleSummary(x=rv$maf_obj), extensions = c("Buttons"), fillContainer = T,
                          rownames = FALSE,
                      options = list(paging = TRUE, scrollX= TRUE, ordering = TRUE, scrollY= TRUE,
                                     dom = 'Bfrtip', pageLength = 5,
                                     buttons = c('copy', 'csv', 'excel')), class = "display")}, server = FALSE)
        #gene summary
        #server = TRUE
        output$gene_summ<-DT::renderDataTable({
            DT::datatable(getGeneSummary(x=rv$maf_obj), extensions = c("Buttons"), fillContainer = T,
                          rownames = FALSE,
                          options = list(paging = TRUE, scrollX=TRUE, ordering = TRUE, pageLength = 5,
                                         scrollY=TRUE, dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')), 
                          class = "display")}, server = FALSE)
        #print(head(rv$maf))
    })
    #print(head(rv$maf))
    # #Gene view features
    observeEvent(input$gene, {
        # current_maf<-reactive({
        #     req(rv)
        #     rv$maf_obj
        #     })
        #current_maf_obj<-reactive({rv$maf_obj()})
        #head(current_maf_obj@data)
        output$gene_lollipop<-renderPlot({lollipopPlot(maf = rv$maf_obj, gene = input$gene)})

        # rv$gene_flt_data<-isolate(rv$maf_obj@data) %>% filter(Hugo_Symbol==input$gene,
        #                                           Variant_Classification %in% input$var_class,
        #                                           Variant_Type %in% input$var_type)
        # #print(head(rv$gene_flt_data))
        # output$gene_data<-DT::renderDataTable({
        #     DT::datatable(x=rv$gene_flt_data, extensions = c("Buttons"), fillContainer = T,
        #                   rownames = FALSE,
        #                   options = list(paging = TRUE, scrollX=TRUE, ordering = TRUE, pageLength = 5,
        #                                  scrollY=TRUE, dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')),
        #                   class = "display")}, server = FALSE)

    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
