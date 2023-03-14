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
library(mclust)
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
            checkboxGroupInput("var_type", label = "Variant type", choices = character())
            , width = 3),
        mainPanel(
            tabsetPanel(
                tabPanel("Batch",
                         br(),
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
                tabPanel("Sample",
                         fluidRow(column(width = 6, h4("VAF plot"), plotOutput("vafplot")),
                                  column(width = 6, h4("Hyper mutated regions (Kataegis)"),
                                         DT::dataTableOutput("raintable")),
                         fluidRow(column(width = 6, 
                                         h4("Tumor heterogeneity"), 
                                         plotOutput("hetplot")),
                                  column(width = 6, 
                                         h4("VAF Cluster summary"), 
                                         DT::dataTableOutput("vaf_cluster_summ"))
                                  ),
                                  )
                         ),
                tabPanel("Gene",
                         fluidRow(
                             column(width = 12, h4("Amino acid changes on protein structure"), 
                                    plotOutput("gene_lollipop"),
                                    h4("Variants in selected gene"),
                                    DT::dataTableOutput("gene_data"))) 
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
        rv$maf<-rv$maf %>% filter(FILTER=="PASS")
        # print(head(rv$maf))
        #create maf object
        rv$maf_obj<-read.maf(rv$maf)
        
        rv$batch_name<-str_remove(input$maf_file$name, pattern = fixed(".maf"))
        
        #get sample names
        updateSelectInput(inputId = "sample", choices = unique(rv$maf_obj@data$Tumor_Sample_Barcode))
        #get gene names
        updateSelectInput(inputId = "gene", choices = unique(rv$maf_obj@data$Hugo_Symbol))
        #get variant classes
        updateCheckboxGroupInput(inputId = "var_class", choices = unique(rv$maf_obj@data$Variant_Classification), 
                                 selected = unique(rv$maf_obj@data$Variant_Classification))
        #get variant type
        updateCheckboxGroupInput(inputId = "var_type", choices = unique(rv$maf_obj@data$Variant_Type), 
                                 selected = unique(rv$maf_obj@data$Variant_Type))
        
        #maf summary plot
        output$batch_maf_summary <- renderPlot({
           plotmafSummary(rv$maf_obj, addStat = 'median', textSize = 1)
        })
        #oncoplot
        output$oncoplot <- renderPlot({
            oncoplot(maf = rv$maf_obj, fontSize = 0.8)
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
    #Gene view features
    observeEvent(input$gene, {
        # current_maf<-reactive({
        #     req(rv)
        #     rv$maf_obj
        #     })
        #current_maf_obj<-reactive({rv$maf_obj()})
        #head(current_maf_obj@data)
        output$gene_lollipop<-renderPlot({lollipopPlot(maf = read.maf(rv$maf %>%
                                                                          filter(Variant_Classification %in% input$var_class, 
                                                                                 Variant_Type %in% input$var_type,
                                                                                 Hugo_Symbol == input$gene)),
                                                       gene = input$gene, showDomainLabel = F, domainLabelSize = 1, labPosSize = 1,
                                                       legendTxtSize = 1)})

        # rv$gene_flt_data<-isolate(rv$maf_obj@data) %>% filter(Hugo_Symbol==input$gene,
        #                                           Variant_Classification %in% input$var_class,
        #                                           Variant_Type %in% input$var_type)
        # #print(head(rv$gene_flt_data))
        #rv$maf_obj@data %>% 
        #filter(Hugo_Symbol==input$gene, Variant_Classification %in% input$var_class, Variant_Type %in% input$var_type)
        maf_data<-reactive({rv$maf_obj@data})
        output$gene_data<-DT::renderDataTable({
            DT::datatable(rv$maf_obj@data %>% mutate(DOMAINS=str_trim(DOMAINS)) %>% select(-DOMAINS) %>%
                              filter(Hugo_Symbol==input$gene, Variant_Classification %in% input$var_class, Variant_Type %in% input$var_type),
                          extensions = c("Buttons"), fillContainer = T,
                          rownames = FALSE,
                          options = list(paging = TRUE, scrollX=TRUE, ordering = TRUE, pageLength = 5,
                                         scrollY=TRUE, dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')),
                          class = "display")}, server = FALSE)

    })
    
    #Sample view features
    observeEvent(input$sample, {
        #rainfall plot
        # output$rainplot<-renderPlot({rainfallPlot(maf = read.maf(rv$maf %>%
        #                                                              filter(Variant_Classification %in% input$var_class, 
        #                                                                     Variant_Type %in% input$var_type,
        #                                                                     Tumor_Sample_Barcode == input$sample)), 
        #                                           detectChangePoints = TRUE,
        #                                           pointSize = 0.4)
        # })
        #VAF plot
        output$vafplot<-renderPlot({plotVaf(maf = read.maf(rv$maf %>%
                                                                     filter(Variant_Classification %in% input$var_class, 
                                                                            Variant_Type %in% input$var_type,
                                                                            Tumor_Sample_Barcode == input$sample)), 
                                                  top = 20)
        })
        #tumor heterogeneity
        # sample_het<-reactive({inferHeterogeneity(maf = read.maf(rv$maf %>%
        #                                                   filter(Variant_Classification %in% input$var_class, 
        #                                                          Variant_Type %in% input$var_type,
        #                                                          Tumor_Sample_Barcode == input$sample)),
        #                                tsb = input$sample)}),
        output$hetplot<-renderPlot({plotClusters(clusters = inferHeterogeneity(maf = read.maf(rv$maf %>%
                                                                                                  filter(Variant_Classification %in% input$var_class, 
                                                                                                         Variant_Type %in% input$var_type,
                                                                                                         Tumor_Sample_Barcode == input$sample)),
                                                                               tsb = input$sample))})
        output$vaf_cluster_summ<-DT::renderDataTable({
            DT::datatable(inferHeterogeneity(maf = read.maf(rv$maf %>% filter(Variant_Classification %in% input$var_class, 
                                                                              Variant_Type %in% input$var_type,
                                                                              Tumor_Sample_Barcode == input$sample)),
                                             tsb = input$sample)$clusterMeans,
                          extensions = c("Buttons"), fillContainer = T,
                          rownames = FALSE,
                          options = list(paging = TRUE, scrollX=TRUE, ordering = TRUE, pageLength = 5,
                                         scrollY=TRUE, dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')),
                          class = "display")}, server = FALSE)
                          #})
    })
                 
                 
    
}

# Run the application 
shinyApp(ui = ui, server = server)
