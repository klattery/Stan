if (exists("env_shiny")) rm(env_shiny)
env_shiny <- new.env(parent = emptyenv())

# Define UI for data upload app ----
env_shiny$ui_1 <- fluidPage(
  # App title ----
  titlePanel("Upload files for HB into R"),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      # Input: Select a file ----
      fileInput("file_conjoint", "File with Choices (csv or rds)",
                multiple = FALSE,
                accept = c(".csv",".rds")),
      fileInput("file_specs","File with Excel Specs",
                multiple = FALSE,
                accept = c(".xls",".xlsx")),
      checkboxInput("use_cov", "Use Covariates (Optional)"),
      conditionalPanel(
        condition = "input.use_cov == true",
        fileInput("file_cov", "Optional: File with Covariates",
                  multiple = FALSE,
                  accept = c(".csv",".rds"))),
      tags$hr(),
      # Checkboxes
      checkboxInput("add_none", "Add None Variable when All Attributes are 0", TRUE),
      checkboxInput("check_collinearity", "Check Collinearity", TRUE),
      checkboxInput("est_aggmodel", "Estimate Aggregate Model", TRUE),
      checkboxInput("auto_stop", "Stop Instance after Running", TRUE),
      checkboxInput("est_EB", "Estimate Empirical Bayes with Draws", FALSE),
      textInput("out_prefix", "Text you want to prefix output", value = "MyCBC", width = NULL, placeholder = NULL),
      #textInput("password", "Enter SKIMVERSE password to automatically STOP instance after running:", value = "", width = NULL, placeholder = NULL),
      tags$hr(),
      actionButton("setup_ready","Save Changes to R & Exit Upload", class = "btn-primary")
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      h4("Choice Data (first 10 rows):"),
      h6(tableOutput("file1")),
      h4("Attribute Coding from Excel:"),
      h6(tableOutput("file2")),
      h4("Optional Covariate Data (first 10 rows):"),
      h6(tableOutput("file3"))
    )
  )
)


env_shiny$server_1 <- function(input, output) {
  options(shiny.maxRequestSize=1000*1024^2)
  data1 <- reactive({
    if (is.null(input$file_conjoint)) {
      result <- NULL
    } else {
      data_file <- input$file_conjoint$datapath
      ftype <- toupper(substr(data_file,nchar(data_file)-2, nchar(data_file)))
      if (ftype %in% c("CSV", "RDS")){
        if (ftype == "CSV") {result <- read.csv(data_file, as.is=TRUE, check.names = FALSE)}
        if (ftype == "RDS") {result <- readRDS(data_file)}
      }  
    } 
    return(result)
  })
  
  data2 <- reactive({
    if (is.null(input$file_specs)) {
      result <- NULL
    } else{
      result <- list()
      result$specs_att_coding <- data.frame(read_xlsx(input$file_specs$datapath,
                                                      sheet = "Att_Coding", col_types = c("text","text","numeric")))
      result$specs_pair_constraints <- data.frame(read_xlsx(input$file_specs$datapath,
                                                            sheet = "Pair_Constraints", col_types = c("text","numeric","numeric")))
      result$specs_cov_coding <- data.frame(read_xlsx(input$file_specs$datapath,
                                                      sheet = "Cov_Coding", col_types = c("text","text")))
    } 
    return(result)
  })
  
  data3 <- reactive({
    if (is.null(input$file_cov)) {
      result <- NULL
    } else {
      data_file <- input$file_cov$datapath
      ftype <- toupper(substr(data_file,nchar(data_file)-2, nchar(data_file)))
      if (ftype %in% c("CSV", "RDS")){
        if (ftype == "CSV") {result <- read.csv(data_file, as.is=TRUE, check.names = FALSE)}
        if (ftype == "RDS") {result <- readRDS(data_file)}
      }
    }   
    return(result)
  })
  
  output$file1 <- renderTable({data1()[1:10,]})
  output$file2 <- renderTable({data2()$specs_att_coding})
  output$file3 <- renderTable({data3()[1:10,]})
  
  observeEvent(input$setup_ready, {
    .GlobalEnv$control_code$add_none <- input$add_none # To be like Sawtooth
    .GlobalEnv$control_code$check_collinearity <- input$check_collinearity # May take a few minutes to check if your coding is deficient
    .GlobalEnv$control_code$est_aggmodel <- input$est_aggmodel # Option to estimate aggregate model
    .GlobalEnv$control_code$est_EB <- input$est_EB # Option to estimate aggregate model
    .GlobalEnv$control_code$auto_stop <- input$auto_stop
    .GlobalEnv$out_prefix <- input$out_prefix
    if (!is.null(data1())) {.GlobalEnv$data_conjoint <- data1()}
    if (!is.null(data2())) {
      .GlobalEnv$specs_att_coding <- data2()$specs_att_coding
      .GlobalEnv$specs_pair_constraints <- data2()$specs_pair_constraints
      .GlobalEnv$specs_cov_coding <- data2()$specs_cov_coding
    }
    if (input$use_cov){
      .GlobalEnv$data_cov <- data3()
    } else .GlobalEnv$data_cov <- NULL
    stopApp()
  })
}
attach(env_shiny)