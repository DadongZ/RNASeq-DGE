ui <- fluidPage(
  tabsetPanel(
    #Input tabPanel
    tabPanel("Input", fluid = TRUE,
             
             # App title ----
             titlePanel("Upload data"),
             
             # Sidebar layout with input and output definitions ----
             sidebarLayout(
               
               # Sidebar panel for inputs ----
               sidebarPanel(
                 
                 uiOutput("mexample"),
                 
                 # Input: Select a file ----
                 fileInput("file1", "Count matrix File (.xlsx)",
                           multiple = TRUE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),

                 uiOutput("pexample"),
                 fileInput("file2", "Manifest File (.xlsx)",
                           multiple = TRUE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 
                 textInput("design", "Column name for analysis", " "),

                 # Horizontal line ----
                 tags$hr(),
                 
                 # Input: Checkbox if file has header ----
                 checkboxInput("header", "Header", TRUE),
                 
                 # Horizontal line ----
                 tags$hr(),
                 
                 # Input: Select number of rows to display ----
                 radioButtons("disp", "Display",
                              choices = c(Head = "head",
                                          All = "all"),
                              selected = "head"),
                 
                 uiOutput("README"),
                 
                 uiOutput("issue")
                 
               ),
               # Main panel for displaying outputs ----
               mainPanel(
                 
                 # Output: Data file ----
                 tableOutput("matrix"),
                 tableOutput("pdat")
               )
             )
    ),
    
    #Results tabPanel
    tabPanel("Results", fluid = TRUE,
             # App title ----
             titlePanel("Download results"),
             
             # Sidebar layout with input and output definitions ----
             sidebarLayout(
               
               # Sidebar panel for inputs ----
               sidebarPanel(
                 
                 # Input: Choose dataset ----
                 selectInput("results", "Choose a dataset:",
                             choices = c("Results", "Normalized matrix")),
                 
                 # Button
                 downloadButton("downloadData", "Download")
                 
               ),
               
               # Main panel for displaying outputs ----
               mainPanel(
                 
                 tableOutput("table")
                 
               )
               
             )             
             
    ),
    #plot tabPanel 
    tabPanel("Plots", fluid = TRUE,
             fluidRow(
               column(width = 8,
                      plotOutput("plot1", height = 800,
                                 # Equivalent to: click = clickOpts(id = "plot_click")
                                 click = "plot1_click",
                                 brush = brushOpts(
                                   id = "plot1_brush"
                                 )
                      )
               ),
               column(width = 4,
                      h4("Brushed points"),
                      verbatimTextOutput("brush_info")
               )
             )
    )
  )
)


