##Fri 17 Jan 2020 10:27:22 AM EST
##Dadong Zhang
##https://github.com/DadongZ/RNASeqDGE

rm(list=ls())
library(tidyverse)
library(readxl)
library(shiny)
library(ggplot2)
library(DESeq2)
library(reshape2)
library(calibrate)
options(stringsAsFactors=FALSE)

countfile<-"example_rnaseq_count_matrix.xlsx"
phenofile<-"pheno_data.xlsx"

##Local func
getres<-function(countfile, phenofile, AdjustedCutoff=0.05, FCCutoff=0.5){
  count <- read_excel(countfile)%>%
    column_to_rownames(var="Gene")
  pheno <- read_excel(phenofile)%>%
    column_to_rownames(var="Samples")
  
  if (!identical(colnames(count), rownames(pheno))) {
    message("Sample names in count matrix must be identical to the sample names in pheno data")
    break
  } 
  
  dds <- DESeqDataSetFromMatrix(countData=count, 
                                colData=pheno, 
                                design=~Condition)
  dds<-dds[rowSums(counts(dds))>0, ]
  normalized_dat<-assay(rlog(dds, blind=TRUE))%>%data.frame%>%rownames_to_column(var="Gene")
  dds<-DESeq(dds)
  res<-results(dds)%>%data.frame%>%
    rownames_to_column(var="Gene")%>%
    filter(complete.cases(.))%>%
    mutate(log10padj=-log10(padj))%>%
    arrange(pvalue)
  res<-res%>%
    mutate(Significance=ifelse(log2FoldChange < (-FCCutoff) & padj < AdjustedCutoff, "FC_FDR_Down",
                               ifelse(log2FoldChange > FCCutoff & padj < AdjustedCutoff, "FC_FDR_Up", "NS")))
  
  p <-ggplot(res, aes(x=log2FoldChange, y=log10padj)) +
    geom_point(aes(color=factor(Significance)), alpha=1/2, size=2) +
    theme_bw(base_size=16) +
    xlab(bquote(~Log[2]~ "fold change")) +
    ylab(bquote(~-Log[10]~adjusted~italic(P))) +
    geom_vline(xintercept=c(-FCCutoff, FCCutoff), linetype="longdash", colour="black", size=0.4) +
    geom_hline(yintercept=-log10(AdjustedCutoff), linetype="longdash", colour="black", size=0.4)
  
  return(list(counts=count%>%rownames_to_column(var="Gene"),
              pdat=pheno%>%rownames_to_column(var="Samples"),
              results=res,
              normal=normalized_dat,
              plot=p))
}

# Define UI for data upload app ----
ui <- fluidPage(
  tabsetPanel(
    tabPanel("Input", fluid = TRUE,
             
             # App title ----
             titlePanel("Uploading Files"),
             
             # Sidebar layout with input and output definitions ----
             sidebarLayout(
               
               # Sidebar panel for inputs ----
               sidebarPanel(
                 
                 # Input: Select a file ----
                 fileInput("file1", "Count matrix File (.xlsx)",
                           multiple = TRUE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 fileInput("file2", "Manifest File (.xlsx)",
                           multiple = TRUE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 
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
                              selected = "head")
                 
               ),
               # Main panel for displaying outputs ----
               mainPanel(
                 
                 # Output: Data file ----
                 tableOutput("contents"),
                 tableOutput("contents2")
               )
             )
    ),
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
# Define server logic to read selected file ----
server <- function(input, output) {
  ##main results output
  resobj <- reactive({
    req(input$file1)
    req(input$file2)
    res<-getres(input$file1$datapath, input$file2$datapath)
    return(list(counts=res[["counts"]], 
                pheno=res[["pdat"]], 
                normal=res[["normal"]],
                results=res[["results"]],
                volcano=res[["plot"]]))
  })
  ##input tab
  
  ### matrix file
  output$contents <- renderTable({
    if(input$disp == "head") {
      return(head(resobj()[["counts"]]))
    }
    else {
      return(resobj()[["counts"]])
    }
  })
  
  ### pheno file 
  output$contents2 <- renderTable({
    
    
    if(input$disp == "head") {
      return(head(resobj()[["pheno"]]))
    }
    else {
      return(resobj()[["pheno"]])
    }
  })
  
  ##results panel
  todowndat <- reactive({
    switch(input$results,
           "Results" = resobj()[["results"]],
           "Normalized matrix" = resobj()[["normal"]]
    )
  })
  
  output$table <- renderTable({
    todowndat()
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$results, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(todowndat(), file, row.names = FALSE)
    }
  ) 
  
  ##plot panel
  output$plot1 <- renderPlot({
    resobj()[["volcano"]]
  })
  
  output$brush_info <- renderPrint({
    showdf<-resobj()[["results"]]%>%dplyr::select(Gene, log2FoldChange, pvalue, padj, log10padj) 
    brushedPoints(showdf, input$plot1_brush)
  })
}
# Run the app ----
shinyApp(ui, server)
