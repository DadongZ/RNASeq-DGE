# Define Server ----
server <- function(input, output) {
  ##main results output
  datobj <- reactive({
    req(input$file1)
    req(input$file2)
    count <- read_excel(input$file1$datapath)
    pheno <- read_excel(input$file2$datapath)
    return(list(counts=count, 
                pheno=pheno))
  })
  
  ##input tab
  ### matrix file

  output$matrix <- renderTable({
    if(input$disp == "head") {
      return(head(datobj()[["counts"]]))
    }
    else {
      return(datobj()[["counts"]])
    }
  })
  
  url1 <- a("Count matrix example", href="https://github.com/DadongZ/RNASeqDGE")
  output$mexample <- renderUI({
    tagList(url1)
  })
  
  ### pheno file 
  url2 <- a("Manifest example", href="https://github.com/DadongZ/RNASeqDGE")
  output$pexample <- renderUI({
    tagList( url2)
  })
  
  readme <- a("README", href="https://github.com/DadongZ/RNASeqDGE/blob/master/README.md")
  output$README <- renderUI({
    tagList("Need help? ", readme)
  })
  
  issue <- a("issues", href="https://github.com/DadongZ/RNASeqDGE/issues")
  output$issue <- renderUI({
    tagList("Please report issues at: ", issue)
  })
  
  output$pdat <- renderTable({
    if(input$disp == "head") {
      return(datobj()[["pheno"]])
    }
    else {
      return(datobj()[["pheno"]])
    }
  })
  

  resobj <- reactive({
    res<-getres(input$file1$datapath, input$file2$datapath, comparison=input$design)
    return(list(normal=res[["normal"]],
                results=res[["results"]],
                volcano=res[["plot"]]))
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


