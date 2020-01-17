##Fri 17 Jan 2020 10:27:22 AM EST
##Dadong Zhang
##https://github.com/DadongZ/RNASeq-DGE

rm(list=ls())
library(tidyverse)
library(readxl)
library(shiny)
library(ggplot2)
library(DESeq2)
library(reshape2)
library(calibrate)
options(stringsAsFactors=FALSE)

##Local func
getres<-function(countfile, phenofile, AdjustedCutoff=0.05, FCCutoff=1){
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
                                  design=~Dose)
    dds<-dds[rowSums(counts(dds))>0, ]
    normalized_dat<-rlog(dds, blind=TRUE)
    dds<-DESeq(dds)
    res<-results(dds)%>%data.frame%>%
        rownames_to_column(var="Gene")%>%
        filter(complete.cases(.))%>%
        mutate(log10padj=-log10(padj))%>%
        arrange(pvalue)
    res$Significance <- "NS"
    res$Significance[(abs(res$log2FoldChange) > FCCutoff)] <- "FC"
    res$Significance[(res$padj<AdjustedCutoff)] <- "FDR"
    res$Significance[(res$padj<AdjustedCutoff) & (abs(res$log2FoldChange)>FCCutoff)] <- "FC_FDR"
    res$Significance <- factor(res$Significance, levels=c("NS", "FC", "FDR", "FC_FDR"))

    p <-ggplot(res, aes(x=log2FoldChange, y=log10padj)) +
        geom_point(aes(color=factor(Significance)), alpha=1/2, size=2) +
        theme_bw(base_size=16) +
        xlab(bquote(~Log[2]~ "fold change")) +
        ylab(bquote(~-Log[10]~adjusted~italic(P))) +
        geom_vline(xintercept=c(-FCCutoff, FCCutoff), linetype="longdash", colour="black", size=0.4) +
        geom_hline(yintercept=-log10(AdjustedCutoff), linetype="longdash", colour="black", size=0.4)

    return(list(Results=res, plot=p))
}

#DEG analysis
countfile <-"example_rnaseq_count_matrix.xlsx"
phenofile <-"pheno_data.xlsx"
reslst<-getres(countfile, phenofile)

#plot
getshiny<-function(){
  ui <- fluidPage(
    fluidRow(
      column(width = 8,
        plotOutput("plot1", height = 900,
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

server <- function(input, output) {
    output$plot1 <- renderPlot({
        reslst[["plot"]]
    })
    output$brush_info <- renderPrint({
        showdf<-reslst[["Results"]]%>%dplyr::select(Gene, log2FoldChange, pvalue, log10padj) 
        brushedPoints(showdf, input$plot1_brush)
    })
}
   shinyApp(ui, server)
}

getshiny()

