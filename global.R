library(tidyverse)
library(readxl)
library(shiny)
library(ggplot2)
library(DESeq2)
library(reshape2)
library(calibrate)
options(stringsAsFactors=FALSE)

getres<-function(countfile, phenofile, comparison, AdjustedCutoff=0.05, FCCutoff=0.5){
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
                                design=as.formula(paste0("~", comparison)))
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
  
  return(list(results=res,
              normal=normalized_dat,
              plot=p))
}

