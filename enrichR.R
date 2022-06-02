
setwd("/stor/scratch/WCAAR/emily_scratch/")
outdir <- "/stor/home/eg33747/NTS_snRNA/"

# Convenience functions
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(outdir, name, ".", type),
        width = width, height = height, units = "in", res = 200)
  } else {
    pdf(paste0(outdir, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}

#Use Enrichr for R 

install.packages("enrichR")
#install_github("wjawaid/enrichR") if the above doesn't work try this
library(devtools)
library("enrichR")



########## following code is to make 2 functions for helping make bargraphs with enrichr output ##############################
########### you only have to run it once and then you'll have the functions in your environment for use #################


enrichrdf <- function(x, n=NULL, database=NULL) {
  cc1 <- x  #rename output dataframe from enrichr function
  colnames(cc1) <- c("Term", "overlap", "P.value", "padj", "old", "old", "zscore", "comb", "genes") #take terms and p-values, then subset only significant 
  df <- as.data.frame(cc1$Term)
  df$P.value <- cc1$P.value
  colnames(df) <- c("Term", "P-value")
  df$`P-value` <- -log10(df$`P-value`) #Take top number and sort by -log10 pvalue 
  df <- df[order(df$`P-value`, decreasing = TRUE),]
  
  if ( is.null(n) )
    n <- 5
  
  df <- as.data.frame(df[1:n,])
  df <- df[order(df$`P-value`, decreasing = FALSE),]
  if (is.null(database))
    database <- ""
  if (database == "GO")
    df$Term <- gsub("\\(GO:.*)", "", df$Term, fixed = F)  #To remove GO ID numbers from the terms (for GO categories)
  if (database == "KEGG")
    df$Term <- gsub("_Homo sapiens_hsa.*", "", df$Term, fixed = F)
  if (database =="wiki")
    df$Term <- gsub("_WP.*", "", df$Term, fixed = F)
  if (database =="reactome")
    df$Term <- gsub("_R-HSA-.*", "", df$Term, fixed = F)
  df$Term <- factor(df$Term, levels=unique(df$Term)) #Turn "Term" column into an ordered factor 
  dfbp <- df
  return(dfbp)
}

library(ggplot2)
plotenrichr <- function(df, title, barcolor) {
  n <- nrow(df)
  ggplot(data = df, aes(x = Term, y = `P-value`)) +
    geom_hline(yintercept = 1.3)+
    geom_bar( stat = 'identity', color = "white", fill = barcolor ) +
    annotate("text", x = 1:n, y = (rep(.1, times =n)), 
             label = df$Term, hjust = 0) + 
    coord_flip() + 
    theme_classic() + 
    scale_y_continuous(expand = c(0,0)) + 
    ylab("-log10(p-value)") +
    ggtitle(paste(title)) +
    theme(axis.text.y = element_blank (), axis.ticks.y = element_blank())
}


########## end of function making ##############################
###########  ################# ###########  ################# 


#Now to use enrichr

#Put together a list of libraries you'd like to look at
#For example - Gene Ontology datasets, KEGG, WikiPathways, etc...
#to see available databases, run the following:
listEnrichrSites()
setEnrichrSite("Enrichr") 
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

#select databases from the dbs table that you are interested in and set as libraries

libraries <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", 
               "KEGG_2021_Human", "CellMarker_Augmented_2021", "UK_Biobank_GWAS_v1", "PheWeb_2019", 
               "BioPlanet_2019", "WikiPathways_2019_Mouse", "PanglaoDB_Augmented_2021")

#generate gene list
cluster_markers <- subset(markers, markers$cluster == "1")
genelist <- cluster_markers$gene

#run enrichment
enrichment <- enrichr(genes = genelist, databases = libraries)

#make plots (I am usually just interested in cell-type and biol process for single cell but you can add as many plots as you want)
 plotdata1 <- enrichrdf(enrichment$PanglaoDB_Augmented_2021, n=6, database = "GO")
 plot1 <- plotenrichr(plotdata1, "Cell-type 1", "#adadad")
 plotdata2 <- enrichrdf(enrichment$GO_Biological_Process_2021, n=6, database = "GO")
 plot2 <- plotenrichr(plotdata2, "Biol Proc 1", "#adadad")
 plot1 + plot2

#function for putting enriched terms into dataframe that you can then export to excel
if (getOption("enrichR.live")) {
  enrichRLive <- TRUE
  dbs <- listEnrichrDbs()
  if(is.null(dbs)) enrichRLive <- FALSE
  dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021",
           "GO_Biological_Process_2021", "PanglaoDB_Augmented_2021")
  enriched <- enrichr(genelist, dbs)
  if (enrichRLive) printEnrich(enriched)
}
write.csv(enriched$PanglaoDB_Augmented_2021, "/stor/home/eg33747/NTS_snRNA/cleaned_cluster_markers_soupX_harmony/cluster1_pangalao.csv")

#plots that show # of genes in each term and use color to show p-value
p1<- plotEnrich(enrichment$PanglaoDB_Augmented_2021,
           showTerms = 20,
           numChar = 40,
           y = "Count",
           orderBy = "P.value")

p2<- plotEnrich(enrichment$GO_Biological_Process_2021,
                showTerms = 20,
                numChar = 40,
                y = "Count",
                orderBy = "P.value")
p3<- p1+p2

SaveFigure(p3, "cluster1", width = 16, height = 4)


