
#Use Enrichr from R 


install.packages("enrichR")
library(devtools)
install_github("wjawaid/enrichR")
install.packages("enrichR")
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

#Put together a list of libraries you'd like to look at - can just be one or two libraries for thesis, whatever you've been using in the web version of enrichr
#For example - Gene Ontology datasets, KEGG, WikiPathways, etc...

libraries <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", "KEGG_2021_Human", "CellMarker_Augmented_2021", "UK_Biobank_GWAS_v1", "PheWeb_2019", "BioPlanet_2019", "WikiPathways_2019_Mouse", "PanglaoDB_Augmented_2021")

#Prepare a list of genes 
#You should have upreg and downnreg lists already ready from your DESeq code, but in case you need to make those dataframes again:
 # upreg_FSS_0 <- subset(sig_FSS_vs_CNTL_0, sig_FSS_vs_CNTL_0$log2FoldChange > 0)
 # downreg_FSS_0 <- subset(sig_FSS_vs_CNTL_0, sig_FSS_vs_CNTL_0$log2FoldChange < 0)
cluster_markers <- subset(nts_seur_markers, nts_seur_markers$cluster == "0")

#Gene list that will be entered into enrichr (just change out for every list and re-run the code after saving the results) 
genelist <- cluster_markers$gene

#function for putting enriched terms into dataframe
if (getOption("enrichR.live")) {
    enrichRLive <- TRUE
    dbs <- listEnrichrDbs()
    if(is.null(dbs)) enrichRLive <- FALSE
    dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021",
              "GO_Biological_Process_2021", "PanglaoDB_Augmented_2021")
     enriched <- enrichr(genelist, dbs)
     if (enrichRLive) printEnrich(enriched)
}

#Create enrichment object out of input gene list and selected enrichr libraries to query 
enrichment <- enrichr(genes = genelist, databases = libraries)

#Use enrichrdf() and plotenrichr() to pull out top n terms from category of your choice by p-value, then plot in bar graph
plotdata <- enrichrdf(enrichment$PanglaoDB_Augmented_2021, n=6, database = "GO")
#If you don't want GO ID numbers to show on the plot, you can put:
#enrichrdf(enrichment$GO_Biological_Process_2018, n=10, database = "GO")
# There are similar options for KEGG ("KEGG"), Wikipathways ("wiki"), and reactome ("reactome") but I made them specifically for older versions so idk how it will fit with the updated libraries
plotenrichr(plotdata, "Cell-type Cluster 25", "#adadad")

#Code for making tables (you can cut this shorter by only keeping the libraries you will use in the thesis)

# X_BioPlanet <- if (websiteLive) enrichment[["BioPlanet_2019"]]
# X_BiolProc <- if (websiteLive) enrichment[["GO_Biological_Process_2021"]]
# X_MolFun <- if (websiteLive) enrichment[["GO_Molecular_Function_2021"]]
# X_CellComp <- if (websiteLive) enrichment[["GO_Cellular_Component_2021"]]
# X_KEGGhum <- if (websiteLive) enrichment[["KEGG_2021_Human"]]
# X_WikiPaths <- if (websiteLive) enrichment[["WikiPathways_2019_Mouse"]]
# X_CellMarker <- if (websiteLive) enrichment[["CellMarker_Augmented_2021"]]
X_Panglao <- if (websiteLive) enrichment[["PanglaoDB_Augmented_2021"]]

#save the table as a .csv to your home directory 
# write.csv(X_BioPlanet, "/stor/home/eg33747/NTS_snRNA/cluster0_BioPlanet.csv")
# write.csv(X_BiolProc, "/stor/home/eg33747/Becker_Tagseq_BNST/Time0/DEG_Enrichment/table_FSS_0_down_BiolProc.csv")
# write.csv(X_MolFun, "/stor/home/eg33747/Becker_Tagseq_BNST/Time0/DEG_Enrichment/table_FSS_0_down_MolFun.csv")
# write.csv(X_CellComp, "/stor/home/eg33747/Becker_Tagseq_BNST/Time0/DEG_Enrichment/table_FSS_0_down_CellComp.csv")
# write.csv(X_KEGGhum, "/stor/home/eg33747/Becker_Tagseq_BNST/Time0/DEG_Enrichment/table_FSS_0_down_KEGGhum.csv")
# write.csv(X_WikiPaths, "/stor/home/eg33747/Becker_Tagseq_BNST/Time0/DEG_Enrichment/table_FSS_0_down_WikiPaths.csv")

#generate gene list according to specific cluster
cluster_markers <- subset(cleaned_markers, cleaned_markers$cluster == "Chond")
genelist <- cluster_markers$gene
enrichment <- enrichr(genes = genelist, databases = libraries)
#write.csv(enrichment$PanglaoDB_Augmented_2021, "/stor/home/eg33747/NTS_snRNA/cleaned_cluster_markers_soupX_harmony/cluster22_pangalao.csv")
plotdata <- enrichrdf(enrichment$PanglaoDB_Augmented_2021, n=6, database = "GO")
plotenrichr(plotdata, "Cell-type Chond", "#adadad")
plotdata <- enrichrdf(enrichment$GO_Biological_Process_2021, n=6, database = "GO")
plotenrichr(plotdata, "Biol Proc Chond", "#adadad")

new.cluster.ids <- c( "Neu1",	"Neu2",	"Soup1",	"Neu3", "Oligo1", "Ast1", "Oligo2", "Soup2", "Endo1", "Ast2", "Inter1", "Endo2", "Inter2", "OPCs", "MG", "Pod", "Oligo3", "Epend", "Soup3", "Inter3", "Chond", "Neu4", "Neu5")
old.cluster.ids <- c( "0",	"1",	"2",	"3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22")
names(new.cluster.ids) <- levels(x_cleaned_all)
x_cleaned_all <- RenameIdents(x_cleaned_all, new.cluster.ids)

#ffab61 CIE
#78ffd0 CIE FSS
#c6aefc FSS




