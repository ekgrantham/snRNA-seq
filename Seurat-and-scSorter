#Analyzing samples 1-3 and 5-7 - naive nts samples for basic single cell profiling of the brain region

#preparing the seurat object ####
#running seurat pipeline on the  aggregate 6 samples feb 16, 2022
setwd("/stor/scratch/WCAAR/emily_scratch/")
outdir <- "/stor/home/eg33747/NTS_snRNA/"
library (Seurat)
library(dplyr)
library(patchwork)
library(reshape2)
library(ggplot2)
library(sctransform)
#the seurat tutorial https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

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

#in command line untar the filtered_feature_bc_matrix.tar.gz
# tar -zxvf filtered_feature_bc_matrix.tar.gz

# Load the NTS nuclei
nts_nuclei_full<- Read10X(data.dir = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S1-3_4-7_Naive_NTS/count/")
# Initialize the Seurat object with the raw (non-normalized data).
nts_seurat_full <- CreateSeuratObject(counts = nts_nuclei_full, project = "NTS_nuclei_6samps", min.cells = 3, min.features = 200)
nts_seurat_full
# An object of class Seurat 
# 24360 features across 20256 samples within 1 assay 
# Active assay: RNA (24360 features, 0 variable features)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats

#taking soupX cleaned Seurat object and continuing with Seurat pipeline - also look at ribosomal RNA content
cleaned_all[["percent.mt"]] <- PercentageFeatureSet(cleaned_all, pattern = "^mt-")
cleaned_all[["percent.rb"]] <- PercentageFeatureSet(cleaned_all, pattern = "^Rp[sl][[:digit:]]")

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(cleaned_all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cleaned_all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(cleaned_all, feature1 = "nCount_RNA", feature2 = "percent.rb", pt.size=0.1)
plot4 <- FeatureScatter(cleaned_all, feature1 = "nFeature_RNA", feature2 = "percent.rb", pt.size=0.1)
plot1 + plot2
plot3 + plot4


# Visualize QC metrics as a violin plots
VlnPlot(cleaned_all, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,pt.size = 0)
VlnPlot(cleaned_all, features = c("percent.mt", "percent.rb"), ncol = 2,pt.size = 0)

#Filter cells with feature counts between 200-7500 and greater than 10% mt counts based on violin plots
#save unfiltered seurat object in case you change your mind about filtering parameters
#unfiltered_seur <- cleaned_all
#unfiltered_seur_man <- cleaned_all

#First calculate some basic statistics on the various QC parameters, which can be helpful for choosing cutoffs. For example:
min <- min(cleaned_all@meta.data$nFeature_RNA);
m <- median(cleaned_all@meta.data$nFeature_RNA)
max <- max(cleaned_all@meta.data$nFeature_RNA)    
s <- sd(cleaned_all@meta.data$nFeature_RNA)
min1 <- min(cleaned_all@meta.data$nCount_RNA)
max1 <- max(cleaned_all@meta.data$nCount_RNA)
m1 <- mean(cleaned_all@meta.data$nCount_RNA)
s1 <- sd(cleaned_all@meta.data$nCount_RNA)
Count93 <- quantile(cleaned_all@meta.data$nCount_RNA, 0.93) # calculate value in the 93rd percentile
print(paste("Feature stats:",min,m,max,s));
print(paste("UMI stats:",min1,m1,max1,s1,Count93));

cleaned_all <- subset(x = cleaned_all, subset = nFeature_RNA > 300  & nCount_RNA < Count93 & percent.mt < 5)

#visualize QC parameters again post-filtering
VlnPlot(cleaned_all, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,pt.size = 0)
VlnPlot(cleaned_all, features = c("percent.mt", "percent.rb"), ncol = 2,pt.size = 0)


      # Original Normalization pipeline (use scTransform now)
      # #seur<- NormalizeData(seur, normalization.method = "LogNormalize", scale.factor = 10000)
      # 
      # #identification of highly variable features
      # #nts_seurat_full<- FindVariableFeatures(nts_seurat_full, selection.method = "vst", nfeatures = 2000)
      # #seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 2000)
      # # Identify the 10 most highly variable genes
      # #top10 <- head(VariableFeatures(seur), 10)
      # #head(VariableFeatures(seur), 30)
      # # [1] "9630013A20Rik" "Gm32647"       "Dnah12"        "Il1rapl2"      "Flt1"          "Zfp804b"       "Slco1a4"      
      # # [8] "Otx2os1"       "Hmcn1"         "Cntn5"         "Adarb2"        "Prkg1"         "Mecom"         "Cpne4"        
      # # [15] "Vtn"           "Itpr1"         "Slit3"         "Unc5cl"        "Gmnc"          "Schip1"        "6330411D24Rik"
      # # [22] "Acta2"         "Atp13a5"       "Cfap299"       "Galntl6"       "Tafa1"         "Gm34544"       "Htr2c"        
      # # [29] "6530403H02Rik" "Fbxl7"
      # 
      # # plot variable features with and without labels
      # VariableFeaturePlot(seur)
      # LabelPoints(plot = plot1, points = top10, repel = TRUE)
      # dev.off()
      # 
      # #scale the data
      # all.genes <- rownames(seur)
      # seur<- ScaleData(seur, features = all.genes)

#instead of normalize, scale, and find variable features, use scTransform

# store mitochondrial percentage in object meta data
#seur <- PercentageFeatureSet(seur, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
#seur <- SCTransform(seur, vars.to.regress = "percent.mt", verbose = FALSE)

#use updated glmGamPoi for scTransform which improves the speed of the learning procedure
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# 
# BiocManager::install("glmGamPoi")
# seur <- SCTransform(seur, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

#perform linear dimensional reduction
# These are now standard steps in the Seurat workflow for visualization and clustering
# seur <- RunPCA(seur, verbose = FALSE)
# seur <- RunUMAP(seur, dims = 1:30, verbose = FALSE)
# 
# seur <- FindNeighbors(seur, dims = 1:30, verbose = FALSE)
# seur <- FindClusters(seur, verbose = FALSE)
# DimPlot(seur, label = TRUE) + NoLegend()

#The following code replicates the full end-to-end workflow for Seurat, in a single command:
cleaned_all <- cleaned_all %>%
   PercentageFeatureSet(pattern = "^mt-", col.name = "percent.mt") %>%
   SCTransform(vars.to.regress = "percent.mt") %>%
   RunPCA() %>%
   FindNeighbors(dims = 1:30) %>%
   RunUMAP(dims = 1:30) %>%
   FindClusters()

#Look at variable genes per PC (first 5)
PClist_1 <- names(sort(Loadings(object=cleaned_all, reduction="pca")[,1], decreasing=TRUE));
PClist_2 <- names(sort(Loadings(object=cleaned_all, reduction="pca")[,2], decreasing=TRUE));
PClist_3 <- names(sort(Loadings(object=cleaned_all, reduction="pca")[,3], decreasing=TRUE));
PClist_4 <- names(sort(Loadings(object=cleaned_all, reduction="pca")[,4], decreasing=TRUE));
PClist_5 <- names(sort(Loadings(object=cleaned_all, reduction="pca")[,5], decreasing=TRUE));

#feature.pal = rev(colorRampPalette(brewer.pal(11,"Spectral"))(50)); # a useful color palette
pdf(sprintf("%s/umap.%d.colorby.UMI.pdf", outdir, 30), width = 10, height = 8);
fp <- FeaturePlot(object = cleaned_all, features = c("percent.rb"), pt.size=0.1, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank()); # the text after the ‘+’ simply removes the axis using ggplot syntax
print(fp);
dev.off();

# Look at cluster IDs of the first 5 cells, and number of clusters identified
head(Idents(cleaned_all), 5)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
orig.idents <- Idents(cleaned_all)

DimPlot(cleaned_all, reduction = "umap", group.by = "orig.ident")

DimPlot(object = cleaned_all, reduction = "pca", pt.size = .1, group.by = "orig.ident")
VlnPlot(object = cleaned_all, features = "PC_1", group.by = "orig.ident", pt.size = .1)

#Clustering is happening by sample, need to remove batch effects - run through Harmony script
install.packages("harmony")
library(harmonny)

#add the sex info
barcode_info <- cleaned_all@assays$SCT@counts@Dimnames[[2]]
barcode_info <- as.character(barcode_info)

group <-ifelse(grepl(pattern = 's1|s2|s3' ,barcode_info)==TRUE,"Male","Female")
group <- as.data.frame(group)
row.names(group) <- cleaned_all@assays$SCT@counts@Dimnames[[2]]
cleaned_all <- AddMetaData(cleaned_all, group)

DimPlot(cleaned_all, reduction = "umap", group.by = "group")

table(cleaned_all@meta.data$group)


#same with sample id
# sample_info<-SampleID
# row.names(sample_info) <- barcode_info
# x_cleaned_all<-AddMetaData(x_cleaned_all, sample_info)
# table(cleaned_all@meta.data$SampleID)
# S1_M S2_M S3_M S5_F S6_F S7_F 
# 1176 5009 2228  698 1228 3390 

cleaned_all <- RunHarmony(cleaned_all, "orig.ident", assay.use = "SCT", project.dim = F)
DimPlot(cleaned_all, reduction = "harmony", group.by = "orig.ident")
harmony_embeddings <- Embeddings(cleaned_all, 'harmony')
harmony_embeddings[1:5, 1:5]

cleaned_all <- cleaned_all %>% 
   RunUMAP(reduction = "harmony", dims = 1:30) %>% 
   FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
   FindClusters() %>% 
   identity()

DimPlot(cleaned_all, reduction = "umap", group.by = "orig.ident")
DimPlot(cleaned_all, reduction = "umap", group.by = "group")
#DimPlot(cleaned_all, reduction = "tsne") #no longer run tsne with SCTransform

#visualize QC parameters again post-filtering by cluster
VlnPlot(cleaned_all, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,pt.size = 0)
VlnPlot(cleaned_all, features = c("percent.mt", "percent.rb"), ncol = 2,pt.size = 0)

#save object
#saveRDS(cleaned_all, file = "/stor/home/eg33747/NTS_snRNA/seurat_after_soupX_manual.rds")
#saveRDS(cleaned_all, file = "/stor/home/eg33747/NTS_snRNA/seurat_after_soupX_auto.rds")
saveRDS(cleaned_all, file = "/stor/home/eg33747/NTS_snRNA/seurat_after_soupX_harmony.rds")


# find markers for every cluster compared to all remaining cells, report only the positive ones
# seur_markers <- FindAllMarkers(seur, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# seur_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
# write.csv(seur_markers,file="/stor/home/eg33747/NTS_snRNA/cluster_markers_snRNA_after_diem.csv")

#run same code on unfiltered data:
cleaned_markers <- FindAllMarkers(cleaned_all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cleaned_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(cleaned_markers,file="/stor/home/eg33747/NTS_snRNA/cleaned_markers_soupX_harmony.csv")

top5 <-cleaned_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
to_plot <- unique(top5$gene)
plot <- DotPlot(x_cleaned_all, features = to_plot) + coord_flip() + theme(axis.text.x = element_text(angle = 45))
plot
SaveFigure(plot, "dplot_top5_postmanX_harm_mtrbrm", width = 9, height = 20)

# find all markers distinguishing cluster 0 from clusters 2 and 6
# cluster0.markers <- FindMarkers(nts_seurat_full, ident.1 = 0, ident.2 = c(8, 11), min.pct = 0.25)
# head(cluster0.markers, n = 5)

#use ROC test to return classification power for individual markers
# cluster2.markers <- FindMarkers(nts_seur, ident.1 = Neu2, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
# head(cluster2.markers, n = 5)

#Select top two markers and show violin plots
#VlnPlot(nts_seurat_full, features = c("Cd14", "Cx3cr1", "Tmem119", "Tnf", "Ccr2", "Top2a", "Pf4", "Mrc1",
#                                      "S100a9"))

VlnPlot(x_cleaned_all, features = c("Nr3c2"), pt.size = 0)

FeaturePlot(cleaned_all, features = c("Slc17a6", "Slc32a1"))

table(cleaned_all@meta.data$seurat_clusters)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22 
# 2200 2071 1666 1311 1288 1148 1065  967  742  470  439  277  270  236  211  136  113   86   85   78   78   51   45 

seur_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_logFC) -> top10
DoHeatmap(seur, features = top10$gene)

#new.cluster.ids <- c( "Neu1",	"Ast1",	"Neu2",	"Neu3",	"Neu4",	"Oligo1",	"Neu5",	"Oligo2",	"Endo1",	"Neu6",	"Ast2",	"Endo2",	"Inter1",	"Oligo3",	"OPCs",	"Pod",	"Neu7",	"Oligo4",	"MG",	"Epend",	"Inter2",	"Neu8")
# names(new.cluster.ids) <- levels(cleaned_all)
# cleaned_all <- RenameIdents(cleaned_all, new.cluster.ids)

DimPlot(x_cleaned_all, reduction = "umap", label = TRUE, pt.size = 0.5)

#Remove uninteresting genes
ribGenes=grep("Rp", x = rownames(x = cleaned_all@assays$SCT@data), value = TRUE)
red_cleaned_all = cleaned_all[!rownames(cleaned_all@assays$SCT@data) %in% ribGenes,]
ribFrac = colSums(cleaned_all[ribGenes,])/colSums(cleaned_all)
red_cleaned_all@meta.data$ribFrac = ribFrac

mtGenes=grep('^mt-', x = rownames(x = red_cleaned_all@assays$SCT@data),value=TRUE)
x_cleaned_all = red_cleaned_all[!rownames(red_cleaned_all@assays$SCT@data) %in% mtGenes,]
mtFrac = colSums(red_cleaned_all[mtGenes,])/colSums(red_cleaned_all)
x_cleaned_all@meta.data$mtFrac = mtFrac

#exclude mt genes from the highly variable genes used for calculating PCs and downstream analyses
# cleaned_all@assays$SCT@var.features = setdiff(cleaned_all@assays$SCT@var.features,mtGenes)
# 
# DimPlot(red_cleaned_all, reduction = "umap", label = TRUE, pt.size = 0.5)

cleaned_markers <- FindAllMarkers(x_cleaned_all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cleaned_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(cleaned_markers,file="/stor/home/eg33747/NTS_snRNA/cleaned_markers_soupX_harmony_rbmtrm.csv")


#scSorter
install.packages('scSorter')
library(scSorter)


if (!require("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("DropletUtils")

install.packages("remotes")
remotes::install_github("MarioniLab/DropletUtils")
library(DropletUtils)

write10xCounts(x = x_cleaned_all@assays$SCT@counts, path = "/stor/home/eg33747/NTS_snRNA/cleaned_cloupe_2")
write.csv(x_cleaned_all@assays$SCT@counts, file = "/stor/home/eg33747/NTS_snRNA/cleaned_cloupe.csv")


###
#save the image
#save.image("snRNA_full_Nov21/aggr_ranger_8samps/processing8samps.RData")
###Now use scSorter 
#use the human cell marker set
#Markers from here :https://www.nature.com/articles/s41598-018-27293-5#Sec2

#gross cell types
#markers here
human_cell_type_marker<-read.csv("/stor/scratch/WCAAR/nihal_scratch/human_cell_type_markers.csv")
#convert to mouse
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# Basic function to convert mouse to human gene names
convertHumanGeneList <- function(x){
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=F)
  return(genesV2)
}  
mouse_gene<-convertHumanGeneList(human_cell_type_marker$gene)
#merge that  with the lineage data
mouse_cell_type_marker<-merge(human_cell_type_marker,mouse_gene,by.x="gene",by.y="HGNC.symbol")
colnames(mouse_cell_type_marker)[5:6]<-c("Type","Marker")
#now scSorter
#now run the scSorter pipeline
start_time <- Sys.time()
scaled_data_seur<-x_cleaned_all@assays$SCT@scale.data
sorted_seur<-scSorter(scaled_data_seur,mouse_cell_type_marker)
table(sorted_seur$Pred_Type)
# ast     end     mic     neu     oli Unknown 
# 2368    4554    2234    4780    2392       3

end_time <- Sys.time() #this took 6 hours
#now add that to the small object
cleaned_all@meta.data$cell_type<-sorted_seur$Pred_Type
table(cleaned_all@meta.data$seurat_clusters,cleaned_all@meta.data$cell_type)


# ast  end  mic  neu  oli Unknown
# 0    29   23   59 2502   35       0
# 1   142 1324  397   33   32       0
# 2   109  904  599   66   43       2
# 3   136  512  525  397   57       1
# 4    18   54   56   14  785       0
# 5    43  170   81  527   17       0
# 6   731   30   35   25    7       0
# 7     0    2    0    0  734       0
# 8    17  136   81    6  450       0
# 9   402   67   31   92   13       0
# 10   23  549   15    5    3       0
# 11    2  504   11    5    3       0
# 12    4    2    4  428    1       0
# 13  386    0    1    0    0       0
# 14   33   99   72   69    4       0
# 15   30   13   15   19  172       0
# 16    0    1    1  223    0       0
# 17    0    3  196    1    0       0
# 18  100   37    8   19    5       0
# 19    2    2    3  143    2       0
# 20    3   16   13   83    1       0
# 21    0    0    1  108    0       0
# 22   78    7    6    9    1       0
# 23   68    0    1    1   23       0
# 24    6   44   20    5    4       0
# 25    6   55    3    0    0       0


##plot seurat clusters by cell type
celltype_by_cluster<-table(cleaned_all@meta.data$seurat_clusters,cleaned_all@meta.data$cell_type)
celltype_by_cluster_ggplot<-as.data.frame(celltype_by_cluster)
ggplot(data=celltype_by_cluster_ggplot,aes(x=Var1,y=Freq,fill=Var2))+geom_bar(stat="identity")+xlab("Seurat Cluster")+ylab("number of cells")
#plot mitochondrial % by sample and cell type
# Idents(seur)<-"SampleID"
# VlnPlot(seur, features = c("percent.mt"), ncol = 1,pt.size = 0)
# Idents(seur)<-"cell_type"
# VlnPlot(seur, features = c("percent.mt"), ncol = 1,pt.size = 0)


#plot cell type by sample
celltype_by_sample<-table(cleaned_all@meta.data$SampleID,cleaned_all@meta.data$cell_type)
#change into percentages by sample
celltype_by_sample<-celltype_by_sample/rowSums(celltype_by_sample)
celltype_by_sample_ggplot<-as.data.frame(celltype_by_sample)
ggplot(data=celltype_by_sample_ggplot,aes(x=Var1,y=Freq,fill=Var2))+geom_bar(stat="identity")+xlab("Sample")+ylab("proportion of cells")+geom_text(aes(label=round(Freq,2)),position = position_fill(vjust = 0.5))
#will do stats on that
celltype_by_sample_ggplot$group<-ifelse(grepl("M",celltype_by_sample_ggplot$Var1)==TRUE,"Male","Female")
ggplot(data=celltype_by_sample_ggplot,aes(x=group,y=Freq,fill=Var2))+geom_bar(stat="identity")+xlab("Sample")+ylab("number of cells")+geom_text(aes(label=round(Freq,2)),position = position_fill(vjust = 0.5))
#do a t-test for each cell type , percentage in control vs percentage in alcohol
neu<-celltype_by_sample_ggplot[which(celltype_by_sample_ggplot$Var2=="neu"),]
t.test(data=neu,Freq~group)
end<-celltype_by_sample_ggplot[which(celltype_by_sample_ggplot$Var2=="end"),]
t.test(data=end,Freq~group)
mic<-celltype_by_sample_ggplot[which(celltype_by_sample_ggplot$Var2=="mic"),]
t.test(data=mic,Freq~group)
oli<-celltype_by_sample_ggplot[which(celltype_by_sample_ggplot$Var2=="oli"),]
t.test(data=oli,Freq~group)


###do this by cluster but put it in a loop for each of the 26 clusters
#cluster by sample
#plot cluster by sample
cluster_by_sample<-table(seur@meta.data$SampleID,seur@meta.data$seurat_clusters)
#change into percentages by sample
cluster_by_sample<-cluster_by_sample/rowSums(cluster_by_sample)
cluster_by_sample_ggplot<-as.data.frame(cluster_by_sample)
cluster_by_sample_ggplot$group<-ifelse(grepl("M",cluster_by_sample_ggplot$Var1)==TRUE,"Male","Female")
ggplot(data=cluster_by_sample_ggplot,aes(x=Var1,y=Freq,fill=Var2))+geom_bar(stat="identity")+xlab("Sample")+ylab("proportion of cells")+geom_text(aes(label=round(Freq,2)),position = position_fill(vjust = 0.5))
##t-test in each cluster
for (i in 0:17){cluster_i<-cluster_by_sample_ggplot[which(cluster_by_sample_ggplot$Var2==i),]
print(i)
print(t.test(data=cluster_i,Freq~group))}

