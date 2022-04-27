#This is a script for using soupX to fitler out ambient rna from snRNA-seq data
#using the vignette available at the soupX github: https://github.com/constantAmateur/SoupX/blob/master/vignettes/pbmcTutorial.Rmd 

#If you have some 10X data which has been mapped with cellranger, the typical SoupX work flow would be:
install.packages('SoupX')
#If you want to use the latest experimental features, you can install the development version from github using the [devtools](https://devtools.r-lib.org/) `install_github` function as follows:
devtools::install_github("constantAmateur/SoupX",ref='devel') 
library(SoupX)

#Load data and estimate soup profile
#clean####
#Loading the data, create soup channel object, and run seurat on each individual sample####

#Sample1
untar(tarfile = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S1/filtered_feature_bc_matrix.tar.gz", exdir = "./NTS_snRNA/txg-linux-v1.2.0/S1/filtered_feature_bc_matrix")
untar(tarfile = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S1/raw_feature_bc_matrix.tar.gz", exdir = "./NTS_snRNA/txg-linux-v1.2.0/S1/raw_feature_bc_matrix")
sc_1 = load10X("/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S1/")

#However, there are some bits of meta data that are so essential that they have their own special loading functions.  
#The most essential is clustering information.  Without it, SoupX will still work, but you won't be able to automatically 
#estimate the contamination fraction and the correction step will be far less effective.  Metadata associated with our PBMC dataset is also bundled with SoupX.  
#We can use it to add clustering data by running,

nts_nuclei_1<- Read10X(data.dir = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S1/filtered_feature_bc_matrix/")
nts_seurat_1 <- CreateSeuratObject(counts = nts_nuclei_1, project = "NTSNuclei_6samps", min.cells = 3, min.features = 200)
nts_seurat_1    <- SCTransform(nts_seurat_1, verbose = F)
nts_seurat_1    <- RunPCA(nts_seurat_1, verbose = F)
nts_seurat_1    <- RunUMAP(nts_seurat_1, dims = 1:30, verbose = F)
nts_seurat_1   <- FindNeighbors(nts_seurat_1, dims = 1:30, verbose = F)
nts_seurat_1    <- FindClusters(nts_seurat_1, verbose = T)
nts_seurat_1 [["percent.mt"]] <- PercentageFeatureSet(nts_seurat_1 , pattern = "^mt-")

#In general you can add any meta data by providing a `data.frame` with row names equal to the column names of 
#the `toc` when building the `SoupChannel` object. 

meta_1 <- nts_seurat_1@meta.data
sc_1 = setClusters(sc_1,setNames(meta_1$seurat_clusters, rownames(meta_1)))
sc_1 = autoEstCont(sc_1)
out_1 = adjustCounts(sc_1)

head(sc_1$soupProfile[order(sc_1$soupProfile$est, decreasing = T), ], n = 20)

#Sample2
untar(tarfile = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S2/filtered_feature_bc_matrix.tar.gz", exdir = "./NTS_snRNA/txg-linux-v1.2.0/S2/filtered_feature_bc_matrix")
untar(tarfile = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S2/raw_feature_bc_matrix.tar.gz", exdir = "./NTS_snRNA/txg-linux-v1.2.0/S2/raw_feature_bc_matrix")
sc_2 = load10X("/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S2/")

nts_nuclei_2<- Read10X(data.dir = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S2/filtered_feature_bc_matrix/")
nts_seurat_2 <- CreateSeuratObject(counts = nts_nuclei_2, project = "NTSNuclei_6samps", min.cells = 3, min.features = 200)
nts_seurat_2    <- SCTransform(nts_seurat_2, verbose = F)
nts_seurat_2    <- RunPCA(nts_seurat_2, verbose = F)
nts_seurat_2    <- RunUMAP(nts_seurat_2, dims = 1:30, verbose = F)
nts_seurat_2   <- FindNeighbors(nts_seurat_2, dims = 1:30, verbose = F)
nts_seurat_2    <- FindClusters(nts_seurat_2, verbose = T)
nts_seurat_2 [["percent.mt"]] <- PercentageFeatureSet(nts_seurat_2 , pattern = "^mt-")

meta_2 <- nts_seurat_2@meta.data
sc_2 = setClusters(sc_2,setNames(meta_2$seurat_clusters, rownames(meta_2)))
sc_2 = autoEstCont(sc_2)
out_2 = adjustCounts(sc_2)

head(sc_2$soupProfile[order(sc_2$soupProfile$est, decreasing = T), ], n = 20)

#Sample3
untar(tarfile = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S3/filtered_feature_bc_matrix.tar.gz", exdir = "./NTS_snRNA/txg-linux-v1.2.0/S3/filtered_feature_bc_matrix")
untar(tarfile = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S3/raw_feature_bc_matrix.tar.gz", exdir = "./NTS_snRNA/txg-linux-v1.2.0/S3/raw_feature_bc_matrix")
sc_3 = load10X("/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S3/")

nts_nuclei_3<- Read10X(data.dir = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S3/filtered_feature_bc_matrix/")
nts_seurat_3 <- CreateSeuratObject(counts = nts_nuclei_3, project = "NTSNuclei_6samps", min.cells = 3, min.features = 200)
nts_seurat_3    <- SCTransform(nts_seurat_3, verbose = F)
nts_seurat_3    <- RunPCA(nts_seurat_3, verbose = F)
nts_seurat_3    <- RunUMAP(nts_seurat_3, dims = 1:30, verbose = F)
nts_seurat_3   <- FindNeighbors(nts_seurat_3, dims = 1:30, verbose = F)
nts_seurat_3    <- FindClusters(nts_seurat_3, verbose = T)
nts_seurat_3 [["percent.mt"]] <- PercentageFeatureSet(nts_seurat_3 , pattern = "^mt-")

meta_3 <- nts_seurat_3@meta.data
sc_3 = setClusters(sc_3,setNames(meta_3$seurat_clusters, rownames(meta_3)))
sc_3 = autoEstCont(sc_3)
out_3 = adjustCounts(sc_3)

head(sc_3$soupProfile[order(sc_3$soupProfile$est, decreasing = T), ], n = 20)

#Sample5
untar(tarfile = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S5/filtered_feature_bc_matrix.tar.gz", exdir = "./NTS_snRNA/txg-linux-v1.2.0/S5/filtered_feature_bc_matrix")
untar(tarfile = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S5/raw_feature_bc_matrix.tar.gz", exdir = "./NTS_snRNA/txg-linux-v1.2.0/S5/raw_feature_bc_matrix")
sc_5 = load10X("/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S5/")

nts_nuclei_5<- Read10X(data.dir = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S5/filtered_feature_bc_matrix/")
nts_seurat_5 <- CreateSeuratObject(counts = nts_nuclei_5, project = "NTSNuclei_6samps", min.cells = 3, min.features = 200)
nts_seurat_5    <- SCTransform(nts_seurat_5, verbose = F)
nts_seurat_5    <- RunPCA(nts_seurat_5, verbose = F)
nts_seurat_5    <- RunUMAP(nts_seurat_5, dims = 1:30, verbose = F)
nts_seurat_5   <- FindNeighbors(nts_seurat_5, dims = 1:30, verbose = F)
nts_seurat_5    <- FindClusters(nts_seurat_5, verbose = T)
nts_seurat_5 [["percent.mt"]] <- PercentageFeatureSet(nts_seurat_5 , pattern = "^mt-")

meta_5 <- nts_seurat_5@meta.data
sc_5 = setClusters(sc_5,setNames(meta_5$seurat_clusters, rownames(meta_5)))
sc_5 = autoEstCont(sc_5)
out_5 = adjustCounts(sc_5)

head(sc_5$soupProfile[order(sc_5$soupProfile$est, decreasing = T), ], n = 20)

#Sample6
untar(tarfile = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S6/filtered_feature_bc_matrix.tar.gz", exdir = "./NTS_snRNA/txg-linux-v1.2.0/S6/filtered_feature_bc_matrix")
untar(tarfile = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S6/raw_feature_bc_matrix.tar.gz", exdir = "./NTS_snRNA/txg-linux-v1.2.0/S6/raw_feature_bc_matrix")
sc_6 = load10X("/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S6/")

nts_nuclei_6 <- Read10X(data.dir = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S6/filtered_feature_bc_matrix/")
nts_seurat_6 <- CreateSeuratObject(counts = nts_nuclei_6, project = "NTSNuclei_6samps", min.cells = 3, min.features = 200)
nts_seurat_6    <- SCTransform(nts_seurat_6, verbose = F)
nts_seurat_6    <- RunPCA(nts_seurat_6, verbose = F)
nts_seurat_6    <- RunUMAP(nts_seurat_6, dims = 1:30, verbose = F)
nts_seurat_6   <- FindNeighbors(nts_seurat_6, dims = 1:30, verbose = F)
nts_seurat_6    <- FindClusters(nts_seurat_6, verbose = T)
nts_seurat_6 [["percent.mt"]] <- PercentageFeatureSet(nts_seurat_6 , pattern = "^mt-")

meta_6 <- nts_seurat_6@meta.data
sc_6 = setClusters(sc_6,setNames(meta_6$seurat_clusters, rownames(meta_6)))
sc_6 = autoEstCont(sc_6)
out_6 = adjustCounts(sc_6)

head(sc_6$soupProfile[order(sc_6$soupProfile$est, decreasing = T), ], n = 20)


#Sample7
untar(tarfile = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S7/filtered_feature_bc_matrix.tar.gz", exdir = "./NTS_snRNA/txg-linux-v1.2.0/S7/filtered_feature_bc_matrix")
untar(tarfile = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S7/raw_feature_bc_matrix.tar.gz", exdir = "./NTS_snRNA/txg-linux-v1.2.0/S7/raw_feature_bc_matrix")
sc_7 = load10X("/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S7/")

nts_nuclei_7 <- Read10X(data.dir = "/stor/scratch/WCAAR/emily_scratch/NTS_snRNA/txg-linux-v1.2.0/S7/filtered_feature_bc_matrix/")
nts_seurat_7 <- CreateSeuratObject(counts = nts_nuclei_7, project = "NTSNuclei_6samps", min.cells = 3, min.features = 200)
nts_seurat_7    <- SCTransform(nts_seurat_7, verbose = F)
nts_seurat_7    <- RunPCA(nts_seurat_7, verbose = F)
nts_seurat_7    <- RunUMAP(nts_seurat_7, dims = 1:30, verbose = F)
nts_seurat_7   <- FindNeighbors(nts_seurat_7, dims = 1:30, verbose = F)
nts_seurat_7    <- FindClusters(nts_seurat_7, verbose = T)
nts_seurat_7 [["percent.mt"]] <- PercentageFeatureSet(nts_seurat_7 , pattern = "^mt-")

meta_7 <- nts_seurat_7@meta.data
sc_7 = setClusters(sc_7,setNames(meta_7$seurat_clusters, rownames(meta_7)))
sc_7 = autoEstCont(sc_7)
out_7 = adjustCounts(sc_7)

head(sc_7$soupProfile[order(sc_7$soupProfile$est, decreasing = T), ], n = 20)

# plotMarkerDistribution(sc_1)
# plotMarkerDistribution(sc_2)
# plotMarkerDistribution(sc_3)
# plotMarkerDistribution(sc_5)
# plotMarkerDistribution(sc_6)
# plotMarkerDistribution(sc_7)

## Profiling the soup

#Usually the only reason to not have `estimateSoup` run automatically is if you want to change the default parameters or 
#have some other way of calculating the soup profile.  One case where you may want to do the latter is if you only have the 
#table of counts available and not the empty droplets.  In this case you can proceed by running

library(Matrix)
toc = sc_1$toc
scNoDrops = SoupChannel(toc,toc,calcSoupProfile=FALSE)
#Calculate soup profile
soupProf = data.frame(row.names = rownames(toc),
                      est = rowSums(toc)/sum(toc),
                      counts = rowSums(toc))
scNoDrops = setSoupProfile(scNoDrops,soupProf)

#In this case the `setSoupProfile` command is used instead of `estimateSoup` and directly adds 
#the custom estimation of the soup profile to the `SoupChannel` object.  Note that we have loaded the `Matrix` library 
#to help us manipulate the sparse matrix `toc`.

#This is usually not needed when using the `load10X` function as the cellranger produced values are automatically loaded.


###Automated method:###
#Estimating the contamination fraction using the automated is as simple as running:
sc_1 = autoEstCont(sc_1)

###Manual method:####
#non expressed genes####
#I tried using ribosomal genes only or mitochondrial only or both , or using the cytoplasmic enriched markers
#at the end decided on using the mitochondrial plus ribosomal

#cell markers
# cell_nuc_markers <- read_csv("cell_nuc_markers.csv")
#add a column saying that positive log fold change means high in cell , and low means high in nucleus and high p value means no change
# library(dplyr)
# cell_nuc_markers <-cell_nuc_markers %>% mutate(enrichment = case_when(
#   logFC>0 & P.Value<0.06  ~ "cell_marker",
#   logFC<0 & P.Value<0.06  ~ "nuclear_marker",
#   
#   TRUE ~ "Not_marker"
# ))
#to get cell markers, extract which enrichment = cell marker

# cell_markers <- subset(cell_nuc_markers, cell_nuc_markers$enrichment == "cell_marker")
# rownames(cell_markers) <- cell_markers$gene
# rownames_to_column(cell_markers, var = "gene")
#cell_markers$Row.names[grep("Rp",cell_markers$Row.names)]

# we selected a short list of gene sets and provided these to SoupX’s estimateNonExpressingCells function 
# (Mature Oligodendrocytes = Plp1, Mog, Pde4b, St18, Slc24a2, Pcdh9; choroid plexus cells = Ttr, Htr2c; 
#  Astrocytes = Apoe, Slc1a2, Slc1a3, Prex2, Cst3, Gabrb1, Gpc6; Neurons = Lsamp, Kcnip4, Tenm2, Nrg3, Celf2, Grin2a). 
# SoupX was used to estimate and remove ambient RNA contamination for every sample individually and the results 
# were verified using SoupX’s visualisation tools.

# nonExpressedGeneList1<-list(ne=c(intersect(rownames(sc$toc),cell_markers$gene[grep("Rp",cell_markers$gene)])))
# 
# nonExpressedGeneListX <- list(Oligo = c("Plp1", "Mog", "Pde4b", "St18", "Slc24a2", "Pcdh9"), CPCs = c("Ttr", "Htr2c"), Ast = c("Apoe", "Slc1a2", "Slc1a3", "Prex2", "Cst3", "Gabrb1", "Gpc6"), Neu = c("Lsamp", "Kcnip4", "Tenm2", "Nrg3", "Celf2", "Grin2a"))
# 
# nonExpressedGeneList2<-list(ne=c(intersect(rownames(sc$toc),cell_markers$gene[grep("Rp",cell_markers$gene)]),rownames(sc$toc)[grep("mt-",rownames(sc$toc))]))
# 
# nonExpressedGeneList3<-list(ne=rownames(sc$toc)[grep("mt-",rownames(sc$toc))])
# 
# nonExpressedGeneListall<-list(ne=intersect(rownames(sc$toc),cell_markers$gene))

# ribo_genes <- rownames(sc_1$toc)[grep("Rp",rownames(sc_1$toc))]
# mt_genes <- rownames(sc_1$toc)[grep("mt-",rownames(sc_1$toc))]
# nonExpressedGeneListall <- list(ribo_genes, mt_genes)

#using mitochondrial and ribosomal genes

# nonExpressedGeneList<-nonExpressedGeneListall


#estimate contamination and adjust counts#### (must run seurat on cleaned_all object FIRST before proceeding with manual - see Seurat script)
#Sample1 adjust counts
# useToEst = estimateNonExpressingCells(sc_1,nonExpressedGeneList=nonExpressedGeneList,clusters=FALSE)
# 
# sc_1 = calculateContaminationFraction(sc_1,nonExpressedGeneList,useToEst)
# 
# out_1 = adjustCounts(sc_1,verbose=2)
# 
# #zeroed cells#### s1
# cntSoggys1_cells= colSums(sc_1$toc>0)
# 
# cntStraineds1_cells= colSums(out_1>0)
# 
# mostZeroeds1_cells =data.frame((cntSoggys1_cells-cntStraineds1_cells)/cntSoggys1_cells)
# colnames(mostZeroeds1_cells)<-"percent_zeroed_s1"
# mostZeroeds1_cells$barcode<-gsub("\\-.*","",rownames(mostZeroeds1_cells ))
# meta_for_merge_s1<-cleaned_all@meta.data[which(cleaned_all@meta.data$orig.ident=="s1"),]
# meta_for_merge_s1$Barcode<-gsub("^.{0,3}","",rownames(meta_for_merge_s1))
# meta_for_merge_s1$Barcode<-gsub(".{2}$","", meta_for_merge_s1$Barcode)
# mostZeroeds1_cells<-merge(mostZeroeds1_cells,meta_for_merge_s1,by.x="barcode",by.y="Barcode")
# ggplot(mostZeroeds1_cells,aes(x=percent.mt,y=percent_zeroed_s1))+geom_point()
# ggplot(mostZeroeds1_cells,aes(x=seurat_clusters,y=percent_zeroed_s1))+geom_boxplot()
# 
# #Sample2 adjust counts
# useToEst = estimateNonExpressingCells(sc_2,nonExpressedGeneList=nonExpressedGeneList,clusters=FALSE)
# 
# sc_2 = calculateContaminationFraction(sc_2,nonExpressedGeneList,useToEst)
# 
# out_2 = adjustCounts(sc_2,verbose=2)
# 
# #zeroed cells#### s2
# cntSoggys2_cells= colSums(sc_2$toc>0)
# 
# cntStraineds2_cells= colSums(out_2>0)
# 
# mostZeroeds2_cells =data.frame((cntSoggys2_cells-cntStraineds2_cells)/cntSoggys2_cells)
# colnames(mostZeroeds2_cells)<-"percent_zeroed_s2"
# mostZeroeds2_cells$barcode<-gsub("\\-.*","",rownames(mostZeroeds2_cells ))
# meta_for_merge_s2<-cleaned_all@meta.data[which(cleaned_all@meta.data$orig.ident=="s2"),]
# meta_for_merge_s2$Barcode<-gsub("^.{0,3}","",rownames(meta_for_merge_s2))
# meta_for_merge_s2$Barcode<-gsub(".{2}$","", meta_for_merge_s2$Barcode)
# mostZeroeds2_cells<-merge(mostZeroeds2_cells,meta_for_merge_s1,by.x="barcode",by.y="Barcode")
# ggplot(mostZeroeds2_cells,aes(x=percent.mt,y=percent_zeroed_s2))+geom_point()
# ggplot(mostZeroeds2_cells,aes(x=seurat_clusters,y=percent_zeroed_s2))+geom_boxplot()
# 
# #Sample3 adjust counts
# useToEst = estimateNonExpressingCells(sc_3,nonExpressedGeneList=nonExpressedGeneList,clusters=FALSE)
# 
# sc_3 = calculateContaminationFraction(sc_3,nonExpressedGeneList,useToEst)
# 
# out_3 = adjustCounts(sc_3,verbose=2)
# 
# #zeroed cells#### s3
# cntSoggys3_cells= colSums(sc_3$toc>0)
# 
# cntStraineds3_cells= colSums(out_3>0)
# 
# mostZeroeds3_cells =data.frame((cntSoggys3_cells-cntStraineds3_cells)/cntSoggys3_cells)
# colnames(mostZeroeds3_cells)<-"percent_zeroed_s3"
# mostZeroeds3_cells$barcode<-gsub("\\-.*","",rownames(mostZeroeds3_cells ))
# meta_for_merge_s3<-cleaned_all@meta.data[which(cleaned_all@meta.data$orig.ident=="s3"),]
# meta_for_merge_s3$Barcode<-gsub("^.{0,3}","",rownames(meta_for_merge_s3))
# meta_for_merge_s3$Barcode<-gsub(".{2}$","", meta_for_merge_s3$Barcode)
# mostZeroeds3_cells<-merge(mostZeroeds3_cells,meta_for_merge_s1,by.x="barcode",by.y="Barcode")
# ggplot(mostZeroeds3_cells,aes(x=percent.mt,y=percent_zeroed_s3))+geom_point()
# ggplot(mostZeroeds3_cells,aes(x=seurat_clusters,y=percent_zeroed_s3))+geom_boxplot()
# 
# #Sample5 adjust counts ##accidentally changed all 5's to 6's so change this portion back
# useToEst = estimateNonExpressingCells(sc_5,nonExpressedGeneList=nonExpressedGeneList,clusters=FALSE)
# 
# sc_5 = calculateContaminationFraction(sc_5,nonExpressedGeneList,useToEst)
# 
# out_5 = adjustCounts(sc_5,verbose=2)
# 
# #zeroed cells#### s5
# cntSoggys5_cells= colSums(sc_5$toc>0)
# 
# cntStraineds5_cells= colSums(out_5>0)
# 
# mostZeroeds5_cells =data.frame((cntSoggys5_cells-cntStraineds5_cells)/cntSoggys5_cells)
# colnames(mostZeroeds5_cells)<-"percent_zeroed_s5"
# mostZeroeds5_cells$barcode<-gsub("\\-.*","",rownames(mostZeroeds5_cells ))
# meta_for_merge_s5<-cleaned_all@meta.data[which(cleaned_all@meta.data$orig.ident=="s5"),]
# meta_for_merge_s5$Barcode<-gsub("^.{0,3}","",rownames(meta_for_merge_s5))
# meta_for_merge_s5$Barcode<-gsub(".{2}$","", meta_for_merge_s5$Barcode)
# mostZeroeds5_cells<-merge(mostZeroeds5_cells,meta_for_merge_s1,by.x="barcode",by.y="Barcode")
# ggplot(mostZeroeds5_cells,aes(x=percent.mt,y=percent_zeroed_s5))+geom_point()
# ggplot(mostZeroeds5_cells,aes(x=seurat_clusters,y=percent_zeroed_s5))+geom_boxplot()
# 
# #Sample6 adjust counts
# useToEst = estimateNonExpressingCells(sc_6,nonExpressedGeneList=nonExpressedGeneList,clusters=FALSE)
# 
# sc_6 = calculateContaminationFraction(sc_6,nonExpressedGeneList,useToEst)
# 
# out_6 = adjustCounts(sc_6,verbose=2)
# 
# #zeroed cells#### s6
# cntSoggys6_cells= colSums(sc_6$toc>0)
# 
# cntStraineds6_cells= colSums(out_6>0)
# 
# mostZeroeds6_cells =data.frame((cntSoggys6_cells-cntStraineds6_cells)/cntSoggys6_cells)
# colnames(mostZeroeds6_cells)<-"percent_zeroed_s6"
# mostZeroeds6_cells$barcode<-gsub("\\-.*","",rownames(mostZeroeds6_cells ))
# meta_for_merge_s6<-cleaned_all@meta.data[which(cleaned_all@meta.data$orig.ident=="s6"),]
# meta_for_merge_s6$Barcode<-gsub("^.{0,3}","",rownames(meta_for_merge_s6))
# meta_for_merge_s6$Barcode<-gsub(".{2}$","", meta_for_merge_s6$Barcode)
# mostZeroeds6_cells<-merge(mostZeroeds6_cells,meta_for_merge_s1,by.x="barcode",by.y="Barcode")
# ggplot(mostZeroeds6_cells,aes(x=percent.mt,y=percent_zeroed_s6))+geom_point()
# ggplot(mostZeroeds6_cells,aes(x=seurat_clusters,y=percent_zeroed_s6))+geom_boxplot()
# 
# #Sample7 adjust counts
# useToEst = estimateNonExpressingCells(sc_7,nonExpressedGeneList=nonExpressedGeneList,clusters=FALSE)
# 
# sc_7 = calculateContaminationFraction(sc_7,nonExpressedGeneList,useToEst)
# 
# out_7 = adjustCounts(sc_7,verbose=2)
# 
# #zeroed cells#### s7
# cntSoggys7_cells= colSums(sc_7$toc>0)
# 
# cntStraineds7_cells= colSums(out_7>0)
# 
# mostZeroeds7_cells =data.frame((cntSoggys7_cells-cntStraineds7_cells)/cntSoggys7_cells)
# colnames(mostZeroeds7_cells)<-"percent_zeroed_s7"
# mostZeroeds7_cells$barcode<-gsub("\\-.*","",rownames(mostZeroeds7_cells ))
# meta_for_merge_s7<-cleaned_all@meta.data[which(cleaned_all@meta.data$orig.ident=="s7"),]
# meta_for_merge_s7$Barcode<-gsub("^.{0,3}","",rownames(meta_for_merge_s7))
# meta_for_merge_s7$Barcode<-gsub(".{2}$","", meta_for_merge_s7$Barcode)
# mostZeroeds7_cells<-merge(mostZeroeds7_cells,meta_for_merge_s1,by.x="barcode",by.y="Barcode")
# ggplot(mostZeroeds7_cells,aes(x=percent.mt,y=percent_zeroed_s7))+geom_point()
# ggplot(mostZeroeds7_cells,aes(x=clusters,y=percent_zeroed_s7))+geom_boxplot()

#merge all and compare the % zeroed of each genes --- i cant find mostZeroed genes objects.. not sure how to create

# merged_zeros<-merge(mostZeroeds1_genes,mostZeroeds1_genes,by="row.names")
# merged_zeros<-merge(merged_zeros,mostZeroeds2_genes,by.x="Row.names",by.y="row.names")
# merged_zeros<-merge(merged_zeros,mostZeroeds3_genes,by.x="Row.names",by.y="row.names")
# merged_zeros<-merge(merged_zeros,mostZeroeds5_genes,by.x="Row.names",by.y="row.names")
# merged_zeros<-merge(merged_zeros,mostZeroeds6_genes,by.x="Row.names",by.y="row.names")
# merged_zeros<-merge(merged_zeros,mostZeroeds7_genes,by.x="Row.names",by.y="row.names")
# colnames(merged_zeros)[6]<-"c21"
# merged_zeros[merged_zeros$Row.names=="Lcn2",]


#merge cluster numbers with zero percentages
# clusters_markers <- read_csv("snRNA_full_Nov21/aggr_13_samps/clusters_markers_snRNA_cie_combined.csv")
# 
# cluster_markers_zero_merges<-merge(clusters_markers,merged_zeros,by="gene")
# ggplot(cluster_markers_zero_merges,aes(x=as.character(cluster),y=c17))+geom_boxplot()
# ggplot(cluster_markers_zero_merges,aes(x=as.character(cluster),y=c19))+geom_boxplot()
# library(reshape)
# cluster_markers_zero_merges_melt<- melt(cluster_markers_zero_merges, id=c("gene","avg_log2FC","cluster"),measure.vars=c("c17","c19","c21","v16"))
# colnames(cluster_markers_zero_merges_melt)[4:6]<-c("sample","percentzeroed")
# pdf("cluster_markers_vs_percentzero.pdf",height=6)
# for (i in 0:27){
#   p<-ggplot(cluster_markers_zero_merges_melt[which(cluster_markers_zero_merges_melt$cluster==i),],aes(x=avg_log2FC,y=percentzeroed))+geom_point()+ylim(c(0,0.9))+xlim(c(0,4))+ylab("percent zeroed")+ggtitle(i)+facet_grid(~sample)
#   print(p)}
# dev.off()

##put all of that in a big seurat
library(Seurat)
scs = list(s1=out_1, s2=out_2, s3=out_3, s5=out_5, s6=out_6, s7=out_7)
srat<-list()
for(nom in names(scs)){
  #Clean channel named 'nom'
  #tmp = adjustCounts(scs[[nom]])
  #already adjusted 
  tmp<-scs[[nom]]
  #Add experiment name to cell barcodes to make them unique
  colnames(tmp) = paste0(nom,'_',colnames(tmp))
  #Store the result
  srat[[nom]] = tmp
}
#Combine all count matricies into one matrix
srat = do.call(cbind,srat)
cleaned_all = CreateSeuratObject(srat)
cleaned_all [["percent.mt"]] <- PercentageFeatureSet(cleaned_all , pattern = "^mt-")

View(cleaned_all@meta.data)
#extract the metadata,add some labels and add it again
# meta_clean<-cleaned_all@meta.data
# meta_clean$group<-ifelse(grepl("M",as.character(meta_clean$orig.ident)),"Male","Female")
# meta_clean$group<-as.factor(meta_clean$group)
meta_clean$Barcode <- as.factor(rownames(meta_clean))
meta_clean$Barcode<-gsub("^.{0,3}","",meta_clean$Barcode)
#add sampleID info as this will have sex
sample_info<-SampleID
#join sampleID with meta_clean by barcode (rownames)

# #add metadata
cleaned_all <- AddMetaData(
   object = cleaned_all,
    metadata = meta_clean, )

View(cleaned_all@meta.data)

View(srat@meta.data)

dim(cleaned_all)
table(cleaned_all@meta.data$orig.ident) 

#save the cleaned object and start a clean seurat
save(cleaned_all,file="NTS_cleaned.RData")
save(scs, sc_1, sc_2, sc_3, sc_5, sc_6, sc_7, file = "soupxoutputs_manual.Rdata")
