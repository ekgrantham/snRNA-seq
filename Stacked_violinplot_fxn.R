library(Seurat)
library(patchwork)
library(ggplot2)


features<- c("Meg3", "Ubb", "Ttr", "mt-Nd2", "Pde4b", "Slc4a4", "Flt1", "Nkain2", "Frmd4a", "Trpm3", "Dlc1", "Prkg1", "Tnr", "Tgfbr1", "Schip1", "Enpp2", "Ntng1", "Cfap299", "Chn2")
features = c("Sox10", "Cldn10", "Chat", "Dynlrb2", "Esam", "C1qc", "Syt1", "Gabra6", "Rbfox3", "Slc17a7")
features = c("Syt1", "Slc17a7", "Gad1", "Vcan", "Aqp4", "Slc4a4", "Mobp", "Htr2c", "Flt1", "Cfap299", "Pdgfrb", "Apbb1ip", "Ptprc")
features = c("Slc17a7", "Satb2", "Camk2a", "Gfap", "Aldh1l1", "Slc1a2", "Slc4a4", "Mbp", "Mobp", "Plp1", "Gad1", "Gad2", "Vcan", "Pdgfra", "Pcdh15", "Csf1r", "Apbb1ip", "P2ry12", "Flt1", "B2m")
features = c("Xist", "Uty")
features = c("Cck", "Gal", "Mt1", "Nr3c2", "Nr3c1", "Chat", "Mt2", "Sst", "Slc30a10", "Ache", "Gabra6", "Gabrd", "Grm4", "Crh")
features = c("Tlr7", "Tlr3", "Irf7", "Irf3", "Il33", "Il1r1", "Nfkb1", "Irak1")

#a <- VlnPlot(nts_seurat_full, features, stack = TRUE, sort = TRUE) +
  #theme(legend.position = "none") + ggtitle("identities on y-axis")

b <- VlnPlot(x_cleaned_all, features, stack = TRUE, sort = TRUE, flip = TRUE, split.by = "group") +
  theme(legend.position = "none") + ggtitle("Cell type markers by sex")

b
 
