require(Seurat)

#require(hdf5r)

require(harmony)

require(ggplot2)

require(patchwork)

require(tidyverse)


setwd(".")
output_dir <- "output"


######################################
# Loading in count data

dirs <- list.files("/mnt/vast/hpc/MenonLab/SenNet/snRNAseq", pattern="-GEX", include.dirs=T)
allcounts=c()
allmeta=c()


for (ii in dirs) {

  filenam=(paste0(ii,"-GEX/",ii,"out2/outs/",ii,"_cellbender.h5")) 
  tempall=Read10X_h5(filenam)
  temp1=tempall[[1]]
  remgenes=grep("^MT-|\\.|-PS|-AS|^RP[0-9]",rownames(temp1))
  temp2=temp1[,which(colSums(temp1[-remgenes,])>=500)]
  temp1=temp2
  cellid=paste0(colnames(temp1),"_",ii)
  meta.data=data.frame(batch=ii,cellid=cellid,bc=colnames(temp1),mtpct=colSums(temp1[grep("^MT-",rownames(temp1)),])/colSums(temp1[-remgenes,]))
  rownames(meta.data)=cellid
  colnames(temp1)=cellid
  keepnuc=rownames(meta.data)[meta.data$mtpct<=0.05]
  allmeta=rbind(allmeta,meta.data[keepnuc,])
  allcounts=cbind(allcounts,temp1[,keepnuc])
  print(c(ii,nrow(meta.data),nrow(allmeta)))

}

rm(temp1,temp2,tempall);gc()

obj1=CreateSeuratObject(counts=allcounts[-remgenes,],meta.data=allmeta)
rm(allcounts,allmeta);gc()


obj1=NormalizeData(obj1)
obj1=ScaleData(obj1)
obj1=FindVariableFeatures(obj1)
obj1=RunPCA(obj1)
obj1=RunHarmony(obj1,group.by="batch",dims.use=1:30,reduction="pca")
obj1=RunUMAP(obj1,reduction="harmony",dim=1:30)
obj1=FindNeighbors(obj1,reduction="harmony",dims.use=1:30)
obj1=FindClusters(obj1,resolution=1)

library(DoubletFinder)
nExp_poi <- round(0.05*nrow(obj1@meta.data))
obj1 <- doubletFinder_v3(obj1, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
colnames(obj1@meta.data)[tail(1:ncol(obj1@meta.data), 2)] <- str_extract(colnames(obj1@meta.data)[tail(1:ncol(obj1@meta.data), 2)] ,
                                                                         "pANN|DF.classifications" )

require(tidyverse)
meta_file <- read_csv("metadata.csv")
obj1@meta.data <- obj1@meta.data %>% left_join(meta_file,by = c("batch"="CU_Core_ID"))
rownames(obj1@meta.data) <- colnames(obj1)


keep_clusters <- obj1@meta.data %>% select(seurat_clusters, DF.classifications) %>%
  group_by(seurat_clusters) %>%
  summarize(prop = sum(DF.classifications=="Singlet")/dplyr::n(), count = n()) %>%
  ungroup() %>% filter(prop>=.8) %>% pull(seurat_clusters)

obj2 <- subset(obj1, subset=seurat_clusters %in% keep_clusters)
obj2 <- subset(obj2, subset=DF.classifications =="Singlet")


obj2=NormalizeData(obj2)
obj2=ScaleData(obj2)
obj2=FindVariableFeatures(obj2)
obj2=RunPCA(obj2)
obj2=RunHarmony(obj2,group.by="batch",dims=1:30,reduction="pca")
obj2=RunUMAP(obj2,reduction="harmony",dim=1:30)
obj2=FindNeighbors(obj2,reduction="harmony",dims=1:30)
obj2=FindClusters(obj2,resolution=1)


dp2 <- DotPlot(obj2,features=c("SNAP25","SLC17A7","GAD1","GAD2","MOG","MAG","PDGFRA","AQP4","FGFR3","TMEM119",
                               "AIF1","CLDN5","PDGFRB"))+RotatedAxis()
ggsave(paste0(output_dir,"DIM_PLOT2_ALL_BATCHES.png"), dp1, height=12, width=14)

png(paste0(output_dir,"DOT_PLOT2_ALL_BATCHES.png"), height=10,width=9,units='in',res=800)
print(dp2)
dev.off()


broadclass=rep("Mixed",nrow(obj2@meta.data))
broadclass[which(obj2$seurat_clusters %in% c(0,2,6,8,13,14,15,18,20,23,26,27,31,32,35))]="NeuronGlut"
broadclass[which(obj2$seurat_clusters %in% c(9,10,11,17,22,24,28,29))]="NeuronGABA"
broadclass[which(obj2$seurat_clusters %in% c(1,5,12,25,30))]="Oligodendrocytes"
broadclass[which(obj2$seurat_clusters %in% c(3,21))]="Astrocytes"
broadclass[which(obj2$seurat_clusters %in% c(7))]="Microglia"
broadclass[which(obj2$seurat_clusters %in% c(4))]="OPCs"
broadclass[which(obj2$seurat_clusters %in% c(16,19))]="Vascular"

obj2$broadclass <- broadclass


######################################################
# The next steps include cleaning the data and sub-clustering

clean_cells_1 <- function(ct){
  print(ct)
  obj3 <- subset(obj2, subset=broadclass %in% ct)
  obj3=NormalizeData(obj3)
  obj3=ScaleData(obj3)
  obj3=FindVariableFeatures(obj3)
  obj3=RunPCA(obj3)
  obj3=RunHarmony(obj3,group.by="batch",dims=1:30,reduction="pca")
  obj3=RunUMAP(obj3,reduction="harmony",dim=1:30)
  obj3=FindNeighbors(obj3,reduction="harmony",dims=1:30)
  obj3=FindClusters(obj3,resolution=1)
  
  dp1 <- (DimPlot(obj3,label=T, group.by="broadclass")+DimPlot(obj3, group.by="Sex", label=T))/
    (DimPlot(obj3, group.by="batch", label=T)+DimPlot(obj3, group.by="seurat_clusters", label=T))
  dp2 <- DotPlot(obj3,features=c("SNAP25","SLC17A7","GAD1","GAD2",
                                 "MOG","MAG","PDGFRA","AQP4","FGFR3",
                                 "TMEM119","AIF1","CLDN5","PDGFRB"))+RotatedAxis()
  
  ggsave(paste0(output_dir,"DIM_PLOT2_ALL_BATCHES_",ct,".png"), dp1, height=12, width=14)
  
  png(paste0(output_dir,"DOT_PLOT2_ALL_BATCHES_",ct,".png"), height=10,width=9,units='in',res=800)
  print(dp2)
  dev.off()
  
  save(obj3, file = paste0(output_dir,ct,"_initial.RData"))
  rm(obj3,dp1,dp2)
  
  return()
}

sapply(unique(obj2$broadclass), clean_cells_1)


####

filter_cell_types <- function(ct, filter_out){
  load(file = paste0(output_dir,ct,"_initial.RData"))
  subset(obj3, subset = seurat_clusters %in% filter_out, invert=T)
  assign(ct, obj3,env=.GlobalEnv)
}

filter_cell_types("Oligodendrocytes",filter_out=c(8,10,13,14))
filter_cell_types("NeuronGlut",filter_out=c(8,15,17,18,19,22,24,25))
filter_cell_types("Astrocytes",filter_out=c(10,11,12))
filter_cell_types("Mixed",filter_out=c(0,1,2,3)) # empty
filter_cell_types("NeuronGABA",filter_out=c(3,9,18,25,27,29,30))
filter_cell_types("Microglia",filter_out=c(5,7,8,9,11,13))
filter_cell_types("OPCs",filter_out=c(5,6,9))
filter_cell_types("Vascular",filter_out=c(5,9,10,11,12))


obj <- merge(Oligodendrocytes, y =c(NeuronGlut,Astrocytes,NeuronGABA,Microglia,OPCs,Vascular))
save(obj, file = paste0(output_dir,"filtered_object.RData"))

####


obj=NormalizeData(obj)
obj=ScaleData(obj)
obj=FindVariableFeatures(obj)
obj=RunPCA(obj)
obj=RunHarmony(obj,group.by="batch",dims=1:30,reduction="pca")
obj=RunUMAP(obj,reduction="harmony",dim=1:30)
obj=FindNeighbors(obj,reduction="harmony",dims=1:30)
obj=FindClusters(obj,resolution=.3)

dp1 <- (DimPlot(obj,label=T, group.by="broadclass")+DimPlot(obj, group.by="Sex", label=T))/
  (DimPlot(obj, group.by="batch", label=T)+DimPlot(obj, group.by="seurat_clusters", label=T))
dp2 <- DotPlot(obj,features=c("SNAP25","SLC17A7","GAD1","GAD2","MOG","MAG","PDGFRA","AQP4","FGFR3","TMEM119","AIF1","CLDN5","PDGFRB"))+RotatedAxis()

ggsave(paste0(output_dir,"DIM_PLOT3_ALL_BATCHES_.png"), dp1, height=12, width=14)

png(paste0(output_dir,"DOT_PLOT3_ALL_BATCHES.png"), height=10,width=9,units='in',res=800)
print(dp2)
dev.off()

save(obj, file = paste0(output_dir,"filtered_object.RData"))

#####


clean_cell_type_2 <- function(ct){
  print(ct)
  obj3 <- subset(obj, subset=broadclass %in% ct)
  
  obj3=NormalizeData(obj3)
  obj3=ScaleData(obj3)
  obj3=FindVariableFeatures(obj3)
  obj3=RunPCA(obj3)
  obj3=RunHarmony(obj3,group.by="batch",dims=1:30,reduction="pca")
  obj3=RunUMAP(obj3,reduction="harmony",dim=1:30)
  obj3=FindNeighbors(obj3,reduction="harmony",dims=1:30)
  obj3=FindClusters(obj3,resolution=.5)
  
  obj3$sub_clusters <- paste0(ct,"_", obj3$seurat_clusters)
  
  dp1 <- (DimPlot(obj3,label=T, group.by="broadclass")+DimPlot(obj3, group.by="Sex", label=T))/
    (DimPlot(obj3, group.by="batch", label=T)+DimPlot(obj3, group.by="sub_clusters", label=T))
  ggsave(paste0(output_dir,"DIM_PLOT2_ALL_BATCHES_",ct,"_2.png"), dp1, height=12, width=14)
  
  Idents(obj3) <- obj3$sub_clusters
  roc_markers <- FindAllMarkers(obj3,test.use="roc")
  roc_markers %>% write_csv(paste0(output_dir,ct,"_roc_markers.csv"))
  
  cluster_ids <- data.frame(cellid = rownames(obj3@meta.data), sub_clusters = obj3$sub_clusters)
  
  save(cluster_ids, file = paste0(output_dir,ct,"_initial_cluster_ids.RData"))
  rm(obj3,dp1)
  
  return()
}

sapply(c("NeuronGABA","NeuronGlut"), clean_cell_type_2)


Vasc <- subset(obj, subset=broadclass %in% "Vascular")
load(paste0(output_dir,"Vascular_initial_cluster_ids.RData"))
Vasc$sub_clusters <- paste0(str_extract(cluster_ids$sub_clusters ,"Astro|Olig|OPC|Micro|Vasc|Glut|GABA"),"_", str_extract(cluster_ids$sub_cluster,"[0-9]+"))
Vasc$sub_cell_types <- "Vasc"
Vasc$sub_cell_types[Vasc$sub_clusters %in% paste0("Vasc",c(5,8,9,10))] <- "Mixed"
Vasc$sub_cell_types[Vasc$sub_clusters %in% paste0("Vasc",c(0,4))] <- "Vasc_PARD3_ZEB1"
Vasc$sub_cell_types[Vasc$sub_clusters %in% paste0("Vasc",c(1,2,3,6))] <- "Vasc_SYNE1"
Vasc$sub_cell_types[Vasc$sub_clusters %in% paste0("Vasc",c(7))] <- "Vasc_TSHZ2_PLAG2G4A"

OPC <- subset(obj, subset=broadclass %in% "OPCs")
load(paste0(output_dir,"OPCs_initial_cluster_ids.RData"))
OPC$sub_clusters <- paste0(str_extract(cluster_ids$sub_clusters ,"Astro|Olig|OPC|Micro|Vasc|Glut|GABA"),"_", str_extract(cluster_ids$sub_cluster,"[0-9]+"))
OPC$sub_cell_types <- "OPC"
OPC$sub_cell_types[OPC$sub_clusters %in% paste0("OPC_",c(7,8,9,12,13))] <- "Mixed"
OPC$sub_cell_types[OPC$sub_clusters %in% paste0("OPC_",c(0,4,1,6,10))] <- "OPC_SGCZ"
OPC$sub_cell_types[OPC$sub_clusters %in% paste0("OPC_",c(3,5))] <- "OPC_BCAS1"

Micro <- subset(obj, subset=broadclass %in% "Microglia")
load(paste0(output_dir,"Microglia_initial_cluster_ids.RData"))
Micro$sub_clusters <- paste0(str_extract(cluster_ids$sub_clusters ,"Astro|Olig|OPC|Micro|Vasc|Glut|GABA"),"_", str_extract(cluster_ids$sub_cluster,"[0-9]+"))
Micro$sub_cell_types <- "Micro"
Micro$sub_cell_types[Micro$sub_clusters %in% paste0("Micro_",c(6,10,12,9))] <- "Mixed"
Micro$sub_cell_types[Micro$sub_clusters %in% paste0("Micro_",c(7))] <- "Micro_CD83"

Astro <- subset(obj, subset=broadclass %in% "Astrocytes")
load(paste0(output_dir,"Astrocytes_initial_cluster_ids.RData"))
Astro$sub_clusters <- paste0(str_extract(cluster_ids$sub_clusters ,"Astro|Olig|OPC|Micro|Vasc|Glut|GABA"),"_", str_extract(cluster_ids$sub_cluster,"[0-9]+"))
Astro$sub_cell_types <- "Astro"
Astro$sub_cell_types[Astro$sub_clusters %in% paste0("Astro_",c(5,8,11))] <- "Mixed"
Astro$sub_cell_types[Astro$sub_clusters %in% paste0("Astro_",c(3,4,10))] <- "Astro_DPP10"
Astro$sub_cell_types[Astro$sub_clusters %in% paste0("Astro_",c(9))] <- "Astro_SYT1"
Astro$sub_cell_types[Astro$sub_clusters %in% paste0("Astro_",c(0,2,7,1,6,12))] <- "Astro_HPSE2"

Olig <- subset(obj, subset=broadclass %in% "Oligodendrocytes")
load(paste0(output_dir,"Oligodendrocytes_initial_cluster_ids.RData"))
Olig$sub_clusters <- paste0(str_extract(cluster_ids$sub_clusters ,"Astro|Olig|OPC|Micro|Vasc|Glut|GABA"),"_", str_extract(cluster_ids$sub_cluster,"[0-9]+"))
Olig$sub_cell_types <- "Olig"
Olig$sub_cell_types[Olig$sub_clusters %in% paste0("Olig_",c(13))] <- "Mixed"
Olig$sub_cell_types[Olig$sub_clusters %in% paste0("Olig_",c(6,8,12))] <- "Olig_MAP1B"
Olig$sub_cell_types[Olig$sub_clusters %in% paste0("Olig_",c(11))] <- "Olig_MDGA2_NLGN1"
Olig$sub_cell_types[Olig$sub_clusters %in% paste0("Olig_",c(1,5))] <- "Olig_TANC2"
Olig$sub_cell_types[Olig$sub_clusters %in% paste0("Olig_",c(9))] <- "Olig_LINC01099"
Olig$sub_cell_types[Olig$sub_clusters %in% paste0("Olig_",c(0,7))] <- "Olig_NLGN1"

GABA <- subset(obj, subset=broadclass %in% "NeuronGABA")
load(paste0(output_dir,"NeuronGABA_initial_cluster_ids.RData"))
GABA$sub_clusters <- paste0(str_extract(cluster_ids$sub_clusters ,"Astro|Olig|OPC|Micro|Vasc|Glut|GABA"),"_", str_extract(cluster_ids$sub_cluster,"[0-9]+"))
GABA$sub_cell_types <- "GABA"
GABA$sub_cell_types[GABA$sub_clusters %in% paste0("GABA_",c(0,16,23))] <- "GABA_TRPC5_TAC1"
GABA$sub_cell_types[GABA$sub_clusters %in% paste0("GABA_",c(8,14))] <- "GABA_SRGAP1"
GABA$sub_cell_types[GABA$sub_clusters %in% paste0("GABA_",c(9,4,15,18))] <- "GABA_TPD52L1"
GABA$sub_cell_types[GABA$sub_clusters %in% paste0("GABA_",c(7,17,24))] <- "Mixed" #NCKAP1
GABA$sub_cell_types[GABA$sub_clusters %in% paste0("GABA_",c(1,6,8,11,19))] <- "GABA_PDE1A"
GABA$sub_cell_types[GABA$sub_clusters %in% paste0("GABA_",c(2,5,20))] <- "GABA_SVEP1"
GABA$sub_cell_types[GABA$sub_clusters %in% paste0("GABA_",c(21,22))] <- "GABA_PTPRK"
GABA$sub_cell_types[GABA$sub_clusters %in% paste0("GABA_",c(10))] <- "GABA_SRGAP1_PTPRK"

Glut <- subset(obj, subset=broadclass %in% "NeuronGlut")
load(paste0(output_dir,"NeuronGlut_initial_cluster_ids.RData"))
Glut$sub_clusters <- paste0(str_extract(cluster_ids$sub_clusters ,"Astro|Olig|OPC|Micro|Vasc|Glut|GABA"),"_", str_extract(cluster_ids$sub_cluster,"[0-9]+"))
Glut$sub_cell_types <- "Glut"
Glut$sub_cell_types[Glut$sub_clusters %in% paste0("Glut_",c(12,13))] <- "Glut_ERBB4"
Glut$sub_cell_types[Glut$sub_clusters %in% paste0("Glut_",c(4,5,10))] <- "Glut_VWC2L"
Glut$sub_cell_types[Glut$sub_clusters %in% paste0("Glut_",c(3,6,7,8,16))] <- "Glut_GRIK3" #NCKAP1
Glut$sub_cell_types[Glut$sub_clusters %in% paste0("Glut_",c(17))] <- "Glut_HTR2C"
Glut$sub_cell_types[Glut$sub_clusters %in% paste0("Glut_",c(9,14))] <- "Mixed" #NR4A2
Glut$sub_cell_types[Glut$sub_clusters %in% paste0("Glut_",c(1,2))] <- "Glut_PRAG1"
Glut$sub_cell_types[Glut$sub_clusters %in% paste0("Glut_",c())] <- "Glut_SRGAP1_PTPRK"


clean_cell_type_3 <- function(obj,ct){
  
  obj3 <- subset(obj, subset=sub_cell_types %in% "Mixed", invert=T)
  
  obj3=NormalizeData(obj3)
  obj3=ScaleData(obj3)
  obj3=FindVariableFeatures(obj3)
  obj3=RunPCA(obj3)
  obj3=RunHarmony(obj3,group.by="batch",dims=1:30,reduction="pca")
  obj3=RunUMAP(obj3,reduction="harmony",dim=1:30)
  obj3=FindNeighbors(obj3,reduction="harmony",dims=1:30)
  obj3=FindClusters(obj3,resolution=.3)
  
  obj3$sub_clusters <- paste0(ct,"_", obj3$seurat_clusters)
  
  dp1 <- (DimPlot(obj3,label=T, group.by="broadclass")+DimPlot(obj3, group.by="Sex", label=T))/
    (DimPlot(obj3, group.by="sub_cell_types", label=T)+DimPlot(obj3, group.by="sub_clusters", label=T))
  ggsave(paste0(output_dir,"DIM_PLOT2_ALL_BATCHES_",ct,"_3.png"), dp1, height=12, width=15)
  
  Idents(obj3) <- obj3$sub_clusters
  roc_markers <- FindAllMarkers(obj3,test.use="roc")
  roc_markers %>% write_csv(paste0(output_dir,ct,"_roc_markers_2.csv"))
  
  roc_markers <- read.csv(paste0(output_dir,ct,"_roc_markers_2.csv"))
  
  features_input <- roc_markers %>% group_by(cluster) %>% filter(abs(myAUC-.5)>=.25) %>% top_n(2) %>% ungroup() %>% arrange(cluster) %>% 
    pull(gene) %>% unique()
  Idents(obj3) <- obj3$sub_clusters
  dp <- DotPlot(obj3, features = features_input)+RotatedAxis()
  png(paste0(output_dir,ct,"_DOT_PLOT_4.png"), height=10,width=12,units='in',res=800)
  print(dp)
  dev.off()
  
  save(obj3, file = paste0(output_dir,ct,"_filtered_2.RData"))
  rm(obj3,dp1)
  
  return()
}

clean_cell_type_3(Vasc,"Vasc")
clean_cell_type_3(Astro,"Astro")
clean_cell_type_3(OPC,"OPC")
clean_cell_type_3(Micro, "Micro")
clean_cell_type_3(Olig, "Olig")
clean_cell_type_3(GABA,"GABA")
clean_cell_type_3(Glut,"Glut")


clean_cell_type_4 <- function(ct, remove){
  load(file = paste0(output_dir,ct,"_filtered_2.RData"))
  obj3 <- subset(obj3, subset=sub_clusters %in% paste0(ct,"_",remove), invert=T)
  
  obj3=NormalizeData(obj3)
  obj3=ScaleData(obj3)
  obj3=FindVariableFeatures(obj3)
  obj3=RunPCA(obj3)
  obj3=RunHarmony(obj3,group.by="batch",dims=1:30,reduction="pca")
  obj3=RunUMAP(obj3,reduction="harmony",dim=1:30)
  obj3=FindNeighbors(obj3,reduction="harmony",dims=1:30)
  obj3=FindClusters(obj3,resolution=ifelse(ct %in% c("Glut","GABA"),.5,.3))
  
  obj3$sub_clusters <- paste0(ct,"_", obj3$seurat_clusters)
  
  dp1 <- (DimPlot(obj3,label=T, group.by="broadclass")+DimPlot(obj3, group.by="Sex", label=T))/
    (DimPlot(obj3, group.by="sub_cell_types", label=T)+DimPlot(obj3, group.by="sub_clusters", label=T))
  ggsave(paste0(output_dir,"DIM_PLOT2_ALL_BATCHES_",ct,"_4.png"), dp1, height=12, width=15)
  
  Idents(obj3) <- obj3$sub_clusters
  roc_markers <- FindAllMarkers(obj3,test.use="roc")
  roc_markers %>% write_csv(paste0(output_dir,ct,"_roc_markers_3.csv"))
  
  features_input <- c("SNAP25","SLC17A7","CUX2","RORB","THEMIS","CTGF")
  Idents(obj3) <- obj3$sub_clusters
  dp <- DotPlot(obj3, features = features_input)+RotatedAxis()
  png(paste0(output_dir,ct,"_DOT_PLOT_5.png"), height=10,width=12,units='in',res=800)
  print(dp)
  dev.off()
  
  save(obj3, file = paste0(output_dir,ct,"_filtered_3.RData"))
  rm(obj3,dp1)
  
  return(obj3)
}

Vasc <- clean_cell_type_4("Vasc", c(2,6,7))
Micro <- clean_cell_type_4("Micro", c(4))
Olig <- clean_cell_type_4("Olig", c(5))
Glut <- clean_cell_type_4("Glut", c(12))
GABA <- clean_cell_type_4("GABA", c(9))

clean_cell_type_5 <- function(ct, remove){
  load(file = paste0(output_dir,ct,"_filtered_3.RData"))
  obj3 <- subset(obj3, subset=sub_clusters %in% paste0(ct,"_",remove), invert=T)
  
  obj3=NormalizeData(obj3)
  obj3=ScaleData(obj3)
  obj3=FindVariableFeatures(obj3)
  obj3=RunPCA(obj3)
  obj3=RunHarmony(obj3,group.by="batch",dims=1:30,reduction="pca")
  obj3=RunUMAP(obj3,reduction="harmony",dim=1:30)
  obj3=FindNeighbors(obj3,reduction="harmony",dims=1:30)
  obj3=FindClusters(obj3,resolution=ifelse(ct %in% c("Glut","GABA"),.5,.3))
  
  obj3$sub_clusters <- paste0(ct,"_", obj3$seurat_clusters)
  
  dp1 <- (DimPlot(obj3,label=T, group.by="broadclass")+DimPlot(obj3, group.by="Sex", label=T))/
    (DimPlot(obj3, group.by="sub_cell_types", label=T)+DimPlot(obj3, group.by="sub_clusters", label=T))
  ggsave(paste0(output_dir,"DIM_PLOT2_ALL_BATCHES_",ct,"_5.png"), dp1, height=12, width=15)
  
  Idents(obj3) <- obj3$sub_clusters
  save(obj3, file = paste0(output_dir,ct,"_filtered_4.RData"))
  
  roc_markers <- FindAllMarkers(obj3,test.use="roc")
  roc_markers %>% write_csv(paste0(output_dir,ct,"_roc_markers_4.csv"))
  
  roc_markers <- read.csv(paste0(output_dir,ct,"_roc_markers_4.csv"))
  features_input <- c("SNAP25","SLC17A7","CUX2","RORB","THEMIS","CTGF")
  dp <- DotPlot(obj3, features = features_input)+RotatedAxis()
  png(paste0(output_dir,ct,"_DOT_PLOT_6.png"), height=10,width=12,units='in',res=800)
  print(dp)
  dev.off()
  
  rm(obj3,dp1)
  
  return(obj3)
}

Glut$clusters <- Glut$sub_clusters
Glut$clusters[Glut$sub_clusters %in% c("Glut_2")] <- "Glut_CUX2"
Glut$clusters[Glut$sub_clusters %in% c("Glut_6")] <- "Glut_THEMIS"
Glut$clusters[Glut$sub_clusters %in% c("Glut_9","Glut_15")] <- "Glut_CTGF"
Glut$clusters[Glut$sub_clusters %in% c("Glut_0","Glut_10","Glut_7","Glut_3","Glut_4","Glut_12","Glut_13")] <- "Glut_RORB"
Glut$clusters[Glut$sub_clusters %in% c('Glut_1',"Glut_5")] <- "Glut_CUX2_RORB"
Glut$clusters[Glut$sub_clusters %in% c("Glut_11")] <- "Glut_THEMIS_RORB"
Glut$clusters[Glut$sub_clusters %in% c("Glut_8")] <- "Glut_TRPM3"
Glut$clusters[Glut$sub_clusters %in% c("Glut_16")] <- "Glut_GRIK2_ROBO2"
Glut$clusters[Glut$sub_clusters %in% c("Glut_18")] <- "Glut_ZNF385D_HTR2C"
Glut$clusters[Glut$sub_clusters %in% c("Glut_14","Glut_17")] <- "Mixed"


Glut <- clean_cell_type_5("Glut", c(14,17))

GABA$clusters <- GABA$sub_clusters
GABA$clusters[GABA$sub_clusters %in% c("GABA_7","GABA_2","GABA_5","GABA_16","GABA_14")] <- "GABA_SST"
GABA$clusters[GABA$sub_clusters %in% c("GABA_8","GABA_4")] <- "GABA_LAMP5"
GABA$clusters[GABA$sub_clusters %in% c("GABA_10","GABA_0","GABA_9","GABA_11","GABA_3","GABA_18")] <- "GABA_PVALB"
GABA$clusters[GABA$sub_clusters %in% c("GABA_6","GABA_1","GABA_13","GABA_12","GABA_15")] <- "GABA_VIP"
GABA$clusters[GABA$sub_clusters %in% c("GABA_17")] <- "GABA_NDNF"
GABA$clusters[GABA$sub_clusters %in% c("GABA_19")] <- "GABA_SST_CHODL"

GABA$cluster_level <- GABA$clusters
GABA$cluster_level[GABA$clusters %in% c("GABA_NDNF")] <- "GABA_L1_LAMP5_NDNF"
GABA$cluster_level[GABA$clusters %in% c("GABA_SST_CHODL")] <- "GABA_L1_SST"
GABA$cluster_level[GABA$clusters %in% "GABA_VIP"] <- "GABA_L1_VIP"

GABA$cluster_level[GABA$sub_clusters %in% "GABA_18"] <- "GABA_L3_L6_KCNIP4"
GABA$cluster_level[GABA$sub_clusters %in% "GABA_3"] <- "GABA_PVALB"

GABA$cluster_level[GABA$sub_clusters %in% "GABA_9"] <- "GABA_L5_L6_PVALB_WIF1"
GABA$cluster_level[GABA$sub_clusters %in% c("GABA_10","GABA_0","GABA_11")] <- "GABA_L1_PVALB_ZNF804A"

GABA$cluster_level[GABA$clusters %in% "GABA_SST"] <- "GABA_L5_SST_GRIK1"
GABA$cluster_level[GABA$clusters %in% c("GABA_LAMP5")] <- "GABA_L1_LAMP5_RAB11FIP1"
GABA <- cell_type_v4("GABA", c(3))



load(file = paste0(output_dir,"Vasc_filtered_3.RData"));Vasc <- obj3
load(file = paste0(output_dir,"Micro_filtered_3.RData"));Micro <- obj3
load(file = paste0(output_dir,"Olig_filtered_3.RData")) ; Olig <- obj3
load(file = paste0(output_dir,"OPC_filtered_2.RData"));OPC <- obj3
load(file = paste0(output_dir,"Astro_filtered_2.RData")) ; Astro <- obj3
rm(obj3)


Vasc$clusters <- "Endo"
Vasc$clusters[Vasc$sub_clusters %in% c("Vasc_1")] <- "Peri"
Vasc$cluster_level <- Vasc$clusters


Micro$clusters <- Micro$sub_clusters
Micro$clusters <- "Micro"
Micro$cluster_level <- Micro$clusters

Olig$clusters <- Olig$sub_clusters
Olig$clusters <- "Olig"

Olig$cluster_level <- Olig$clusters

OPC$clusters <- OPC$sub_clusters
OPC$clusters <- "OPC"
OPC$cluster_level <- OPC$clusters

Astro$clusters <- Astro$sub_clusters
Astro$clusters <- "Astro"
Astro$cluster_level <- Astro$clusters

Idents(GABA) <- GABA$sub_clusters
dp <- DotPlot(GABA, features = c("GAD1","GAD2","PVALB","SST","CHODL","VIP","LAMP5","NDNF"))
png(paste0(output_dir,"DOTPLOT4_GABA_6.png"), height=10,width=9,units='in',res=800)
print(dp)
dev.off()

load(file = paste0(output_dir,"GABA_filtered_4.RData"));GABA <- obj3


GABA$clusters <- GABA$sub_clusters
GABA$clusters[GABA$sub_clusters %in% c("GABA_7","GABA_2","GABA_11","GABA_3","GABA_15","GABA_20")] <- "GABA_SST"
GABA$clusters[GABA$sub_clusters %in% c("GABA_8","GABA_4")] <- "GABA_LAMP5"
GABA$clusters[GABA$sub_clusters %in% c("GABA_9","GABA_0","GABA_6")] <- "GABA_PVALB"
GABA$clusters[GABA$sub_clusters %in% c("GABA_5","GABA_1","GABA_12","GABA_13","GABA_14")] <- "GABA_VIP"
GABA$clusters[GABA$sub_clusters %in% c("GABA_16")] <- "GABA_NDNF"
GABA$clusters[GABA$sub_clusters %in% c("GABA_19")] <- "GABA_SST_CHODL"
GABA$clusters[GABA$sub_clusters %in% c("GABA_17","GABA_18")] <- "GABA_PTPRK"
GABA$clusters[GABA$sub_clusters %in% c("GABA_10")] <- "GABA_DPP10"


GABA$cluster_level <- GABA$clusters
GABA$cluster_level[GABA$sub_clusters %in% c("GABA_17","GABA_18","GABA_12","GABA_10")] <- "Mixed"
GABA$cluster_level[GABA$sub_clusters %in% "GABA_18"] <- "GABA_L1_VIP"
GABA$cluster_level[GABA$sub_clusters %in% "GABA_19"] <- "GABA_L1_L6_NOS1_SST"
GABA$cluster_level[GABA$sub_clusters %in% "GABA_16"] <- "GABA_L1_CXCL14_NDNF"
GABA$cluster_level[GABA$sub_clusters %in% c("GABA_1")] <- "GABA_L1_L2_ASIC2_VIP"
GABA$cluster_level[GABA$sub_clusters %in% c("GABA_13")] <- "GABA_L3_L5_NCAM2_VIP"
GABA$cluster_level[GABA$sub_clusters %in% c("GABA_15")] <- "GABA_L1_L2_FBN2_SST"
GABA$cluster_level[GABA$sub_clusters %in% c("GABA_6")] <- "GABA_L3_L5_SULF1_PVALB"
GABA$cluster_level[GABA$sub_clusters %in% c("GABA_5")] <- "GABA_L1_L2_KCNT2_PVALB"
GABA$cluster_level[GABA$sub_clusters %in% c("GABA_0","GABA_9")] <- "GABA_L2_L5_ADAMTS17_PVALB"
GABA$cluster_level[GABA$sub_clusters %in% c("GABA_4")] <- "GABA_L1_L6_PRELID2_LAMP5"
GABA$cluster_level[GABA$sub_clusters %in% c("GABA_8")] <- "GABA_L1_L6_C8orf34_LAMP5"
GABA$cluster_level[GABA$sub_clusters %in% c("GABA_2","GABA_3")] <- "GABA_L5_GRIK1_SST"
GABA$cluster_level[GABA$sub_clusters %in% c("GABA_11")] <- "GABA_L1_L3_CDH12_SST"
GABA$cluster_level[GABA$sub_clusters %in% c("GABA_7")] <- "GABA_L5_L6_CDH9_SST"
GABA$cluster_level[GABA$sub_clusters %in% c("GABA_15")] <- "GABA_L5_PAWR_SST"
GABA$cluster_level[GABA$sub_clusters %in% c("GABA_20")] <- "GABA_L1_L6_NPY_SST"
GABA$cluster_level[GABA$sub_clusters %in% c("GABA_14")] <- "GABA_L1_L3_ZBTB20_VIP"

GABA <- subset(GABA,subset = cluster_level =="Mixed", invert=T)

load(file = paste0(output_dir,"Glut_filtered_4.RData")) ; Glut <- obj3


Glut$clusters <- Glut$sub_clusters
Glut$clusters[Glut$sub_clusters %in% c("Glut_1","Glut_5","Glut_2")] <- "Glut_CUX2"
Glut$clusters[Glut$sub_clusters %in% c("Glut_7")] <- "Glut_THEMIS"
Glut$clusters[Glut$sub_clusters %in% c("Glut_10","Glut_12")] <- "Glut_CTGF"
Glut$clusters[Glut$sub_clusters %in% c("Glut_0","Glut_3","Glut_4","Glut_8","Glut_11","Glut_14")] <- "Glut_RORB"
Glut$clusters[Glut$sub_clusters %in% c("Glut_6")] <- "Glut_CUX2_RORB"
Glut$clusters[Glut$sub_clusters %in% c("Glut_13")] <- "Glut_THEMIS_RORB"
Glut$clusters[Glut$sub_clusters %in% c("Glut_15")] <- "Glut_GRIK2_ROBO2"
Glut$clusters[Glut$sub_clusters %in% c("Glut_9")] <- "Glut_HS3ST4_SEMA3E"
Glut$clusters[Glut$sub_clusters %in% c("Glut_16")] <- "Glut_ZNF385D_TSHZ2"

Glut$cluster_level <- Glut$clusters
Glut$cluster_level[Glut$sub_clusters %in% c("Glut_2","Glut_6")] <- "Glut_L3_L5_CUX2_RORB"
Glut$cluster_level[Glut$sub_clusters %in% c("Glut_5","Glut_1")] <- "Glut_L2_L3_IT_CUX2"
Glut$cluster_level[Glut$sub_clusters %in% c("Glut_7")] <- "Glut_L5_THEMIS"
Glut$cluster_level[Glut$sub_clusters %in% c("Glut_12")] <- "Mixed"
Glut$cluster_level[Glut$sub_clusters %in% c("Glut_10")] <- "Glut_L6b_TLE4_CTGF"
Glut$cluster_level[Glut$sub_clusters %in% c("Glut_3")] <- "Glut_L3_L5_THEMIS_RORB"
Glut$cluster_level[Glut$sub_clusters %in% c("Glut_15")] <- "Glut_L5_ET_NRG1"
Glut$cluster_level[Glut$sub_clusters %in% c("Glut_9")] <- "Glut_L6_CT_SEMA3E"
Glut$cluster_level[Glut$sub_clusters %in% c("Glut_16")] <- "Glut_L5/6_NP_ASIC2"
Glut$cluster_level[Glut$sub_clusters %in% c("Glut_0","Glut_8","Glut_14")] <- "Glut_L3_L5_SPARCL1_RORB"
Glut$cluster_level[Glut$sub_clusters %in% c("Glut_11")] <- "Glut_L3_DPP10_RORB"
Glut$cluster_level[Glut$sub_clusters %in% c("Glut_4","Glut_13")] <- "Glut_L5_IT_TOX_RORB"
Glut$cluster_level[Glut$sub_clusters %in% c("Glut_3")] <- "Glut_L3_L5_LRRK1_RORB"
Glut <- subset(Glut, subset = cluster_level=="Mixed", invert=T)


###############################################
# Mering cleaned objects into final Seurat object

obj <- merge(Olig, y =c(Glut,Astro,GABA,Micro,OPC,Vasc))
save(obj, file = paste0(output_dir,"filtered_final_object_2.RData"))


obj=NormalizeData(obj)
obj=ScaleData(obj)
obj=FindVariableFeatures(obj)
obj=RunPCA(obj)
obj=RunHarmony(obj,group.by="batch",dims=1:30,reduction="pca")
obj=RunUMAP(obj,reduction="harmony",dim=1:30)
obj=FindNeighbors(obj,reduction="harmony",dims=1:30)
obj=FindClusters(obj,resolution=1.3)

dp1 <- ((DimPlot(obj,label=T, group.by="broadclass")+NoLegend())+(DimPlot(obj, group.by="seurat_clusters", label=T)+NoLegend()))/
  ((DimPlot(obj, group.by="cluster_level", repel=T,label=T)+NoLegend()))
ggsave(paste0(output_dir,"DIM_PLOT5_ALL_BATCHES_.png"), dp1, height=14, width=14)

Idents(obj) <- obj$seurat_clusters
Idents(obj) <- obj$cluster_level

dp2 <- DotPlot(obj,features=c("SNAP25","SLC17A7","RORB","THEMIS","CUX2","CTGF",
                              "GAD1","GAD2","VIP","SST","NPY","PVALB","NDNF",
                              "MOG","MAG","PDGFRA","AQP4","FGFR3","TMEM119","AIF1","CLDN5","PDGFRB",
                              "CDKN1A","CDKN2A"))+RotatedAxis()


png(paste0(output_dir,"DOT_PLOT6_ALL_BATCHES.png"), height=10,width=11,units='in',res=800)
print(dp2)
dev.off()

save(obj, file = paste0(output_dir,"filtered_final_object_v3.RData"))

