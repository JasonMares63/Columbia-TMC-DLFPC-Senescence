library(tidyverse)
library(gtools)
library(Seurat)

######################

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

mapping_res <- read.csv("/mnt/vast/hpc/MenonLab/SenNet/Jason_work/p400_to_sennet_mapping_5sets.csv")
mapping_res$mode_cluster <- sapply(1:nrow(mapping_res), function(x) mapping_res[x,grepl("state_name",colnames(mapping_res))] %>%
                                     unlist() %>% Mode()
)

mapping_res <- mapping_res %>% filter(cellid %in% obj$cellid)
obj$mapped_cluster <- mapping_res$mode_cluster[order(match(mapping_res$cellid,obj$cellid))]

dp <- DimPlot(obj, group.by = "mapped_cluster",label=T)+NoLegend()
ggsave("Mapped_clusters_UMAP.png",dp,height=9,width=9)

mdf <- data.frame(mapped = obj$mapped_cluster) %>%
  mutate(broad = str_extract(mapped,"Ast|Exc|Inh|Mic|Oli|OPC|Peri|End")) %>%
  group_by(mapped) %>%
  mutate(mapped_n = dplyr::n()) %>%
  ungroup() %>%
  group_by(broad) %>%
  mutate(broad_n = dplyr::n()) 

mdf$broad[is.na(mdf$broad)] <- "Misc."
table(mdf$mapped[grepl("Inh",mdf$mapped)])

sub_obj <- subset(obj, subset = mapped_cluster %in% c("Inh.12","Inh.4"))
Idents(sub_obj) <- sub_obj$mapped_cluster
mark <- FindMarkers(sub_obj, test.use="roc", ident.1 = "Inh.4")

sub_obj <- subset(obj, subset = mapped_cluster %in% c("Inh.10","Inh.9"))
Idents(sub_obj) <- sub_obj$mapped_cluster
mark1 <- FindMarkers(sub_obj, test.use="roc",ident.1 = "Inh.9")

sub_obj <- subset(obj, subset = mapped_cluster %in% c("Inh.11","Inh.9"))
Idents(sub_obj) <- sub_obj$mapped_cluster
mark2 <- FindMarkers(sub_obj, test.use="roc",ident.1 = "Inh.9")

### Merge smaller clusters

mdf$mapped2 <- mdf$mapped
mdf$mapped2[mdf$mapped %in% c("Inh.1","Inh.6")] <- "Inh.1"
mdf$mapped2[mdf$mapped %in% c("Inh.13","Inh.16")] <- "Inh.13"
mdf$mapped2[mdf$mapped %in% c("Inh.14","Inh.15")] <- "Inh.14"
mdf$mapped2[mdf$mapped %in% c("Inh.5","Inh.7")] <- "Inh.7"
mdf$mapped2[mdf$mapped %in% c("Inh.9","Inh.11")] <- "Inh.9"
mdf$mapped2[mdf$mapped %in% c("Inh.2","Inh.4")] <- "Inh.2"

mdf$mapped2[mdf$mapped %in% c("Inh.8")] <- NA
mdf$mapped2[mdf$mapped %in% paste0("End.",c(1,3,4,5))] <- "End.1"
mdf$mapped2[mdf$mapped %in% c("Arteriole","Venule")] <- "End.1"
mdf$mapped2[mdf$mapped %in% c("Peri.1","Peri.2")] <- "Peri.1"

mdf$mapped2[!grepl("Ast|End|Exc|Inh|Mic|Oli|OPC|Peri",mdf$mapped2)] <- NA

mdf$mapped2[mdf$mapped %in% paste0("Ast.",c(1,10,3,5,6,7,8,9))] <- "Ast.1"
mdf$mapped2[mdf$mapped %in% paste0("Exc.",c(12,15,16))] <- "Exc.12"
mdf$mapped2[mdf$mapped %in% paste0("Mic.",c(12:16,2,5,7,8,9))] <- "Mic.2"
mdf$mapped2[mdf$mapped %in% paste0("OPC.",c(1,3))] <- "OPC.1"
mdf$mapped2[mdf$mapped %in% paste0("Oli.",c(1:3,6:12))] <- "Oli.1"

# Rename clusters by cluster size
mdf2 <- mdf %>%
  mutate(bc = str_extract(mapped2,"Ast|End|Exc|Inh|Mic|Oli|OPC|Peri")) %>%
  group_by(bc, mapped2) %>%
  summarize(mapped2_n = dplyr::n()) %>%
  ungroup() %>%
  arrange(bc, desc(mapped2_n)) %>%
  group_by(bc) %>%
  mutate(rank = 1:dplyr::n()) %>%
  mutate(mapped3 = paste0(bc,".",rank)) %>% ungroup() %>%
  select(mapped2, mapped3)

mdf <- mdf %>%
  left_join(mdf2, by= "mapped2")

mdf$mapped3[grepl("NA", mdf$mapped3)] <- NA

obj$mapped_cluster3 <- mdf$mapped3
dp <- DimPlot(obj, group.by = "mapped_cluster3",label=T)

ggsave("Mapped_clusters_UMAP_merging.png",dp,height=9,width=9)


mapped_clusters_df <- data.frame(cellid = colnames(obj), mapped_cluster = obj$mapped_cluster3)
mapped_clusters_df %>% write.csv("mapped_clusters.csv",row.names=F)


convert <-  table(obj$mapped_cluster,obj$mapped_cluster3) %>%
  as.data.frame %>%
  filter(Freq!=0) %>%
  arrange(Var2)
colnames(convert) <- c("supertype_name","internal_name","N")
convert %>% write.table("Supplementary_Table_11.csv",quote=F,row.names=F)

######################
# Plot top markers
obj$broadclass2 <- str_extract(obj$mapped_cluster3,"Ast|End|Exc|Inh|Mic|Oli|OPC|Peri")

roc_markers <- map_dfr(na.omit(unique(obj$broadclass2)), function(x){
  
  obj1 <- subset(obj, subset = broadclass2  %in% x)
  Idents(obj1) <- obj1$mapped_cluster3
  mark <- FindAllMarkers(obj1,test.use = "roc",pos.only=T)
  return(mark %>% mutate(bc = x))
}
)
Idents(obj1) <- obj1$mapped_cluster3
mark <- FindAllMarkers(obj1,test.use = "roc",pos.only=T)


marks <- roc_markers %>% filter(myAUC >= .5) %>% group_by(cluster) %>% top_n(10, myAUC) %>%
  mutate(cluster = factor(cluster , levels = gtools::mixedsort(unique(obj$mapped_cluster3)))) %>%
  arrange(cluster) 

marks %>% write.csv("ROC_markers_mapped_sub_clusters.csv")

obj1 <- subset(obj, subset = mapped_cluster3 %in% na.omit(unique(obj$mapped_cluster3)))
Idents(obj1) <- factor(obj1$mapped_cluster3, gtools::mixedsort(unique(obj1$mapped_cluster3)))


marks <- roc_markers %>% filter(myAUC >= .8, !(bc %in% c("Exc","Inh"))) %>% group_by(cluster) %>% 
  top_n(2, myAUC) %>%
  mutate(cluster = factor(cluster , levels =c(paste0("Oli.",1:3),paste0("Opc.",1:2),
                                              paste0("Ast.",1:3), paste0("Mic.",1:3),
                                              paste0("End.",1:2)))) %>%
  arrange(cluster) %>%
  pull(gene) %>% unique()

levs <- c(paste0("Exc.",1:14),paste0("Inh.",1:9),paste0("Oli.",1:3),paste0("Opc.",1:2),
          paste0("Ast.",1:3), paste0("Mic.",1:3), paste0("End.",1:2),"Peri.1")

obj1 <- subset(obj, subset = mapped_cluster3 %in% grep("Oli|Opc|Ast|Mic|End|Inh|Exc|Per",obj$mapped_cluster3, value=T))
Idents(obj1) <- factor(obj1$mapped_cluster3, levels =levs)

png("ROC_markers_mapped_clusters_all.png", height=11,width=12,units="in",res=1000)
DotPlot(obj1,features = marks )+coord_flip()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,face="bold"),
        axis.text.y = element_text(face="bold"), axis.title=element_blank())+
  scale_colour_gradient2(low="#FAFAFA",mid="#CECECE",high="firebrick1")+
  scale_size_continuous()+
  guides(colour=guide_legend(title="Average\nExpression",guide = "colorbar"),
         size = guide_legend(title="Percent\nExpressed"))
dev.off()

png("ROC_markers_mapped_clusters_Exc_Inh_classic_markers.png", height=6,width=12,units="in",res=1000)
DotPlot(obj1,features = c("PVALB","SST","CHODL","VIP","LAMP5","NDNF","CUX2","RORB","THEMIS","CTGF") )+coord_flip()+theme(axis.text.x = element_text(angle=45,vjust=0,hjust=0))
dev.off()


##################################
# MASC analysis

library(MASC)

dat <- obj@meta.data %>% as.data.frame() %>%
  dplyr::select(Batch,Sex, age_group4, mapped_cluster3, broadclass)

masc_fn <- function(var, var1,var2){
  
  dat1 <- dat[dat[[var]] %in% c(var1,var2),]
  dat1[[var]] <- factor(dat1[[var]], levels = c(var1,var2))
  masc <- MASC(dat1, cluster = dat1$broadclass, contrast = var,
               random_effects = c("Batch"),
               fixed_effects = c("Sex"))
  
  masc  %>% write.csv(paste0("MASC_output_bc_",var,"_",var1,"_",var2,".csv"))
  return()
}

masc_fn("age_group4","Young (<=44)","Old (>55)")
masc_fn("age_group4","Young (<=44)","Middle (45-55)")

masc_fn <- function(var, var1,var2){
  
  dat1 <- dat[dat[[var]] %in% c(var1,var2),]
  dat1[[var]] <- factor(dat1[[var]], levels = c(var1,var2))
  masc <- MASC(dat1, cluster = dat1$mapped_cluster3, contrast = var,
               random_effects = c("Batch"),
               fixed_effects = c("Sex"))
  
  masc  %>% write.csv(paste0("MASC_output_mc_",var,"_",var1,"_",var2,".csv"))
  return()
}

masc_fn("age_group4","Young (<=44)","Old (>55)")
masc_fn("age_group4","Young (<=44)","Middle (45-55)")

df <- map_dfr(list.files(".", pattern = "MASC_output_mc"), function(x){
  y <- read.csv(x)
  colnames(y)[5:7] <- c("Odds_ratio","lower_bound","upper_bound")
  return(y)
})

df$group <- c(rep("Middle",38),rep("Old",38))

levs <- c(paste0("Exc.",1:14),paste0("Inh.",1:9),paste0("Oli.",1:3),paste0("Opc.",1:2),
          paste0("Ast.",1:3), paste0("Mic.",1:3), paste0("End.",1:2))


df <- df %>%
  mutate(significant = ifelse(model.pvalue <= .05,"n.s.","sig"),
         cluster = str_remove(cluster,"cluster")) %>%
  filter(!grepl("Per|Misc",cluster)) %>%
  filter(upper_bound <= 2^100) %>%
  mutate(bc =str_extract(cluster,"[a-zA-Z]+"))%>%
  mutate(cluster1 = factor(str_to_title(cluster),levels = rev(levs)))

for(ag in c("Middle","Old")){
  plt <- df %>%
    filter(group ==ag) %>%
    ggplot(aes(x = cluster1, y= log2(Odds_ratio),color = cluster1))+
    geom_hline(yintercept = 0, linetype="dashed",color="black",alpha=.7)+
    theme(axis.text =element_text(color = "black"))+
    geom_errorbar(aes(ymin =log2(lower_bound),
                      ymax=log2(upper_bound),linetype=significant),width=.2, )+
    geom_point(size=3)+
    ylab("log2(Odds ratio)")+
    geom_point(aes(alpha=significant), color="white",
               size=2)+
    scale_alpha_manual(values=c("sig"=1,"n.s."=0))+
    guides(linetype="none",color="none")+
    theme_bw()+xlab("Cell type")+coord_flip()+
    theme(axis.title.y = element_blank(),axis.text = element_text(color="black",size=13),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),legend.position="none")
  ggsave(paste0("MASC-mapped-class_OR_",ag,".png"),plt, height=9,width=5,dpi=900)
  
}















