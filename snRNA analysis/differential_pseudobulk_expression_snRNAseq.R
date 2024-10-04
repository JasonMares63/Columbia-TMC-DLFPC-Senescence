
suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(Matrix)
  library(Matrix.utils)
  library(edgeR)
  library(limma)
  library(RColorBrewer)
  library(cowplot)
  library(gridExtra)
})


setwd(".")
output_dir <- "output"

obj <- readRDS(file = paste0(output_dir,"filtered_final_object_v3.RData"))

###########################################
# Pseudobulk differential expression

convert_data_to_pseudobulk <- function(seurat_object, sample_variable, #pseudo=individual level
                                       cluster_variable, # clusters
                                       subset_variable=NULL, subset_include=NULL){
  if(!is.null(subset_variable) && !is.null(subset_include)){
    eval(parse(text = paste0("seurat <- subset(seurat_object, subset = ",
                             subset_variable," %in% c('",
                             paste0(subset_include,collapse="','"),
                             "'))" )))
  } else{
    seurat <- seurat_object
  }
  expr <- FetchData(seurat, vars = cluster_variable)
  seurat <- seurat[, which(!is.na(expr[[cluster_variable]]))]
  counts <- seurat[["RNA"]]$counts
  metadata <- seurat@meta.data
  Idents(object = seurat) <- cluster_variable
  if(is.factor(Idents(seurat))){
    metadata$cluster_id <- seurat@active.ident
    
  } else{
    metadata$cluster_id <- factor(seurat@active.ident)
  }
  sce <- SingleCellExperiment(assays = list(counts = counts),colData = metadata)
  sce <- sce[rowSums(counts(sce) > 1) >= 30, ]
  
  
  kids <- purrr::set_names(levels(sce$cluster_id))
  nk <- length(kids)
  colData(sce)[[sample_variable]] <- paste0(colData(sce)[[sample_variable]],"_",sce$cluster_id)
  sids <- purrr::set_names(levels(as.factor(colData(sce)[[sample_variable]])))
  sce$sampleInfo <- factor(colData(sce)[[sample_variable]],
                           labels = unique(colData(sce)[[sample_variable]]),
                           levels= unique(colData(sce)[[sample_variable]]))
  ns <- length(sids)
  
  n_cells <- table(sce$sampleInfo) %>%  as.vector()
  names(n_cells) <- names(table(sce$sampleInfo))
  
  m <- match(names(n_cells), sce$sampleInfo)
  ei <- data.frame(colData(sce)[m, ], 
                   n_cells, row.names = NULL) %>% 
    dplyr::select("sampleInfo", "cluster_id", "n_cells")

  groups <- colData(sce)[, c("cluster_id", "sampleInfo")]
  groups$sampleInfo <- factor(groups$sampleInfo)
  
  pb <- aggregate.Matrix(t(counts(sce)), 
                         groupings = groups, fun = "sum") 
  
  splitf <- sapply(stringr::str_split(rownames(pb), 
                                      pattern = "_",n = 2), `[`, 1)
  ei$sampleInfo <- str_remove(ei$sampleInfo,"_.+$")
  pb <- split.data.frame(pb,factor(splitf)) %>%
    lapply(function(u) 
      magrittr::set_colnames(t(u), str_extract(rownames(u),paste0(unique(ei$sampleInfo),collapse="|"))))
  pb <- lapply(pb, function(u)
    u[rowSums(u > 1) >= 4, ]
    
  ) 
  output <- list(pb = pb, ei = ei)
  return(output)
}



limma_fun <- function(pb, cluster_class_specific, var_of_interest= "log2Age"){
  
  d0 <- DGEList(pb$pb[[cluster_class_specific]])
  print("NormFact")
  d0 <- calcNormFactors(d0)
  cutoff <- 1
  d <- d0
  meta <- pb$ei[pb$ei$sampleInfo %in% colnames(pb$pb[[cluster_class_specific]]),]
  model_input <- as.formula(paste0("~1+meta$",var_of_interest,"+meta$Sex"))
  mm <- model.matrix(model_input)
  y <- voom(d, mm, plot = F)
  print("lmfit")
  fit <- lmFit(y, mm,method = "robust")
  head(coef(fit))
  tmp <- eBayes(fit)
  top.table <- topTable(tmp, sort.by = "none", n = Inf, coef = paste0("meta$",var_of_interest)) %>%
    mutate(cluster = cluster_class_specific)
  return(top.table)
}



pb <- convert_data_to_pseudobulk(obj,sample_variable = "batch" , 
                                 cluster_variable= "broadclass",
                                 subset_variable="batch",
                                 subset_include=c("HN003","HN004"))
meta_file <- read_csv("metadata.csv")


pb$ei <- pb$ei %>% left_join(meta_file %>% select(CU_Core_ID,Age_at_Death,Sex),
                             by = c("sampleInfo"="CU_Core_ID"))
pb$ei <- pb$ei %>%
  distinct(sampleInfo,.keep_all=T) %>%
  mutate(Sex = factor(Sex, levels = c("M","F"), labels = c(0,1)))
pb$ei$log2Age <- log2(pb$ei$Age_at_Death)

results <- deg_analysis(pb,var_of_interest="log2Age", ctrl_vars="Sex")
results %>% write_csv(paste0(output_dir,"broadclass_log2Age.csv"))



pb <- convert_data_to_pseudobulk(subset(obj,subset = mapped_cluster3 %in% na.omit(unique(obj$mapped_cluster3))),
                                 sample_variable = "batch" , 
                                 cluster_variable= "mapped_cluster3",
                                 subset_variable="batch", subset_include=c("HN003","HN004"))
meta_file <- read_csv("metadata.csv")


pb$ei <- pb$ei %>% left_join(meta_file %>% select(CU_Core_ID,Age_at_Death,Sex), by = c("sampleInfo"="CU_Core_ID"))
pb$ei <- pb$ei %>%
  distinct(sampleInfo,.keep_all=T) %>% mutate(Sex = factor(Sex, levels = c("M","F"), labels = c(0,1)))
pb$ei$log2Age <- log2(pb$ei$Age_at_Death)

results <- deg_analysis(pb,var_of_interest="Age_at_Death", ctrl_vars="Sex")
results %>% write_csv(paste0(output_dir,"mapped_cluster3_log2Age.csv"))


###################################
# Plot expression of some markers


gene_list <- rep("Astrocytes",9)
gene_list <- c("Microglia","Microglia","NeuronGlut","NeuronGlut","NeuronGABA","NeuronGABA","NeuronGlut","OPCs","OPCs")

names(gene_list) <- c("FASTKD3","HIF1A","PLK2","RORA","ABCA1","GREB1L","TBC1D9","DOCK5")#,"TRPS1")



plt_list <- lapply(1:length(gene_list), function(x){
  print(x)
  gene_input <- names(gene_list)[x]
  bc_input <- gene_list[[x]]
  bc_name <- switch(bc_input, "Oligodendrocytes"="Oli","Microglia"="Mic","Astrocytes"="Ast","NeuronGABA"="Inh",
                    "NeuronGlut"="Exc","OPCs"="Opc")
  
  df <- data.frame(gene =1e6*pb$pb[[bc_input]][gene_input,]/apply(pb$pb[[bc_input]],2, sum),
                   age =meta_file$Age_at_Death[-c(3,4,25,26)]) %>% arrange(age) 
  plt <- df %>% ggplot(aes(x = age, y= gene))+geom_point()+#xlab("Age at death")+
    geom_smooth(method="lm",se=F)+theme_bw()+#ylab("Transcript abundance (cpm)")+
    theme(axis.text = element_text(color = "black", size=6),axis.title = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.title = element_text(size=9,color="black",face="bold"))+
    labs(title= paste0(gene_input," (",bc_name,")"))+
    scale_y_continuous(#labels = ~ sprintf(fmt = "%0.1e", .)
      labels= ~ sub('\\+0', '\\+',  sprintf(fmt="%0.1e", .)))
  
  if(x <= 5){
    plt <- plt+
      theme(axis.text.y = element_text(color = "black",size=6),axis.title = element_blank(),
            axis.text.x = element_blank(), axis.ticks.x=element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            plot.title = element_text(size=9,color="black",face="bold"))
  }
  
  return(plt)
})


plts <- grid.arrange(grobs = plt_list, nrow = 2, left = "Transcript abundance (cpm)",
                     bottom = "Age at death",as.table = FALSE) ## display plot
ggsave(paste0(output_dir,"plot_list_selected_markers.png"), plts, height=8,width=10) 





