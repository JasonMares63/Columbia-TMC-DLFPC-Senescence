

# Calculting module scores for senescence flavors across cell types

flavors <- read.csv("flavors_of_sen_0911_SYMBOL.csv") %>%
  dplyr::select(-San.Diego.TMC)
modules <- colnames(flavors)

flavors <- lapply(1:ncol(flavors), function(col) {
  vals <- str_to_upper(flavors[flavors[,col]!="",col])
  genes1 <- vals[vals %in% rownames(obj)]
  genes2 <- vals[str_to_upper(vals) %in% str_to_upper(rownames(obj))]
  genes2 <- genes2[!(genes2 %in% genes1)]
  if(length(genes2)>0){
    genes2 <- sapply(genes2, function(x) rownames(obj)[which(grepl(paste0("^",x,"$"), str_to_upper(rownames(obj))))])
  }
  genes <- unique(c(genes1,genes2))
  return(genes)
}
)

names(flavors) <- modules
set.seed(100)
pool.features <- sample(setdiff(rownames(obj),unlist(flavors)),2400,replace=F)
for(i in 1:length(flavors)){
  pool.features1 <- c(pool.features, unlist(flavors[[i]]) %>% unique())
  obj1 <- AddModuleScore(obj, features=flavors[i], pool=pool.features1,name = names(flavors)[i])
  colnames(obj@meta.data)[ncol(obj@meta.data)] <- names(flavors)[i]
}

# obj2 <- subset(obj, subset = mapped_cluster3 %in% grep("Exc",levels(obj1$mapped_cluster3),value=T))
# obj2@meta.data[,67:81] <- NULL
# for(i in 1:length(flavors)){
#   pool.features1 <- c(pool.features, unlist(flavors[[i]]) %>% unique())
#   obj2 <- AddModuleScore(obj2, features=flavors[i], pool=pool.features1,name = names(flavors)[i])
#   colnames(obj2@meta.data)[ncol(obj2@meta.data)] <- names(flavors)[i]
# }

meta <- obj@meta.data[,c("cellid","broadclass","Batch","Age_at_Death","mapped_cluster3",
                         "Sex","age_group4",names(flavors))]
meta  %>%write.csv("senescence_module_scores.csv", quote=F, row.names=F)


vals <- c()
for(i in 1:length(flavors)){
  pres <- apply((obj@assays$RNA@counts[flavors[[i]],] > 0),1,sum)/ncol(obj) 
  vals <- c(vals,median(pres,na.rm=T))
}

names(vals) <- names(flavors)
vals <- data.frame(vals) %>%
  rownames_to_column("module")


thresh <- sapply(vals$module, function(x){
  quantile(meta[,x], 1-vals[vals$module==x,"vals"])
  
})
names(thresh) <- names(flavors)


flavors1 <- flavors[c(11,12,10,13,1,2,8,9)]
names(flavors1)[6] <- "SASP.1"


cols <- c("Ast"="#87CEFA",
          "End" = "#DAA520",
          "Exc"="#B22222",
          "Inh"="#00008b",
          "Mic"="#FFC0CB",
          "Opc"="#D3D3D3",
          "Oli"="#808080",
          "Peri"="#8B4513")

plts <- lapply(names(flavors1), function(x){
  vln <- VlnPlot(obj2, features = x,ncol=1,pt.size=0,group.by="broadclass")+NoLegend()+
    theme(axis.text.x = element_text(size=8,hjust=1,vjust=1,angle=45,face="bold"),
          axis.title.y = element_blank(),axis.text.y = element_text(size=6),
          axis.title.x = element_blank(),
          plot.title=element_blank())+
    geom_hline(yintercept = thresh[x], linetype="dotted",color ="black")+
    scale_fill_manual(values = cols)
  
  return(vln)
})


vln <- grid.arrange(grobs = plts, ncol = 2) ## display plot
ggsave(paste0("VlnPlots_all_flavors_broadclass.png"), vln,height=13,width=6, dpi=900)

###############
# 3+ Hallmark SnCs
meta <- meta %>%
  pivot_longer(-any_of(c("cellid","broadclass","bc","Batch","Age_at_Death",
                         "age_group4","batch","Sex","mapped_cluster3")),
               names_to="module",values_to = "exp" ) %>%
  left_join(thresh, by = "module") 
meta$sen <- sapply(1:nrow(meta), function(x) ifelse(meta$exp[x] >= meta$thresh[x],"sen","no"))

meta <- meta %>%
  filter(module %in%  keep) %>%
  mutate(module = factor(module, levels = c(keep,"3+ Hallmarks")))

multiple <- meta %>%
  group_by(cellid) %>%
  mutate(n_sen = sum(sen=="sen")) %>%
  ungroup() %>%
  filter(n_sen >= 3) %>%
  pull(cellid) %>% unique()


sub_meta <- meta %>% distinct(cellid,.keep_all=T)
sub_meta$module <- "3+ Hallmarks"
sub_meta$thresh <- NA
sub_meta$expre <- NA
sub_meta$sen <- sapply(1:nrow(sub_meta), function(x) ifelse(sub_meta$cellid[x] %in% multiple,"sen","no"))

meta <- rbind(meta,sub_meta)

################
# Enrichment of SnCs across cell types


fish_res_bc <- map_dfr(unique(meta$module), function(m){
  meta1 <- meta %>% filter(module ==m)
  
  res <- map_dfr(unique(meta$broadclass), function(mo){
    tab <- table(meta1$broadclass==mo,meta1$sen=="sen")
    tab[tab==0] <- .5
    result <- fisher_test(tab, detailed=T,alternative="greater") %>%
      dplyr::select(-alternative, method) %>%
      mutate(module = m,
             bc =mo) %>%
      relocate(module,bc)
    return(result)
  })
  return(res)
})

fish_res_bc$p_adj <- p.adjust(fish_res_bc$p,"BH")

file <- map_dfr(unique(fish_res_bc$module), function(i){
  file <- meta %>%
    filter(module ==i) %>%
    group_by(module,broadclass) %>%
    summarise(Normal = sum(sen=="no"),Senescent = sum(sen=="sen")) %>%
    ungroup() %>%
    left_join(fish_res_bc %>% filter(module ==i) %>%
                dplyr::select(broadclass,p_adj), by = "bc") %>%
    mutate(p_adj = round(-log10(p_adj),2)) %>%
    arrange(module,bc)
  colnames(file)[5] <- "-log10(adj P value)"
  return(file)
})

file %>% write.table(paste0("Eight_modules_enrichment_results.csv"), quote=F,row.names=F)

################
# UpSet plots across senescence hallmarks
library(ComplexHeatmap)

for(i in c("Exc","Mic","Ast","Peri","End","Oli","Opc","Inh")){
  print(i)
  rid_list <- lapply(unique(meta$module), function(x){
    if(x != "3+ Hallmarks"){
      meta %>% filter(module==x, sen=="sen",bc ==i) %>% 
        pull(cellid)%>% return()
    }
  }) 
  
  min_num <- round(0.2*min(sapply(rid_list,length)))
  
  names(rid_list) <- unique(meta$module)
  rids_mat <- make_comb_mat(list_to_matrix(rid_list), mode = "intersect")
  rids_mat <- rids_mat[comb_degree(rids_mat)>=3]
  plt <- UpSet(rids_mat[comb_size(rids_mat) >= min_num],row_names_gp = grid::gpar(fontsize = 7,fontface="bold"),
               comb_order = order(comb_size(rids_mat[comb_size(rids_mat) >= min_num]),decreasing = T))
  png(paste0("Upset_Plot_broad_",i,"_v2.png"), height=3,width=6, units = "in",res=900)
  print(plt)
  dev.off()
}

#######################
### Old/Middle vs Young proportion of SnCs across cell types

prop_df <- meta %>%
  drop_na(age_group4)%>%
  dplyr::group_by(age_group4,batch, bc,module) %>%
  #dplyr::mutate(bc_count = dplyr::n()) %>%
  ungroup() %>%
  dplyr::group_by(age_group4,batch, bc,module,thresh) %>%
  dplyr::summarise(prop = 100*sum(sen=="sen")/dplyr::n())  %>%
  ungroup()


props <- prop_df %>%
  mutate(age_group4 = factor(age_group4,levels =c("Young (<=44)","Middle (45-55)","Old (>55)")))%>%
  group_by(module,bc) %>%
  t_test(prop~age_group4,ref.group="Young (<=44)",detailed=T) %>%
  ungroup() 


props$symbol <- ""
props$symbol[props$p <= .1] <- "."
props$symbol[props$p <= .05] <- "*"
props$symbol[props$p <= .01] <- "**"
props$symbol[props$p <= .001] <- "***"
props$symbol[props$p <= .0001] <- "****"

prop_heat <- prop_df %>%
  group_by(module) %>%
  summarise(prop = quantile(prop,1)*.95) %>%
  ungroup() %>%
  right_join(props, by = c("module"))



group <- "Old"
for(group in c("Old","Middle")){
heat <- props %>%
  filter(grepl(group,group2)) %>%
  mutate(statistic = -statistic) %>%
  dplyr::select(bc,module,statistic,symbol)

heat_stat <- heat %>% dplyr::select(-symbol) %>%
  pivot_wider(names_from = "broadclass",values_from="statistic")

heat_p <- heat %>% dplyr::select(-statistic) %>%
  pivot_wider(names_from = "broadclass",values_from="symbol")

heat_stat[is.na(heat_stat)] <- 0
heat_stat <- as.data.frame(heat_stat)
rownames(heat_stat) <- heat_stat$module
heat_stat$module <- NULL
clust1 <- hclust(dist(as.matrix(heat_stat)))
clust2 <- hclust(dist(t(as.matrix(heat_stat))))

clust1_levels <- rownames(heat_stat)[clust1$order]
clust2_levels <- colnames(heat_stat)[clust2$order]

clust1_levels <- c(setdiff(clust1_levels, "3+ Hallmarks"),"3+ Hallmarks")

age_plt <- props %>%
  filter(grepl(group,group2)) %>%
  mutate(mean_beta =-statistic) %>%
  mutate(mean_beta = round(mean_beta,2)) %>%
  mutate(mean_beta = ifelse(mean_beta >=2,2, mean_beta)) %>%
  mutate(mean_beta = ifelse(mean_beta <=-3,-3, mean_beta)) %>%
  mutate(module = factor(module, levels = clust1_levels,
                         labels= str_wrap(clust1_levels,10))) %>%
  mutate(bc = factor(bc, levels = clust2_levels)) %>%
  ggplot(aes(y= module, x= bc, fill = mean_beta,label=symbol)) +
  geom_tile()+
  scale_fill_gradient2(name="t-stat",low="#008ECE",mid="white",high="firebrick1",limits=c(-3,2))+
  geom_text(size=4)+
  theme_bw()+
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  theme( panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),panel.border = element_blank(),axis.ticks= element_blank(),
         axis.line = element_blank(), axis.text =element_text(size=5.5,color = "black",face="bold"),
         axis.title = element_blank(),legend.position="none")

ggsave(paste0("All_hallmarks_box_heat_",group,".png"), age_plt, height=3,width=5,dpi=750)

}


