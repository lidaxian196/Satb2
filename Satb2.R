rm(list = ls())
library(SeuratDisk)
library(patchwork)
library(Seurat)
library(ggsignif)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggsci)
library(dplyr)
library(monocle)
library(paletteer)
library(ggplot2)
library(gghalves)
library(tidyverse)
library(ggpubr)
library(VennDiagram)
library(RColorBrewer)



#Whole brain data read
scRNA <- readRDS("/public/home/liyu/stereo_seq/output_7.0/houxufenxi_xiayouwenjian/sterepy_1_qc_norm_log1p/output/seurat_1_final.rds")
scRNA.metadata <- scRNA@meta.data

#Distinguish the wild-type from the mutant according to the coordinate position
scRNA.metadata[which(scRNA.metadata$y > 12500),'group'] <- 'wt'
scRNA.metadata[which(scRNA.metadata$y < 13000),'group'] <- 'mutant'

#mutant
satb2.mutant <- subset(scRNA,group == "mutant")
satb2.mutant.metadata <- satb2.mutant@meta.data

#wt
satb2.wt <- subset(scRNA,group == "wt")
satb2.wt.metadata <- satb2.wt@meta.data

#Fig 2D
p <- ggplot(data = satb2_mutant.metadata,aes(x=x,y=y))+
      geom_point(aes(color = spatial_leiden),size = 1.5,alpha = 0.8) +
      theme_classic()+
      scale_color_manual(values=mycol)

#Fig 2E
p <- DimPlot(satb2.wt, reduction = "umap", group.by = c("spatial_leiden", "group"))



#cortex data read
cortex <- readRDS("/public/home/liyu/stereo_seq/output_7.0/houxufenxi_xiayouwenjian/sterepy_1_qc_norm_log1p/output/cortex.rds")

#wt
cortex_wt <- subset(cortex,group == "wt")
cortex_wt.metadata <- cortex_wt@meta.data

#mutant
cortex_mutant <- subset(cortex,group == "mutant")
cortex_mutant.metadata <- cortex_mutant@meta.data

#Fig 2F
p <- ggplot(data = cortex_mutant.metadata,aes(x=x,y=y))+
      geom_point(aes(color = spatial_leiden),size = 3) +
      theme_classic()+
      scale_color_manual(values=c("#FFFF99","#666666","#7570B3"))

cortex.marker <- c("Tbr1","Bcl11b","Sox5","Eomes","Pax6","Sox2")
p <- DotPlot(cortex, group.by = 'spatial_leiden',features = unique(cortex.marker)) + RotatedAxis()
exp <- p$data
library(forcats)
exp$features.plot <- as.factor(exp$features.plot) 
exp$features.plot <- fct_inorder(exp$features.plot)
p1 <- ggplot(exp,aes(x=id,y=features.plot))+
  geom_point(aes(size=`pct.exp`,color=`avg.exp.scaled`))+
  theme_classic()+ 
  coord_flip()+
  theme(panel.grid = element_blank(),
        text = element_text(size=20))+
  scale_color_gradient2(low = "blue", mid = "lightgrey", high =  "red")+labs(x=NULL,y=NULL)


#Fig 2G
p <- SpatialFeaturePlot(satb2.mutant.metadata, features = "Tbr1",alpha = c(0.1, 1),pt.size.factor = 100) + scale_fill_gradient(low="lightgrey",high="#DE1F1F")
p1 <- SpatialFeaturePlot(satb2.wt.metadata, features = "Tbr1",alpha = c(0.1, 1),pt.size.factor = 100) + scale_fill_gradient(low="lightgrey",high="#DE1F1F")


#Fig 3A
data <- as(as.matrix(cortex@assays$Spatial@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = cortex@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
print(head(fData(monocle_cds)))
expressed_genes <- row.names(subset(fData(monocle_cds),num_cells_expressed >= 10))
diff <- differentialGeneTest(monocle_cds[expressed_genes,],fullModelFormulaStr="~spatial_leiden",cores=1) 
deg <- subset(diff, qval < 0.01) 
deg <- deg[order(deg$qval,decreasing=F),]
ordergene <- rownames(deg) 
monocle_cds <- setOrderingFilter(monocle_cds, ordergene)  
ordergene <- row.names(deg)[order(deg$qval)][1:2000]
monocle_cds <- reduceDimension(monocle_cds, 
                               max_components = 2,
                               reduction_method = 'DDRTree',
                               residualModelFormulaStr = "~group")
monocle_cds <- orderCells(monocle_cds)

pdf("./train.monocle.celltype.cortex.pdf",width = 7,height = 7)
plot_cell_trajectory(monocle_cds,color_by="spatial_leiden", size=1,show_backbone=TRUE) + scale_color_manual(values=c("#FFFF99","#666666","#7570B3"))
dev.off()

rb.genes <- rownames(deg)[grep("^Rp[sl]",rownames(deg))]
deg <- deg[-which(rownames(deg) %in% rb.genes),]
deg[c("Sox4","Sox11","Sox9","Bcl11b","Rnd2","Neurog2","Mapt","Map2","Tbr1","Eomes"),]
ordergene <- row.names(deg)[order(deg$qval)][1:2500]

#Fig 3B
Time_diff <- differentialGeneTest(monocle_cds[ordergene,], cores = 5, fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] 
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
p <- plot_pseudotime_heatmap(monocle_cds[Time_genes,], 
                             num_clusters=3, 
                             show_rownames=T, return_heatmap=T,
                             hmcols = colorRampPalette(c("navy","white","firebrick3"))(1000),
                             use_gene_short_name = T)


#Fig 3C
clusters <- cutree(p$tree_row, k = 5)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering) 
gene <- as.data.frame(rownames(subset(clustering,Gene_Clusters == "3")))
colnames(gene) <- "gene"
gene <- as.list(gene)
cortex <- AddModuleScore(object = cortex,features = gene,ctrl = 100,name = 'cluster3_Features')
df$ident <- factor(df$ident,levels = c('16','11','5'))
p <- ggplot()+
    geom_half_violin(
    data = df %>% filter(split == "wt"),
    aes(x = ident,y = cluster2_Features1),colour="white",fill="#ce0665",side = "l")+
    geom_half_violin(
    data = df %>% filter(split == "mutant"),
    aes(x = ident,y = cluster2_Features1),colour="white",fill="#008000",side = "r")+
  theme_classic()+
  xlab("")+
  ylab("Gene set score")+
  geom_point(data = df, aes(x = ident,y = cluster2_Features1, fill = split),
             stat = 'summary', fun=mean,
             position = position_dodge(width = 0.2))+
  stat_summary(data = df, aes(x = ident,y = cluster2_Features1, fill = split),
               fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = 'errorbar', color='black',
               width=0.01,size=0.5,
               position = position_dodge(width = 0.2))+
  stat_compare_means(data = df, aes(x = ident,y = cluster2_Features1, fill = split),
                     # 修改显著性标注：
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "-")),
                     label = "p.signif",
                     label.y = max(df$cluster2_Features1),
                     hide.ns = F)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "top",
        legend.justification = "right")+
        ylim(-0.05,0.2)
  
#Fig 3D
CP <- read.csv("./5_deg.CSV",row.names = 1)
rb.genes <- rownames(CP)[grep("^Rp[sl]",rownames(CP))]
CP <- CP[-which(rownames(CP) %in% rb.genes),]
CP[which(CP$p_val<0.05 & CP$avg_log2FC < 0),"deg"] <- "down" 
CP[which(CP$p_val<0.05 & CP$avg_log2FC > 0),"deg"] <- "up" 
table(CP$deg)
CP$symbol <- rownames(CP)

IZ <- read.csv("./11_deg.CSV",row.names = 1)
rb.genes <- rownames(IZ)[grep("^Rp[sl]",rownames(IZ))]
IZ <- IZ[-which(rownames(IZ) %in% rb.genes),]
IZ[which(IZ$p_val<0.05 & IZ$avg_log2FC < 0),"deg"] <- "down" 
IZ[which(IZ$p_val<0.05 & IZ$avg_log2FC > 0),"deg"] <- "up" 
table(IZ$deg)
IZ$symbol <- rownames(IZ)

SVZ <- read.csv("./16_deg.CSV",row.names = 1)
rb.genes <- rownames(SVZ)[grep("^Rp[sl]",rownames(SVZ))]
SVZ <- SVZ[-which(rownames(SVZ) %in% rb.genes),]
SVZ[which(SVZ$p_val<0.05 & SVZ$avg_log2FC < 0),"deg"] <- "down" 
SVZ[which(SVZ$p_val<0.05 & SVZ$avg_log2FC > 0),"deg"] <- "up" 
table(SVZ$deg)
SVZ$symbol <- rownames(SVZ)

library(VennDiagram)
library(RColorBrewer)
col <- brewer.pal(7, "Set1")[1:2]
venn <- venn.diagram(
  x = list(c1 = rownames(subset(clustering,Gene_Clusters == "1")), 
           SVZ = rownames(subset(SVZ,deg == "up" | deg == "down"))),  
  filename = NULL,      
  fill = col,
  scaled =FALSE,
  cex = 1,
  col = NA
)

pdf("./Time_heatmapAll_cortex_c2_intersection_CP_up.pdf")
grid.draw(venn)
dev.off()

inter <- get.venn.partitions(x)
interset_mRNA_IZ<-as.data.frame(inter$..values..[1])
ego <- enrichGO(gene = interset_mRNA_IZ$X1, 
                        OrgDb = org.Mm.eg.db,
                        ont = "ALL", 
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        keyType ='SYMBOL')
ego_result <- ego@result
write.csv(ego_result,"./CP_deg_up_monocle2_c2.csv")

#Fig 3E
p <- ggplot(df, aes(Pseudotime, a,color=group)) + 
            geom_smooth()+ 
            geom_line(linewidth = NA)+
            theme(text = element_text(size=20),
            axis.ticks.length.y=unit(0.5, "cm"),
            axis.ticks.length.x=unit(0.5, "cm"),
            legend.position="none",
            axis.title.x = element_blank(),
            axis.title.y = element_blank())+ 
            theme_classic()+
            scale_color_manual(values=c("#008000","#ce0665"))













