#"Orange2","DarkOrchid"

#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())
#加载R包
#.libPaths("/home/data/t210333/R/x86_64-pc-linux-gnu-library/4.2")
#.libPaths("/usr/local/lib/R/library")
#.libPaths("/home/data/refdir/Rlib")
#.libPaths()

library(Seurat)
library(data.table)
library(dplyr)
library(patchwork)
library(Matrix)
library(readr)
library(tibble)
library(ggplot2)
library(SingleR)
library(tidyverse) 
library(future)
library(harmony)
library(pheatmap)
library(DoubletFinder)
library(devtools)
library(copykat)
library(GSVA)
library(AUCell)
library(msigdbr)
library(clusterProfiler)
library(monocle)
library(CellChat)
library(SCENIC)
library(AUCell)
library(ggpubr)

#设置工作路径
setwd("~/3/2action")
#加载功能包
source("sc_function.R")

#读入数据
if(f){
  pbmc <-fread("ResData/count.txt", sep="\t")
  #将第一列转变为行名
  pbmc <- pbmc %>% column_to_rownames("V1")
  dim(pbmc)
  #创建seurat对象
  pbmc <- CreateSeuratObject(pbmc, min.cells = 3, min.features = 200)
  dim(pbmc)
  
  #meta.data添加信息
  #读入meta信息
  meta.data <-fread("ResData/metadata.txt", sep="\t")
  #创建data.frame数据框
  proj_name <- data.frame(meta.data)
  #替换行名
  rownames(proj_name) <- row.names(pbmc@meta.data)
  #删除第一列
  proj_name <- select(proj_name,c(-1,-2))
  #增加meta注释
  pbmc <- AddMetaData(pbmc, proj_name)
  #查看
  table(pbmc@meta.data$orig.ident)
  table(pbmc@meta.data$site)
  #删除不要的样本
  pbmc <- subset(pbmc, orig.ident != "HCC07P")
  pbmc@meta.data$orig.ident <- as.character(pbmc@meta.data$orig.ident)
  pbmc <- subset(pbmc, orig.ident != "HCC08P")
  pbmc@meta.data$orig.ident <- as.character(pbmc@meta.data$orig.ident)
  pbmc <- subset(pbmc, orig.ident != "HCC10")
  pbmc@meta.data$orig.ident <- as.character(pbmc@meta.data$orig.ident)
  pbmc <- subset(pbmc, orig.ident != "HCC10L")
  pbmc@meta.data$orig.ident <- as.character(pbmc@meta.data$orig.ident)
  #查看
  table(pbmc@meta.data$orig.ident)
  table(pbmc@meta.data$site)
  
  #检测双细胞,暂未运行
  #pbmc <- findDoublets(data = pbmc, SingleRref = "Doublecell/ref_Human_all.RData")
  #table(pbmc@meta.data$DF.classifications)
  
  #保存seurat
  saveRDS(pbmc, file = "ResData/1pbmc_Seurat.rds")
}

#细胞质控
if(f){
  #读取seurat
  pbmc <- readRDS("ResData/1pbmc_Seurat.rds")
  #创建QC文件夹
  dir.create("1QC")
  #计算质控指标
  #计算细胞中线粒体基因比例
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  #计算细胞中核糖体基因比例
  pbmc[["percent.rb"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[LS]")
  #计算红细胞比例
  HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB.genes <- CaseMatch(HB.genes, rownames(pbmc))
  pbmc[["percent.HB"]]<-PercentageFeatureSet(pbmc, features=HB.genes) 
  
  #查看质控指标
  #设置绘图元素
  theme.set2 = theme(axis.title.x=element_blank())
  plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
  group = "orig.ident"
  #质控前小提琴图
  plots = list()
  for(i in seq_along(plot.featrures)){
    plots[[i]] = VlnPlot(pbmc, group.by=group, pt.size = 0,
                         features = plot.featrures[i]) + theme.set2 + NoLegend()}
  violin <- wrap_plots(plots = plots, nrow=2)    
  ggsave("1QC/1vlnplot_before_qc.pdf", plot = violin, width = 14, height = 10) 
  
  #设置质控指标
  quantile(pbmc$nFeature_RNA, seq(0.01, 0.1, 0.01))
  quantile(pbmc$nFeature_RNA, seq(0.9, 1, 0.01))
  plots[[1]] + geom_hline(yintercept = 500) + geom_hline(yintercept = 4500)
  quantile(pbmc$nCount_RNA, seq(0.9, 1, 0.01))
  plots[[2]] + geom_hline(yintercept = 22000)
  quantile(pbmc$percent.mt, seq(0.9, 1, 0.01))
  plots[[3]] + geom_hline(yintercept = 20)
  quantile(pbmc$percent.HB, seq(0.9, 1, 0.01))
  plots[[4]] + geom_hline(yintercept = 1)
  
  #质控
  #设置质控标准
  minGene=600
  maxGene=6000
  maxUMI=60000
  pctMT=10
  pctHB=1
  
  #数据质控并绘制小提琴图
  pbmc <- subset(pbmc, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & nCount_RNA < maxUMI &
                   percent.mt < pctMT & percent.HB < pctHB)
  plots = list()
  for(i in seq_along(plot.featrures)){
    plots[[i]] = VlnPlot(pbmc, group.by=group, pt.size = 0,
                         features = plot.featrures[i]) + theme.set2 + NoLegend()}
  violin <- wrap_plots(plots = plots, nrow=2)    
  ggsave("1QC/2vlnplot_after_qc.pdf", plot = violin, width = 14, height = 10) 
  dim(pbmc)
  
  #细胞周期评分
  pbmc <- NormalizeData(pbmc)  #解决每个细胞测序深度不同的问题
  g2m_genes <- cc.genes$g2m.genes
  g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(pbmc))
  s_genes <- cc.genes$s.genes    
  s_genes <- CaseMatch(search=s_genes, match=rownames(pbmc))
  pbmc <- CellCycleScoring(pbmc, g2m.features=g2m_genes, s.features=s_genes)
  colnames(pbmc@meta.data)
  table(pbmc$Phase)
  
  #保存
  saveRDS(pbmc, file = "ResData/2pbmc_qc.rds")
}

#批次矫正
if(f){
  #读取seurat
  pbmc <- readRDS("ResData/2pbmc_qc.rds")
  
  #标准流程
  #SCT
  pbmc <- SCTransform(pbmc, vars.to.regress = c("S.Score", "G2M.Score"),variable.features.n = 3000)
  #PCA
  pbmc <- RunPCA(pbmc, verbose = F)
  ElbowPlot(pbmc, ndims = 50)
  pcs = 1:30
  
  table(pbmc@meta.data$orig.ident)
  #去除批次
  dir.create("2Harmony")
  #harmony
  pbmc <- RunHarmony(pbmc, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20)
  
  #降维聚类
  pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = pcs) %>% FindClusters(resolution = 0.5)
  pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = pcs) %>% RunTSNE(reduction = "harmony", dims = pcs)
  #画图
  p <- DimPlot(pbmc, reduction = "umap", label = T)
  ggsave("2Harmony/1Cluster_reduction.png", p, width = 7, height = 6)
  p = DimPlot(pbmc,reduction = "umap",label = F,group.by = "orig.ident")
  ggsave("2Harmony/2Cluster_reduction_sample.png", p, width = 7, height = 6)
  p = DimPlot(pbmc,reduction = "umap",label = F,group.by = "site")
  ggsave("2Harmony/3Cluster_reduction_site.png", p, width = 7, height = 6)
  #保存
  saveRDS(pbmc, file = "ResData/3pbmc_Harmony.rds")
  #读取
  pbmc <- readRDS("ResData/3pbmc_Harmony.rds")
  
}

#细胞注释
if(f){
  dir.create("3Celltype")
  #Marker基因标记
  markers <- c("COL1A2","COL1A1","ACTA2", #Fib
               "FCGR2B","PECAM1","VWF", #Endo
               "CD79A","JCHAIN","IGHG1", #B
               "KIT","CLEC4C","CD1C","CD163","CD14","CD68",#Myeloid
               "NKG7","CD3E","CD3D",#T/NK
               "EPCAM","HNF4A","SERPINA1","ALB") #Hepatocyte
  #featureplot
  #将表达矩阵SCT转为RNA
  DefaultAssay(pbmc) <- "RNA"
  p <- FeaturePlot(pbmc, features = markers, ncol = 3)
  ggsave("3Celltype/1Markers_featureplot.pdf", p, width = 15, height = 32)
  #dotplot
  p <- DotPlot(pbmc, features = markers) + RotatedAxis()
  ggsave("3Celltype/2Markers_dotplot.pdf", p, width = 14, height = 6)
  #vlnplot
  p <- VlnPlot(pbmc, features = markers, stack = T, flip = T) + NoLegend()
  ggsave("3Celltype/3Markers_vlnplot.pdf", p, width = 14, height = 6)
  
  #T/NK B Myeloid Endothelial Fibroblast Hepatocyte
  pbmc$celltype.main <- recode(pbmc@meta.data$seurat_clusters,
                               "0" = "Myeloid",
                               "1" = "T/NK",
                               "2" = "Hepatocyte",
                               "3" = "T/NK",
                               "4" = "T/NK",
                               "5" = "T/NK",
                               "6" = "T/NK",
                               "7" = "Myeloid",
                               "8" = "Myeloid",
                               "9" = "Endothelial",
                               "10" = "Fibroblast",
                               "11" = "Hepatocyte",
                               "12" = "Hepatocyte",
                               "13" = "Endothelial",
                               "14" = "B",
                               "15" = "B",
                               "16" = "Hepatocyte",
                               "17" = "B",
                               "18" = "T/NK",
                               "19" = "Hepatocyte",
                               "20" = "Hepatocyte",
                               "21" = "Myeloid",
                               "22" = "22")
  #查看各个细胞群的亚型
  table(pbmc@meta.data$celltype.main)
  #画图
  DimPlot(pbmc, reduction = "umap", label = T,group.by = "celltype.main")
  #剔除23群细胞
  pbmc <- subset(pbmc, celltype.main != "22")
  pbmc$celltype.main <- as.character(pbmc$celltype.main)
  table(pbmc@meta.data$celltype.main)
  #画图
  p = DimPlot(pbmc, reduction = "umap", label = T,group.by = "celltype.main")
  ggsave("3Celltype/4Cluster_umap.png", p, width = 7, height = 6)
  p = DimPlot(pbmc, reduction = "tsne", label = T,group.by = "celltype.main")
  ggsave("3Celltype/5Cluster_tsne.png", p, width = 7, height = 6)
  #保存
  saveRDS(pbmc, file = "ResData/4pbmc_cluster.rds")
}

#肿瘤细胞鉴定
if(f){
  #读取
  pbmc <- readRDS("ResData/4pbmc_cluster.rds")
  table(pbmc@meta.data$celltype.main)
  
  #寻找肿瘤细胞
  dir.create("4Copykat")
  #提取肿瘤来源细胞和部分正常参考细胞
  pbmc1 <- subset(pbmc, celltype.main %in% 'Hepatocyte')
  pbmc2 <- subset(pbmc, celltype.main %in% "B")
  #提取用于分析的表达矩阵
  pbmc <- merge(pbmc1, pbmc2)
  counts <- as.matrix(pbmc@assays$RNA@counts)
  #设置正常参考细胞
  ref <- colnames(pbmc2)
  #运行
  res <- copykat(rawmat=counts, ngene.chr=5, norm.cell.names=ref, sam.name="all",n.cores=40)
  #保存结果
  saveRDS(res, file = "copykat.res.rds")
  
  #copykat预测结果导入seurat对象
  pbmc <- readRDS("~/3/2action/ResData/4pbmc_cluster.rds")
  malignant <- read.delim("~/3/2action/4Copykat/all_copykat_prediction.txt")
  malignant <- malignant[!duplicated(malignant$cell.names),]
  malignant <- data.frame(copykat.pred = malignant$copykat.pred, row.names = malignant$cell.names)
  pbmc <- AddMetaData(pbmc, metadata = malignant)
  
  #结果展示
  table(pbmc$copykat.pred)
  pbmc$copykat <- recode(pbmc@meta.data$copykat.pred,
                         "not.defined" = "diploid",
                         "NA" = "diploid")
  table(pbmc$copykat)
  p1 <- DimPlot(pbmc, group.by = "copykat", cols = c("red", "gray50", "gray50"))
  ggsave("4Copykat/2copykat_res.png", p1, width = 14, height = 7)
  #保存
  saveRDS(pbmc, file = "ResData/5pbmc_copykat.rds")
  
}

#单类型细胞差异
if(f){
  #读取
  pbmc <- readRDS("ResData/5pbmc_copykat.rds")
  #提取全部肝细胞
  table(pbmc@meta.data$copykat.pred)
  table(pbmc@meta.data$celltype.main)
  #修改名字
  pbmc@meta.data$copykat.pred <- recode(pbmc@meta.data$copykat.pred,
                                        "not.defined" = "Normal_Hepatocyte",
                                        "diploid" = "Normal_Hepatocyte",
                                        "aneuploid" = "Tumor_Hepatocyte")
  table(pbmc@meta.data$copykat.pred)
  table(pbmc@meta.data$celltype.main)
  #提取
  pbmc1 <- subset(pbmc, celltype.main %in% "Hepatocyte")
  table(pbmc1@meta.data$copykat.pred)
  table(pbmc1@meta.data$celltype.main)
  #保存
  saveRDS(pbmc1, file = "ResData/6pbmc_Hepatocyte.rds")
  
  #读取
  pbmc <- readRDS("ResData/6pbmc_Hepatocyte.rds")
  
  #差异基因及富集分析
  dir.create("5Diff")
  table(pbmc@meta.data$copykat.pred)
  
  #不同细胞差异基因
  options(future.globals.maxSize = 20 * 1024^3)
  plan(multisession, workers = 40)
  #分别修改矩阵及组别
  pbmc@active.assay = "RNA"
  Idents(pbmc) <- "copykat.pred"
  markers<-FindAllMarkers(pbmc,only.pos = T,min.pct = 0,logfc.threshold = 0)
  write.csv(markers,file="5Diff/markers.csv")
  
}

#拟时序分析
if(f){
  
  #剩下部分必须在R4.1.3中运行
  rm(list=ls())
  setwd("~/3/2action")
  dir.create("6Monocle")
  
  #创建monocle分析对象
  sco <- readRDS("6pbmc_Hepatocyte.rds")
  DimPlot(sco, group.by = "celltype_rename", label = T)
  DimPlot(sco, group.by = "seurat_clusters", label = T)
  sco$celltype <- sco$copykat.pred
  colnames(sco@meta.data)
  
  #monocle不推荐使用slot="data"
  data <- GetAssayData(sco, assay = "RNA", slot = "counts")
  pd <- new('AnnotatedDataFrame', data = sco@meta.data[,c(1,5,18,21)])
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  mycds <- newCellDataSet(data,
                          phenoData = pd,
                          featureData = fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
  
  #数据预处理
  mycds <- estimateSizeFactors(mycds)
  mycds <- estimateDispersions(mycds, cores=8)
  
  #选择排序基因
  disp_table <- dispersionTable(mycds)
  order.genes <- subset(disp_table, mean_expression >= 0.005 & dispersion_empirical >= 
                          1 * dispersion_fit) %>% pull(gene_id) %>% as.character()
  #增加或剔除特定基因，默认不运行
  if(F){
    order.genes.adj <- order.genes
    #增加基因
    add.genes <- c('CCR7','LEF1','TCF7','SELL')
    order.genes.adj <- unique(c(order.genes.adj, add.genes))
    #减少基因
    del.genes <- unique(c(grep("^MT-", rownames(sco), v=T), "NEAT1","TMSB4X","TMSB10"))
    order.genes.adj <- setdiff(order.genes.adj, del.genes)
    #改变高变基因
    order.genes <- order.genes.adj
  }
  mycds <- setOrderingFilter(mycds, order.genes)
  #设置排序基因
  p <- plot_ordering_genes(mycds)
  ggsave("OrderGenes.pdf", p, width = 8, height = 6)
  
  #降维排序
  mycds <- reduceDimension(mycds, max_components = 2, reduction_method = 'DDRTree', 
                           residualModelFormulaStr = "~orig.ident")
  mycds <- orderCells(mycds)
  
  #结果可视化
  #State轨迹分布图
  p1 <- plot_cell_trajectory(mycds, color_by = "State")
  ggsave("Trajectory_State.pdf", plot = p1, width = 10, height = 6.5)
  #Pseudotime轨迹图
  p2 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
  ggsave("Trajectory_Pseudotime.pdf", plot = p2, width = 10, height = 6.5)
  #Celltype轨迹分布图
  p3 <- plot_cell_trajectory(mycds, color_by = "celltype")
  ggsave("Trajectory_Celltype.pdf", plot = p3, width = 10, height = 6.5)
  #Sample轨迹分布图
  p4 <- plot_cell_trajectory(mycds, color_by = "orig.ident")
  ggsave("Trajectory_Sample.pdf", plot = p4, width = 10, height = 6.5)
  
  #调整排序重新画图（看情况而定，不是每次都要调整排序）
  if(F){
    mycds <- orderCells(mycds, root_state = 5)
    p1 <- plot_cell_trajectory(mycds, color_by = "State")
    ggsave("Trajectory_State5.pdf", plot = p1, width = 10, height = 6.5)
    #Pseudotime轨迹图
    p2 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
    ggsave("Trajectory_Pseudotime5.pdf", plot = p2, width = 10, height = 6.5)
    #Celltype轨迹分布图
    p3 <- plot_cell_trajectory(mycds, color_by = "celltype")
    ggsave("Trajectory_Celltype5.pdf", plot = p3, width = 10, height = 6.5)
    #Sample轨迹分布图
    p4 <- plot_cell_trajectory(mycds, color_by = "orig.ident")
    ggsave("Trajectory_Sample5.pdf", plot = p4, width = 10, height = 6.5)
  }
  
  #提取拟时信息给seurat对象
  pdata <- Biobase::pData(mycds)
  sco <- AddMetaData(sco, metadata = pdata[,c("Pseudotime","State")])
  
  #寻找拟时差异基因一
  #拟时差异基因分析
  Time_diff <- differentialGeneTest(mycds, cores = 10, fullModelFormulaStr = "~sm.ns(Pseudotime)")
  write.csv(Time_diff, "Time_diff_all.csv", row.names = F)
  #显著差异基因的作图
  Time_genes <- Time_diff[order(Time_diff$qval), "gene_short_name"][1:100] #提取qval最小的100个基因
  p = plot_pseudotime_heatmap(mycds[Time_genes,], num_clusters=3, show_rownames=T, return_heatmap=T)
  ggsave("Time_heatmap.pdf", p, width = 5, height = 10)
  #显著差异基因按热图结果排序并保存
  hp.genes <- p$tree_row$labels[p$tree_row$order]
  Time_diff_sig <- Time_diff[hp.genes, c("gene_short_name", "pval", "qval")]
  write.csv(Time_diff_sig, "Time_diff_sig.csv", row.names = F)
  
  #寻找拟时差异基因二
  Idents(sco) <- "State"
  deg <- FindAllMarkers(sco, logfc.threshold = 0.5, only.pos = T)
  State_genes <- group_by(deg, cluster) %>% top_n(25, avg_log2FC) %>% pull(gene) %>% unique()
  p = plot_pseudotime_heatmap(mycds[State_genes,], num_clusters=3, show_rownames=T, return_heatmap=T)
  ggsave("State_heatmap.pdf", p, width = 5, height = 10)
  #显著差异基因按热图结果排序并保存
  hp.genes <- p$tree_row$labels[p$tree_row$order]
  State_diff_sig <- deg[hp.genes, ]
  write.csv(State_diff_sig, "State_diff_sig.csv", row.names = F)
  
  #指定基因的可视化
  s.genes <- c("IFNG","FOXP3","CCR7")
  p1 <- plot_genes_jitter(mycds[s.genes,], grouping = "State", color_by = "State")
  p2 <- plot_genes_violin(mycds[s.genes,], grouping = "State", color_by = "State")
  p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "State")
  plotc <- p1|p2|p3
  ggsave("Genes_Jitterplot.pdf", plot = plotc, width = 16, height = 8)
  #模型基因的可视化
  s.genes <- c("S100A9","CFH","PPP1R16A","POP4","PRDX6","GAGE2A",
               "UGT2B15","LDHA","TMEM106C","ALAS1","RTN3","DNAJB4")
  plot_genes_in_pseudotime(mycds[s.genes,], color_by = "Pseudotime")
  ggsave("Model_gene5.pdf", plot = plotc, width = 8, height = 16)
  #保存
  saveRDS(sco, file = "7pbmc.pseudotime.rds")
  
}

#细胞通讯分析
if(F){
  
  options(stringsAsFactors = FALSE)
  dir.create("8CellChat")
  #创建CellChat对象
  sco.brca.9s <- readRDS("ResData/5pbmc_copykat.rds")
  #s.cells <- c('CAFs MSC iCAF-like','CAFs myCAF-like','Macrophage','Monocyte','NK cells','NKT cells','T cells CD4+','T cells CD8+')
  #sco <- subset(sco.brca.9s, subtype == "ER+"&celltype_minor %in% s.cells)
  #sco <- subset(sco.brca.9s, subtype == "TNBC"&celltype_minor %in% s.cells)
  table(sco.brca.9s@meta.data[["celltype.main"]])
  table(sco.brca.9s@meta.data[["copykat.pred"]])
  #导出meta.data信息
  write.table(sco.brca.9s@meta.data,"8CellChat/metadata_change.xlsx",sep="\t",quote=F,col.names=T)
  #自行修改
  #增加meta.data的信息
  metadata_celltype_finnal <-fread("8CellChat/metadata_change2.xlsx", sep="\t",header = T)
  #将第一列转变为行名
  metadata_celltype_finnal <- metadata_celltype_finnal %>% column_to_rownames("V1")
  dim(metadata_celltype_finnal)
  #添加
  sco.brca.9s <- AddMetaData(sco.brca.9s, metadata = metadata_celltype_finnal)
  table(sco.brca.9s@meta.data$celltype.final)
  #转换
  sco <- sco.brca.9s
  sco <- Seurat::NormalizeData(sco)
  sco$celltype <- sco$celltype.final
  cellchat <- createCellChat(sco@assays$RNA@data, meta = sco@meta.data, group.by = "celltype")
  table(sco$celltype,sco$orig.ident)
  #分析初始化
  levels(cellchat@idents)
  groupSize <- as.numeric(table(cellchat@idents))
  CellChatDB <- CellChatDB.human    #小鼠用CellChatDB.mouse
  showDatabaseCategory(CellChatDB)
  #可以选择数据库子集用于分析
  #CellChatDB <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat@DB <- CellChatDB
  
  #数据预处理
  #提取细胞通讯信号基因
  cellchat <- subsetData(cellchat)
  #识别过表达的配体、受体基因
  future::plan("multisession", workers = 40)
  cellchat <- identifyOverExpressedGenes(cellchat)
  #识别过表达的配体-受体对
  cellchat <- identifyOverExpressedInteractions(cellchat)
  future::plan("sequential", workers = 6)
  #数据校正（可选）
  cellchat <- projectData(cellchat, PPI.human)
  
  #计算互作概率
  #配体受体层面计算互作概率
  cellchat <- computeCommunProb(cellchat, raw.use = F) #如果用非校正数据，raw.use = T
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  #信号通路层面计算互作概率
  cellchat <- computeCommunProbPathway(cellchat)
  #提取细胞互作结果
  df.net <- subsetCommunication(cellchat)
  write.csv(df.net, "8CellChat/CommunProb.csv", row.names = F)
  df.netP <- subsetCommunication(cellchat, slot.name = "netP")
  write.csv(df.netP, "8CellChat/CommunProbPathway.csv", row.names = F)
  #聚合通讯概率或数量（细胞间会多个配体-受体互作关系，这一步是统计数量或聚合概率）
  cellchat <- aggregateNet(cellchat)
  
  #细胞互作概览
  pdf("8CellChat/fig1a_netVisual_overview_all.pdf", width = 8, height = 6)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                   label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                   label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
  
  pdf("8CellChat/fig1b_netVisual_overview_split.pdf", width = 6, height = 5)
  mat <- cellchat@net$weight
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  dev.off()
  
  #按信号通路探索细胞互作
  mypathways <- cellchat@netP$pathways
  mypathways <- mypathways[1:10] #用10个信号通路演示，实际分析时不运行这行代码！
  
  #细胞间互作概率展示
  #网络图展示
  pdf("8CellChat/fig3a_netVisual_pathways_circle.pdf", width = 6, height = 5)
  for(pathways.show in mypathways){
    netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  }
  dev.off()
  #和弦图展示互作概览
  pdf("8CellChat/fig3b_netVisual_pathways_chord.pdf", width = 10, height = 8)
  for(pathways.show in mypathways){
    netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
  }
  dev.off()
  #热图展示互作概览
  pdf("8CellChat/fig3c_netVisual_pathways_heatmap.pdf", width = 10, height = 6.5)
  for(pathways.show in mypathways){
    p <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
    plot(p)
  }
  dev.off()
  
  #信号通路内的配体-受体
  if(!file.exists("8CellChat/pathways_detail")) dir.create("8CellChat/pathways_detail")
  for(pathways.show in mypathways){
    pdf(paste0("8CellChat/pathways_detail/", pathways.show, ".pdf"), width = 8, height = 6.5)
    #配体-受体贡献度展示
    netAnalysis_contribution(cellchat, signaling = pathways.show)
    pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)$interaction_name
    for(LR.show in pairLR){
      #网络图展示细胞间的配体-受体互作
      netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
    }
    for(LR.show in pairLR){
      #和弦图展示细胞间的配体-受体互作
      netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
    }
    dev.off()
  }
  
  #按配体-受体探索细胞互作
  levels(cellchat@idents)
  
  #显示所有的细胞间配体-受体互作
  p <- netVisual_bubble(cellchat, sources.use = 1:8, targets.use = 1:8, remove.isolate = FALSE)
  ggsave("8CellChat/fig4a_CCI_all.pdf", p, width = 15, height = 50, limitsize = F)
  
  
  
  
  
  
  
  
  #指定细胞间的配体-受体互作
  
  pairLR.use <- c("FGF1_FGFR1","FGF1_FGFR2","FGF1_FGFR3","FGF1_FGFR4","FGF10_FGFR1","FGF10_FGFR2","FGF2_FGFR1","FGF2_FGFR2","FGF2_FGFR3","FGF2_FGFR4","FGF4_FGFR1","FGF4_FGFR2",
                  "FGF4_FGFR3","FGF4_FGFR4","FGF5_FGFR1","FGF5_FGFR2","FGF5_FGFR3","FGF5_FGFR4","FGF7_FGFR1","FGF7_FGFR2","FGF9_FGFR1","FGF9_FGFR2","FGF9_FGFR3","FGF9_FGFR4",
                  "FGF21_FGFR1","FGF21_FGFR3","AREG_(EGFR+ERBB2)","AREG_EGFR","EGF_(EGFR+ERBB2)","EGF_EGFR","HBEGF_(EGFR+ERBB2)","HBEGF_EGFR","PGF_VEGFR1","TGFA_(EGFR+ERBB2)",
                  "TGFA_EGFR","VEGFA_VEGFR1","VEGFA_VEGFR1R2","VEGFA_VEGFR2","VEGFB_VEGFR1","VEGFC_VEGFR2","VEGFC_VEGFR2R3","VEGFC_VEGFR3","HGF_MET")
  pairLR.use <- data.frame(interaction_name = pairLR.use)
  p <- netVisual_bubble(cellchat, sources.use = 1:6, targets.use = 7, pairLR.use =  pairLR.use, remove.isolate = FALSE)
  ggsave("8CellChat/Tumor/targets.pdf", p, width = 6, height = 10, limitsize = F)
  p <- netVisual_bubble(cellchat, sources.use = 7, targets.use = 1:6, pairLR.use =  pairLR.use, remove.isolate = FALSE)
  ggsave("8CellChat/Tumor/sources.pdf", p, width = 6, height = 10, limitsize = F)
  
  
  #指定细胞和信号通路的配体-受体互作
  p <- netVisual_bubble(cellchat, sources.use = 1:4, targets.use = 7:8, signaling = c("CCL","CXCL"), remove.isolate = FALSE)
  ggsave("8CellChat/fig4c_CCI_subcell_subpathway.pdf", p, width = 6, height = 6, limitsize = F)
  
  #指定细胞和配体-受体
  #pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CXCL"))
  pairLR.use <- c('CXCL2_CXCR2','CXCL12_CXCR4','CCL3_CCR5','CCL4_CCR5' )
  pairLR.use <- data.frame(interaction_name = pairLR.use)
  p <- netVisual_bubble(cellchat, sources.use = 1:4, targets.use = 7, pairLR.use =  pairLR.use, remove.isolate = FALSE)
  ggsave("8CellChat/fig4d_CCI_subcell_subLR.pdf", p, width = 6, height = 4, limitsize = F)
  
  #配体、受体基因表达
  p <- plotGeneExpression(cellchat, signaling = "CXCL")
  ggsave("8CellChat/fig5a_GeneExpression_violin_sig.pdf", p, width = 10, height = 9, limitsize = F)
  p <- plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)
  ggsave("8CellChat/fig5b_GeneExpression_violin_all.pdf", p, width = 10, height = 9, limitsize = F)
  
  #信号网络生态分析
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  pdf("8CellChat/fig6_signalingRole.pdf", width = 6, height = 4.5)
  for(pathways.show in mypathways){
    netAnalysis_signalingRole_network(cellchat, signaling=pathways.show, width=8, height=2.5, font.size=10)
  }
  dev.off()
  
  #保存结果
  saveRDS(cellchat, file = "ResData/8pbmc_CellChat.rds")
  
  #读取
  cellchat <- readRDS("ResData/8pbmc_CellChat.rds")
  
}

#01TCGA
if(f){
  
  #设置路径
  setwd("F:/1/3TCGAGEO\2action\01.TCGA")
  library(tidyverse)#读入meta.data文件
  json <- jsonlite::fromJSON("metadata.cart.2023-01-13.json")
  #查看
  #View(json)
  #获取需要的数据
  sample_id <- sapply(json$associated_entities,function(x){x[,1]})
  file_sample <- data.frame(sample_id,file_name=json$file_name)  
  #读入counts文件
  count_file <- list.files('gdc_download_20230113_075345.926514/',
                           pattern = '*.tsv',recursive = TRUE)
  #整理
  count_file_name <- strsplit(count_file,split='/')
  count_file_name <- sapply(count_file_name,function(x){x[2]})
  #必要时修改下面的基因数
  matrix = data.frame(matrix(nrow=60660,ncol=0))
  #下面的修改样本例数
  for (i in 1:407){
    path = paste0('gdc_download_20230113_075345.926514//',count_file[i])   #Counts文件夹名
    data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
    colnames(data)<-data[2,]
    data <-data[-c(1:6),]
    #数据类型，选择其中之一 3：unstranded；4：stranded_first；5：stranded_second；6：tpm_unstranded；7：fpkm_unstranded；8：fpkm_uq_unstranded
    #data <- data[3]
    data <- data[6]
    #data <- data[7]
    colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
    matrix <- cbind(matrix,data)
  }
  #选择3，即为counts
  #write.csv(matrix,'1TCGA/LUAD_Counts.csv',row.names = TRUE)
  #选择6，即为TPM
  #write.csv(matrix,'TPM.csv',row.names = TRUE)
  #选择7，即为FPKM
  #write.csv(matrix,'1TCGA/LUAD_FPKM.csv',row.names = TRUE)
  
  #转化为gene_symbol
  path = paste0('gdc_download_20230113_075345.926514//',count_file[1])
  data<- as.matrix(read.delim(path,fill = TRUE,header = FALSE,row.names = 1))
  gene_name <-data[-c(1:6),1]
  matrix0 <- cbind(gene_name,matrix)
  gene_type <- data[-c(1:6),2]
  matrix0 <- cbind(gene_type,matrix0)
  #将gene_name列去除重复的基因，保留每个基因最大表达量结果
  matrix0 <- aggregate( . ~ gene_name,data=matrix0, max)    
  #保留mRNA
  matrix0 <- subset(x = matrix0, gene_type == "protein_coding")
  #将gene_name列设为行名
  rownames(matrix0) <- matrix0[,1]
  matrix0 <- matrix0[,-c(1,2)]
  #write.csv(matrix0,'1TCGA/LUAD_Counts_GeneSymbol.csv',row.names = TRUE)
  #write.csv(matrix0,'1TCGA/LUAD_TPM_GeneSymbol.txt',row.names = TRUE)
  #write.csv(matrix0,'1TCGA/LUAD_FPKM_GeneSymbol.csv',row.names = TRUE)
  #write.table(matrix0,'TCGA.txt', sep="\t", quote=F, row.names = TRUE)
  
  matrix1 = data.frame(ID=rownames(matrix0),matrix0)
  #colnames(matrix1)
  name123 = gsub('[.]', '-', colnames(matrix1))
  colnames(matrix1) = name123
  #write.table(matrix1,'TCGA_TPM.txt', sep="\t", quote=F, row.names = F)
  write.table(matrix1,'TCGA.txt', sep="\t", quote=F, row.names = F)
  
}

#02GGEO
if(f){
  
  setwd("F:/1/3TCGAGEO\2action\02.GEO")
  #下载矩阵,destdir填写后会自动识别
  gset <- getGEO("GSE54129",destdir = "F:/1/3TCGAGEO\2action\02.GEO",AnnotGPL = F,getGPL = F) 
  gset
  #提取gset中第一个元素（包含基因表达矩阵和注释信息），并赋值给变量a
  a=gset[[1]]
  #提取a中的表达矩阵并赋值给dat，这时的表达矩阵的基因名字为探针ID
  dat=exprs(a)
  #展示表达矩阵的前6行，看下数据是否经过log转换，一般数据在20以内，是经过log转换的，若有成百上千的数据，表示这个数据还需要进一步log处理。
  head(dat)
  #也可以输入上面这段代码，如果数据已经经log后，显示log2 transform not needed；如果数据没有尽行log，需要log程序按照代码尽行log转换，完成后显示log2 transform finished。
  ex <- dat
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) { ex[which(ex <= 0)] <- NA
  dat <- log2(ex)
  print("log2 transform finished")}else{print("log2 transform not needed")}
  #查看a的临床信息，为后面选择用于分组的变量做准备
  pd=pData(a)
  #查看pd
  #View(pd)
  #输出临床信息
  write.csv(pd,'clinical_GSE54129.csv',row.names = TRUE)
  
  
  #另外一种方式
  #下载注释文件
  #gpl <- getGEO('GPL10558', destdir="")
  #df2=read.table("GPL25318_family.soft",sep = "\t")
  GPL=getGEO(filename = 'GPL570_family.soft.gz',destdir = "F:/1/3TCGAGEO\2action\02.GEO") 
  gpl=GPL@dataTable@table
  colnames(gpl)
  
  #我们将平台文件转为list格式，并赋值给gpl1，将gpl1保存为R文件，方便以后调用。
  #gpl1<-Table(df2)
  #查看平台文件的列名，我们看到有ID和gene symbol,记住gene symbol在第几列，我们这里在第11列。
  #colnames(Table(gpl))
  #再了解一下平台文件的数据，当然这里大家可以直接选择gene symbol字段也行，主要了解一下symbol中基因symbol值是否只有一个。
  #View(gpl1)
  #我们这里的gene symbol字段中的symbol，有的基因就不止一个名称，
  #后面有重名，我们需要第一个名字，所以需要用字符分割包再处理一下,原理同上面处理title一样。
  #ids=gpl[,c(1,7)]
  #提取平台文件gpl1中的ID和gene symbol，并赋值给probe2symbol_df
  #gpl$`Symbol`=str_split(gpl[,7],'///',simplify = T)[,1]
  probe2symbol_df<-gpl[,c(1,5)]
  ids <- probe2symbol_df
  #probe2symbol_df <- ids
  #将列名改为probe_id和symbol
  #这两步是因为我懒，不想再调整代码，为了方便后面代码运行，我改的，大家也可以不改。
  colnames(probe2symbol_df)=c('probe_id','symbol')
  #查看symbol为唯一值的个数
  length(unique(probe2symbol_df$symbol))
  #查看每个基因出现n次的个数,我们可以看到，symbol出现一次的基因有7050个，出现2次有1493个。。。
  table(sort(table(probe2symbol_df$symbol)))
  #去掉没有注释symbol的探针（其实这里没有注释的探针数量即为上面出现次数最多的基因440，也就是说有440个探针没有symbol）
  ids=probe2symbol_df[probe2symbol_df$symbol != '',]
  #%in%用于判断是否匹配，
  #注意这里dat是gset表达矩阵的数据，这一步就是把平台文件的ID和矩阵中的ID匹配。
  ids=probe2symbol_df[probe2symbol_df$probe_id %in%rownames(dat),]
  #取表达矩阵中可以与探针名匹配的那些，去掉无法匹配的表达数据
  dat=dat[ids$probe_id,]
  #ids新建mean这一列，列名为mean，同时对dat这个矩阵按行操作，取每一行（即每个样本在这个探针上的）的均值，将结果给到mean这一列的每一行，这里也可以改为中位值，median.
  ids$mean=apply(dat,1,mean)
  #即先按symbol排序，相同的symbol再按照均值从大到小排列，方便后续保留第一个值。
  ids=ids[order(ids$symbol,ids$mean,decreasing = T),]
  #将symbol这一列取出重复项，'!'为否，即取出不重复的项，去除重复的gene
  #取median的话，这样就保留每个基因最大表达量结果.最后得到n个基因。
  ids=ids[!duplicated(ids$symbol),]
  #新的ids取出probe_id这一列，将dat按照取出的这一列，用他们的每一行组成一个新的dat
  dat=dat[ids$probe_id,]
  #把ids的symbol这一列中的每一行给dat作为dat的行名
  rownames(dat)=ids$symbol
  View(dat)
  #但是我们看到矩阵最后一行，是没有symbol名字的，我们把他去掉，数字自己更改。
  dim(dat)
  #dat<-dat[-13437,]
  #最后我们把结果保存。
  write.table(data.frame(ID=rownames(dat),dat),file="GSE54129.txt", sep="\t", quote=F, row.names = F)
  
}

#03TRG
if(f){
  
  library(limma)            #引用包
  expFile="symbol.txt"      #表达数据文件
  geneFile="gene.txt"       #基因列表文件
  setwd("E:/6/1data/3TCGAGEO/2action/03.TRG")     #设置工作目录
  
  #读取输入文件，并对数据进行处理
  rt=read.table(expFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  data=data[rowMeans(data)>0,]
  
  #读取基因列表文件,提取铜死亡相关基因的表达量
  gene=read.table(geneFile, header=F, sep="\t", check.names=F)
  sameGene=intersect(as.vector(gene[,1]), rownames(data))
  geneExp=data[sameGene,]
  
  #输出结果
  outTab=rbind(ID=colnames(geneExp),geneExp)
  write.table(outTab, file="cuproptosisExp.txt", sep="\t", quote=F, col.names=F)
  
}

#04intersect
if(f){
  
  #引用包
  library(limma)
  library(sva)
  tcgaExpFile="cuproptosisExp.txt"         #TCGA表达数据文件
  geoExpFile="geneMatrix.txt"       #GEO表达数据文件
  hubGeneFile="gene.txt"        #PPI核心基因
  setwd("E:/6/1data/3TCGAGEO/2action/04.intersect")   #设置工作目录
  
  #读取TCGA基因表达文件,并对数据进行处理
  rt = read.table(tcgaExpFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  tcga=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  tcga=avereps(tcga)
  tcga=log2(tcga+1)
  
  #读取geo基因表达文件,并对数据进行处理
  rt = read.table(geoExpFile,header=T,sep="\t",check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  geo=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  geo=avereps(geo)
  geo=log2(geo+1)      #需要修改，如果数值很大，去掉前面的#，如果数值很小，保留#
  geo=normalizeBetweenArrays(geo)
  
  #对基因取交集,分别输出交集基因在TCGA矩阵和GEO矩阵的表达量
  sameGene=intersect(row.names(tcga),row.names(geo))
  tcgaOut=tcga[sameGene,]
  geoOut=geo[sameGene,]
  
  #批次矫正
  all=cbind(tcgaOut,geoOut)
  batchType=c(rep(1,ncol(tcgaOut)),rep(2,ncol(geoOut)))
  outTab=ComBat(all, batchType,par.prior=TRUE)
  tcgaOut=outTab[,colnames(tcgaOut)]
  tcgaOut[tcgaOut<0]=0
  geoOut=outTab[,colnames(geoOut)]
  geoOut[geoOut<0]=0
  
  #读取ppi网络中的核心基因
  hubGene=read.table(hubGeneFile,header=F,sep="\t",check.names=F)
  sameHubGene=intersect(row.names(tcgaOut), as.vector(hubGene[,1]))
  tcgaHubExp=tcgaOut[sameHubGene,]
  geoHubExp=geoOut[sameHubGene,]
  
  #输出缺氧基因的表达量
  tcgaHubExp=rbind(ID=colnames(tcgaHubExp),tcgaHubExp)
  write.table(tcgaHubExp,file="tcgaHypoxia.share.txt",sep="\t",quote=F,col.names=F)
  geoHubExp=rbind(ID=colnames(geoHubExp),geoHubExp)
  write.table(geoHubExp,file="geoHypoxia.share.txt",sep="\t",quote=F,col.names=F)
  
  #输出矫正后的TCGA和GEO数据
  tcgaOut=rbind(ID=colnames(tcgaOut),tcgaOut)
  write.table(tcgaOut,file="tcga.normalize.txt",sep="\t",quote=F,col.names=F)
  geoOut=rbind(ID=colnames(geoOut),geoOut)
  write.table(geoOut,file="geo.normalize.txt",sep="\t",quote=F,col.names=F)
  
  
}

#05tcgaMergeTime
if(f){
  
  library(limma)                       #引用包
  expFile="tcgaHypoxia.share.txt"      #表达数据文件
  cliFile="time.txt"                   #临床数据
  setwd("E:/6/1data/3TCGAGEO/2action/05.tcgaMergeTime")    #设置工作目录
  
  #读取表达文件，并对输入文件整理
  rt=read.table(expFile,sep="\t",header=T,check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  data=data[rowMeans(data)>0,]
  
  #删掉正常样品
  group=sapply(strsplit(colnames(data),"\\-"),"[",4)
  group=sapply(strsplit(group,""),"[",1)
  group=gsub("2","1",group)
  data=data[,group==0]
  colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(data))
  data=t(data)
  data=avereps(data)
  
  #读取生存数据
  cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)     #读取临床文件
  
  #数据合并并输出结果
  sameSample=intersect(row.names(data),row.names(cli))
  data=data[sameSample,]
  cli=cli[sameSample,]
  out=cbind(cli,data)
  out=cbind(id=row.names(out),out)
  write.table(out,file="tcga.expTime.txt",sep="\t",row.names=F,quote=F)
  
}

#06geoMergeTime
if(f){
  
  library(limma)                       #引用包
  expFile="geoHypoxia.share.txt"       #表达数据文件
  cliFile="time.txt"                   #临床数据
  setwd("E:/6/1data/3TCGAGEO/2action/06.geoMergeTime")   #设置工作目录
  
  #读取表达文件，并对输入文件整理
  rt=read.table(expFile,sep="\t",header=T,check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  data=data[rowMeans(data)>0,]
  data=t(data)
  
  #读取生存数据
  cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)     #读取临床文件
  
  #数据合并并输出结果
  sameSample=intersect(row.names(data),row.names(cli))
  data=data[sameSample,]
  cli=cli[sameSample,]
  out=cbind(cli,data)
  out=cbind(id=row.names(out),out)
  write.table(out,file="geo.expTime.txt",sep="\t",row.names=F,quote=F)
  
}

#08model
if(f){
  
  #引用包
  library(survival)
  library(caret)
  library(glmnet)
  library(survminer)
  library(timeROC)
  
  coxPfilter=0.05        #单因素cox方法显著性的过滤标准
  setwd("E:/6/1data/3TCGAGEO/2action/08.model")      #设置工作目录
  
  #读取输入文件
  rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)
  rt$futime[rt$futime<=0]=1
  rt$futime=rt$futime/365
  rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)
  
  ############定义森林图函数############
  bioForest=function(coxFile=null,forestFile=null,forestCol=null){
    #读取输入文件
    rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
    gene <- rownames(rt)
    hr <- sprintf("%.3f",rt$"HR")
    hrLow  <- sprintf("%.3f",rt$"HR.95L")
    hrHigh <- sprintf("%.3f",rt$"HR.95H")
    Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
    pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
    
    #输出图形
    pdf(file=forestFile, width=7, height=6)
    n <- nrow(rt)
    nRow <- n+1
    ylim <- c(1,nRow)
    layout(matrix(c(1,2),nc=2),width=c(3,2.5))
    
    #绘制森林图左边的临床信息
    xlim = c(0,3)
    par(mar=c(4,2.5,2,1))
    plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
    text.cex=0.8
    text(0,n:1,gene,adj=0,cex=text.cex)
    text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,adj=1)
    text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,adj=1,)
    
    #绘制森林图
    par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
    LOGindex = 10 
    hrLow = log(as.numeric(hrLow),LOGindex)
    hrHigh = log(as.numeric(hrHigh),LOGindex)
    hr = log(as.numeric(hr),LOGindex)
    xlim = c(floor(min(hrLow,hrHigh)),ceiling(max(hrLow,hrHigh)))
    plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
    arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
    abline(v=log(1,LOGindex),col="black",lty=2,lwd=2)
    boxcolor = ifelse(as.numeric(hr) > log(1,LOGindex), forestCol[1],forestCol[2])
    points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
    a1 = axis(1,labels=F,tick=F)
    axis(1,a1,10^a1)
    dev.off()
  }
  ############绘制森林图函数############
  
  #对数据进行分组，构建模型
  n=1      #分组的数目
  for(i in 1:n){
    #############对数据进行分组#############
    inTrain<-createDataPartition(y=rt[,2], p=0.99, list=F)
    train<-rt[inTrain,]
    test<-rt[-inTrain,]
    trainOut=cbind(id=row.names(train),train)
    testOut=cbind(id=row.names(test),test)
    
    #单因素cox分析
    outUniTab=data.frame()
    sigGenes=c("futime","fustat")
    for(i in colnames(train[,3:ncol(train)])){
      #cox分析
      cox <- coxph(Surv(futime, fustat) ~ train[,i], data = train)
      coxSummary = summary(cox)
      coxP=coxSummary$coefficients[,"Pr(>|z|)"]
      
      #保留显著性基因
      if(coxP<coxPfilter){
        sigGenes=c(sigGenes,i)
        outUniTab=rbind(outUniTab,
                        cbind(id=i,
                              HR=coxSummary$conf.int[,"exp(coef)"],
                              HR.95L=coxSummary$conf.int[,"lower .95"],
                              HR.95H=coxSummary$conf.int[,"upper .95"],
                              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
        )
      }
    }
    uniSigExp=train[,sigGenes]
    uniSigExpOut=cbind(id=row.names(uniSigExp),uniSigExp)
    if(ncol(uniSigExp)<6){next}
    
    #lasso回归
    x=as.matrix(uniSigExp[,c(3:ncol(uniSigExp))])
    y=data.matrix(Surv(uniSigExp$futime,uniSigExp$fustat))
    fit <- glmnet(x, y, family = "cox", maxit = 1000)
    cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
    coef <- coef(fit, s = cvfit$lambda.min)
    index <- which(coef != 0)
    actCoef <- coef[index]
    lassoGene=row.names(coef)[index]
    lassoSigExp=uniSigExp[,c("futime", "fustat", lassoGene)]
    lassoSigExpOut=cbind(id=row.names(lassoSigExp), lassoSigExp)
    geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
    if(nrow(geneCoef)<2){next}
    
    #############构建COX模型#############
    multiCox <- coxph(Surv(futime, fustat) ~ ., data = lassoSigExp)
    multiCox=step(multiCox, direction = "both")
    multiCoxSum=summary(multiCox)
    
    #输出模型的公式
    outMultiTab=data.frame()
    outMultiTab=cbind(
      coef=multiCoxSum$coefficients[,"coef"],
      HR=multiCoxSum$conf.int[,"exp(coef)"],
      HR.95L=multiCoxSum$conf.int[,"lower .95"],
      HR.95H=multiCoxSum$conf.int[,"upper .95"],
      pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
    outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
    outMultiTab=outMultiTab[,1:2]
    
    #输出train组风险文件
    riskScore=predict(multiCox,type="risk",newdata=train)      #利用train得到模型预测train样品风险
    coxGene=rownames(multiCoxSum$coefficients)
    coxGene=gsub("`","",coxGene)
    outCol=c("futime","fustat",coxGene)
    medianTrainRisk=median(riskScore)
    risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
    trainRiskOut=cbind(id=rownames(cbind(train[,outCol],riskScore,risk)),cbind(train[,outCol],riskScore,risk))
    
    #输出test组风险文件
    riskScoreTest=predict(multiCox,type="risk",newdata=test)     #利用train得到模型预测test样品风险
    riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
    testRiskOut=cbind(id=rownames(cbind(test[,outCol],riskScoreTest,riskTest)),cbind(test[,outCol],riskScore=riskScoreTest,risk=riskTest))
    
    #比较高低风险组的生存差异，得到差异的pvalue	
    diff=survdiff(Surv(futime, fustat) ~risk,data = train)
    pValue=1-pchisq(diff$chisq, df=1)
    diffTest=survdiff(Surv(futime, fustat) ~riskTest,data = test)
    pValueTest=1-pchisq(diffTest$chisq, df=1)
    
    
    #ROC曲线下面积
    predictTime=1    #预测时间
    roc=timeROC(T=train$futime, delta=train$fustat,
                marker=riskScore, cause=1,
                times=c(predictTime), ROC=TRUE)
    rocTest=timeROC(T=test$futime, delta=test$fustat,
                    marker=riskScoreTest, cause=1,
                    times=c(predictTime), ROC=TRUE)	
    
    if((pValue<1)){
      #输出分组结果
      write.table(trainOut,file="data.train.txt",sep="\t",quote=F,row.names=F)
      write.table(testOut,file="data.test.txt",sep="\t",quote=F,row.names=F)
      #输出单因素结果
      write.table(outUniTab,file="uni.trainCox.txt",sep="\t",row.names=F,quote=F)
      write.table(uniSigExpOut,file="uni.SigExp.txt",sep="\t",row.names=F,quote=F)
      bioForest(coxFile="uni.trainCox.txt",forestFile="uni.foreast.pdf",forestCol=c("red","green"))
      #lasso结果
      write.table(lassoSigExpOut,file="lasso.SigExp.txt",sep="\t",row.names=F,quote=F)
      pdf("lasso.lambda.pdf")
      plot(fit, xvar = "lambda", label = TRUE)
      dev.off()
      pdf("lasso.cvfit.pdf")
      plot(cvfit)
      abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
      dev.off()
      #输出多因素结果
      write.table(outMultiTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)
      write.table(trainRiskOut,file="risk.train.txt",sep="\t",quote=F,row.names=F)
      write.table(testRiskOut,file="risk.test.txt",sep="\t",quote=F,row.names=F)
      #所有样品的风险值
      allRiskOut=rbind(trainRiskOut, testRiskOut)
      write.table(allRiskOut,file="risk.all.txt",sep="\t",quote=F,row.names=F)
      break
    }
  }
  
}

#09model
if(f){
  
  #引用包
  library(glmnet)
  library(survival)
  trainFile="lasso.SigExp.txt"      #train组输入文件
  testFile="geo.expTime.txt"          #test组输入文件
  setwd("E:/6/1data/3TCGAGEO/2action/09.model")                     #设置工作目录
  rt=read.table(trainFile, header=T, sep="\t", row.names=1,check.names=F)    #读取train组输入文件
  
  #COX模型构建
  multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
  multiCox=step(multiCox,direction = "both")
  multiCoxSum=summary(multiCox)
  
  #输出模型相关信息
  outMultiTab=data.frame()
  outMultiTab=cbind(
    coef=multiCoxSum$coefficients[,"coef"],
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
  write.table(outMultiTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)
  
  #输出train组风险值
  trainScore=predict(multiCox,type="risk",newdata=rt)
  coxGene=rownames(multiCoxSum$coefficients)
  coxGene=gsub("`","",coxGene)
  outCol=c("futime","fustat",coxGene)
  risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
  outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
  write.table(cbind(id=rownames(outTab),outTab),file="trainRisk.txt",sep="\t",quote=F,row.names=F)
  
  #输出test组风险值
  rt=read.table(testFile, header=T, sep="\t", row.names=1, check.names=F)
  rt$futime=rt$futime/365
  testFinalGeneExp=rt[,coxGene]
  testScore=predict(multiCox,type="risk",newdata=rt)
  outCol=c("futime","fustat",coxGene)
  risk=as.vector(ifelse(testScore>median(testScore),"high","low"))
  outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
  write.table(cbind(id=rownames(outTab),outTab),file="testRisk.txt",sep="\t",quote=F,row.names=F)
  
  
  ############绘制森林图函数############
  bioForest=function(coxFile=null,forestFile=null,forestCol=null){
    #读取输入文件
    rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
    gene <- rownames(rt)
    hr <- sprintf("%.3f",rt$"HR")
    hrLow  <- sprintf("%.3f",rt$"HR.95L")
    hrHigh <- sprintf("%.3f",rt$"HR.95H")
    Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
    pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
    
    #输出图形
    pdf(file=forestFile, width = 8.3,height = 4.5)
    n <- nrow(rt)
    nRow <- n+1
    ylim <- c(1,nRow)
    layout(matrix(c(1,2),nc=2),width=c(3,2.5))
    
    #绘制森林图左边的临床信息
    xlim = c(0,3)
    par(mar=c(4,2.5,2,1))
    plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
    text.cex=0.8
    text(0,n:1,gene,adj=0,cex=text.cex)
    text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
    text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
    
    #绘制森林图
    par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
    xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
    plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
    arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
    abline(v=1,col="black",lty=2,lwd=2)
    boxcolor = ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])
    points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
    axis(1)
    dev.off()
  }
  ############绘制森林图函数############
  
  bioForest(coxFile="multiCox.txt",forestFile="multiCox.pdf",forestCol=c("red","green"))
  
}

#09unicox
if(f){
  
  inputFile="input.txt"        
  outFile="forest.pdf"         
  etwd("E:/6/1data/3TCGAGEO/2action/09.model/cox/unicox")  
  
  rt=read.table(inputFile,header=T,sep="\t",row.names=1,check.names=F)
  gene=rownames(rt)
  hr=sprintf("%.3f",rt$"HR")
  hrLow=sprintf("%.3f",rt$"HR.95L")
  hrHigh=sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio=paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal=ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #出图格式
  pdf(file=outFile, width = 7, height =15)
  n=nrow(rt)
  nRow=n+1
  ylim=c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2))
  
  #森林图左边的基因信息
  xlim = c(0,3)
  par(mar=c(4,2,1.5,1.5))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,1,1.5,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.03,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, 'DarkOrchid', 'Orange2')
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
  axis(1)
  dev.off()
  
}

#09muiticox
if(f){
  
  inputFile="input.txt"        
  outFile="forest.pdf"         
  setwd("E:/6/1data/3TCGAGEO/2action/09.model/cox/multiCox")
  
  rt=read.table(inputFile,header=T,sep="\t",row.names=1,check.names=F)
  gene=rownames(rt)
  hr=sprintf("%.3f",rt$"HR")
  hrLow=sprintf("%.3f",rt$"HR.95L")
  hrHigh=sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio=paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal=ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #出图格式
  pdf(file=outFile, width = 7, height =6.5)
  n=nrow(rt)
  nRow=n+1
  ylim=c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2))
  
  #森林图左边的基因信息
  xlim = c(0,3)
  par(mar=c(4,2,1.5,1.5))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,1,1.5,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.03,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, 'DarkOrchid', 'Orange2')
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
  axis(1)
  dev.off()
  
}

#10survival
if(f){
  
  #引用包
  library(survival)
  library(survminer)
  setwd("E:/6/1data/3TCGAGEO/2action/10.survival")       #设置工作目录
  
  #绘制生存曲线函数
  bioSurvival=function(inputFile=null,outFile=null){
    #读取输入文件
    rt=read.table(inputFile,header=T,sep="\t")
    #比较高低风险组生存差异，得到显著性p值
    diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
    pValue=1-pchisq(diff$chisq,df=1)
    if(pValue<0.001){
      pValue="p<0.001"
    }else{
      pValue=paste0("p=",sprintf("%.03f",pValue))
    }
    fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
    
    #绘制生存曲线
    surPlot=ggsurvplot(fit, 
                       data=rt,
                       conf.int=T,
                       pval=pValue,
                       pval.size=6,
                       risk.table=TRUE,
                       legend.labs=c("High risk", "Low risk"),
                       legend.title="Risk",
                       xlab="Time(years)",
                       break.time.by = 1,
                       risk.table.title="",
                       palette=c("DarkOrchid","Orange2"),
                       risk.table.height=.25)
    pdf(file=outFile,onefile = FALSE,width = 6.5,height =6)
    print(surPlot)
    dev.off()
  }
  bioSurvival(inputFile="trainRisk.txt",outFile="trainSurv.pdf")
  bioSurvival(inputFile="testRisk.txt",outFile="testSurv.pdf")
  
}

#11ROC
if(f){
  
  #引用包
  library(survival)
  library(survminer)
  library(timeROC)
  
  riskFile="trainRisk.txt"     #风险文件
  cliFile="clinical.txt"      #临床数据文件
  setwd("E:/6/1data/3TCGAGEO/2action/11.ROC")    #修改工作目录
  
  #读取风险输入文件
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
  risk=risk[,c("futime", "fustat", "riskScore")]
  
  #读取临床数据文件
  cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
  
  #合并数据
  samSample=intersect(row.names(risk), row.names(cli))
  risk1=risk[samSample,,drop=F]
  cli=cli[samSample,,drop=F]
  rt=cbind(risk1, cli)
  
  #定义颜色
  bioCol=c("DarkOrchid","Orange2","NavyBlue","MediumSeaGreen","Firebrick3")
  
  ######绘制1 3 5年的ROC曲线######
  ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
                 marker=risk$riskScore,cause=1,
                 weighting='aalen',
                 times=c(1,3,5),ROC=TRUE)
  pdf(file="ROC.pdf", width=5, height=5)
  plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=4)
  plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=4)
  plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=4)
  legend('bottomright',
         c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
         col=bioCol[1:3], lwd=4, bty = 'n')
  dev.off()
  
  
  ######绘制临床的ROC曲线######
  predictTime=5     #定义预测年限
  aucText=c()
  pdf(file="cliROC.pdf", width=5.5, height=5.5)
  #绘制风险得分的ROC曲线
  i=5
  ROC_rt=timeROC(T=risk$futime,
                 delta=risk$fustat,
                 marker=risk$riskScore, cause=1,
                 weighting='aalen',
                 times=c(predictTime),ROC=TRUE)
  plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=4)
  aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
  abline(0,1)
  #对临床数据进行循环，绘制临床数据的ROC曲线
  for(i in 4:ncol(rt)){
    ROC_rt=timeROC(T=rt$futime,
                   delta=rt$fustat,
                   marker=rt[,i], cause=1,
                   weighting='aalen',
                   times=c(predictTime),ROC=TRUE)
    plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=4, add=TRUE)
    aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
  }
  #绘制图例，得到ROC曲线下的面积
  legend("bottomright", aucText,lwd=4,bty="n",col=bioCol[1:(ncol(rt)-1)])
  dev.off()
  
}

#12cliGroupSur
if(f){
  
  #引用包
  library(survival)
  library(survminer)
  
  riskFile="trainRisk.txt"      #风险文件
  cliFile="clinical.txt"       #临床数据文件
  setwd("E:/6/1data/3TCGAGEO/2action/12.cliGroupSur")       #设置工作目录
  
  #读取风险文件
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
  
  #读取临床文件
  cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
  cliName=colnames(cli)[1]
  
  #数据合并
  sameSample=intersect(row.names(cli), row.names(risk))
  risk=risk[sameSample,,drop=F]
  cli=cli[sameSample,,drop=F]
  rt=cbind(futime=risk[,1], fustat=risk[,2], cli, risk[,"risk",drop=F])
  colnames(rt)=c("futime", "fustat", "clinical", "Risk")
  tab=table(rt[,"clinical"])
  tab=tab[tab!=0]
  
  #对每个临床信息里面的每个分组进行循环
  for(j in names(tab)){
    rt1=rt[(rt[,"clinical"]==j),]
    tab1=table(rt1[,"Risk"])
    tab1=tab1[tab1!=0]
    labels=names(tab1)
    if(length(labels)!=2){next}
    if((cliName=="age") | (cliName=="Age") | (cliName=="AGE")){
      titleName=paste0("age",j)
    }
    
    #计算高低风险组差异pvalue
    diff=survdiff(Surv(futime, fustat) ~Risk,data = rt1)
    pValue=1-pchisq(diff$chisq,df=1)
    if(pValue<0.001){
      pValue="p<0.001"
    }else{
      pValue=paste0("p=",sprintf("%.03f",pValue))
    }
    
    #绘制生存曲线
    fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt1)
    surPlot=ggsurvplot(fit, 
                       data=rt1,
                       conf.int=F,
                       pval=pValue,
                       pval.size=6,
                       title=paste0("Patients with ",j),
                       legend.title="Risk",
                       legend.labs=labels,
                       font.legend=12,
                       xlab="Time(years)",
                       break.time.by = 1,
                       palette=c("DarkOrchid", "Orange2"),
                       risk.table=F,
                       risk.table.title="",
                       risk.table.col = "strata",
                       risk.table.height=.25)
    
    #输出图片
    j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
    pdf(file=paste0("survival.",cliName,"_",j,".pdf"), onefile = FALSE,
        width = 6,        #图片的宽度
        height =5)        #图片的高度
    print(surPlot)
    dev.off()
  }
  
}

#13indep
if(f){
  
  library(survival)      #引用包
  setwd("E:/6/1data/3TCGAGEO/2action/13.indep")      #设置工作目录
  
  ############绘制森林图函数############
  bioForest=function(coxFile=null, forestFile=null, forestCol=null){
    #读取输入文件
    rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
    gene <- rownames(rt)
    hr <- sprintf("%.3f",rt$"HR")
    hrLow  <- sprintf("%.3f",rt$"HR.95L")
    hrHigh <- sprintf("%.3f",rt$"HR.95H")
    Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
    pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
    
    #输出图形
    pdf(file=forestFile, width=6.6, height=4.5)
    n <- nrow(rt)
    nRow <- n+1
    ylim <- c(1,nRow)
    layout(matrix(c(1,2),nc=2),width=c(3,2.5))
    
    #绘制森林图左边的临床信息
    xlim = c(0,3)
    par(mar=c(4,2.5,2,1))
    plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
    text.cex=0.8
    text(0,n:1,gene,adj=0,cex=text.cex)
    text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
    text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
    
    #绘制右边的森林图
    par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
    xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
    plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
    arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
    abline(v=1,col="black",lty=2,lwd=2)
    boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
    points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
    axis(1)
    dev.off()
  }
  ############绘制森林图函数############
  
  ############独立预后分析函数#############
  indep=function(riskFile=null, cliFile=null, project=null){
    risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #读取风险文件
    cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      #读取临床文件
    
    #数据合并
    sameSample=intersect(row.names(cli),row.names(risk))
    risk=risk[sameSample,]
    cli=cli[sameSample,]
    rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
    
    #单因素独立预后分析
    uniCoxFile=paste0(project,".uniCox.txt")
    uniCoxPdf=paste0(project,".uniCox.pdf")
    uniTab=data.frame()
    for(i in colnames(rt[,3:ncol(rt)])){
      cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
      coxSummary = summary(cox)
      uniTab=rbind(uniTab,
                   cbind(id=i,
                         HR=coxSummary$conf.int[,"exp(coef)"],
                         HR.95L=coxSummary$conf.int[,"lower .95"],
                         HR.95H=coxSummary$conf.int[,"upper .95"],
                         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
      )
    }
    write.table(uniTab,file=uniCoxFile,sep="\t",row.names=F,quote=F)
    bioForest(coxFile=uniCoxFile, forestFile=uniCoxPdf, forestCol="Orange2")
    
    
    #多因素独立预后分析
    multiCoxFile=paste0(project,".multiCox.txt")
    multiCoxPdf=paste0(project,".multiCox.pdf")
    uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
    rt1=rt[,c("futime","fustat",as.vector(uniTab[,"id"]))]
    multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
    multiCoxSum=summary(multiCox)
    multiTab=data.frame()
    multiTab=cbind(
      HR=multiCoxSum$conf.int[,"exp(coef)"],
      HR.95L=multiCoxSum$conf.int[,"lower .95"],
      HR.95H=multiCoxSum$conf.int[,"upper .95"],
      pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
    multiTab=cbind(id=row.names(multiTab),multiTab)
    write.table(multiTab, file=multiCoxFile, sep="\t", row.names=F, quote=F)
    bioForest(coxFile=multiCoxFile, forestFile=multiCoxPdf, forestCol="DarkOrchid")
  }
  ############独立预后分析函数#############
  
  #独立预后分析
  indep(riskFile="trainRisk.txt", cliFile="clinical.txt", project="all")
  
}

#14C-index
if(f){
  
  #引用包
  library(dplyr)
  library(survival)
  library(rms)
  library(pec)
  
  riskFile="trainRisk.txt"     #风险文件
  cliFile="clinical.txt"      #临床数据文件
  setwd("E:/6/1data/3TCGAGEO/2action/14.C-index")    #修改工作目录
  
  #读取风险输入文件
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
  risk=risk[,c("futime", "fustat", "riskScore")]
  
  #读取临床数据文件
  cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
  
  #合并数据
  samSample=intersect(row.names(risk), row.names(cli))
  risk1=risk[samSample,,drop=F]
  cli=cli[samSample,,drop=F]
  rt=cbind(risk1, cli)
  
  #定义颜色
  bioCol=c("DarkOrchid","Orange2","NavyBlue","MediumSeaGreen","Firebrick3")
  
  #C-index值计算
  riskScore=cph(Surv(futime,fustat)~riskScore, data=rt, surv=TRUE)
  Age=cph(Surv(futime,fustat)~Age, data=rt, surv=TRUE)
  Grade=cph(Surv(futime,fustat)~Grade, data=rt, surv=TRUE)
  Gender=cph(Surv(futime,fustat)~Gender, data=rt, surv=TRUE)
  Stage=cph(Surv(futime,fustat)~Stage, data=rt, surv=TRUE)
  c_index  <- cindex(list("Risk score"=riskScore, 
                          "Age"=Age,
                          "Gender"=Gender,
                          "Grade"=Grade,
                          "Stage"=Stage),
                     formula=Surv(futime,fustat)~ .,
                     data=rt,
                     eval.times=seq(0,12,1),
                     splitMethod="bootcv",
                     B=1000
  )
  #输出图形
  pdf(file="C-index.pdf", width=6, height=6)
  plot(c_index, xlim=c(0,12), ylim=c(0.4,0.8), col=bioCol, legend.x=6, legend.y=0.82, legend.cex=1)
  dev.off()
  
}

#15DCA
if(f){
  
  #引用包
  library(survival)
  library(ggDCA)
  riskFile="trainRisk.txt"         #风险输入文件
  cliFile="clinical.txt"      #临床数据文件
  setwd("E:/6/1data/3TCGAGEO/2action/15.DCA")    #修改工作目录
  
  #读取风险输入文件
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
  risk=risk[,c("futime", "fustat", "risk")]
  
  #读取临床数据文件
  cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
  
  #合并数据
  samSample=intersect(row.names(risk), row.names(cli))
  risk1=risk[samSample,,drop=F]
  cli=cli[samSample,,drop=F]
  rt=cbind(risk1, cli)
  rt[,"Age"]=ifelse(rt[,"Age"]>65, 1, 0)
  
  #DCA分析
  predictTime=1    #预测时间
  Risk<-coxph(Surv(futime,fustat)~risk,rt)
  Age<-coxph(Surv(futime,fustat)~Age,rt)
  Gender<-coxph(Surv(futime,fustat)~Gender,rt)
  Grade<-coxph(Surv(futime,fustat)~Grade,rt)
  Stage<-coxph(Surv(futime,fustat)~Stage,rt)
  
  #绘制决策曲线
  pdf(file="DCA.pdf", width=6.5, height=6)
  d_train=dca(Risk,Age,Gender,Grade,Stage, times=predictTime)
  ggplot(d_train, linetype=1,col=c("DarkOrchid","Orange2","NavyBlue","MediumSeaGreen","Firebrick3","DeepPink1","Sienna1"))
  dev.off()
  
}

#16Nomo
if(f){
  
  #引用包
  library(survival)
  library(regplot)
  library(rms)
  
  riskFile="trainRisk.txt"      #风险输入文件
  cliFile="clinical.txt"       #临床数据文件
  setwd("E:/6/1data/3TCGAGEO/2action/16.Nomo")    #修改工作目录
  
  #读取风险输入文件
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
  
  #读取临床数据文件
  cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
  cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
  cli$Age=as.numeric(cli$Age)
  
  #合并数据
  samSample=intersect(row.names(risk), row.names(cli))
  risk1=risk[samSample,,drop=F]
  cli=cli[samSample,,drop=F]
  rt=cbind(risk1[,c("futime", "fustat", "risk")], cli)
  
  #绘制列线图
  res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
  nom1=regplot(res.cox,
               plots = c("density", "boxes"),
               clickable=F,
               title="",
               points=TRUE,
               droplines=TRUE,
               observation=rt[50,],
               rank="sd",
               failtime = c(1,3,5),
               prfail = F)
  
  #列线图风险打分
  nomoRisk=predict(res.cox, data=rt, type="risk")
  rt=cbind(risk1, Nomogram=nomoRisk)
  outTab=rbind(ID=colnames(rt), rt)
  write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)
  
  #校准曲线
  pdf(file="calibration.pdf", width=6, height=6)
  #1年校准曲线
  f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1)
  cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
  plot(cal, xlim=c(0,1), ylim=c(0,1),
       xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=3, col="DarkOrchid", sub=F)
  #3年校准曲线
  f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=3)
  cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
  plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=3, col="Orange2", sub=F, add=T)
  #5年校准曲线
  f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=5)
  cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
  plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=3, col="NavyBlue", sub=F, add=T)
  legend('bottomright', c('1-year', '3-year', '5-year'),
         col=c("DarkOrchid","Orange2","NavyBlue"), lwd=3, bty = 'n')
  dev.off()
  
}

#17riskDiff
if(f){
  
  #引用包
  library(limma)
  expFile="mRNA.txt"          #表达数据文件
  riskFile="trainRisk.txt"       #风险文件
  logFCfilter=1                 #logFC过滤条件
  fdrFilter=0.05                #fdr过滤条件
  setwd("E:/6/1data/3TCGAGEO/2action/17.riskDiff")   #设置工作目录
  
  #读取表达数据文件,并对输入文件整理
  rt=read.table(expFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp), colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
  data=avereps(data)
  
  ##去除正常样品
  group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
  group=sapply(strsplit(group,""), "[", 1)
  group=gsub("2", "1", group)
  data=data[,group==0]
  data=t(data)
  rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
  data=avereps(data)
  data=t(data)
  
  #读取risk文件
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
  sameSample=intersect(colnames(data), row.names(risk))
  data=data[,sameSample]
  risk=risk[sameSample,]
  
  #提取low risk和high risk样品
  riskLow=risk[risk$risk=="low",]
  riskHigh=risk[risk$risk=="high",]
  dataLow=data[,row.names(riskLow)]
  dataHigh=data[,row.names(riskHigh)]
  data=cbind(dataLow,dataHigh)
  data=data[rowMeans(data)>1,]
  conNum=ncol(dataLow)
  treatNum=ncol(dataHigh)
  Type=c(rep(1,conNum), rep(2,treatNum))
  
  #差异分析
  outTab=data.frame()
  for(i in row.names(data)){
    rt=data.frame(expression=data[i,], Type=Type)
    wilcoxTest=wilcox.test(expression ~ Type, data=rt)
    pvalue=wilcoxTest$p.value
    conGeneMeans=mean(data[i,1:conNum])
    treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
    logFC=log2(treatGeneMeans)-log2(conGeneMeans)
    conMed=median(data[i,1:conNum])
    treatMed=median(data[i,(conNum+1):ncol(data)])
    diffMed=treatMed-conMed
    if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
      outTab=rbind(outTab,cbind(gene=i,lowMean=conGeneMeans,highMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
    }
  }
  pValue=outTab[,"pValue"]
  fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
  outTab=cbind(outTab, fdr=fdr)
  
  #输出差异表格
  outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
  write.table(outDiff, file="riskDiff.txt", sep="\t", row.names=F, quote=F)
  
}

#18GO
if(f){
  
  #引用包
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(circlize)
  library(RColorBrewer)
  library(dplyr)
  library(ComplexHeatmap)
  
  pvalueFilter=0.05       #p值过滤条件
  qvalueFilter=1       #矫正后的p值过滤条件
  
  #定义颜色
  colorSel="qvalue"
  if(qvalueFilter>0.05){
    colorSel="pvalue"
  }
  ontology.col=c("#00AFBB", "#E7B800", "#90EE90")
  
  setwd("E:/6/1data/3TCGAGEO/2action/18.GO/22.GO")      #设置工作目录
  rt=read.table("interGene.txt", header=F, sep="\t", check.names=F)     #读取输入文件
  
  #基因名字转换为基因id
  genes=unique(as.vector(rt[,1]))
  entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
  entrezIDs=as.character(entrezIDs)
  gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
  #gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)
  
  #GO富集分析
  kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
  GO=as.data.frame(kk)
  GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
  #保存富集结果
  write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)
  
  #定义显示GO的数目
  showNum=10
  if(nrow(GO)<30){
    showNum=nrow(GO)
  }
  
  #柱状图
  pdf(file="barplot.pdf", width=10, height=7)
  bar=barplot(kk, drop=TRUE, showCategory=showNum, label_format=30, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
  print(bar)
  dev.off()
  
  #气泡图
  pdf(file="bubble.pdf", width=10, height=7)
  bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=30, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
  print(bub)
  dev.off()
  
  ###########绘制GO圈图###########
  data=GO[order(GO$pvalue),]
  datasig=data[data$pvalue<0.05,,drop=F]
  BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
  CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
  MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]
  BP = head(BP,6)
  CC = head(CC,6)
  MF = head(MF,6)
  data = rbind(BP,CC,MF)
  main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]
  
  #整理圈图数据
  BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
  Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
  ratio = Gene/BgGene
  logpvalue = -log(data$pvalue,10)
  logpvalue.col = brewer.pal(n = 8, name = "Reds")
  f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
  BgGene.col = f(logpvalue)
  df = data.frame(GO=data$ID,start=1,end=max(BgGene))
  rownames(df) = df$GO
  bed2 = data.frame(GO=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
  bed3 = data.frame(GO=data$ID,start=1,end=Gene,BgGene=Gene)
  bed4 = data.frame(GO=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
  bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5
  
  #绘制圈图主体部分
  pdf("GO.circlize.pdf",width=10,height=10)
  par(omi=c(0.1,0.1,0.1,1.5))
  circos.par(track.margin=c(0.01,0.01))
  circos.genomicInitialize(df,plotType="none")
  circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
  }, track.height = 0.08, bg.border = NA,bg.col = main.col)
  
  for(si in get.all.sector.index()) {
    circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
                major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
  }
  f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
  circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                      panel.fun = function(region, value, ...) {
                        i = getI(...)
                        circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                           border = NA, ...)
                        circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                      })
  circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                      panel.fun = function(region, value, ...) {
                        i = getI(...)
                        circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                           border = NA, ...)
                        circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                      })
  circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                      panel.fun = function(region, value, ...) {
                        cell.xlim = get.cell.meta.data("cell.xlim")
                        cell.ylim = get.cell.meta.data("cell.ylim")
                        for(j in 1:9) {
                          y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                          circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                        }
                        circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                           border = NA, ...)
                        #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                      })
  circos.clear()
  #绘制圈图中间的图例
  middle.legend = Legend(
    labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
    type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
    title="",nrow=3,size= unit(3, "mm")
  )
  circle_size = unit(1, "snpc")
  draw(middle.legend,x=circle_size*0.42)
  #绘制GO分类的图例
  main.legend = Legend(
    labels = c("Biological Process", "Molecular Function","Cellular Component"),  type="points",pch=15,
    legend_gp = gpar(col=ontology.col), title_position = "topcenter",
    title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
    grid_width = unit(5, "mm")
  )
  #绘制pvalue的图例
  logp.legend = Legend(
    labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
    type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
    title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
    size = unit(3, "mm")
  )
  lgd = packLegend(main.legend,logp.legend)
  circle_size = unit(1, "snpc")
  print(circle_size)
  draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
  dev.off()
  
}

#19KEGG
if(f){
  
  #引用包
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(circlize)
  library(RColorBrewer)
  library(dplyr)
  library(ComplexHeatmap)
  
  pvalueFilter=0.05      #p值过滤条件
  qvalueFilter=1      #矫正后的p值过滤条件
  
  #定义颜色
  colorSel="qvalue"
  if(qvalueFilter>0.05){
    colorSel="pvalue"
  }
  
  setwd("E:/6/1data/3TCGAGEO/2action/19.KEGG")          #设置工作目录
  rt=read.table("riskDiff.txt", header=T, sep="\t", check.names=F)     #读取输入文件
  
  #提取差异基因的名称,将基因名字转换为基因id
  genes=unique(as.vector(rt[,1]))
  entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
  entrezIDs=as.character(entrezIDs)
  rt=data.frame(genes, entrezID=entrezIDs)
  gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
  #gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)
  
  #KEGG富集分析
  R.utils::setOption("clusterProfiler.download.method","auto")
  kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
  KEGG=as.data.frame(kk)
  KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
  KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
  #保存显著富集的结果
  write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)
  
  #定义显示通路的数目
  showNum=30
  if(nrow(KEGG)<showNum){
    showNum=nrow(KEGG)
  }
  
  #柱状图
  pdf(file="barplot.pdf", width=9, height=7)
  barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
  dev.off()
  
  #气泡图
  pdf(file="bubble.pdf", width = 9, height = 7)
  dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
  dev.off()
  
  ###########绘制KEGG圈图###########
  Pathway.col=c("#90EE90", "#E7B800", "#00AFBB")
  showNum=18
  data=KEGG[order(KEGG$p.adjust),]
  if(nrow(KEGG)>showNum){
    data=data[1:showNum,]
  }
  data$Pathway="KEGG"
  main.col = Pathway.col[as.numeric(as.factor(data$Pathway))]
  
  #整理圈图数据
  BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
  Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
  ratio = Gene/BgGene
  logpvalue = -log(data$pvalue,10)
  logpvalue.col = brewer.pal(n = 8, name = "Reds")
  f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
  BgGene.col = f(logpvalue)
  df = data.frame(KEGG=data$ID,start=1,end=max(BgGene))
  rownames(df) = df$KEGG
  bed2 = data.frame(KEGG=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
  bed3 = data.frame(KEGG=data$ID,start=1,end=Gene,BgGene=Gene)
  bed4 = data.frame(KEGG=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
  bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5
  
  #绘制圈图主体部分
  pdf(file="KEGG.circlize.pdf",width=10,height=10)
  par(omi=c(0.1,0.1,0.1,1.5))
  circos.par(track.margin=c(0.01,0.01))
  circos.genomicInitialize(df,plotType="none")
  circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
  }, track.height = 0.08, bg.border = NA,bg.col = main.col)
  
  for(si in get.all.sector.index()) {
    circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
                major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
  }
  f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
  circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                      panel.fun = function(region, value, ...) {
                        i = getI(...)
                        circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                           border = NA, ...)
                        circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                      })
  circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                      panel.fun = function(region, value, ...) {
                        i = getI(...)
                        circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                           border = NA, ...)
                        circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                      })
  circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                      panel.fun = function(region, value, ...) {
                        cell.xlim = get.cell.meta.data("cell.xlim")
                        cell.ylim = get.cell.meta.data("cell.ylim")
                        for(j in 1:9) {
                          y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                          circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                        }
                        circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                           border = NA, ...)
                        #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                      })
  circos.clear()
  #绘制圈图中间的图例
  middle.legend = Legend(
    labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
    type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',Pathway.col[1])),
    title="",nrow=3,size= unit(3, "mm")
  )
  circle_size = unit(1, "snpc")
  draw(middle.legend,x=circle_size*0.42)
  #绘制KEGG分类的图例
  main.legend = Legend(
    labels = c("KEGG"),  type="points",pch=15,
    legend_gp = gpar(col=Pathway.col), title_position = "topcenter",
    title = "Pathway", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
    grid_width = unit(5, "mm")
  )
  #绘制pvalue的图例
  logp.legend = Legend(
    labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
    type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
    title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
    size = unit(3, "mm")
  )
  lgd = packLegend(main.legend,logp.legend)
  circle_size = unit(1, "snpc")
  print(circle_size)
  draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
  dev.off()
  
}

#20preMaftools

#21maftools
if(f){
  
  library(maftools)       #引用包
  setwd("E:/6/1data/3TCGAGEO/2action/21.maftools")    #设置工作目录
  
  #读取风险文件
  risk=read.table("trainRisk.txt", header=T, sep="\t", check.names=F)
  outTab=risk[,c(1, ncol(risk))]
  colnames(outTab)=c("Tumor_Sample_Barcode", "Risk")
  write.table(outTab, file="ann.txt", sep="\t", quote=F, row.names=F)
  
  #读取基因突变的文件
  geneNum=20
  geneMut=read.table("geneMut.txt", header=T, sep="\t", check.names=F, row.names=1)
  gene=row.names(geneMut)[1:geneNum]
  
  #注释的颜色
  ann_colors=list()
  col=c("Orange2", "DarkOrchid")
  names(col)=c("low", "high")
  ann_colors[["Risk"]]=col
  
  #绘制低风险组瀑布图
  pdf(file="low.pdf", width=6, height=6)
  maf=read.maf(maf="low.maf", clinicalData="ann.txt")
  oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
  dev.off()
  
  #绘制高风险组瀑布图
  pdf(file="high.pdf", width=6, height=6)
  maf=read.maf(maf="high.maf", clinicalData="ann.txt")
  oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
  dev.off()
  
}

#22TIDE
if(f){
  
  #引用包
  library(limma)
  library(ggpubr)
  tideFile="TIDE.txt"          #TIDE的打分文件
  riskFile="trainRisk.txt"      #风险文件
  setwd("E:/6/1data/3TCGAGEO/2action/22.TIDE")     #设置工作目录
  
  #读取TIDE数据
  tide=read.table(tideFile, header=T, sep="\t", check.names=F, row.names=1)
  group=sapply(strsplit(row.names(tide),"\\-"), "[", 4)
  group=sapply(strsplit(group,""), "[", 1)
  group=gsub("2", "1", group)
  tide=tide[group==0,,drop=F]
  row.names(tide)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(tide))
  tide=avereps(tide)
  
  #读取风险数据文件
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
  
  #合并数据
  sameSample=intersect(row.names(tide), row.names(risk))
  tide=tide[sameSample, , drop=F]
  risk=risk[sameSample, "risk", drop=F]
  data=cbind(tide, risk)
  
  #设置比较组
  data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
  group=levels(factor(data$risk))
  data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
  group=levels(factor(data$risk))
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  #绘制小提琴图
  gg1=ggviolin(data, x="risk", y="TIDE", fill = "risk", 
               xlab="", ylab="TIDE",
               palette=c("Orange2","DarkOrchid"),
               legend.title="Risk",
               add = "boxplot", add.params = list(fill="white"))+ 
    stat_compare_means(comparisons = my_comparisons)
  
  #输出图形	
  pdf(file="TIDE.pdf", width=5, height=4.5)
  print(gg1)
  dev.off()
  
}

#23riskTMB
if(f){
  
  #引用包
  library(limma)
  library(ggpubr)
  setwd("E:/6/1data/3TCGAGEO/2action/23.riskTMB")      #设置工作目录
  
  #读取肿瘤突变负荷文件
  tmb=read.table("TMB.txt", header=T, sep="\t", check.names=F, row.names=1)
  
  #读取风险数据文件
  risk=read.table("trainRisk.txt", header=T, sep="\t", check.names=F, row.names=1)
  
  #合并数据
  sameSample=intersect(row.names(tmb), row.names(risk))
  tmb=tmb[sameSample,,drop=F]
  risk=risk[sameSample,,drop=F]
  data=cbind(tmb, risk)
  data$TMB=log2(data$TMB+1)
  
  #设置比较组
  data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
  group=levels(factor(data$risk))
  data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  #绘制小提琴图
  boxplot=ggviolin(data, x="risk", y="TMB", fill="risk",
                   xlab="",
                   ylab="Tumor tmbation burden (log2)",
                   legend.title="",
                   palette = c("Orange2","DarkOrchid"),
                   add = "boxplot", add.params = list(fill="white"))+ 
    stat_compare_means(comparisons = my_comparisons)
  
  #输出图片
  pdf(file="riskTMB.pdf", width=5, height=4.5)
  print(boxplot)
  dev.off()
  
}

#24tmbSur
if(f){
  #引用包
  library(survival)
  library(survminer)
  tmbFile="TMB.txt"            #肿瘤突变负荷文件
  riskFile="trainRisk.txt"      #风险文件
  setwd("E:/6/1data/3TCGAGEO/2action/24.tmbSur")    #修改工作目录
  
  #读取输入文件
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #读取风险文件
  tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)      #读取TMB数据文件
  
  #合并数据
  sameSample=intersect(row.names(tmb), row.names(risk))
  tmb=tmb[sameSample,,drop=F]
  risk=risk[sameSample,,drop=F]
  data=cbind(risk, tmb)
  
  #获取肿瘤突变负荷最优cutoff
  res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("TMB"))
  cutoff=as.numeric(res.cut$cutpoint[1])
  tmbType=ifelse(data[,"TMB"]<=cutoff, "L-TMB", "H-TMB")
  scoreType=ifelse(data$risk=="low", "low risk", "high risk")
  mergeType=paste0(tmbType, "+", scoreType)
  
  #定义生存分析函数
  bioSurvival=function(surData=null, outFile=null){
    diff=survdiff(Surv(futime, fustat) ~ group, data=surData)
    length=length(levels(factor(surData[,"group"])))
    pValue=1-pchisq(diff$chisq, df=length-1)
    if(pValue<0.001){
      pValue="p<0.001"
    }else{
      pValue=paste0("p=",sprintf("%.03f",pValue))
    }
    fit <- survfit(Surv(futime, fustat) ~ group, data = surData)
    #print(surv_median(fit))
    
    #绘制生存曲线
    bioCol=c("DarkOrchid","Orange2","NavyBlue","MediumSeaGreen","Firebrick3")
    bioCol=bioCol[1:length]
    surPlot=ggsurvplot(fit, 
                       data=surData,
                       conf.int=F,
                       pval=pValue,
                       pval.size=6,
                       legend.title="",
                       legend.labs=levels(factor(surData[,"group"])),
                       font.legend=10,
                       legend = c(0.8, 0.8),
                       xlab="Time(years)",
                       break.time.by = 1,
                       palette = bioCol,
                       surv.median.line = "hv",
                       risk.table=F,
                       cumevents=F,
                       risk.table.height=.25)
    #输出图形
    pdf(file=outFile, onefile = FALSE, width=5.5, height=4.8)
    print(surPlot)
    dev.off()
  }
  
  #绘制肿瘤突变负荷的生存曲线
  data$group=tmbType
  bioSurvival(surData=data, outFile="TMB.survival.pdf")
  
  #绘制肿瘤突变负荷联合高低风险的生存曲线
  data$group=mergeType
  bioSurvival(surData=data, outFile="TMB-risk.survival.pdf")
  
}

#25immuneCor
if(f){
  
  #引用包
  library(limma)
  library(scales)
  library(ggplot2)
  library(ggtext)
  riskFile="trainRisk.txt"      #风险输入文件
  immFile="infiltration_estimation_for_tcga.csv"     #免疫细胞浸润文件
  setwd("E:/6/1data/3TCGAGEO/2action/25.immuneCor")     #设置工作目录
  
  #读取风险输入文件
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
  
  #读取免疫细胞浸润文件
  immune=read.csv(immFile, header=T, sep=",", check.names=F, row.names=1)
  immune=as.matrix(immune)
  rownames(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(immune))
  immune=avereps(immune)
  
  #对风险文件和免疫细胞浸润文件取交集，得到交集样品
  sameSample=intersect(row.names(risk), row.names(immune))
  risk=risk[sameSample, "riskScore"]
  immune=immune[sameSample,]
  
  #对风险打分和免疫细胞进行相关性分析
  x=as.numeric(risk)
  outTab=data.frame()
  for(i in colnames(immune)){
    y=as.numeric(immune[,i])
    corT=cor.test(x, y, method="spearman")
    cor=corT$estimate
    pvalue=corT$p.value
    if(pvalue<0.05){
      outTab=rbind(outTab,cbind(immune=i, cor, pvalue))
    }
  }
  #输出相关性结果
  write.table(file="corResult.txt", outTab, sep="\t", quote=F, row.names=F)
  
  #绘制气泡图
  corResult=read.table("corResult.txt", head=T, sep="\t")
  corResult$Software=sapply(strsplit(corResult[,1],"_"), '[', 2)
  corResult$Software=factor(corResult$Software,level=as.character(unique(corResult$Software[rev(order(as.character(corResult$Software)))])))
  b=corResult[order(corResult$Software),]
  b$immune=factor(b$immune,levels=rev(as.character(b$immune)))
  colslabels=rep(hue_pal()(length(levels(b$Software))),table(b$Software))     #定义颜色
  pdf(file="cor.pdf", width=10, height=6.5)       #保存图片
  ggplot(data=b, aes(x=cor, y=immune, color=Software))+
    labs(x="Correlation coefficient",y="Immune cell")+
    geom_point(size=4.1)+
    theme(panel.background=element_rect(fill="white",size=1,color="black"),
          panel.grid=element_line(color="grey75",size=0.5),
          axis.ticks = element_line(size=0.5),
          axis.text.y = ggtext::element_markdown(colour=rev(colslabels)))
  dev.off()
  
}

#26immFunction
if(f){
  
  #引用包
  library(limma)
  library(GSVA)
  library(GSEABase)
  library(ggpubr)
  library(reshape2)
  expFile="symbol.txt"             #表达输入文件
  gmtFile="immune.gmt"             #免疫功能数据集文件
  riskFile="trainRisk.txt"              #风险文件
  socreFile="immFunScore.txt"      #免疫功能打分的输出文件
  setwd("E:/6/1data/3TCGAGEO/2action/26.immFunction")      #设置工作目录
  
  #读取表达输入文件，并对输入文件处理
  rt=read.table(expFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  mat=avereps(mat)
  mat=mat[rowMeans(mat)>0,]
  
  #读取数据集文件
  geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())
  
  #ssgsea分析
  ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
  #定义ssGSEA score矫正函数
  normalize=function(x){
    return((x-min(x))/(max(x)-min(x)))}
  #对ssGSEA score进行矫正
  data=normalize(ssgseaScore)
  ssgseaOut=rbind(id=colnames(data), data)
  write.table(ssgseaOut, file=socreFile, sep="\t", quote=F, col.names=F)
  
  #去除正常样品
  group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
  group=sapply(strsplit(group,""), "[", 1)
  group=gsub("2", "1", group)
  data=t(data[,group==0])
  rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
  data=avereps(data)
  
  #读取风险文件
  risk=read.table(riskFile,header=T,sep="\t",row.names=1,check.names=F)
  
  #合并数据
  sameSample=intersect(row.names(data),row.names(risk))
  data=data[sameSample,,drop=F]
  risk=risk[sameSample,"risk",drop=F]
  rt1=cbind(data, risk)
  
  #对免疫相关功能绘制箱线图
  data=melt(rt1,id.vars=c("risk"))
  colnames(data)=c("Risk","Type","Score")
  data$Risk=factor(data$Risk, levels=c("low","high"))
  p=ggboxplot(data, x="Type", y="Score", color = "Risk",
              ylab="Score",add = "none",xlab="",palette = c("Orange2","DarkOrchid") )
  p=p+rotate_x_text(50)
  p=p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
  
  #输出图片文件
  pdf(file="immFunction.pdf", width=10, height=5)
  print(p)
  dev.off()
  
}

#27checkpoint
if(f){
  
  #引用包
  library(limma)
  library(reshape2)
  library(ggplot2)
  library(ggpubr)
  expFile="symbol.txt"      #表达输入文件
  riskFile="trainRisk.txt"       #风险输入文件
  geneFile="gene.txt"       #免疫检查点的基因文件
  setwd("E:/6/1data/3TCGAGEO/2action/27.checkpoint")    #设置工作目录
  
  #读取基因表达文件,并对数据进行处理
  rt=read.table(expFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  
  #读取基因文件
  gene=read.table(geneFile, header=F, sep="\t", check.names=F)
  sameGene=intersect(row.names(data),as.vector(gene[,1]))
  data=t(data[sameGene,])
  data=log2(data+1)
  
  #删除正常样品
  group=sapply(strsplit(row.names(data),"\\-"),"[",4)
  group=sapply(strsplit(group,""),"[",1)
  group=gsub("2","1",group)
  data=data[group==0,]
  row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(data))
  data=avereps(data)
  
  #合并数据
  risk=read.table(riskFile, sep="\t", header=T, check.names=F, row.names=1)
  sameSample=intersect(row.names(data),row.names(risk))
  rt1=cbind(data[sameSample,],risk[sameSample,])
  rt1=rt1[,c(sameGene,"risk")]
  
  #提取显著差异的基因
  sigGene=c()
  for(i in colnames(rt1)[1:(ncol(rt1)-1)]){
    if(sd(rt1[,i])<0.001){next}
    wilcoxTest=wilcox.test(rt1[,i] ~ rt1[,"risk"])
    pvalue=wilcoxTest$p.value
    if(wilcoxTest$p.value<0.05){
      sigGene=c(sigGene, i)
    }
  }
  sigGene=c(sigGene, "risk")
  rt1=rt1[,sigGene]
  
  #把数据转换成ggplot2输入文件
  rt1=melt(rt1,id.vars=c("risk"))
  colnames(rt1)=c("risk","Gene","Expression")
  
  #设置比较组
  group=levels(factor(rt1$risk))
  rt1$risk=factor(rt1$risk, levels=c("low","high"))
  comp=combn(group,2)
  my_comparisons=list()
  for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
  
  #绘制箱线图
  boxplot=ggboxplot(rt1, x="Gene", y="Expression", fill="risk",
                    xlab="",
                    ylab="Gene expression",
                    legend.title="Risk",
                    width=0.8,
                    palette = c("Orange2","DarkOrchid") )+
    rotate_x_text(50)+
    stat_compare_means(aes(group=risk),
                       method="wilcox.test",
                       symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
  
  #输出图片
  pdf(file="checkpoint.diff.pdf", width=10, height=4.5)
  print(boxplot)
  dev.off()
  
}

#28pRRophetic
if(f){
  
  #引用包
  library(limma)
  library(ggpubr)
  library(pRRophetic)
  library(ggplot2)
  set.seed(12345)
  
  pFilter=0.001                 #pvalue的过滤条件
  expFile="symbol.txt"         #表达数据文件
  riskFile="trainRisk.txt"      #风险文件
  setwd("E:/6/1data/3TCGAGEO/2action/28.pRRophetic")     #设置工作目录
  allDrugs=c("A.443654", "A.770041", "ABT.263", "ABT.888", "AG.014699", "AICAR", "AKT.inhibitor.VIII", "AMG.706", "AP.24534", "AS601245", "ATRA", "AUY922", "Axitinib", "AZ628", "AZD.0530", "AZD.2281", "AZD6244", "AZD6482", "AZD7762", "AZD8055", "BAY.61.3606", "Bexarotene", "BI.2536", "BIBW2992", "Bicalutamide", "BI.D1870", "BIRB.0796", "Bleomycin", "BMS.509744", "BMS.536924", "BMS.708163", "BMS.754807", "Bortezomib", "Bosutinib", "Bryostatin.1", "BX.795", "Camptothecin", "CCT007093", "CCT018159", "CEP.701", "CGP.082996", "CGP.60474", "CHIR.99021", "CI.1040", "Cisplatin", "CMK", "Cyclopamine", "Cytarabine", "Dasatinib", "DMOG", "Docetaxel", "Doxorubicin", "EHT.1864", "Elesclomol", "Embelin", "Epothilone.B", "Erlotinib", "Etoposide", "FH535", "FTI.277", "GDC.0449", "GDC0941", "Gefitinib", "Gemcitabine", "GNF.2", "GSK269962A", "GSK.650394", "GW.441756", "GW843682X", "Imatinib", "IPA.3", "JNJ.26854165", "JNK.9L", "JNK.Inhibitor.VIII", "JW.7.52.1", "KIN001.135", "KU.55933", "Lapatinib", "Lenalidomide", "LFM.A13", "Metformin", "Methotrexate", "MG.132", "Midostaurin", "Mitomycin.C", "MK.2206", "MS.275", "Nilotinib", "NSC.87877", "NU.7441", "Nutlin.3a", "NVP.BEZ235", "NVP.TAE684", "Obatoclax.Mesylate", "OSI.906", "PAC.1", "Paclitaxel", "Parthenolide", "Pazopanib", "PD.0325901", "PD.0332991", "PD.173074", "PF.02341066", "PF.4708671", "PF.562271", "PHA.665752", "PLX4720", "Pyrimethamine", "QS11", "Rapamycin", "RDEA119", "RO.3306", "Roscovitine", "Salubrinal", "SB.216763", "SB590885", "Shikonin", "SL.0101.1", "Sorafenib", "S.Trityl.L.cysteine", "Sunitinib", "Temsirolimus", "Thapsigargin", "Tipifarnib", "TW.37", "Vinblastine", "Vinorelbine", "Vorinostat", "VX.680", "VX.702", "WH.4.023", "WO2009093972", "WZ.1.84", "X17.AAG", "X681640", "XMD8.85", "Z.LLNle.CHO", "ZM.447439")
  
  #读取表达输入文件,并对数据进行处理
  rt = read.table(expFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  data=data[rowMeans(data)>0.5,]
  
  #删掉正常样品
  group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
  group=sapply(strsplit(group,""), "[", 1)
  group=gsub("2","1",group)
  data=data[,group==0]
  data=t(data)
  rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(data))
  data=avereps(data)
  data=t(data)
  
  #读取风险输入文件
  riskRT=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
  
  for(drug in allDrugs){
    #预测药物敏感性
    sensitivity=pRRopheticPredict(data, drug, selection=1)
    sensitivity=sensitivity[sensitivity!="NaN"]
    #sensitivity[sensitivity>quantile(sensitivity,0.99)]=quantile(sensitivity,0.99)
    
    #风险文件和药物敏感性结果合并
    sameSample=intersect(row.names(riskRT), names(sensitivity))
    risk=riskRT[sameSample, "risk",drop=F]
    sensitivity=sensitivity[sameSample]
    rt=cbind(risk, sensitivity)
    
    #设置比较组
    rt$risk=factor(rt$risk, levels=c("low", "high"))
    type=levels(factor(rt[,"risk"]))
    comp=combn(type, 2)
    my_comparisons=list()
    for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
    
    #获取高低风险组差异pvalue
    test=wilcox.test(sensitivity~risk, data=rt)
    
    if(test$p.value<pFilter){
      #绘制箱线图
      boxplot=ggboxplot(rt, x="risk", y="sensitivity", fill="risk",
                        xlab="Risk",
                        ylab=paste0(drug, " sensitivity (IC50)"),
                        legend.title="Risk",
                        palette=c("Orange2","DarkOrchid")
      )+ 
        stat_compare_means(comparisons=my_comparisons)
      pdf(file=paste0("durgsensitivity.", drug, ".pdf"), width=2.3, height=4.3)
      print(boxplot)
      dev.off()
    }
  }
  
}

#29pancancer
if(f){
  
  #使用国内镜像安装包
  options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
  Sys.setenv(LANGUAGE = "en") #显示英文报错信息
  options(stringsAsFactors = FALSE) #禁止chr转成factor
  #加载包
  library(ggplot2)
  library(data.table)
  library(cowplot)
  library(ggpubr)
  library(GSVA)
  library(SimDesign)
  library(tidyr)
  setwd("E:/已整完/21肝癌/1data/3TCGAGEO/2action/29.pancancer")
  #自定义函数将gmt文件读取为list
  gmt2list <- function(annofile){
    if (!file.exists(annofile)) {
      stop("There is no such gmt file.")
    }
    
    if (tools::file_ext(annofile) == "xz") {
      annofile <- xzfile(annofile)
      x <- scan(annofile, what="", sep="\n", quiet=TRUE)
      close(annofile)
    } else if (tools::file_ext(annofile) == "gmt") {
      x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    } else {
      stop ("Only gmt and gmt.xz are accepted for gmt2list")
    }
    
    y <- strsplit(x, "\t")
    names(y) <- sapply(y, `[[`, 1)
    
    annoList <- lapply(y, `[`, c(-1,-2))
  }
  
  # 读取风险基因以及对应系数（来自原文补充材料表格Table S6）
  risk.coeff <- read.table("multiCox.txt",sep = "\t", row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
  # 读取肿瘤注释文件
  rawAnno <- read.delim("merged_sample_quality_annotations.tsv",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
  rawAnno$simple_barcode <- substr(rawAnno$aliquot_barcode,1,15)
  samAnno <- rawAnno[!duplicated(rawAnno$simple_barcode),c("cancer type", "simple_barcode")]
  samAnno <- samAnno[which(samAnno$`cancer type` != ""),]
  write.table(samAnno,"output_simple_sample_annotation.txt",sep = "\t",row.names = F,col.names = T,quote = F)
  # 快速读取表达谱数据并做数据预处理
  expr <- fread("EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",sep = "\t",stringsAsFactors = F,check.names = F,header = T)
  expr <- as.data.frame(expr); rownames(expr) <- expr[,1]; expr <- expr[,-1]
  gene <- sapply(strsplit(rownames(expr),"|",fixed = T), "[",1)
  expr$gene <- gene
  expr <- expr[!duplicated(expr$gene),]
  rownames(expr) <- expr$gene; expr <- expr[,-ncol(expr)]
  expr[expr < 0] <- 0 # 对于这份泛癌数据，将略小于0的数值拉到0，否则不能取log（其他途径下载的泛癌数据可能不需要此操作）
  colnames(expr) <- substr(colnames(expr),1,15)
  gc()
  # 去掉对于风险基因存在NA值的样本
  expr.sub <- expr[risk.coeff$Gene, ] # 提取仅有风险相关基因的表达谱子集
  expr.sub <- as.data.frame(t(na.omit(t(expr.sub)))) # 对列做去空值，而非对行做
  keepSam <- colnames(expr.sub) # 提取被保留的样本
  expr <- expr[,keepSam] # 重构表达谱
  
  # 读取生存数据(虽然在本代码中没有用到，但是原文使用的样本是具有生存数据的)
  surv <- read.delim("Survival_SupplementalTable_S1_20171025_xena_sp", sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) 
  # 确定肿瘤样本以及对应肿瘤类型
  sam <- samAnno[which(samAnno$`cancer type` != "LAML"),"simple_barcode"] # 去掉白血病样本
  comsam <- intersect(intersect(colnames(expr), sam), rownames(surv)) # 得到与表达谱以及生存的共有样本
  tumsam <- comsam[substr(comsam,14,14) == "0"] # 仅提取肿瘤样本
  tumAnno <- samAnno[which(samAnno$simple_barcode %in% tumsam),] # 获取这些肿瘤样本的注释信息
  tumAnno <- tumAnno[order(tumAnno$`cancer type`),] # 根据肿瘤类型排序
  tumors <- unique(tumAnno$`cancer type`) # 得到32个肿瘤
  
  # 在所有样本中计算TRGs得分(在本代码中仅仅是为了确定根据cox-based TRGs score确定肿瘤的level)
  TRGs.score <- list() # 初始化列表
  TRGs.mean <- c() # 初始化得分均值向量
  outTab <- NULL
  for (i in tumors) {
    sam <- tumAnno[which(tumAnno$`cancer type` == i),"simple_barcode"] # 提取当前肿瘤类型的肿瘤样本
    expr.sub <- log2(expr[risk.coeff$Gene,sam] + 1) # 提取表达谱子集并对数化
    TRGs <- scale(apply(expr.sub,2,function(x) {x %*% risk.coeff$Coefficient})) # 计算经过z-score的TRGs得分
    TRGs.score[[i]] <- TRGs
    TRGs.mean <- c(TRGs.mean, mean(TRGs))
    outTab <- rbind.data.frame(outTab, # 保存得分的计算结果
                               data.frame(tumor = i, # 肿瘤类型
                                          TRGs = as.numeric(TRGs), # 当前得分
                                          row.names = sam,
                                          stringsAsFactors = F),
                               stringsAsFactors = F)
  }
  sapply(TRGs.score, range) # 不存在空值
  write.table(outTab, file = "output_TRGs score of all tumor sample across 32 tumor types.txt",sep = "\t",row.names = T,col.names = F,quote = F)
  names(TRGs.mean) <- tumors
  TRGs.mean <- sort(TRGs.mean, decreasing = T) # 根据均值对肿瘤进行排序
  tumor.level <- names(TRGs.mean) # 将排序结果作为肿瘤因子的等级
  
  # 在所有样本中通过z-score计算致癌通路以及TRGs得分（注意此时TRGs得分不再是由cox系数计算，而是由zscore算法下的单样本富集得到）
  oncosig <- gmt2list("oncogenic.gmt") # 将原文补充材料S4以及S6的基因制作成gmt文件，并将gmt文件读取为list
  oncosig[[4]][1:13] = risk.coeff$Gene
  oncosig[[4]][14:19] = ""
  oncosig <- sapply(oncosig, function(x) setdiff(x,"")) # 去掉list中的空值
  zscore.list <- list()
  outSig <- NULL
  for (i in tumors) {
    message(i)
    sam <- tumAnno[which(tumAnno$`cancer type` == i),"simple_barcode"] # 提取当前肿瘤类型的肿瘤样本
    expr.sub <- log2(expr[,sam] + 1) # 提取表达谱子集并对数化
    zscore.list[[i]] <- quiet(gsva(as.matrix(expr.sub), gset.idx.list = oncosig, method = "zscore")) # 方法选择zscore
    outSig <- rbind.data.frame(outSig, # 保存得分的计算结果
                               cbind.data.frame(tumor = i,
                                                as.data.frame(t(zscore.list[[i]]))),
                               stringsAsFactors = F)
  }
  write.table(outSig, file = "output_oncogenic and TRGs score of all tumor sample across 32 tumor types.txt",sep = "\t",row.names = T,col.names = F,quote = F)
  
  # 设置颜色
  mycol <- c("#A6CEE3",
             "#1F78B4",
             "#B2DF8A",
             "#33A02C",
             "#FB9A99",
             "#E31A1C",
             "#FDBF6F",
             "#FF7F00",
             "#CAB2D6",
             "#6A3D9A",
             "#B15928",
             "#8DD3C7",
             "#BEBADA",
             "#FB8072",
             "#80B1D3",
             "#FDB462",
             "#B3DE69",
             "#FCCDE5",
             "#D9D9D9",
             "#BC80BD",
             "#CCEBC5",
             "#FFED6F",
             "#8C510A",
             "#BF812D",
             "#DFC27D",
             "#F6E8C3",
             "#80CDC1",
             "#35978F",
             "#01665E",
             "#003C30",
             "#8E0152",
             "#C51B7D")
  
  # 制作绘图数据并绘图
  plotdata <- outSig
  plotdata <- gather(plotdata, oncogenic, zscore, Angiogenesis:`Cell cycle`, factor_key=TRUE)
  plotdata$tumor <- factor(plotdata$tumor, levels = tumor.level)
  
  p1 <- ggplot(data = plotdata, aes(x = zscore, y = TRGs)) + 
    geom_point(aes(color=tumor),size=1.5,alpha = 0.5) +
    scale_color_manual(values = mycol) + 
    geom_smooth(method = "lm", se = FALSE) +
    geom_ribbon(stat = "smooth",method = "lm",se = TRUE,alpha = 0,linetype = "dashed") + 
    xlab("Oncogenic (z-score)") + ylab("TRGs (z-score)") + 
    stat_cor(method = "pearson", label.x = -40, label.y = 10) + 
    facet_wrap(.~oncogenic, nrow = 1) + 
    theme_bw() + 
    theme(axis.text.x = element_text(vjust = 0.5, size = 12, color = "black"),
          axis.ticks = element_line(size=0.2, color="black"),
          axis.ticks.length = unit(0.2, "cm"),
          legend.position = "right",
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12))
  ## Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
  ## ℹ Please use the `linewidth` argument instead.
  ggsave(filename = "correlation scatter plot of zscored oncogenic and TRGs in pancancer.pdf", width = 15,height = 5)
  ## `geom_smooth()` using formula = 'y ~ x'
  ## `geom_smooth()` using formula = 'y ~ x'
  
  
  tmp1 <- plotdata[which(plotdata$oncogenic == "Angiogenesis"),]
  p2 <- ggplot(data = tmp1, aes(x = zscore, y = TRGs)) + 
    geom_point(aes(color=tumor),size=1.2,alpha = 0.5) +
    scale_color_manual(values = mycol) + 
    geom_smooth(method = "lm", se = FALSE) +
    geom_ribbon(stat = "smooth",method = "lm",se = TRUE,alpha = 0,linetype = "dashed") + 
    stat_cor(method = "pearson", label.x = -40, label.y = 10) + 
    xlab("Angiogenesis (z-score)") + ylab("TRGs (z-score)") + 
    facet_wrap(.~tumor, ncol = 8) + 
    theme_bw() + 
    theme(axis.text.x = element_text(vjust = 0.5, size = 12, color = "black"),
          axis.ticks = element_line(size=0.2, color="black"),
          axis.ticks.length = unit(0.2, "cm"),
          legend.position = "none",
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12))
  ggsave(filename = "correlation scatter plot of zscored Angiogenesis and TRGs in pancancer.pdf", width = 15,height = 8)
  ## `geom_smooth()` using formula = 'y ~ x'
  ## `geom_smooth()` using formula = 'y ~ x'
  
  
  tmp2 <- plotdata[which(plotdata$oncogenic == "EMT"),]
  p3 <- ggplot(data = tmp2, aes(x = zscore, y = TRGs)) + 
    geom_point(aes(color=tumor),size=1.2,alpha = 0.5) +
    scale_color_manual(values = mycol) + 
    geom_smooth(method = "lm", se = FALSE) +
    geom_ribbon(stat = "smooth",method = "lm",se = TRUE,alpha = 0,linetype = "dashed") + 
    stat_cor(method = "pearson", label.x = -40, label.y = 10) + 
    xlab("EMT (z-score)") + ylab("TRGs (z-score)") + 
    facet_wrap(.~tumor, ncol = 8) + 
    theme_bw() + 
    theme(axis.text.x = element_text(vjust = 0.5, size = 12, color = "black"),
          axis.ticks = element_line(size=0.2, color="black"),
          axis.ticks.length = unit(0.2, "cm"),
          legend.position = "none",
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12))
  ggsave(filename = "correlation scatter plot of zscored EMT and TRGs in pancancer.pdf", width = 15,height = 8)
  ## `geom_smooth()` using formula = 'y ~ x'
  ## `geom_smooth()` using formula = 'y ~ x'
  p <- plot_grid(p1,p2,p3, align = "v", ncol = 1)
  ## `geom_smooth()` using formula = 'y ~ x'
  ## `geom_smooth()` using formula = 'y ~ x'
  ## `geom_smooth()` using formula = 'y ~ x'
  ## `geom_smooth()` using formula = 'y ~ x'
  ## `geom_smooth()` using formula = 'y ~ x'
  ## `geom_smooth()` using formula = 'y ~ x'
  ## Warning: Graphs cannot be vertically aligned unless the axis parameter is set.
  ## Placing graphs unaligned.
  ggsave(filename = "combined correlation scatter plot of zscored oncogenic and TRGs in pancancer.pdf", width = 18,height = 23)
  
}

