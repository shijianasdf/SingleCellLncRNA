#' @description  分型的评估和生物学功能刻画
#' @author shi jian

## GSE108474识别亚型特异基因
{ 
  load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/ConsensusClusterResult.rda")
  load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/subtype.probe.annotation.rda")
  GSE108474_GPL570_exprSet <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE108474/GSE108474-GPL570-expression.MATRIX.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
  GSE108474_GPL570_pData <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE108474/GSE108474-GPL570-pData.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
  GSE108474_GPL570_pData1 <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE108474/GSE108474_REMBRANDT_clinical.data.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
  GSE108474_GPL570_probeAnnotation <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE108474/GSE108474-GPL570-ProbeAnnotation.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
  identical( colnames(GSE108474_GPL570_exprSet), GSE108474_GPL570_pData$geo_accession)
  
  GSE108474_GPL570_pData$characteristics_ch1.1[grep("glioblastoma multiforme",GSE108474_GPL570_pData$characteristics_ch1.1)] <- "gbm"
  GSE108474_GPL570_pData$characteristics_ch1.1[grep("normal",GSE108474_GPL570_pData$characteristics_ch1.1)] <- "normal" #四个正常样本
  GSE108474_GPL570_pData <- GSE108474_GPL570_pData[which(GSE108474_GPL570_pData$characteristics_ch1.1 == "gbm" | GSE108474_GPL570_pData$characteristics_ch1.1 == "normal"),]
  GSE108474_GPL570_pData <- GSE108474_GPL570_pData[GSE108474_GPL570_pData$geo_accession %in% setdiff(colnames(GSE108474_GPL570_exprSet),c("GSM2899720","GSM2899721","GSM2899722","GSM2899723","GSM2899724","GSM2899725","GSM2899726")),]
  GSE108474.anno.res$exprSet <- GSE108474.anno.res$exprSet[ ,GSE108474_GPL570_pData$characteristics_ch1.1 == "gbm"]
  identical( colnames(GSE108474.anno.res$exprSet), names(results[[5]]$consensusClass) )
  
  design1 <- cbind.data.frame( sample = names(results[[5]]$consensusClass), label = ifelse( results[[5]]$consensusClass == 1,"case" ,"control" ) )
  design2 <- cbind.data.frame( sample = names(results[[5]]$consensusClass), label = ifelse( results[[5]]$consensusClass == 2,"case" ,"control" ) )
  design3 <- cbind.data.frame( sample = names(results[[5]]$consensusClass), label = ifelse( results[[5]]$consensusClass == 3,"case" ,"control" ) )
  design4 <- cbind.data.frame( sample = names(results[[5]]$consensusClass), label = ifelse( results[[5]]$consensusClass == 4,"case" ,"control" ) )
  design5 <- cbind.data.frame( sample = names(results[[5]]$consensusClass), label = ifelse( results[[5]]$consensusClass == 5,"case" ,"control" ) )
 
  marker1 <- limma1(counts = GSE108474.anno.res$exprSet,
         design = design1,
         contrast.col = "label",
         method = c("voom","trend")[2],
         contrast.level = c("case","control"), # case相对control计算差异
         cutoff.lFC = 1,
         cutoff.padj = 0.05,
         save.file = T,
         report  = T,
         names = "marker1")
  
  marker2 <- limma1(counts = GSE108474.anno.res$exprSet,
                    design = design2,
                    contrast.col = "label",
                    method = c("voom","trend")[2],
                    contrast.level = c("case","control"), # case相对control计算差异,NO是control组
                    cutoff.lFC = 1,
                    cutoff.padj = 0.05,
                    save.file = T,
                    report  = T,
                    names = "marker2")
  
  marker3 <- limma1(counts = GSE108474.anno.res$exprSet,
                    design = design3,
                    contrast.col = "label",
                    method = c("voom","trend")[2],
                    contrast.level = c("case","control"), # case相对control计算差异,NO是control组
                    cutoff.lFC = 1,
                    cutoff.padj = 0.05,
                    save.file = T,
                    report  = T,
                    names = "marker3")
  
  marker4 <- limma1(counts = GSE108474.anno.res$exprSet,
                    design = design4,
                    contrast.col = "label",
                    method = c("voom","trend")[2],
                    contrast.level = c("case","control"), # case相对control计算差异,NO是control组
                    cutoff.lFC = 1,
                    cutoff.padj = 0.05,
                    save.file = T,
                    report  = T,
                    names = "marker4")
  
  marker5 <- limma1(counts = GSE108474.anno.res$exprSet,
                    design = design5,
                    contrast.col = "label",
                    method = c("voom","trend")[2],
                    contrast.level = c("case","control"), # case相对control计算差异,NO是control组
                    cutoff.lFC = 1,
                    cutoff.padj = 0.05,
                    save.file = T,
                    report  = T,
                    names = "marker5")
  
  save(marker1,
       marker2,
       marker3,
       marker4,
       marker5,file="D:/Rsources/Project/SingleCelllncRNA/Results/6.evaluateSubtype/subtype.markers.rda")
  
  
}

## TCGA官方亚型差异基因
{
  load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA.ConsensusClusterResult.rda")
  gbm.clinic <- get(load("D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/gbm.clinic.rda"))
  exprSet <- get(load("D:/Rsources/Project/SingleCelllncRNA/Data/TCGA-GBM/Transcriptome Profiling.Gene Expression Quantification.HTSeq - FPKM.rda"))
  exprSet <- assay(exprSet) 
  table(substring(colnames(exprSet),14,15) < 10) # 5个正常样本，169个疾病样本
  exprSet <- exprSet[,substring(colnames(exprSet),14,15) < 10]
  colnames(exprSet) <- substring( colnames(exprSet), 1, 12 )
  #提取表达谱样本中病人id和临床病人id overlap部分
  common_patients <- intersect( colnames(exprSet), gbm.clinic$Patient.ID ) 
  exprSet <- exprSet[,colnames(exprSet) %in% common_patients]
  gbm.clinic <- gbm.clinic[gbm.clinic$Patient.ID %in% common_patients,]
  table(as.character(gbm.clinic$Original.Subtype),useNA="ifany")
  
  ## 对表达谱病人以及临床病人去重复
  gbm.clinic <- gbm.clinic[!duplicated(gbm.clinic$Patient.ID),]
  exprSet <- exprSet[,!duplicated(colnames(exprSet))]
  gbm.clinic$label1 <- ifelse(gbm.clinic$Original.Subtype != "Classical" | is.na(gbm.clinic$Original.Subtype),"control","case")
  gbm.clinic$label2 <- ifelse(gbm.clinic$Original.Subtype == "G-CIMP" | is.na(gbm.clinic$Original.Subtype),"control","case")
  gbm.clinic$label3 <- ifelse(gbm.clinic$Original.Subtype == "Mesenchymal" | is.na(gbm.clinic$Original.Subtype),"control","case")
  gbm.clinic$label4 <- ifelse(gbm.clinic$Original.Subtype == "Neural" | is.na(gbm.clinic$Original.Subtype),"control","case")
  gbm.clinic$label5 <- ifelse(gbm.clinic$Original.Subtype == "Proneural" | is.na(gbm.clinic$Original.Subtype),"control","case")
  rownames(gbm.clinic) <- gbm.clinic$Patient.ID
   
  identical( names(TCGA.results[[5]]$consensusClass), gbm.clinic$Patient.ID )
  identical( colnames(exprSet), gbm.clinic$Patient.ID )
  TCGA.results[[5]]$consensusClass <- TCGA.results[[5]]$consensusClass[ match( gbm.clinic$Patient.ID, names(TCGA.results[[5]]$consensusClass) ) ]
  gbm.clinic$label6 <- TCGA.results[[5]]$consensusClass 
  gbm.clinic$label7 <- ifelse(gbm.clinic$label6 == 1,"case","control")
  gbm.clinic$label8 <- ifelse(gbm.clinic$label6 == 2,"case","control")
  gbm.clinic$label9 <- ifelse(gbm.clinic$label6 == 3,"case","control")
  gbm.clinic$label10 <- ifelse(gbm.clinic$label6 == 4,"case","control")
  gbm.clinic$label11 <- ifelse(gbm.clinic$label6 == 5,"case","control")
  
  TCGA.LC1.marker <- ScreenGenes(exprSet,
                                       gbm.clinic,
                                       contrast.col = "label7",
                                       contrast.level = c("case","control"),
                                       contrast.control = c("control"),
                                       method = c("t.test","wilcox.test")[1],
                                       cutoff.p=0.05,
                                       cutoff.q=0.05,
                                       cutoff.logFC=1,
                                       pmethod = "BH",
                                       verbose = F)
  
  TCGA.LC2.marker <- ScreenGenes(exprSet,
                                 gbm.clinic,
                                 contrast.col = "label8",
                                 contrast.level = c("case","control"),
                                 contrast.control = c("control"),
                                 method = c("t.test","wilcox.test")[1],
                                 cutoff.p=0.05,
                                 cutoff.q=0.05,
                                 cutoff.logFC=1,
                                 pmethod = "BH",
                                 verbose = F)
  TCGA.LC3.marker <- ScreenGenes(exprSet,
                                 gbm.clinic,
                                 contrast.col = "label9",
                                 contrast.level = c("case","control"),
                                 contrast.control = c("control"),
                                 method = c("t.test","wilcox.test")[1],
                                 cutoff.p=0.05,
                                 cutoff.q=0.05,
                                 cutoff.logFC=1,
                                 pmethod = "BH",
                                 verbose = F)
  TCGA.LC4.marker <- ScreenGenes(exprSet,
                                 gbm.clinic,
                                 contrast.col = "label10",
                                 contrast.level = c("case","control"),
                                 contrast.control = c("control"),
                                 method = c("t.test","wilcox.test")[1],
                                 cutoff.p=0.05,
                                 cutoff.q=0.05,
                                 cutoff.logFC=1,
                                 pmethod = "BH",
                                 verbose = F)
  TCGA.LC5.marker <- ScreenGenes(exprSet,
                                 gbm.clinic,
                                 contrast.col = "label11",
                                 contrast.level = c("case","control"),
                                 contrast.control = c("control"),
                                 method = c("t.test","wilcox.test")[1],
                                 cutoff.p=0.05,
                                 cutoff.q=0.05,
                                 cutoff.logFC=1,
                                 pmethod = "BH",
                                 verbose = F)
  save(TCGA.LC1.marker,
       TCGA.LC2.marker,
       TCGA.LC3.marker,
       TCGA.LC4.marker,
       TCGA.LC5.marker,file="D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA.LC.subtype.marker.rda")
  
  TCGA.Classical.marker <- ScreenGenes(exprSet,
                                        gbm.clinic,
                                        contrast.col = "label1",
                                        contrast.level = c("case","control"),
                                        contrast.control = c("control"),
                                        method = c("t.test","wilcox.test")[1],
                                        cutoff.p=0.05,
                                        cutoff.q=0.05,
                                        cutoff.logFC=1,
                                        pmethod = "BH",
                                        verbose = F)
  
  TCGA.G_CIMP.marker <- ScreenGenes(exprSet,
                                       gbm.clinic,
                                       contrast.col = "label2",
                                       contrast.level = c("case","control"),
                                       contrast.control = c("control"),
                                       method = c("t.test","wilcox.test")[1],
                                       cutoff.p=0.05,
                                       cutoff.q=0.05,
                                       cutoff.logFC=1,
                                       pmethod = "BH",
                                       verbose = F)
  
  TCGA.Mesenchymal.marker <- ScreenGenes(exprSet,
                                    gbm.clinic,
                                    contrast.col = "label3",
                                    contrast.level = c("case","control"),
                                    contrast.control = c("control"),
                                    method = c("t.test","wilcox.test")[1],
                                    cutoff.p=0.05,
                                    cutoff.q=0.05,
                                    cutoff.logFC=1,
                                    pmethod = "BH",
                                    verbose = F)
  
  TCGA.Neural.marker <- ScreenGenes(exprSet,
                                         gbm.clinic,
                                         contrast.col = "label4",
                                         contrast.level = c("case","control"),
                                         contrast.control = c("control"),
                                         method = c("t.test","wilcox.test")[1],
                                         cutoff.p=0.05,
                                         cutoff.q=0.05,
                                         cutoff.logFC=1,
                                         pmethod = "BH",
                                         verbose = F)
  
  TCGA.Proneural.marker <- ScreenGenes(exprSet,
                                    gbm.clinic,
                                    contrast.col = "label5",
                                    contrast.level = c("case","control"),
                                    contrast.control = c("control"),
                                    method = c("t.test","wilcox.test")[1],
                                    cutoff.p=0.05,
                                    cutoff.q=0.05,
                                    cutoff.logFC=1,
                                    pmethod = "BH",
                                    verbose = F)
  save( TCGA.Classical.marker,
        TCGA.G_CIMP.marker,
        TCGA.Mesenchymal.marker,
        TCGA.Neural.marker,
        TCGA.Proneural.marker,file="D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA.subtype.marker.rda" )

  
  
  
 ## 从文献中挑出来各个亚型特异性表达基因并且还要在各个亚型中相对于其他亚型差异表达，然后绘制热图
 load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA.subtype.marker.rda")
 load("D:/Rsources/Project/SingleCelllncRNA/Results/6.evaluateSubtype/subtype.markers.rda")
 load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA.subtype.marker.rda")
 load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA.LC.subtype.marker.rda")
 load("D:/Rsources/Project/SingleCelllncRNA/Results/6.evaluateSubtype/subtype.markers.rda")
 
 t1 <- convert( intersect(rownames(marker1$diffSig[marker1$diffSig$lab == "UP",]), rownames(na.omit(TCGA.LC1.marker[TCGA.LC1.marker$lab == "UP",]))) ) 
 t2 <- convert( intersect(rownames(marker2$diffSig[marker2$diffSig$lab == "UP",]), rownames(na.omit(TCGA.LC2.marker[TCGA.LC2.marker$lab == "UP",]))) ) 
 t3 <- convert( intersect(rownames(marker3$diffSig[marker3$diffSig$lab == "DOWN",]), rownames(na.omit(TCGA.LC3.marker[TCGA.LC3.marker$lab == "DOWN",]))) )
 t4 <- convert( intersect(rownames(marker4$diffSig[marker4$diffSig$lab == "DOWN",]), rownames(na.omit(TCGA.LC4.marker[TCGA.LC4.marker$lab == "DOWN",]))) )
 t5 <- convert( intersect(rownames(marker5$diffSig[marker5$diffSig$lab == "DOWN",]), rownames(na.omit(TCGA.LC5.marker[TCGA.LC5.marker$lab == "DOWN",]))) )

 marker1$diffSig$lab1 <- ifelse(abs(marker1$diffSig$logFC) > 0.5 & marker1$diffSig$adj.P.Val < 0.05,ifelse(marker1$diffSig$logFC > 0.5, "UP","DOWN"),"NO")
 TCGA.LC1.marker$lab1 <- ifelse(abs(TCGA.LC1.marker$logFC) > 0.5 & TCGA.LC1.marker$qvalue < 0.05,ifelse(TCGA.LC1.marker$logFC > 0.5, "UP","DOWN"),"NO")
 marker2$diffSig$lab1 <- ifelse(abs(marker2$diffSig$logFC) > 0 & marker2$diffSig$adj.P.Val < 0.05,ifelse(marker2$diffSig$logFC > 0, "UP","DOWN"),"NO")
 TCGA.LC2.marker$lab1 <- ifelse(abs(TCGA.LC2.marker$logFC) > 0 & TCGA.LC2.marker$qvalue < 0.05,ifelse(TCGA.LC2.marker$logFC > 0, "UP","DOWN"),"NO")
 t1 <- convert( intersect(rownames(marker1$diffSig[marker1$diffSig$lab1 == "UP",]), rownames(na.omit(TCGA.LC1.marker[TCGA.LC1.marker$lab1 == "UP",]))) ) 
 t2 <- convert( intersect(rownames(marker2$diffSig[marker2$diffSig$lab1 == "DOWN",]), rownames(na.omit(TCGA.LC2.marker[TCGA.LC2.marker$lab1 == "DOWN",]))) ) 
 t1 <- c("REV3L","ZNF587","SOX9","KMT2C","SACS")
 t2 <- c("MUC1","TMCO2","HNRNPC","RBPMS","CDH4")
 t3 <- c("STH","PXDN","HIST3H3","STPG3","FLNC")
 t4 <- c("CMPK2","SIDT1","CATSPERZ","LINC02683")
 t5 <- c("MAGEB6","LRRTM1","PAK3","PTPRT","ADGRB1")
 setdiff( t1, c(t2,t3,t4,t5) )
 setdiff( t2, c(t1,t3,t4,t5) )
 setdiff( t3, c(t1,t2,t4,t5) )
 setdiff( t4, c(t1,t2,t3,t5) ) 
 setdiff( t5, c(t1,t2,t3,t4) )
 
 ## 对TCGA分型结果进行热图可视化
 library(ComplexHeatmap)
 #identical( colnames(exprSet), gbm.clinic$Patient.ID)
 exprSet <- exprSet[ , names(sort(TCGA.results[[5]]$consensusClass)) ]
 gbm.clinic <- gbm.clinic[match( colnames(exprSet), gbm.clinic$Patient.ID ),]
 identical( colnames(exprSet), names(sort(TCGA.results[[5]]$consensusClass)) ) 
 identical( colnames(exprSet), gbm.clinic$Patient.ID )
 Classical.genes <- setdiff(t1,c(t2,t3,t4,t5)) 
 proneural.genes <- setdiff(t2,c(t1,t3,t4,t5))
 mesenchymal.genes <- setdiff(t3,c(t1,t2,t4,t5))
 neural.genes <- setdiff(t4,c(t1,t2,t3,t5))
 G_CIMP.genes <- setdiff(t5,c(t1,t2,t3,t4))
 mat1 <- exprSet[convert(Classical.genes,fromtype = "SYMBOL",totype = "ENSEMBL"),]
 mat2 <- exprSet[convert(proneural.genes,fromtype = "SYMBOL",totype = "ENSEMBL"),]
 mat3 <- exprSet[convert(mesenchymal.genes,fromtype = "SYMBOL",totype = "ENSEMBL"),]
 mat4 <- exprSet[convert(neural.genes,fromtype = "SYMBOL",totype = "ENSEMBL"),]
 mat5 <- exprSet[convert(G_CIMP.genes,fromtype = "SYMBOL",totype = "ENSEMBL"),]
 mat1 <- t( scale( t( mat1 ) ) )
 mat2 <- t( scale( t( mat2 ) ) )
 mat3 <- t( scale( t( mat3 ) ) )
 mat4 <- t( scale( t( mat4 ) ) )
 mat5 <- t( scale( t( mat5 ) ) )
 rownames(mat1) <- convert(rownames(mat1))
 rownames(mat2) <- convert(rownames(mat2))
 rownames(mat3) <- convert(rownames(mat3))
 rownames(mat4) <- convert(rownames(mat4))
 rownames(mat5) <- convert(rownames(mat5))
 library(circlize)
 col_fun <- colorRamp2(c(-2, 0, 2),colors = c("blue", "#EEEEEE", "red"),space = "RGB")
 lgd1 <- Legend(col_fun = col_fun, title = "Z-score")
 lgd2 <- Legend(at=c("LncRNA LC1","LncRNA LC2","LncRNA LC3","LncRNA LC4","LncRNA LC5"),legend_gp = gpar(fill = c("#8B181B","#DC9D27","#186633","#8252A1","#869D79")),title = "LncRNA Consensus Cluster")
 lgd3 <- Legend(at=c("Classical","Mesenchymal","Proneural","Neural","G-CIMP"),legend_gp = gpar(fill = c("#2CAF51","#25AAE1","#F19F32","#2D317C","#8252A1")),title = "TCGA Subtype")
 pd <-  packLegend(lgd2, lgd3, lgd1)
 ha <-  HeatmapAnnotation( df = data.frame("LncRNA Consensus Cluster" = as.character(sort(TCGA.results[[5]]$consensusClass)),
                                                "TCGA Subtype" = as.character(gbm.clinic$Original.Subtype)),
                           col = list("LncRNA Consensus Cluster" = c("1" = "#8B181B","2"="#DC9D27","3"="#186633","4"="#8252A1","5"="#869D79"),
                                     "TCGA Subtype" = c("Classical" = "#2CAF51","Mesenchymal" = "#25AAE1","Proneural" = "#F19F32","Neural" = "#2D317C","G-CIMP" = "#8252A1")
                          ),
                          na_col = "gray", show_legend = FALSE, annotation_name_gp = gpar(fontsize = 6) )
 ht1 <- Heatmap(mat1, row_names_gp = gpar(fontsize = 6), show_column_names = FALSE, row_names_max_width = unit(3, "cm"), col = col_fun, show_heatmap_legend = FALSE,  cluster_rows=F, cluster_columns=F, row_title_gp = gpar(fill = "#2CAF51", col = "black", border = "black"), row_title = "Classical" )
 ht2 <- Heatmap(mat2, row_names_gp = gpar(fontsize = 6), show_column_names = FALSE, row_names_max_width = unit(3, "cm"), col = col_fun, show_heatmap_legend = FALSE,  cluster_rows=F, cluster_columns=F, row_title_gp = gpar(fill = "#F19F32", col = "black", border = "black"), row_title = "Proneural" )
 ht3 <- Heatmap(mat3, row_names_gp = gpar(fontsize = 6), show_column_names = FALSE, row_names_max_width = unit(3, "cm"), col = col_fun, show_heatmap_legend = FALSE,  cluster_rows=F, cluster_columns=F, row_title_gp = gpar(fill = "#25AAE1", col = "black", border = "black"), row_title = "Mesenchymal" )
 ht4 <- Heatmap(mat4, row_names_gp = gpar(fontsize = 6), show_column_names = FALSE, row_names_max_width = unit(3, "cm"), col = col_fun, show_heatmap_legend = FALSE,  cluster_rows=F, cluster_columns=F, row_title_gp = gpar(fill = "#2D317C", col = "black", border = "black"), row_title = "Neural" )
 ht5 <- Heatmap(mat5, row_names_gp = gpar(fontsize = 6), show_column_names = FALSE, row_names_max_width = unit(3, "cm"), col = col_fun, show_heatmap_legend = FALSE,  cluster_rows=F, cluster_columns=F, row_title_gp = gpar(fill = "#8252A1", col = "black", border = "black"), row_title = "G-CIMP" )
 ht_list <- ha %v% ht1 %v% ht2 %v% ht3 %v% ht4 %v% ht5 
 pdf("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/tcga.subtype.heatmap.pdf");draw(ht_list);dev.off()
 pdf("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/tcga.subtype.heatmap1.pdf");draw(pd);dev.off()
 
 
 ## GSE108474表达谱分型结果可视化
 identical( colnames(GSE108474.anno.res$exprSet),names(sort(results[[5]]$consensusClass)) )
 GSE108474.anno.res$exprSet <- GSE108474.anno.res$exprSet[ ,match( names(sort(results[[5]]$consensusClass)), colnames(GSE108474.anno.res$exprSet)  ) ]
 GSE108474.mat1 <- GSE108474.anno.res$exprSet[convert(Classical.genes,fromtype = "SYMBOL",totype = "ENSEMBL"),]
 GSE108474.mat2 <- GSE108474.anno.res$exprSet[convert(proneural.genes,fromtype = "SYMBOL",totype = "ENSEMBL"),]
 GSE108474.mat3 <- GSE108474.anno.res$exprSet[convert(mesenchymal.genes,fromtype = "SYMBOL",totype = "ENSEMBL"),]
 GSE108474.mat4 <- GSE108474.anno.res$exprSet[convert(neural.genes,fromtype = "SYMBOL",totype = "ENSEMBL"),]
 GSE108474.mat5 <- GSE108474.anno.res$exprSet[convert(G_CIMP.genes,fromtype = "SYMBOL",totype = "ENSEMBL"),]
 GSE108474.mat1 <- t( scale( t( GSE108474.mat1 ) ) )
 GSE108474.mat2 <- t( scale( t( GSE108474.mat2 ) ) )
 GSE108474.mat3 <- t( scale( t( GSE108474.mat3 ) ) )
 GSE108474.mat4 <- t( scale( t( GSE108474.mat4 ) ) )
 GSE108474.mat5 <- t( scale( t( GSE108474.mat5 ) ) )
 rownames(GSE108474.mat1) <- convert(rownames(GSE108474.mat1))
 rownames(GSE108474.mat2) <- convert(rownames(GSE108474.mat2))
 rownames(GSE108474.mat3) <- convert(rownames(GSE108474.mat3))
 rownames(GSE108474.mat4) <- convert(rownames(GSE108474.mat4))
 rownames(GSE108474.mat5) <- convert(rownames(GSE108474.mat5))
 library(circlize)
 col_fun <- colorRamp2(c(-2, 0, 2),colors = c("blue", "#EEEEEE", "red"),space = "RGB")
 lgd1 <- Legend(col_fun = col_fun, title = "Z-score")
 lgd2 <- Legend(at=c("LncRNA LC1","LncRNA LC2","LncRNA LC3","LncRNA LC4","LncRNA LC5"),legend_gp = gpar(fill = c("#8B181B","#DC9D27","#186633","#8252A1","#869D79")),title = "LncRNA Consensus Cluster")
 pd <-  packLegend(lgd2, lgd1)
 GSE108474.ha <- HeatmapAnnotation( df = data.frame("LncRNA Consensus Cluster" = as.character(sort(results[[5]]$consensusClass))),
                                    col = list("LncRNA Consensus Cluster" = c("1" = "#8B181B","2"="#DC9D27","3"="#186633","4"="#8252A1","5"="#869D79")),
                                    na_col = "gray",show_legend = FALSE, annotation_name_gp = gpar(fontsize = 6) )
 GSE108474.ht1 <- Heatmap(GSE108474.mat1, col = col_fun, row_names_gp = gpar(fontsize = 6), show_column_names = FALSE, show_heatmap_legend = FALSE, row_names_max_width = unit(3, "cm"), cluster_rows=F, cluster_columns=F, row_title_gp = gpar(fill = "#2CAF51", col = "black", border = "black"), row_title = "Classical" )
 GSE108474.ht2 <- Heatmap(GSE108474.mat2, col = col_fun, row_names_gp = gpar(fontsize = 6), show_column_names = FALSE,show_heatmap_legend = FALSE, row_names_max_width = unit(3, "cm"), cluster_rows=F, cluster_columns=F, row_title_gp = gpar(fill = "#F19F32", col = "black", border = "black"), row_title = "Proneural" )
 GSE108474.ht3 <- Heatmap(GSE108474.mat3, col = col_fun, row_names_gp = gpar(fontsize = 6), show_column_names = FALSE,show_heatmap_legend = FALSE, row_names_max_width = unit(3, "cm"), cluster_rows=F, cluster_columns=F, row_title_gp = gpar(fill = "#25AAE1", col = "black", border = "black"), row_title = "Mesenchymal" )
 GSE108474.ht4 <- Heatmap(GSE108474.mat4, col = col_fun, row_names_gp = gpar(fontsize = 6), show_column_names = FALSE,show_heatmap_legend = FALSE, row_names_max_width = unit(3, "cm"), cluster_rows=F, cluster_columns=F, row_title_gp = gpar(fill = "#2D317C", col = "black", border = "black"), row_title = "Neural" )
 GSE108474.ht5 <- Heatmap(GSE108474.mat5, col = col_fun, row_names_gp = gpar(fontsize = 6), show_column_names = FALSE,show_heatmap_legend = FALSE, row_names_max_width = unit(3, "cm"), cluster_rows=F, cluster_columns=F, row_title_gp = gpar(fill = "#8252A1", col = "black", border = "black"), row_title = "G-CIMP" )
 GSE108474.ht_list <- GSE108474.ha %v% GSE108474.ht1 %v% GSE108474.ht2 %v% GSE108474.ht3 %v% GSE108474.ht4 %v% GSE108474.ht5
 pdf("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE108474.subtype.heatmap.pdf");draw(GSE108474.ht_list);dev.off()
 pdf("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE108474.subtype.heatmap1.pdf");draw(pd);dev.off()
 
}

## RTN转录因子调控网络分析
{
  load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/ConsensusClusterResult.rda")
  load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA.ConsensusClusterResult.rda")
  GSE108474.anno.res$exprSet
  GSE108474.anno.res$ID2gene <- GSE108474.anno.res$ID2gene[,c(4,1,2,3,5)]
  GSE108474.anno.res$ID2gene <- GSE108474.anno.res$ID2gene[!duplicated(GSE108474.anno.res$ID2gene$ensemble),]
  rownames(GSE108474.anno.res$ID2gene) <- GSE108474.anno.res$ID2gene$ensemble
  GSE108474_GPL570_pData <- GSE108474_GPL570_pData[,c(2,1,3:42)]
  rownames(GSE108474_GPL570_pData) <- GSE108474_GPL570_pData$geo_accession
  GSE108474_GPL570_pData <- GSE108474_GPL570_pData[GSE108474_GPL570_pData$characteristics_ch1.1 == "gbm",]
  GSE108474_GPL570_pData <- cbind.data.frame(GSE108474_GPL570_pData,label2=results[[2]]$consensusClass,
                   label3=results[[3]]$consensusClass,label4=results[[4]]$consensusClass,label5=results[[5]]$consensusClass)
  identical(colnames(GSE108474.anno.res$exprSet),names(results[[5]]$consensusClass))
  identical(GSE108474_GPL570_pData$geo_accession,names(results[[5]]$consensusClass))
  identical(colnames(GSE108474.anno.res$exprSet),GSE108474_GPL570_pData$geo_accession )
  rownames(GSE108474.anno.res$ID2gene[!duplicated(GSE108474.anno.res$ID2gene$ensemble),]) <- GSE108474.anno.res$ID2gene[!duplicated(GSE108474.anno.res$ID2gene$ensemble),]$ensemble
  
  ## 基于表达谱构建转录因子调控网络
  library(RTN)
  library(RTNsurvival)
  tfs <- c("FOXM1","E2F2","E2F3","RUNX2","PTTG1","RB1","FGFR3","GATA3","STAT3","ESR2","EGFR","EGBB2","GATA6","HIF1A","ERBB3","ERBB2"
           ,"FOXA1","RXRA") # 转录因子继续挑选
  rtni <- tni.constructor(expData = as.matrix(GSE108474.anno.res$exprSet), 
                          regulatoryElements = tfs, 
                          rowAnnotation = GSE108474.anno.res$ID2gene[!duplicated(GSE108474.anno.res$ID2gene$ensemble),], 
                          colAnnotation = GSE108474_GPL570_pData)
  rtni <- tni.permutation(rtni, nPermutations = 1000)
  rtni <- tni.bootstrap(rtni)
  rtni <- tni.dpi.filter(rtni) 
  regulons <- tni.get(rtni, what = "regulons.and.mode") # 获取target信息
  metabric_annot <- tni.get(rtni, "colAnnotation")  # 获取样本指数信息
  # 计算转录因子调控活性得分
  rtni1st <- tni.gsea2(rtni, regulatoryElements = tfs) 
  # 提取调控活性得分
  metabric_regact <- tni.get(rtni1st, what = "regulonActivity") 
  
  # Replace TCGA samples
  rtni1st_tcgasamples <- tni.replace.samples(rtni1st, exprSet)
  # Compute regulon activity for the new samples
  rtni1st_tcgasamples <- tni.gsea2(rtni1st_tcgasamples, regulatoryElements = tfs)
  tcga_regact <- tni.get(rtni1st_tcgasamples, what = "regulonActivity")
  (tcga_regact$positive + tcga_regact$negative)/2
  
  GSE108474.regulon.activity <- (metabric_regact$positive + metabric_regact$negative)/2 # TF正向靶和负向靶的活性得分均值
  identical( names(sort(results[[5]]$consensusClass)), rownames(GSE108474.regulon.activity) )

  regulon.activity.list <- split(as.data.frame(GSE108474.regulon.activity),unname(sort(results[[5]]$consensusClass)))
  regulon.activity.list <- lapply(regulon.activity.list,function(x){ colMeans(x) })
  heatmap.matrix <- do.call(rbind,regulon.activity.list)
  
  fisher.data <- cbind.data.frame(ifelse( GSE108474.regulon.activity > 0,"high","low" ),label=results[[5]]$consensusClass)
  fisher.data$label1 <- ifelse(fisher.data$label == 1,"cluster1","other")
  fisher.data$label2 <- ifelse(fisher.data$label == 2,"cluster2","other")
  fisher.data$label3 <- ifelse(fisher.data$label == 3,"cluster3","other")
  fisher.data$label4 <- ifelse(fisher.data$label == 4,"cluster4","other")
  fisher.data$label5 <- ifelse(fisher.data$label == 5,"cluster5","other")
  temp <- data.frame(Cluster1=character(0), Cluster2=character(0), Cluster3=character(0), Cluster4=character(0), Cluster5=character(0))
  for(i in colnames(fisher.data)[1:15]){
    t1 <- FastStatisticsTest(fisher.data,
                       variable = i,by="label1",test = "chisq.test",
                       id = NULL, type = "categorical")
    t2 <- FastStatisticsTest(fisher.data,
                       variable = i,by="label2",test = "chisq.test",
                       id = NULL, type = "categorical")
    t3 <- FastStatisticsTest(fisher.data,
                       variable = i,by="label3",test = "chisq.test",
                       id = NULL, type = "categorical")
    t4 <- FastStatisticsTest(fisher.data,
                       variable = i,by="label4",test = "chisq.test",
                       id = NULL, type = "categorical")
    t5 <- FastStatisticsTest(fisher.data,
                       variable = i,by="label5",test = "chisq.test",
                       id = NULL, type = "categorical")
    temp <- rbind.data.frame(temp,cbind.data.frame(Cluster1=t1,Cluster2=t2,Cluster3=t3,Cluster4=t4,Cluster5=t5))
  }

 library(pheatmap)
 pdf("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE108474.regulon.activity.pdf")
 annotation_col <- data.frame("LncRNA Cluster"=paste("Cluster",1:5))
 pheatmap( t(heatmap.matrix),border_color = "gray",cluster_cols = F,
           cluster_rows = T, annotation_col = annotation_col, show_colnames = F, 
           color =colorRampPalette(c("#2061a6", "white", "#bc2c35"))(50),
           display_numbers =  ifelse(temp < 0.5,"*",""),fontsize_number = 20 )
 dev.off()
 
 tcga.regulon.activity <- (tcga_regact$positive + tcga_regact$negative)/2 # TF正向靶和负向靶的活性得分均值
 identical( names(sort(TCGA.results[[5]]$consensusClass)), rownames(tcga.regulon.activity) )
 tcga.regulon.activity <- tcga.regulon.activity[match(names(sort(TCGA.results[[5]]$consensusClass)),rownames(tcga.regulon.activity)),]
 
 tcga.regulon.activity.list <- split(as.data.frame(tcga.regulon.activity),unname(sort(TCGA.results[[5]]$consensusClass)))
 tcga.regulon.activity.list <- lapply(tcga.regulon.activity.list,function(x){ colMeans(x) })
 tcga.heatmap.matrix <- do.call(rbind,tcga.regulon.activity.list)
 
 tcga.fisher.data <- cbind.data.frame(ifelse( tcga.regulon.activity > 0,"high","low" ),label=TCGA.results[[5]]$consensusClass)
 tcga.fisher.data$label1 <- ifelse(tcga.fisher.data$label == 1,"cluster1","other")
 tcga.fisher.data$label2 <- ifelse(tcga.fisher.data$label == 2,"cluster2","other")
 tcga.fisher.data$label3 <- ifelse(tcga.fisher.data$label == 3,"cluster3","other")
 tcga.fisher.data$label4 <- ifelse(tcga.fisher.data$label == 4,"cluster4","other")
 tcga.fisher.data$label5 <- ifelse(tcga.fisher.data$label == 5,"cluster5","other")
 tcga.temp <- data.frame(Cluster1=character(0), Cluster2=character(0), Cluster3=character(0), Cluster4=character(0), Cluster5=character(0))
 for(i in colnames(tcga.fisher.data)[1:15]){
   t1 <- FastStatisticsTest(tcga.fisher.data,
                            variable = i,by="label1",test = "chisq.test",
                            id = NULL, type = "categorical")
   t2 <- FastStatisticsTest(tcga.fisher.data,
                            variable = i,by="label2",test = "chisq.test",
                            id = NULL, type = "categorical")
   t3 <- FastStatisticsTest(tcga.fisher.data,
                            variable = i,by="label3",test = "chisq.test",
                            id = NULL, type = "categorical")
   t4 <- FastStatisticsTest(tcga.fisher.data,
                            variable = i,by="label4",test = "chisq.test",
                            id = NULL, type = "categorical")
   t5 <- FastStatisticsTest(tcga.fisher.data,
                            variable = i,by="label5",test = "chisq.test",
                            id = NULL, type = "categorical")
   tcga.temp <- rbind.data.frame(tcga.temp,cbind.data.frame(Cluster1=t1,Cluster2=t2,Cluster3=t3,Cluster4=t4,Cluster5=t5))
 }
 
 library(pheatmap)
 pdf("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/tcga.regulon.activity.pdf")
 annotation_col <- data.frame("LncRNA Cluster"=paste("Cluster",1:5))
 pheatmap( t(tcga.heatmap.matrix),border_color = "gray",cluster_cols = F,
           cluster_rows = T, annotation_col = annotation_col, show_colnames = F, 
           color =colorRampPalette(c("#2061a6", "white", "#bc2c35"))(50),
           display_numbers =  ifelse(tcga.temp < 0.5,"*",""),fontsize_number = 20 )
 dev.off()
 
   
  
}

## GSE108474以及TCGA亚型差异刻画(搜集TCGA亚型marker基因，各种信号通路以及signature集合，肿瘤纯度信息)
{
  load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/ConsensusClusterResult.rda")
  load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA.ConsensusClusterResult.rda")
  identical( colnames(GSE108474.anno.res$exprSet), names(results[[5]]$consensusClass) )
  identical( colnames(gsva.hallmark.score.matrix), names(results[[5]]$consensusClass) )
  library(GSVA)
  library(GSEABase)
  rownames(GSE108474.anno.res$exprSet) <- GSE108474.anno.res$ID2gene$symbol[match(rownames(GSE108474.anno.res$exprSet),GSE108474.anno.res$ID2gene$ensemble)]
  hallmark.genesets <- getGmt("D:/Rsources/Project/SingleCelllncRNA/Data/signatureGeneSets/h.all.v7.2.symbols.gmt")
  immunologic.signature.genesets <- getGmt("D:/Rsources/Project/SingleCelllncRNA/Data/signatureGeneSets/c7.all.v7.2.symbols.gmt")
  oncogenic.signature.genesets <- getGmt("D:/Rsources/Project/SingleCelllncRNA/Data/signatureGeneSets/c6.all.v7.2.symbols.gmt")
  biocarta.genesets <- getGmt("D:/Rsources/Project/SingleCelllncRNA/Data/signatureGeneSets/c2.cp.biocarta.v7.2.symbols.gmt")
  kegg.genesets <- getGmt("D:/Rsources/Project/SingleCelllncRNA/Data/signatureGeneSets/c2.cp.kegg.v7.2.symbols.gmt")
  reactome.genesets <- getGmt("D:/Rsources/Project/SingleCelllncRNA/Data/signatureGeneSets/c2.cp.reactome.v7.2.symbols.gmt")
  gsva.hallmark.score.matrix <- gsva(expr = as.matrix(GSE108474.anno.res$exprSet),
                                      gset.idx.list = hallmark.genesets,
                                      method="gsva",
                                      kcdf="Gaussian")
  gsva.immunologic.score.matrix <- gsva(expr = as.matrix(GSE108474.anno.res$exprSet),
                                     gset.idx.list = immunologic.signature.genesets,
                                     method="gsva",
                                     kcdf="Gaussian")
  gsva.oncogenic.score.matrix <- gsva(expr = as.matrix(GSE108474.anno.res$exprSet),
                                        gset.idx.list = immunologic.signature.genesets,
                                        method="gsva",
                                        kcdf="Gaussian")
  gsva.biocarta.score.matrix <- gsva(expr = as.matrix(GSE108474.anno.res$exprSet),
                                        gset.idx.list = biocarta.genesets,
                                        method="gsva",
                                        kcdf="Gaussian")
  gsva.kegg.score.matrix <- gsva(expr = as.matrix(GSE108474.anno.res$exprSet),
                                        gset.idx.list = kegg.genesets,
                                        method="gsva",
                                        kcdf="Gaussian")
  gsva.reactome.score.matrix <- gsva(expr = as.matrix(GSE108474.anno.res$exprSet),
                                        gset.idx.list = reactome.genesets,
                                        method="gsva",
                                        kcdf="Gaussian")
  hallmark.box <- cbind.data.frame(t(gsva.hallmark.score.matrix),group=results[[5]]$consensusClass)
  immunologic.box <- cbind.data.frame(t(gsva.immunologic.score.matrix),group=results[[5]]$consensusClass)
  oncogenic.box <- cbind.data.frame(t(gsva.oncogenic.score.matrix),group=results[[5]]$consensusClass)
  biocarta.box <- cbind.data.frame(t(gsva.biocarta.score.matrix),group=results[[5]]$consensusClass)
  kegg.box <- cbind.data.frame(t(gsva.kegg.score.matrix),group=results[[5]]$consensusClass)
  reactome.box <- cbind.data.frame(t(gsva.reactome.score.matrix),group=results[[5]]$consensusClass)
  
  # Visualize: Specify the comparisons you want
  library(ggpubr)
  my_comparisons <- list( c("1", "2"),c("1", "3"),c("1", "4"),c("1","5"),c("2","3"),c("2","4"),c("2","5"),c("3","4"),c("3","5"),c("4","5") )
  pdf("D:/Rsources/Project/SingleCelllncRNA/Results/6.evaluateSubtype/hallmark.diff.pdf")
  for(i in colnames(hallmark.box)[-ncol(hallmark.box)]){
    p <- ggboxplot( hallmark.box,
                     x = "group",
                     #combine = TRUE,
                     y = i,
                     ylab = i,
                     width = 0.9, 
                     color = "group", 
                     palette = "jco", add = "jitter", size = 0.2,
                     font.label = list(size = 1) ) + 
      stat_compare_means( comparisons = my_comparisons, size = 2 ) # Add pairwise comparisons p-value
    print(p)
  }
  dev.off()
  
  pdf("D:/Rsources/Project/SingleCelllncRNA/Results/6.evaluateSubtype/immunologic.diff.pdf")
  for(i in colnames(immunologic.box)[-ncol(immunologic.box)]){
    p <- ggboxplot( immunologic.box,
                    x = "group",
                    #combine = TRUE,
                    y = i,
                    ylab = i,
                    width = 0.9, 
                    color = "group", 
                    palette = "jco", add = "jitter", size = 0.2,
                    font.label = list(size = 1) ) + 
      stat_compare_means( comparisons = my_comparisons, size = 2 ) # Add pairwise comparisons p-value
    print(p)
  }
  dev.off()
  
  pdf("D:/Rsources/Project/SingleCelllncRNA/Results/6.evaluateSubtype/oncogenic.diff.pdf")
  for(i in colnames(oncogenic.box)[-ncol(oncogenic.box)]){
    p <- ggboxplot( oncogenic.box,
                    x = "group",
                    #combine = TRUE,
                    y = i,
                    ylab = i,
                    width = 0.9, 
                    color = "group", 
                    palette = "jco", add = "jitter", size = 0.2,
                    font.label = list(size = 1) ) + 
      stat_compare_means( comparisons = my_comparisons, size = 2 ) # Add pairwise comparisons p-value
    print(p)
  }
  dev.off()
  
  pdf("D:/Rsources/Project/SingleCelllncRNA/Results/6.evaluateSubtype/biocarta.diff.pdf")
  for(i in colnames(biocarta.box)[-ncol(biocarta.box)]){
    p <- ggboxplot( biocarta.box,
                    x = "group",
                    #combine = TRUE,
                    y = i,
                    ylab = i,
                    width = 0.9, 
                    color = "group", 
                    palette = "jco", add = "jitter", size = 0.2,
                    font.label = list(size = 1) ) + 
      stat_compare_means( comparisons = my_comparisons, size = 2 ) # Add pairwise comparisons p-value
    print(p)
  }
  dev.off()
  
  pdf("D:/Rsources/Project/SingleCelllncRNA/Results/6.evaluateSubtype/kegg.diff.pdf")
  for(i in colnames(kegg.box)[-ncol(kegg.box)]){
    p <- ggboxplot( kegg.box,
                    x = "group",
                    #combine = TRUE,
                    y = i,
                    ylab = i,
                    width = 0.9, 
                    color = "group", 
                    palette = "jco", add = "jitter", size = 0.2,
                    font.label = list(size = 1) ) + 
      stat_compare_means( comparisons = my_comparisons, size = 2 ) # Add pairwise comparisons p-value
    print(p)
  }
  dev.off()
  
  pdf("D:/Rsources/Project/SingleCelllncRNA/Results/6.evaluateSubtype/reactome.diff.pdf")
  for(i in colnames(reactome.box)[-ncol(reactome.box)]){
    p <- ggboxplot( reactome.box,
                    x = "group",
                    #combine = TRUE,
                    y = i,
                    ylab = i,
                    width = 0.9, 
                    color = "group", 
                    palette = "jco", add = "jitter", size = 0.2,
                    font.label = list(size = 1) ) + 
      stat_compare_means( comparisons = my_comparisons, size = 2 ) # Add pairwise comparisons p-value
    print(p)
  }
  dev.off()
  
  
  ## TCGA GSVA活性得分
  identical( colnames(exprSet), gbm.clinic$Patient.ID )
  rownames(exprSet) <- convert( rownames(exprSet) )
  tcga.gsva.hallmark.score.matrix <- gsva(expr = as.matrix(exprSet),
                                     gset.idx.list = hallmark.genesets,
                                     method="gsva",
                                     kcdf="Gaussian")
  tcga.gsva.immunologic.score.matrix <- gsva(expr = as.matrix(exprSet),
                                        gset.idx.list = immunologic.signature.genesets,
                                        method="gsva",
                                        kcdf="Gaussian")
  tcga.gsva.oncogenic.score.matrix <- gsva(expr = as.matrix(exprSet),
                                      gset.idx.list = immunologic.signature.genesets,
                                      method="gsva",
                                      kcdf="Gaussian")
  tcga.gsva.biocarta.score.matrix <- gsva(expr = as.matrix(exprSet),
                                     gset.idx.list = biocarta.genesets,
                                     method="gsva",
                                     kcdf="Gaussian")
  tcga.gsva.kegg.score.matrix <- gsva(expr = as.matrix(exprSet),
                                 gset.idx.list = kegg.genesets,
                                 method="gsva",
                                 kcdf="Gaussian")
  tcga.gsva.reactome.score.matrix <- gsva(expr = as.matrix(exprSet),
                                     gset.idx.list = reactome.genesets,
                                     method="gsva",
                                     kcdf="Gaussian")
  identical(colnames(tcga.gsva.hallmark.score.matrix),names(TCGA.results[[5]]$consensusClass))
  tcga.hallmark.box <- cbind.data.frame( t(tcga.gsva.hallmark.score.matrix), group = TCGA.results[[5]]$consensusClass )
  tcga.immunologic.box <- cbind.data.frame(t(tcga.gsva.immunologic.score.matrix),group=TCGA.results[[5]]$consensusClass)
  tcga.oncogenic.box <- cbind.data.frame(t(tcga.gsva.oncogenic.score.matrix),group=TCGA.results[[5]]$consensusClass)
  tcga.biocarta.box <- cbind.data.frame(t(tcga.gsva.biocarta.score.matrix),group=TCGA.results[[5]]$consensusClass)
  tcga.kegg.box <- cbind.data.frame(t(tcga.gsva.kegg.score.matrix),group=TCGA.results[[5]]$consensusClass)
  tcga.reactome.box <- cbind.data.frame(t(tcga.gsva.reactome.score.matrix),group=TCGA.results[[5]]$consensusClass)
  
  # Visualize: Specify the comparisons you want
  library(ggpubr)
  my_comparisons <- list( c("1", "2"),c("1", "3"),c("1", "4"),c("1","5"),c("2","3"),c("2","4"),c("2","5"),c("3","4"),c("3","5"),c("4","5") )
  pdf("D:/Rsources/Project/SingleCelllncRNA/Results/6.evaluateSubtype/tcga.hallmark.diff.pdf")
  for(i in colnames(tcga.hallmark.box)[-ncol(tcga.hallmark.box)]){
    p <- ggboxplot( tcga.hallmark.box,
                    x = "group",
                    #combine = TRUE,
                    y = i,
                    ylab = i,
                    width = 0.9, 
                    color = "group", 
                    palette = "jco", add = "jitter", size = 0.2,
                    font.label = list(size = 1) ) + 
      stat_compare_means( comparisons = my_comparisons, size = 2 ) # Add pairwise comparisons p-value
    print(p)
  }
  dev.off()
  
  pdf("D:/Rsources/Project/SingleCelllncRNA/Results/6.evaluateSubtype/tcga.immunologic.diff.pdf")
  for(i in colnames(tcga.immunologic.box)[-ncol(tcga.immunologic.box)]){
    p <- ggboxplot( tcga.immunologic.box,
                    x = "group",
                    #combine = TRUE,
                    y = i,
                    ylab = i,
                    width = 0.9, 
                    color = "group", 
                    palette = "jco", add = "jitter", size = 0.2,
                    font.label = list(size = 1) ) + 
      stat_compare_means( comparisons = my_comparisons, size = 2 ) # Add pairwise comparisons p-value
    print(p)
  }
  dev.off()
  
  pdf("D:/Rsources/Project/SingleCelllncRNA/Results/6.evaluateSubtype/tcga.oncogenic.diff.pdf")
  for(i in colnames(tcga.oncogenic.box)[-ncol(tcga.oncogenic.box)]){
    p <- ggboxplot( tcga.oncogenic.box,
                    x = "group",
                    #combine = TRUE,
                    y = i,
                    ylab = i,
                    width = 0.9, 
                    color = "group", 
                    palette = "jco", add = "jitter", size = 0.2,
                    font.label = list(size = 1) ) + 
      stat_compare_means( comparisons = my_comparisons, size = 2 ) # Add pairwise comparisons p-value
    print(p)
  }
  dev.off()
  
  pdf("D:/Rsources/Project/SingleCelllncRNA/Results/6.evaluateSubtype/tcga.biocarta.diff.pdf")
  for(i in colnames(tcga.biocarta.box)[-ncol(tcga.biocarta.box)]){
    p <- ggboxplot( tcga.biocarta.box,
                    x = "group",
                    #combine = TRUE,
                    y = i,
                    ylab = i,
                    width = 0.9, 
                    color = "group", 
                    palette = "jco", add = "jitter", size = 0.2,
                    font.label = list(size = 1) ) + 
      stat_compare_means( comparisons = my_comparisons, size = 2 ) # Add pairwise comparisons p-value
    print(p)
  }
  dev.off()
  
  pdf("D:/Rsources/Project/SingleCelllncRNA/Results/6.evaluateSubtype/tcga.kegg.diff.pdf")
  for(i in colnames(tcga.kegg.box)[-ncol(tcga.kegg.box)]){
    p <- ggboxplot( tcga.kegg.box,
                    x = "group",
                    #combine = TRUE,
                    y = i,
                    ylab = i,
                    width = 0.9, 
                    color = "group", 
                    palette = "jco", add = "jitter", size = 0.2,
                    font.label = list(size = 1) ) + 
      stat_compare_means( comparisons = my_comparisons, size = 2 ) # Add pairwise comparisons p-value
    print(p)
  }
  dev.off()
  
  pdf("D:/Rsources/Project/SingleCelllncRNA/Results/6.evaluateSubtype/tcga.reactome.diff.pdf")
  for(i in colnames(tcga.reactome.box)[-ncol(tcga.reactome.box)]){
    p <- ggboxplot( tcga.reactome.box,
                    x = "group",
                    #combine = TRUE,
                    y = i,
                    ylab = i,
                    width = 0.9, 
                    color = "group", 
                    palette = "jco", add = "jitter", size = 0.2,
                    font.label = list(size = 1) ) + 
      stat_compare_means( comparisons = my_comparisons, size = 2 ) # Add pairwise comparisons p-value
    print(p)
  }
  dev.off()
  
  
  ## 肿瘤纯度差异刻画
  identical( names(TCGA.results[[5]]$consensusClass),gbm.clinic$Patient.ID )
  gbm.clinic1 <- gbm.clinic[match(names(TCGA.results[[5]]$consensusClass),gbm.clinic$Patient.ID),]
  identical( names(TCGA.results[[5]]$consensusClass),gbm.clinic1$Patient.ID )
  gbm.clinic1$concensusSubtype <-  TCGA.results[[5]]$consensusClass
  my_comparisons <- list( c("1", "2"),c("1", "3"),c("1", "4"),c("1","5"),c("2","3"),c("2","4"),c("2","5"),c("3","4"),c("3","5"),c("4","5") )
  pdf("D:/Rsources/Project/SingleCelllncRNA/Results/6.evaluateSubtype/tcga.tumour.purity.diff.pdf")
  p <- ggboxplot( gbm.clinic1,
                  x = "concensusSubtype",
                  #combine = TRUE,
                  y = "ABSOLUTE.purity",
                  ylab = "Tumour purity",
                  width = 0.9, 
                  color = "concensusSubtype", 
                  palette = "jco", add = "jitter", size = 0.2,
                  font.label = list(size = 1) ) + 
    stat_compare_means( comparisons = my_comparisons, size = 2 ) # Add pairwise comparisons p-value
  print(p)
  dev.off()
  
  
  
  
}
