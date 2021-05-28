#' @title lncRNA signature分型分析
#' @description  TCGA GSE7696 GSE42669 GSE108474 GSE72951
#' @author  shi jian

## 导入GSE108474 GSE7696 GSE42669  GSE72951数据集
{
  GSE7696_GPL570_exprSet <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE7696/GSE7696-GPL570-expression.MATRIX.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
  GSE7696_GPL570_pData <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE7696/GSE7696-GPL570-pData.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
  GSE7696_GPL570_probeAnnotation <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE7696/GSE7696-GPL570-ProbeAnnotation.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
  
  GSE108474_GPL570_exprSet <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE108474/GSE108474-GPL570-expression.MATRIX.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
  GSE108474_GPL570_pData <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE108474/GSE108474-GPL570-pData.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
  GSE108474_GPL570_pData1 <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE108474/GSE108474_REMBRANDT_clinical.data.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
  GSE108474_GPL570_probeAnnotation <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE108474/GSE108474-GPL570-ProbeAnnotation.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass

  
  GSE42669_GPL6244_exprSet <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE42669/GSE42669-GPL6244-expression.MATRIX.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
  GSE42669_GPL6244_pData <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE42669/GSE42669-GPL6244-pData.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
  GSE42669_GPL6244_probeAnnotation <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE42669/GSE42669-GPL6244-ProbeAnnotation.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
  
  GSE72951_GPL14951_exprSet <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE72951/GSE72951-GPL14951-expression.MATRIX.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
  GSE72951_GPL14951_pData <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE72951/GSE72951-GPL14951-pData.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
  GSE72951_GPL14951_probeAnnotation <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE72951/GSE72951-GPL14951-ProbeAnnotation.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
}

## 探针注释
{
  ## GSE7696表达谱探针注释
  GSE7696.anno.res <- probeAnnotation(exprSet = GSE7696_GPL570_exprSet,
                                      probeAnnotation = GSE7696_GPL570_probeAnnotation,
                                      id = "ID",
                                      symbol = "Gene.Symbol",
                                      entrezid = "ENTREZ_GENE_ID",
                                      index = "ENTREZ_GENE_ID",
                                      fromtype = "ENTREZID",
                                      totype = c("entrez","symbol","ensemble")[3])
  
  ## GSE108474表达谱探针注释
  GSE108474_GPL570_pData$characteristics_ch1.1[grep("glioblastoma multiforme",GSE108474_GPL570_pData$characteristics_ch1.1)] <- "gbm"
  GSE108474_GPL570_pData$characteristics_ch1.1[grep("normal",GSE108474_GPL570_pData$characteristics_ch1.1)] <- "normal" #四个正常样本
  GSE108474_GPL570_exprSet <- GSE108474_GPL570_exprSet[,which(GSE108474_GPL570_pData$characteristics_ch1.1 == "gbm" | GSE108474_GPL570_pData$characteristics_ch1.1 == "normal")]
  GSE108474_GPL570_exprSet <- GSE108474_GPL570_exprSet[,setdiff(colnames(GSE108474_GPL570_exprSet),c("GSM2899720","GSM2899721","GSM2899722","GSM2899723","GSM2899724","GSM2899725","GSM2899726"))]
  GSE108474_GPL570_pData <- GSE108474_GPL570_pData[which(GSE108474_GPL570_pData$characteristics_ch1.1 == "gbm" | GSE108474_GPL570_pData$characteristics_ch1.1 == "normal"),]
  GSE108474_GPL570_pData <- GSE108474_GPL570_pData[GSE108474_GPL570_pData$geo_accession %in% setdiff(colnames(GSE108474_GPL570_exprSet),c("GSM2899720","GSM2899721","GSM2899722","GSM2899723","GSM2899724","GSM2899725","GSM2899726")),]
  
  GSE108474.anno.res <- probeAnnotation(exprSet = GSE108474_GPL570_exprSet,
                                        probeAnnotation = GSE108474_GPL570_probeAnnotation,
                                        id = "ID",
                                        symbol = "Gene.Symbol",
                                        entrezid = "ENTREZ_GENE_ID",
                                        index = "ENTREZ_GENE_ID",
                                        fromtype = "ENTREZID",
                                        totype = c("entrez","symbol","ensemble")[3])

  
  ## GSE42669表达谱探针注释 GPL6244特殊，特异性处理
  ids <- GSE42669_GPL6244_probeAnnotation[,c("ID","gene_assignment")]
  ids <- ids[which(!(ids[,"gene_assignment"] == "")),]  #过滤掉没有基因注释的探针
  a <- strsplit(as.character(ids[,"gene_assignment"]), " // ")
  a <- lapply(a,function(x){ if( length(x) >= 2 ){x[2]}else{x} })
  ID2gene <- cbind.data.frame(ids[,1], unlist(a))
  colnames(ID2gene) <- c("ID","symbol")
  ID2gene <- ID2gene[ID2gene$symbol != "---", ]
  GSE42669_GPL6244_exprSet <- GSE42669_GPL6244_exprSet[rownames(GSE42669_GPL6244_exprSet) %in% ID2gene[ , "ID"], ]
  ID2gene <- ID2gene[ match(rownames(GSE42669_GPL6244_exprSet), ID2gene[ , "ID"]), ]
  ID2gene <- cbind.data.frame( ID2gene[,1],
                               convert(ID2gene[,2],fromtype="SYMBOL",totype ="ENTREZID"),
                               ID2gene[,2],
                               convert(ID2gene[,2],fromtype="SYMBOL",totype ="ENSEMBL"),
                               convert(ID2gene[,2],fromtype="SYMBOL",totype ="gene_type") )
  colnames(ID2gene) <- c("id","entrez","symbol","ensemble","gene_type")
  table(ID2gene$gene_type) # 418个LncRNA
  #ID2gene <- ID2gene[!is.na(ID2gene$ensemble),]
  MAX <- by(GSE42669_GPL6244_exprSet, ID2gene[ ,"ensemble"],function(x) rownames(x)[which.max(rowMeans(x))])
  MAX <- as.character( MAX )
  GSE42669_GPL6244_exprSet <- GSE42669_GPL6244_exprSet[rownames(GSE42669_GPL6244_exprSet) %in% MAX, ]
  rownames(GSE42669_GPL6244_exprSet) <- ID2gene[match(rownames(GSE42669_GPL6244_exprSet),ID2gene[, 1]), "ensemble"]
  GPL6244_ID2gene <- ID2gene
  
  ## GSE72951表达谱探针注释
  GSE72951.anno.res <- probeAnnotation(exprSet = GSE72951_GPL14951_exprSet,
                                       probeAnnotation = GSE72951_GPL14951_probeAnnotation,
                                       id = "ID",
                                       symbol = "Symbol",
                                       entrezid = "Entrez_Gene_ID",
                                       index = "Entrez_Gene_ID",
                                       fromtype = "ENTREZID",
                                       totype = c("entrez","symbol","ensemble")[3])
  
  GSE42669.anno.res <- list(GSE42669_GPL6244_exprSet = GSE42669_GPL6244_exprSet,GPL6244_ID2gene = GPL6244_ID2gene)
  save(GSE7696.anno.res,
       GSE108474.anno.res,
       GSE42669.anno.res = GSE42669.anno.res,
       GSE72951.anno.res,file = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/subtype.probe.annotation.rda")
}

## 一致性聚类应用于GSE108474
{
  load("D:/Rsources/Project/SingleCelllncRNA/Results/1.singleCellSignature/lnc_sig.rda")
  load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/subtype.probe.annotation.rda")
  lnc_sig <- unique(c("ENSG00000176124","ENSG00000231889","ENSG00000232677","ENSG00000245694","ENSG00000224032","ENSG00000226950","ENSG00000231419","ENSG00000231889","ENSG00000232677","ENSG00000245694","ENSG00000176124","ENSG00000231889","ENSG00000232677","ENSG00000245694"))
  GSE108474.exprSet <- GSE108474.anno.res$exprSet[rownames(GSE108474.anno.res$exprSet) %in% lnc_sig,]
  GSE108474.exprSet <- GSE108474.exprSet[ ,GSE108474_GPL570_pData$characteristics_ch1.1 == "gbm"]
  
  library(ConsensusClusterPlus)
  results <- ConsensusClusterPlus(as.matrix(GSE108474.exprSet),
                                   maxK=8,
                                   reps=500,
                                   pItem=0.8,
                                   pFeature=1,
                                   title="GSE108474",
                                   clusterAlg="pam",
                                   distance="pearson",
                                   seed=1262118388.71279,
                                   writeTable = T,
                                   verbose = T,
                                   plot="pdf")
  optimalCutoff <- function(maxK,Consensus.result){
    # 用PAC的方法确定最佳聚类数,面积最小值对应K为最佳K
    Kvec <- 2 : maxK
    x1 <- 0.1;x2 <- 0.9        # threshold defining the intermediate sub-interval
    PAC <- rep(NA,length(Kvec)) 
    names(PAC) <- paste( "K=", Kvec, sep = "" )  # from 2 to maxK
    for(i in Kvec){
      M <- Consensus.result[[i]]$consensusMatrix
      Fn <- ecdf(M[lower.tri(M)])          # M 为计算出共识矩阵
      PAC[i-1] <- Fn(x2) - Fn(x1)
    } 
    optK <- Kvec[which.min(PAC)]  # 理想的K值
    return(optK)
  }
  optimalCutoff(8,results)  # 8个亚型是最好的
  save( results,file="D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/ConsensusClusterResult.rda" )

}

## GSE108474一致性聚类结果进行生存分析
{
  load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/ConsensusClusterResult.rda")
  # 整理GSE108474生存数据
  GSE108474_GPL570_pData <- GSE108474_GPL570_pData[GSE108474_GPL570_pData$characteristics_ch1.1 == "gbm",]
  GSE108474.clinical <- merge(GSE108474_GPL570_pData,GSE108474_GPL570_pData1,by.x="title",by.y="SUBJECT_ID",all.x = T)
  identical(GSE108474.clinical$geo_accession,names(results[[5]]$consensusClass))
  GSE108474.clinical$label2 <- results[[2]]$consensusClass
  GSE108474.clinical$label3 <- results[[3]]$consensusClass
  GSE108474.clinical$label4 <- results[[4]]$consensusClass
  GSE108474.clinical$label5 <- results[[5]]$consensusClass
  GSE108474.clinical$label6 <- results[[6]]$consensusClass
  GSE108474.clinical$label7 <- results[[7]]$consensusClass
  GSE108474.clinical$label8 <- results[[8]]$consensusClass
  
  head(GSE108474.clinical)
  FastSurvKM(GSE108474.clinical,
             time = "OVERALL_SURVIVAL_MONTHS",
             status = "EVENT_OS",
             marker = "label2",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE108474/")
  
  FastSurvKM(GSE108474.clinical,
             time = "OVERALL_SURVIVAL_MONTHS",
             status = "EVENT_OS",
             marker = "label3",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE108474/")
  
  FastSurvKM(GSE108474.clinical,
             time = "OVERALL_SURVIVAL_MONTHS",
             status = "EVENT_OS",
             marker = "label4",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE108474/")
  
  FastSurvKM(GSE108474.clinical,
             time = "OVERALL_SURVIVAL_MONTHS",
             status = "EVENT_OS",
             marker = "label5",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE108474/")
  
  FastSurvKM(GSE108474.clinical,
             time = "OVERALL_SURVIVAL_MONTHS",
             status = "EVENT_OS",
             marker = "label6",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE108474/")
  
  FastSurvKM(GSE108474.clinical,
             time = "OVERALL_SURVIVAL_MONTHS",
             status = "EVENT_OS",
             marker = "label7",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE108474/")
  
  FastSurvKM(GSE108474.clinical,
             time = "OVERALL_SURVIVAL_MONTHS",
             status = "EVENT_OS",
             marker = "label8",
             color = NULL,
             upper.time = NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE108474/")
  
}

## 一致性聚类应用于GSE7696
{
  load("D:/Rsources/Project/SingleCelllncRNA/Results/1.singleCellSignature/lnc_sig.rda")
  load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/subtype.probe.annotation.rda")
  GSE7696.exprSet <- GSE7696.anno.res$exprSet[ rownames(GSE7696.anno.res$exprSet) %in% lnc_sig, ]
  GSE7696.exprSet <- GSE7696.exprSet[ , grep("GBM",GSE7696_GPL570_pData$characteristics_ch1.1) ]
  GSE7696_GPL570_pData <- GSE7696_GPL570_pData[ grep("GBM",GSE7696_GPL570_pData$characteristics_ch1.1), ]
  identical( colnames(GSE7696.exprSet), GSE7696_GPL570_pData$geo_accession )
  
  
  library(ConsensusClusterPlus)
  GSE7696.results <- ConsensusClusterPlus(as.matrix(GSE7696.exprSet),
                                  maxK=8,
                                  reps=500,
                                  pItem=0.8,
                                  pFeature=1,
                                  title = "GSE7696",
                                  clusterAlg = "pam",
                                  distance = "pearson",
                                  seed = 1262118388.71279,
                                  writeTable = T,
                                  verbose = T,
                                  plot = "pdf")
  optimalCutoff <- function(maxK,
                            Consensus.result){
    # 用PAC的方法确定最佳聚类数,面积最小值对应K为最佳K
    Kvec <- 2 : maxK
    x1 <- 0.1;x2 <- 0.9        # threshold defining the intermediate sub-interval
    PAC <- rep(NA,length(Kvec)) 
    names(PAC) <- paste( "K=", Kvec, sep = "" )  # from 2 to maxK
    for(i in Kvec){
      M <- Consensus.result[[i]]$consensusMatrix
      Fn <- ecdf(M[lower.tri(M)])          # M 为计算出共识矩阵
      PAC[i-1] <- Fn(x2) - Fn(x1)
    } 
    optK <- Kvec[which.min(PAC)]  # 理想的K值
    return(optK)
  }
  optimalCutoff(8,GSE7696.results)  # 8个亚型是最好的
  save( GSE7696.results,file="D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE7696.ConsensusClusterResult.rda" )
  
}

## GSE7696一致性聚类结果进行生存分析
{
  load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE7696.ConsensusClusterResult.rda")
  # 整理GSE7696生存数据
  GSE7696_GPL570_pData <- GSE7696_GPL570_pData[grep("GBM",GSE7696_GPL570_pData$characteristics_ch1.1),]
  identical(GSE7696_GPL570_pData$geo_accession,names(GSE7696.results[[5]]$consensusClass))
  GSE7696_GPL570_pData$label2 <- GSE7696.results[[2]]$consensusClass
  GSE7696_GPL570_pData$label3 <- GSE7696.results[[3]]$consensusClass
  GSE7696_GPL570_pData$label4 <- GSE7696.results[[4]]$consensusClass
  GSE7696_GPL570_pData$label5 <- GSE7696.results[[5]]$consensusClass
  GSE7696_GPL570_pData$label6 <- GSE7696.results[[6]]$consensusClass
  GSE7696_GPL570_pData$label7 <- GSE7696.results[[7]]$consensusClass
  GSE7696_GPL570_pData$label8 <- GSE7696.results[[8]]$consensusClass
  
  head(GSE7696_GPL570_pData)
  library(stringr)
  GSE7696_GPL570_pData$characteristics_ch1.7 <- as.numeric(str_split(GSE7696_GPL570_pData$characteristics_ch1.7,":",simplify=T)[,2]) 
  GSE7696_GPL570_pData$characteristics_ch1.6 <- as.numeric(str_split(GSE7696_GPL570_pData$characteristics_ch1.6,":",simplify=T)[,2]) 
  FastSurvKM(GSE7696_GPL570_pData,
             time = "characteristics_ch1.7",
             status = "characteristics_ch1.6",
             marker = "label2",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE7696/")
  
  FastSurvKM(GSE7696_GPL570_pData,
             time = "characteristics_ch1.7",
             status = "characteristics_ch1.6",
             marker = "label3",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE7696/")
  
  FastSurvKM(GSE7696_GPL570_pData,
             time = "characteristics_ch1.7",
             status = "characteristics_ch1.6",
             marker = "label4",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE7696/")
  
  FastSurvKM(GSE7696_GPL570_pData,
             time = "characteristics_ch1.7",
             status = "characteristics_ch1.6",
             marker = "label5",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE7696/")
  
  FastSurvKM(GSE7696_GPL570_pData,
             time = "characteristics_ch1.7",
             status = "characteristics_ch1.6",
             marker = "label6",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE7696/")
  
  FastSurvKM(GSE7696_GPL570_pData,
             time = "characteristics_ch1.7",
             status = "characteristics_ch1.6",
             marker = "label7",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE7696/")
  
  FastSurvKM(GSE7696_GPL570_pData,
             time = "characteristics_ch1.7",
             status = "characteristics_ch1.6",
             marker = "label8",
             color = NULL,
             upper.time = NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE7696/")
}

## 一致性聚类应用于GSE42669
{
  load("D:/Rsources/Project/SingleCelllncRNA/Results/1.singleCellSignature/lnc_sig.rda")
  load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/subtype.probe.annotation.rda")
  GSE42669.exprSet <- GSE42669.anno.res$GSE42669_GPL6244_exprSet[rownames(GSE42669.anno.res$GSE42669_GPL6244_exprSet) %in% lnc_sig, ]
  GSE42669.exprSet <- GSE42669.exprSet[ ,grep("glioblastoma",GSE42669_GPL6244_pData$characteristics_ch1.1)]
  GSE42669_GPL6244_pData <- GSE42669_GPL6244_pData[ grep("glioblastoma",GSE42669_GPL6244_pData$characteristics_ch1.1), ]
  identical( colnames(GSE42669.exprSet), GSE42669_GPL6244_pData$geo_accession )
  names(GSE42669.anno.res)
  
  library(ConsensusClusterPlus)
  GSE42669.results <- ConsensusClusterPlus(as.matrix(GSE42669.exprSet),
                                          maxK=8,
                                          reps=500,
                                          pItem=0.8,
                                          pFeature=1,
                                          title = "GSE42669",
                                          clusterAlg = "pam",
                                          distance = "pearson",
                                          seed = 1262118388.71279,
                                          writeTable = T,
                                          verbose = T,
                                          plot = "pdf")
  optimalCutoff <- function(maxK,
                            Consensus.result){
    # 用PAC的方法确定最佳聚类数,面积最小值对应K为最佳K
    Kvec <- 2 : maxK
    x1 <- 0.1;x2 <- 0.9        # threshold defining the intermediate sub-interval
    PAC <- rep(NA,length(Kvec)) 
    names(PAC) <- paste( "K=", Kvec, sep = "" )  # from 2 to maxK
    for(i in Kvec){
      M <- Consensus.result[[i]]$consensusMatrix
      Fn <- ecdf(M[lower.tri(M)])          # M 为计算出共识矩阵
      PAC[i-1] <- Fn(x2) - Fn(x1)
    } 
    optK <- Kvec[which.min(PAC)]  # 理想的K值
    return(optK)
  }
  optimalCutoff(8,GSE42669.results)  # 8个亚型是最好的
  save( GSE42669.results,file="D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE42669.ConsensusClusterResult.rda" )
}

## 一致性聚类应用于GSE72951
{
  load("D:/Rsources/Project/SingleCelllncRNA/Results/1.singleCellSignature/lnc_sig.rda")
  load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/subtype.probe.annotation.rda")
  GSE72951.exprSet <- GSE72951.anno.res$exprSet[rownames(GSE72951.anno.res$exprSet) %in% c(gain_driver_lncRNA,loss_driver_lncRNA), ]
  GSE72951.exprSet <- GSE72951.exprSet[, grep("GBM",GSE72951_GPL14951_pData$characteristics_ch1)]
  GSE72951_GPL14951_pData <- GSE72951_GPL14951_pData[ grep("GBM",GSE72951_GPL14951_pData$characteristics_ch1), ]
  identical( colnames(GSE72951.exprSet), GSE72951_GPL14951_pData$geo_accession )
  names(GSE72951.anno.res)
  
  library(ConsensusClusterPlus)
  GSE72951.results <- ConsensusClusterPlus(as.matrix(GSE72951.exprSet),
                                           maxK=8,
                                           reps=500,
                                           pItem=0.8,
                                           pFeature=1,
                                           title = "GSE72951",
                                           clusterAlg = "pam",
                                           distance = "pearson",
                                           seed = 1262118388.71279,
                                           writeTable = T,
                                           verbose = T,
                                           plot = "pdf")
  optimalCutoff <- function(maxK,
                            Consensus.result){
    # 用PAC的方法确定最佳聚类数,面积最小值对应K为最佳K
    Kvec <- 2 : maxK
    x1 <- 0.1;x2 <- 0.9        # threshold defining the intermediate sub-interval
    PAC <- rep(NA,length(Kvec)) 
    names(PAC) <- paste( "K=", Kvec, sep = "" )  # from 2 to maxK
    for(i in Kvec){
      M <- Consensus.result[[i]]$consensusMatrix
      Fn <- ecdf(M[lower.tri(M)])          # M 为计算出共识矩阵
      PAC[i-1] <- Fn(x2) - Fn(x1)
    } 
    optK <- Kvec[which.min(PAC)]  # 理想的K值
    return(optK)
  }
  optimalCutoff(8,GSE72951.results)  # 8个亚型是最好的
  save( GSE72951.results,file="D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/GSE72951.ConsensusClusterResult.rda" )
}

## 基于TCGA数据构建单样本分类器
{
  ## 导入TCGA表达谱数据和表型数据
  exprSet <- get(load("D:/Rsources/Project/SingleCelllncRNA/Data/TCGA-GBM/Transcriptome Profiling.Gene Expression Quantification.HTSeq - FPKM.rda"))
  gbm.clinic <- get(load("D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/gbm.clinic.rda"))
  exprSet <- assay(exprSet) 
  head(exprSet)
  table(substring(colnames(exprSet),14,15) < 10) # 5个正常样本，169个疾病样本
  exprSet <- exprSet[lnc_sig[lnc_sig %in% rownames(exprSet)], ]
  exprSet <- exprSet[ ,substring(colnames(exprSet),14,15) < 10]
  colnames(exprSet) <- substring( colnames(exprSet) , 1, 12)
  #提取表达谱样本中病人id和临床病人id overlap部分
  common_patients <- intersect( colnames(exprSet), gbm.clinic$Patient.ID ) 
  exprSet <- exprSet[ , colnames(exprSet) %in% common_patients ]
  gbm.clinic <- gbm.clinic[ gbm.clinic$Patient.ID %in% common_patients, ]
  
  train_data <- cbind.data.frame(t(exprSet),gbm.clinic[match(colnames(exprSet), gbm.clinic$Patient.ID), "Original.Subtype" ])
  colnames(train_data)[ncol(train_data)] <- "subtype"
  train_data$subtype <- as.character(train_data$subtype)
  train_data <- train_data[,colSums(train_data[,-ncol(train_data)]) > 0]
  train_data[,1:48] <- log2(train_data[,1:48]+1)
  table( train_data$subtype, useNA="ifany" )
  # Classical      G-CIMP Mesenchymal      Neural   Proneural        <NA> 
  #   42           8          56          28          29           4
  
  Classical_train_data <- train_data
  Classical_train_data$label <- ifelse(Classical_train_data$subtype != "Classical" | is.na(Classical_train_data$subtype),0,1)
  Classical_train_data$label1 <- ifelse(Classical_train_data$subtype != "Classical" | is.na(Classical_train_data$subtype),"other","Classical")
  
  Mesenchymal_train_data <- train_data
  Mesenchymal_train_data$label <- ifelse(Mesenchymal_train_data$subtype != "Mesenchymal" | is.na(Mesenchymal_train_data$subtype),0,1)
  Mesenchymal_train_data$label1 <- ifelse(Mesenchymal_train_data$subtype != "Mesenchymal" | is.na(Mesenchymal_train_data$subtype),"other","Mesenchymal")
  
  Neural_train_data <- train_data
  Neural_train_data$label <- ifelse(Neural_train_data$subtype != "Neural" | is.na(Neural_train_data$subtype),0,1)
  Neural_train_data$label1 <- ifelse(Neural_train_data$subtype != "Neural" | is.na(Neural_train_data$subtype),"other","Neural")
  
  Proneural_train_data <- train_data
  Proneural_train_data$label <- ifelse(Proneural_train_data$subtype != "Proneural" | is.na(Proneural_train_data$subtype),0,1)
  Proneural_train_data$label1 <- ifelse(Proneural_train_data$subtype != "Proneural" | is.na(Proneural_train_data$subtype),"other","Proneural")
  
  G_CIMP_train_data <- train_data
  G_CIMP_train_data$label <- ifelse(G_CIMP_train_data$subtype != "G-CIMP" | is.na(G_CIMP_train_data$subtype),0,1)
  G_CIMP_train_data$label1 <- ifelse(G_CIMP_train_data$subtype != "G-CIMP" | is.na(G_CIMP_train_data$subtype),"other","G-CIMP")
  
  ## 构建检测数据集GSE108474 GSE7696 GSE42669  GSE72951
  GSE108474_common_lnc <- intersect( rownames(GSE108474.anno.res$exprSet), lnc_sig )
  GSE108474.exprSet <- GSE108474.anno.res$exprSet[GSE108474_common_lnc,]
  GSE108474.test.data <- t(GSE108474.exprSet)
  
  GSE7696_common_lnc <- intersect( rownames(GSE7696.anno.res$exprSet), lnc_sig )
  GSE7696.exprSet <- GSE7696.anno.res$exprSet[GSE7696_common_lnc,]
  GSE7696.test.data <- t(GSE7696.exprSet)
  
  GSE42669_common_lnc <- intersect( rownames(GSE42669.anno.res$GSE42669_GPL6244_exprSet), lnc_sig )
  GSE42669.exprSet <- GSE42669.anno.res$GSE42669_GPL6244_exprSet[GSE42669_common_lnc,]
  GSE42669.test.data <- t(GSE42669.exprSet)
  
  GSE72951_common_lnc <- intersect( rownames(GSE72951.anno.res$exprSet), lnc_sig )
  GSE72951.exprSet <- GSE72951.anno.res$exprSet[GSE72951_common_lnc,]
  GSE72951.test.data <- t(GSE72951.exprSet)
  
  ## lasso logistic回归分类器
  GSE108474_Classical.res <- FastGlmnet(data = Classical_train_data,
                                        test.data = GSE108474.test.data,
                                        response = "label",
                                        cluster = GSE108474_common_lnc,
                                        col = "label1",
                                        control = "other",
                                        alpha = c(1,0)[2],
                                        k = 10,
                                        family = c("binomial","cox")[1],
                                        optimize.method = c("min","lse")[1],
                                        seed = 123,
                                        plot.label = T,
                                        plot.lwd = 2,
                                        line.lty = 3,
                                        line.lwd = 2,
                                        line.col = "blue",
                                        verbose = F,
                                        name = "GSE108474_Classical")
  table(GSE108474_Classical.res$Data$predicted.classes)
  
  GSE108474_Mesenchymal.res <- FastGlmnet(data = Mesenchymal_train_data,
                                          test.data = GSE108474.test.data,
                                          response = "label",
                                          cluster = GSE108474_common_lnc,
                                          col = "label1",
                                          control = "other",
                                          alpha = c(1,0)[2],
                                          k = 10,
                                          family = c("binomial","cox")[1],
                                          optimize.method = c("min","lse")[1],
                                          seed = 123,
                                          plot.label = T,
                                          plot.lwd = 2,
                                          line.lty = 3,
                                          line.lwd = 2,
                                          line.col = "blue",
                                          verbose = F,
                                          name = "GSE108474_Mesenchymal")
  table(GSE108474_Mesenchymal.res$Data$predicted.classes)
  
  GSE108474_Neural.res <- FastGlmnet(data = Neural_train_data,
                                     test.data = GSE108474.test.data,
                                     response = "label",
                                     cluster = GSE108474_common_lnc,
                                     col = "label1",
                                     control = "other",
                                     alpha = c(1,0)[2],
                                     k = 10,
                                     family = c("binomial","cox")[1],
                                     optimize.method = c("min","lse")[1],
                                     seed = 123,
                                     plot.label = T,
                                     plot.lwd = 2,
                                     line.lty = 3,
                                     line.lwd = 2,
                                     line.col = "blue",
                                     verbose = F,
                                     name = "GSE108474_Neural")
  table(GSE108474_Neural.res$Data$predicted.classes)
  
  GSE108474_Proneural.res <- FastGlmnet(data = Proneural_train_data,
                                        test.data = GSE108474.test.data,
                                        response = "label",
                                        cluster = GSE108474_common_lnc,
                                        col = "label1",
                                        control = "other",
                                        alpha = c(1,0)[2],
                                        k = 10,
                                        family = c("binomial","cox")[1],
                                        optimize.method = c("min","lse")[1],
                                        seed = 123,
                                        plot.label = T,
                                        plot.lwd = 2,
                                        line.lty = 3,
                                        line.lwd = 2,
                                        line.col = "blue",
                                        verbose = F,
                                        name = "GSE108474_Proneural")
  table(GSE108474_Proneural.res$Data$predicted.classes)
  
  GSE108474_G_CIMP.res <- FastGlmnet(data = G_CIMP_train_data,
                                        test.data = GSE108474.test.data,
                                        response = "label",
                                        cluster = GSE108474_common_lnc,
                                        col = "label1",
                                        control = "other",
                                        alpha = c(1,0)[2],
                                        k = 10,
                                        family = c("binomial","cox")[1],
                                        optimize.method = c("min","lse")[1],
                                        seed = 123,
                                        plot.label = T,
                                        plot.lwd = 2,
                                        line.lty = 3,
                                        line.lwd = 2,
                                        line.col = "blue",
                                        verbose = F,
                                        name = "GSE108474_G_CIMP")
  table(GSE108474_G_CIMP.res$Data$predicted.classes)
  
  GSE7696_Classical.res <- FastGlmnet(data = Classical_train_data,
                                        test.data = GSE7696.test.data,
                                        response = "label",
                                        cluster = GSE7696_common_lnc,
                                        col = "label1",
                                        control = "other",
                                        alpha = c(1,0)[2],
                                        k = 10,
                                        family = c("binomial","cox")[1],
                                        optimize.method = c("min","lse")[1],
                                        seed = 123,
                                        plot.label = T,
                                        plot.lwd = 2,
                                        line.lty = 3,
                                        line.lwd = 2,
                                        line.col = "blue",
                                        verbose = F,
                                        name = "GSE7696_Classical")
  table(GSE7696_Classical.res$Data$predicted.classes)
  
  GSE7696_Mesenchymal.res <- FastGlmnet(data = Mesenchymal_train_data,
                                          test.data = GSE7696.test.data,
                                          response = "label",
                                          cluster = GSE7696_common_lnc,
                                          col = "label1",
                                          control = "other",
                                          alpha = c(1,0)[2],
                                          k = 10,
                                          family = c("binomial","cox")[1],
                                          optimize.method = c("min","lse")[1],
                                          seed = 123,
                                          plot.label = T,
                                          plot.lwd = 2,
                                          line.lty = 3,
                                          line.lwd = 2,
                                          line.col = "blue",
                                          verbose = F,
                                          name = "GSE7696_Mesenchymal")
  table(GSE7696_Mesenchymal.res$Data$predicted.classes)
  
  GSE7696_Neural.res <- FastGlmnet(data = Neural_train_data,
                                     test.data = GSE7696.test.data,
                                     response = "label",
                                     cluster = GSE7696_common_lnc,
                                     col = "label1",
                                     control = "other",
                                     alpha = c(1,0)[2],
                                     k = 10,
                                     family = c("binomial","cox")[1],
                                     optimize.method = c("min","lse")[1],
                                     seed = 123,
                                     plot.label = T,
                                     plot.lwd = 2,
                                     line.lty = 3,
                                     line.lwd = 2,
                                     line.col = "blue",
                                     verbose = F,
                                     name = "GSE7696_Neural")
  table(GSE7696_Neural.res$Data$predicted.classes)
  
  GSE7696_Proneural.res <- FastGlmnet(data = Proneural_train_data,
                                        test.data = GSE7696.test.data,
                                        response = "label",
                                        cluster = GSE7696_common_lnc,
                                        col = "label1",
                                        control = "other",
                                        alpha = c(1,0)[2],
                                        k = 10,
                                        family = c("binomial","cox")[1],
                                        optimize.method = c("min","lse")[1],
                                        seed = 123,
                                        plot.label = T,
                                        plot.lwd = 2,
                                        line.lty = 3,
                                        line.lwd = 2,
                                        line.col = "blue",
                                        verbose = F,
                                        name = "GSE7696_Proneural")
  table(GSE7696_Proneural.res$Data$predicted.classes)
  
  GSE7696_G_CIMP.res <- FastGlmnet(data = G_CIMP_train_data,
                                     test.data = GSE7696.test.data,
                                     response = "label",
                                     cluster = GSE7696_common_lnc,
                                     col = "label1",
                                     control = "other",
                                     alpha = c(1,0)[2],
                                     k = 10,
                                     family = c("binomial","cox")[1],
                                     optimize.method = c("min","lse")[1],
                                     seed = 123,
                                     plot.label = T,
                                     plot.lwd = 2,
                                     line.lty = 3,
                                     line.lwd = 2,
                                     line.col = "blue",
                                     verbose = F,
                                     name = "GSE7696_G_CIMP")
  table(GSE7696_G_CIMP.res$Data$predicted.classes)
  
  GSE42669_Classical.res <- FastGlmnet(data = Classical_train_data,
                                      test.data = GSE42669.test.data,
                                      response = "label",
                                      cluster = GSE42669_common_lnc,
                                      col = "label1",
                                      control = "other",
                                      alpha = c(1,0)[2],
                                      k = 10,
                                      family = c("binomial","cox")[1],
                                      optimize.method = c("min","lse")[1],
                                      seed = 123,
                                      plot.label = T,
                                      plot.lwd = 2,
                                      line.lty = 3,
                                      line.lwd = 2,
                                      line.col = "blue",
                                      verbose = F,
                                      name = "GSE42669_Classical")
  table(GSE42669_Classical.res$Data$predicted.classes)
  
  GSE42669_Mesenchymal.res <- FastGlmnet(data = Mesenchymal_train_data,
                                        test.data = GSE42669.test.data,
                                        response = "label",
                                        cluster = GSE42669_common_lnc,
                                        col = "label1",
                                        control = "other",
                                        alpha = c(1,0)[2],
                                        k = 10,
                                        family = c("binomial","cox")[1],
                                        optimize.method = c("min","lse")[1],
                                        seed = 123,
                                        plot.label = T,
                                        plot.lwd = 2,
                                        line.lty = 3,
                                        line.lwd = 2,
                                        line.col = "blue",
                                        verbose = F,
                                        name = "GSE42669_Mesenchymal")
  table(GSE42669_Mesenchymal.res$Data$predicted.classes)
  
  GSE42669_Neural.res <- FastGlmnet(data = Neural_train_data,
                                   test.data = GSE42669.test.data,
                                   response = "label",
                                   cluster = GSE42669_common_lnc,
                                   col = "label1",
                                   control = "other",
                                   alpha = c(1,0)[2],
                                   k = 10,
                                   family = c("binomial","cox")[1],
                                   optimize.method = c("min","lse")[1],
                                   seed = 123,
                                   plot.label = T,
                                   plot.lwd = 2,
                                   line.lty = 3,
                                   line.lwd = 2,
                                   line.col = "blue",
                                   verbose = F,
                                   name = "GSE42669_Neural")
  table(GSE42669_Neural.res$Data$predicted.classes)
  
  GSE42669_Proneural.res <- FastGlmnet(data = Proneural_train_data,
                                      test.data = GSE42669.test.data,
                                      response = "label",
                                      cluster = GSE42669_common_lnc,
                                      col = "label1",
                                      control = "other",
                                      alpha = c(1,0)[2],
                                      k = 10,
                                      family = c("binomial","cox")[1],
                                      optimize.method = c("min","lse")[1],
                                      seed = 123,
                                      plot.label = T,
                                      plot.lwd = 2,
                                      line.lty = 3,
                                      line.lwd = 2,
                                      line.col = "blue",
                                      verbose = F,
                                      name = "GSE42669_Proneural")
  table(GSE42669_Proneural.res$Data$predicted.classes)
  
  GSE42669_G_CIMP.res <- FastGlmnet(data = G_CIMP_train_data,
                                   test.data = GSE42669.test.data,
                                   response = "label",
                                   cluster = GSE42669_common_lnc,
                                   col = "label1",
                                   control = "other",
                                   alpha = c(1,0)[2],
                                   k = 10,
                                   family = c("binomial","cox")[1],
                                   optimize.method = c("min","lse")[1],
                                   seed = 123,
                                   plot.label = T,
                                   plot.lwd = 2,
                                   line.lty = 3,
                                   line.lwd = 2,
                                   line.col = "blue",
                                   verbose = F,
                                   name = "GSE42669_G_CIMP")
  table(GSE42669_G_CIMP.res$Data$predicted.classes)
  
  
  GSE72951_Classical.res <- FastGlmnet(data = Classical_train_data,
                                       test.data = GSE72951.test.data,
                                       response = "label",
                                       cluster = GSE72951_common_lnc,
                                       col = "label1",
                                       control = "other",
                                       alpha = c(1,0)[2],
                                       k = 10,
                                       family = c("binomial","cox")[1],
                                       optimize.method = c("min","lse")[1],
                                       seed = 123,
                                       plot.label = T,
                                       plot.lwd = 2,
                                       line.lty = 3,
                                       line.lwd = 2,
                                       line.col = "blue",
                                       verbose = F,
                                       name = "GSE72951_Classical")
  table(GSE72951_Classical.res$Data$predicted.classes)
  
  GSE72951_Mesenchymal.res <- FastGlmnet(data = Mesenchymal_train_data,
                                         test.data = GSE72951.test.data,
                                         response = "label",
                                         cluster = GSE72951_common_lnc,
                                         col = "label1",
                                         control = "other",
                                         alpha = c(1,0)[2],
                                         k = 10,
                                         family = c("binomial","cox")[1],
                                         optimize.method = c("min","lse")[1],
                                         seed = 123,
                                         plot.label = T,
                                         plot.lwd = 2,
                                         line.lty = 3,
                                         line.lwd = 2,
                                         line.col = "blue",
                                         verbose = F,
                                         name = "GSE72951_Mesenchymal")
  table(GSE72951_Mesenchymal.res$Data$predicted.classes)
  
  GSE72951_Neural.res <- FastGlmnet(data = Neural_train_data,
                                    test.data = GSE72951.test.data,
                                    response = "label",
                                    cluster = GSE72951_common_lnc,
                                    col = "label1",
                                    control = "other",
                                    alpha = c(1,0)[2],
                                    k = 10,
                                    family = c("binomial","cox")[1],
                                    optimize.method = c("min","lse")[1],
                                    seed = 123,
                                    plot.label = T,
                                    plot.lwd = 2,
                                    line.lty = 3,
                                    line.lwd = 2,
                                    line.col = "blue",
                                    verbose = F,
                                    name = "GSE72951_Neural")
  table(GSE72951_Neural.res$Data$predicted.classes)
  
  GSE72951_Proneural.res <- FastGlmnet(data = Proneural_train_data,
                                       test.data = GSE72951.test.data,
                                       response = "label",
                                       cluster = GSE72951_common_lnc,
                                       col = "label1",
                                       control = "other",
                                       alpha = c(1,0)[2],
                                       k = 10,
                                       family = c("binomial","cox")[1],
                                       optimize.method = c("min","lse")[1],
                                       seed = 123,
                                       plot.label = T,
                                       plot.lwd = 2,
                                       line.lty = 3,
                                       line.lwd = 2,
                                       line.col = "blue",
                                       verbose = F,
                                       name = "GSE72951_Proneural")
  table(GSE72951_Proneural.res$Data$predicted.classes)
  
  GSE72951_G_CIMP.res <- FastGlmnet(data = G_CIMP_train_data,
                                    test.data = GSE72951.test.data,
                                    response = "label",
                                    cluster = GSE72951_common_lnc,
                                    col = "label1",
                                    control = "other",
                                    alpha = c(1,0)[2],
                                    k = 10,
                                    family = c("binomial","cox")[1],
                                    optimize.method = c("min","lse")[1],
                                    seed = 123,
                                    plot.label = T,
                                    plot.lwd = 2,
                                    line.lty = 3,
                                    line.lwd = 2,
                                    line.col = "blue",
                                    verbose = F,
                                    name = "GSE72951_G_CIMP")
  table(GSE72951_G_CIMP.res$Data$predicted.classes)
  
  
  ## 线性判别分析分类器
  {
    GSE108474.lda.res <- FastLDA(train.data = train_data,
                                 test.data = as.data.frame(GSE108474.test.data),
                                 response = "subtype",
                                 cluster = GSE108474_common_lnc,
                                 method = c("lda","qda","mda","fda","NaiveBayes")[1],
                                 name = "love")
    
    GSE7696.lda.res <- FastLDA(train.data = train_data,
                               test.data = as.data.frame(GSE7696.test.data),
                               response = "subtype",
                               cluster = GSE7696_common_lnc,
                               method = c("lda","qda","mda","fda","NaiveBayes")[1],
                               name = "love")
    
    GSE42669.lda.res <- FastLDA(train.data = train_data,
                                test.data = as.data.frame(GSE42669.test.data),
                                response = "subtype",
                                cluster = GSE42669_common_lnc,
                                method = c("lda","qda","mda","fda","NaiveBayes")[1],
                                name = "love")
    
    GSE72951.lda.res <- FastLDA(train.data = train_data,
                                test.data = as.data.frame(GSE72951.test.data),
                                response = "subtype",
                                cluster = GSE72951_common_lnc,
                                method = c("lda","qda","mda","fda","NaiveBayes")[1],
                                name = "love")
  }
  
  
}

## TCGA一致性聚类分析
{
  ## 导入TCGA表达谱数据和表型数据
  exprSet <- get(load("D:/Rsources/Project/SingleCelllncRNA/Data/TCGA-GBM/Transcriptome Profiling.Gene Expression Quantification.HTSeq - FPKM.rda"))
  gbm.clinic <- get(load("D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/gbm.clinic.rda"))
  exprSet <- assay(exprSet) 
  head(exprSet)
  table(substring(colnames(exprSet),14,15) < 10) # 5个正常样本，169个疾病样本
  exprSet <- exprSet[lnc_sig[lnc_sig %in% rownames(exprSet)], ]
  exprSet <- exprSet[,substring(colnames(exprSet),14,15) < 10]
  colnames(exprSet) <- substring( colnames(exprSet), 1, 12 )
  #提取表达谱样本中病人id和临床病人id overlap部分
  common_patients <- intersect( colnames(exprSet), gbm.clinic$Patient.ID ) 
  exprSet <- exprSet[ , colnames(exprSet) %in% common_patients ]
  exprSet <- exprSet[,!duplicated(colnames(exprSet))]
  dim(exprSet)
  
  library(ConsensusClusterPlus)
  TCGA.results <- ConsensusClusterPlus(as.matrix(exprSet),
                                           maxK=8,
                                           reps=500,
                                           pItem=0.8,
                                           pFeature=1,
                                           title = "TCGA",
                                           clusterAlg = "pam",
                                           distance = "pearson",
                                           seed = 1262118388.71279,
                                           writeTable = T,
                                           verbose = T,
                                           plot = "pdf")
  optimalCutoff <- function(maxK,
                            Consensus.result){
    # 用PAC的方法确定最佳聚类数,面积最小值对应K为最佳K
    Kvec <- 2 : maxK
    x1 <- 0.1;x2 <- 0.9        # threshold defining the intermediate sub-interval
    PAC <- rep(NA,length(Kvec)) 
    names(PAC) <- paste( "K=", Kvec, sep = "" )  # from 2 to maxK
    for(i in Kvec){
      M <- Consensus.result[[i]]$consensusMatrix
      Fn <- ecdf(M[lower.tri(M)])          # M 为计算出共识矩阵
      PAC[i-1] <- Fn(x2) - Fn(x1)
    } 
    optK <- Kvec[which.min(PAC)]  # 理想的K值
    return(optK)
  }
  optimalCutoff(8,TCGA.results)  # 2个亚型是最好的
  save( TCGA.results,file="D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA.ConsensusClusterResult.rda" )
 
}

## TCGA一致性聚类结果生存分析
{
  load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA.ConsensusClusterResult.rda")
  # 整理TCGA生存数据
  gbm.clinic <- get(load("D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/gbm.clinic.rda"))
  gbm.clinic <- gbm.clinic[!duplicated(gbm.clinic$Patient.ID),]
  gbm.clinic <- gbm.clinic[match(names(TCGA.results[[2]]$consensusClass),gbm.clinic$Patient.ID),]
  identical(gbm.clinic$Patient.ID,names(TCGA.results[[5]]$consensusClass))
  gbm.clinic$label2 <- TCGA.results[[2]]$consensusClass
  gbm.clinic$label3 <- TCGA.results[[3]]$consensusClass
  gbm.clinic$label4 <- TCGA.results[[4]]$consensusClass
  gbm.clinic$label5 <- TCGA.results[[5]]$consensusClass
  gbm.clinic$label6 <- TCGA.results[[6]]$consensusClass
  gbm.clinic$label7 <- TCGA.results[[7]]$consensusClass
  gbm.clinic$label8 <- TCGA.results[[8]]$consensusClass
  
  head(gbm.clinic)
  FastSurvKM(gbm.clinic,
             time = "Overall.Survival..Months.",
             status = "Overall.Survival.Status",
             marker = "label2",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_os")
  
  FastSurvKM(gbm.clinic,
             time = "Months.of.disease.specific.survival",
             status = "Disease.specific.Survival.status",
             marker = "label2",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event= "DSS",#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_dss")
  
  FastSurvKM(gbm.clinic,
             time = "Progress.Free.Survival..Months.",
             status = "Progression.Free.Status",
             marker = "label2",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[2],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_dfs")
  
  FastSurvKM(gbm.clinic,
             time = "Overall.Survival..Months.",
             status = "Overall.Survival.Status",
             marker = "label3",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_os")
  
  FastSurvKM(gbm.clinic,
             time = "Months.of.disease.specific.survival",
             status = "Disease.specific.Survival.status",
             marker = "label3",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event= "DSS",#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_dss")
  
  FastSurvKM(gbm.clinic,
             time = "Progress.Free.Survival..Months.",
             status = "Progression.Free.Status",
             marker = "label3",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[2],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_dfs")
  
  FastSurvKM(gbm.clinic,
             time = "Overall.Survival..Months.",
             status = "Overall.Survival.Status",
             marker = "label4",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_os")
  
  FastSurvKM(gbm.clinic,
             time = "Months.of.disease.specific.survival",
             status = "Disease.specific.Survival.status",
             marker = "label4",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event= "DSS",#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_dss")
  
  FastSurvKM(gbm.clinic,
             time = "Progress.Free.Survival..Months.",
             status = "Progression.Free.Status",
             marker = "label4",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[2],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_dfs")
  
  FastSurvKM(gbm.clinic,
             time = "Overall.Survival..Months.",
             status = "Overall.Survival.Status",
             marker = "label5",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_os")
  
  FastSurvKM(gbm.clinic,
             time = "Months.of.disease.specific.survival",
             status = "Disease.specific.Survival.status",
             marker = "label5",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event= "DSS",#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_dss")
  
  FastSurvKM(gbm.clinic,
             time = "Progress.Free.Survival..Months.",
             status = "Progression.Free.Status",
             marker = "label5",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[2],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_dfs")
  
  FastSurvKM(gbm.clinic,
             time = "Overall.Survival..Months.",
             status = "Overall.Survival.Status",
             marker = "label6",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_os")
  
  FastSurvKM(gbm.clinic,
             time = "Months.of.disease.specific.survival",
             status = "Disease.specific.Survival.status",
             marker = "label6",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event= "DSS",#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_dss")
  
  FastSurvKM(gbm.clinic,
             time = "Progress.Free.Survival..Months.",
             status = "Progression.Free.Status",
             marker = "label6",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[2],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_dfs")
  
  FastSurvKM(gbm.clinic,
             time = "Overall.Survival..Months.",
             status = "Overall.Survival.Status",
             marker = "label7",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_os")
  
  FastSurvKM(gbm.clinic,
             time = "Months.of.disease.specific.survival",
             status = "Disease.specific.Survival.status",
             marker = "label7",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event= "DSS",#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_dss")
  
  FastSurvKM(gbm.clinic,
             time = "Progress.Free.Survival..Months.",
             status = "Progression.Free.Status",
             marker = "label7",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[2],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_dfs")
  
  FastSurvKM(gbm.clinic,
             time = "Overall.Survival..Months.",
             status = "Overall.Survival.Status",
             marker = "label8",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[1],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_os")
  
  FastSurvKM(gbm.clinic,
             time = "Months.of.disease.specific.survival",
             status = "Disease.specific.Survival.status",
             marker = "label8",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event= "DSS",#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_dss")
  
  FastSurvKM(gbm.clinic,
             time = "Progress.Free.Survival..Months.",
             status = "Progression.Free.Status",
             marker = "label8",
             upper.time=NULL,
             xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"), #用于转换x轴坐标刻度
             unit.xlabel = c("year", "month", "week", "day")[2], #所用生存时间的单位,
             surv.median.line = c("none", "hv", "h", "v")[2], #是否画出中位生存时间，默认不给出
             risk.table = c(TRUE, FALSE)[1], #是否显示risk table
             pval = c(TRUE, FALSE)[1], #是否给出log-rank的p值
             conf.int = c(FALSE, TRUE)[1], #是否画出置信区间
             main="", #主标题名字 
             survival.event=c("Overall Survival","Progress Free Survival")[2],#事件类型
             inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA/TCGA_dfs")
}

## TCGA一致性聚类结果和已有分型结果比较
{
  load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/TCGA.ConsensusClusterResult.rda")
  gbm.clinic <- get(load("D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/gbm.clinic.rda"))
  gbm.clinic <- gbm.clinic[!duplicated(gbm.clinic$Patient.ID),]
  gbm.clinic <- gbm.clinic[match(names(TCGA.results[[2]]$consensusClass),gbm.clinic$Patient.ID),]
  ## 克罗恩卡帕指数评估聚类结果的一致性
  library(fmsb)
  kt <- Kappa.test(gbm.clinic$Original.Subtype, gbm.clinic$label5)
  
  ## 卡方检验评估亚型一致性
  gbm.clinic$Original.Subtype <- as.character(gbm.clinic$Original.Subtype)
  label8.p.val <- chisq.test(gbm.clinic[["Original.Subtype"]], gbm.clinic[["label8"]])$p.value
  label7.p.val <- chisq.test(gbm.clinic[["Original.Subtype"]], gbm.clinic[["label7"]])$p.value
  label6.p.val <- chisq.test(gbm.clinic[["Original.Subtype"]], gbm.clinic[["label6"]])$p.value
  label5.p.val <- chisq.test(gbm.clinic[["Original.Subtype"]], gbm.clinic[["label5"]])$p.value
  label4.p.val <- chisq.test(gbm.clinic[["Original.Subtype"]], gbm.clinic[["label4"]])$p.value
  label3.p.val <- chisq.test(gbm.clinic[["Original.Subtype"]], gbm.clinic[["label3"]])$p.value
  label2.p.val <- chisq.test(gbm.clinic[["Original.Subtype"]], gbm.clinic[["label2"]])$p.value
  ## fisher精确检验评估亚型一致性
  flabel8.p.val <- fisher.test(gbm.clinic[["Original.Subtype"]], gbm.clinic[["label8"]],simulate.p.value=TRUE)$p.value
  flabel7.p.val <- fisher.test(gbm.clinic[["Original.Subtype"]], gbm.clinic[["label7"]],simulate.p.value=TRUE)$p.value
  flabel6.p.val <- fisher.test(gbm.clinic[["Original.Subtype"]], gbm.clinic[["label6"]],simulate.p.value=TRUE)$p.value
  flabel5.p.val <- fisher.test(gbm.clinic[["Original.Subtype"]], gbm.clinic[["label5"]],simulate.p.value=TRUE)$p.value
  flabel4.p.val <- fisher.test(gbm.clinic[["Original.Subtype"]], gbm.clinic[["label4"]],simulate.p.value=TRUE)$p.value
  flabel3.p.val <- fisher.test(gbm.clinic[["Original.Subtype"]], gbm.clinic[["label3"]],simulate.p.value=TRUE)$p.value
  flabel2.p.val <- fisher.test(gbm.clinic[["Original.Subtype"]], gbm.clinic[["label2"]],simulate.p.value=TRUE)$p.value

  ## rand measure评估聚类结果一致性
  install.packages("fossil")
  library(fossil)
  gbm.clinic[["Original.Subtype"]] <- factor(gbm.clinic[["Original.Subtype"]],levels = c("Classical","G-CIMP","Mesenchymal","Neural","Proneural"))
  gbm.clinic <- gbm.clinic[!is.na(gbm.clinic[["Original.Subtype"]]),]
  library(plyr)
  subtypes <- mapvalues(gbm.clinic[["Original.Subtype"]], from = c("Classical", "G-CIMP","Mesenchymal","Neural","Proneural"), 
                        to = c(1,2,3,4,5))
  l2ri <- rand.index(as.numeric(subtypes) , gbm.clinic[["label2"]])
  l3ri <- rand.index(as.numeric(subtypes) , gbm.clinic[["label3"]])
  l4ri <- rand.index(as.numeric(subtypes) , gbm.clinic[["label4"]])
  l5ri <- rand.index(as.numeric(subtypes) , gbm.clinic[["label5"]])
  l6ri <- rand.index(as.numeric(subtypes) , gbm.clinic[["label6"]])
  l7ri <- rand.index(as.numeric(subtypes) , gbm.clinic[["label7"]])
  l8ri <- rand.index(as.numeric(subtypes) , gbm.clinic[["label8"]])
  save(kt,label8.p.val,label7.p.val,label6.p.val,label5.p.val,label4.p.val,label3.p.val,label2.p.val,flabel8.p.val,
       flabel7.p.val,flabel6.p.val,flabel5.p.val,flabel4.p.val,flabel3.p.val,flabel2.p.val,l2ri,l3ri,l4ri,l5ri,l6ri,
       l7ri,l8ri,file="D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/kappa.randmeasure.rda")
  load("D:/Rsources/Project/SingleCelllncRNA/Results/5.subtype/kappa.randmeasure.rda")
  
  
}




