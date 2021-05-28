#' @description TCGA 检测lncRNA signature是否倾向于driver lncRNA
#' @author shi jian

## 导入TCGA gbm lncRNA拷贝数segment数据以及signature lncRNA，得到拷贝数变异的signature lncRNA
{
  load("D:/Rsources/Project/SingleCelllncRNA/Data/TCGA-GBM/Copy Number Variation.Masked Copy Number Segment.rda")
  load("D:/Rsources/Project/SingleCelllncRNA/Data/AnnotationData/gtf.annot.rda")
  load("D:/Rsources/Project/SingleCelllncRNA/Results/1.singleCellSignature/lnc_sig.rda")
  cpn <- as.data.frame(data)
  head(cpn)
  cpn$Chromosome <- paste0("chr",cpn$Chromosome)
  cpn <- cpn[,c(7,2,3,4,5,6)]
  cpn$status <- ifelse(abs(cpn$Segment_Mean) > 1.5,ifelse(cpn$Segment_Mean > 1.5,"GAIN","LOSS"),"NO")
  cpn <- cpn[cpn$status != "NO",]
  dim(cpn)
  table(cpn$status)
  cpn$strand <- rep("*",length(cpn$Start))
  cpn <- cpn[,c(2,3,4,8,1,6,7)]
  colnames(cpn) <- c("chr","start","end","strand","sample","segment_mean","status")
  
  load("D:/Rsources/Project/SingleCelllncRNA/Data/TCGA-GBM/Transcriptome Profiling.Gene Expression Quantification.HTSeq - FPKM.rda")
  library(SummarizedExperiment)
  exprSet <- assay(data)
  exprSet[1:4,1:4]
  # 这里需要解析TCGA数据库的ID规律，来判断样本归类问题,去除正常样本
  group_list <- ifelse(as.numeric(substr(colnames(exprSet),14,15)) < 10,'tumor','normal')
  exprSet <- exprSet[ , group_list == "tumor"] ## 去除5个正常样本
  ## 拷贝数变异谱中只保留表达谱中有的样本
  cpn <- cpn[substring(cpn$sample,1,12) %in% substring(colnames(exprSet),1,12),]

  # 提取signature lncRNA区间信息
  #lncRNA <- rownames(exprSet)[convert(rownames(assay(data)),totype="gene_type") == "lncRNA"]
  lncRNA <- unique(lnc_sig)
  lncRNA <- lncRNA[lncRNA %in% rownames(exprSet)]
  lncRNA <- gtf.annotation[match(lncRNA,gtf.annotation$ENSEMBL),]
  lncRNA_region <- cbind.data.frame(lncRNA$chr,lncRNA$start,lncRNA$end,lncRNA$strand,lncRNA$ENSEMBL)
  colnames(lncRNA_region) <- c( "chr","start","end","strand","ensemble" )
  lncRNA_region$strand <- rep( "*", length(lncRNA_region$chr) )
  
  # signature lncRNA和拷贝数变异取overlap
  lncRNA_cnv_overlap <- RegionOverlapping.gr(lncRNA_region,cpn)
  lncRNA_cnv_overlap
  # seqnames     start       end  width strand               X space Qindex Sindex    Qstart      Qend    Sstart      Send
  # 1      chr7 117604791 117647415  42625      * ENSG00000083622  chr7      2   6428 117604791 117647415 116643247 117867411
  # 2     chr11   2140501   2148666   8166      * ENSG00000099869 chr11      4   5406   2140501   2148666   2134983   2202130
  # 3     chr12  23181334  23251499  70166      * ENSG00000115934 chr12      6   3344  23181334  23251499  23235706  23235751
  # 4     chr20  10996293  11029455  33163      * ENSG00000125899 chr20     16  17346  10996293  11029455  10531589  12457656
  # OLstart     OLend OLlength      OLpercQ     OLpercS    OLtype                       sample segment_mean status
  # 1  117604791 117647415    42625 100.00000000   3.4819653    inside TCGA-06-0413-01A-01D-0275-01       3.0594   GAIN
  # 2    2140501   2148666     8166 100.00000000  12.1611962    inside TCGA-06-0129-01A-01G-0289-01      -3.0546   LOSS
  # 3   23235706  23235751       46   0.06555882 100.0000000 contained TCGA-12-1597-01B-01D-0911-01      -3.9716   LOSS
  # 4   10996293  11029455    33163 100.00000000   1.7217980    inside TCGA-28-5214-01A-01D-1479-01       2.4267   GAIN
  
  ## 基于lncRNA扩增缺失将样本分开
  lncRNA_cnv_overlap$sample <- substring(lncRNA_cnv_overlap$sample,1,12)
  lncRNA_gaincnv_overlap <- lncRNA_cnv_overlap[lncRNA_cnv_overlap$status == "GAIN",]
  lncRNA_losscnv_overlap <- lncRNA_cnv_overlap[lncRNA_cnv_overlap$status == "LOSS",]
  
  ## 基于lncRNA ENSEMBLE名字将列表分开,生成每个lncRNA的list
  lncRNA_gaincnv_overlap_list <- split(lncRNA_gaincnv_overlap,lncRNA_gaincnv_overlap$X) 
  lncRNA_losscnv_overlap_list <- split(lncRNA_losscnv_overlap,lncRNA_losscnv_overlap$X) 
  save(lncRNA_gaincnv_overlap_list,lncRNA_losscnv_overlap_list,
       file="D:/Rsources/Project/SingleCelllncRNA/Results/4.1DriverLncRNA/lnc_cnv_overlap.rda")
  ## 剔除只有一个样本的列表元素(因为要做wilcox rank-sum test)，分别生成扩增和缺失的wilcox rank-sum检验的数据
  # lgol <- NULL
  # for(i in 1:length(lncRNA_gaincnv_overlap_list)){
  #   if(length(unique(lncRNA_gaincnv_overlap_list[[i]]$sample)) == 1){
  #     
  #   }else{
  #     lgol <- c(lgol,lncRNA_gaincnv_overlap_list[i])
  #     names(lgol[i]) <- names(lncRNA_gaincnv_overlap_list)[i]
  #   }
  # }
  # lgol <- lgol[-which(sapply(lgol, is.null))]
  # length(lgol) #115个基因
  # 
  # llol <- NULL
  # for(i in 1:length(lncRNA_losscnv_overlap_list)){
  #   if(length(unique(lncRNA_losscnv_overlap_list[[i]]$sample)) == 1){
  #     
  #   }else{
  #     llol <- c(llol,lncRNA_losscnv_overlap_list[i])
  #     names(llol[i]) <- names(lncRNA_losscnv_overlap_list)[i]
  #   }
  # }
  # llol <- llol[-which(sapply(llol, is.null))]
  # length(llol) #143个基因
  
  ## 对扩增lncRNA进行差异分析，得到上调的lncRNA
  colnames(exprSet) <- substring(colnames(exprSet),1,12)
  exprSet <- exprSet[,colnames(exprSet) %in% substring(cpn$sample,1,12)]
  gr <- NULL
  gain_sample <- NULL
  for(i in 1:length(lncRNA_gaincnv_overlap_list)){
    #p.value <- wilcox.test(temp[,"exp"]~temp[,"group_list"],alternative = "greater")$p.value
    if(length(unique(lncRNA_gaincnv_overlap_list[[i]]$sample)) == 1){
      p.value <- wilcox.test(exprSet[unique(lncRNA_gaincnv_overlap_list[[i]]$X),setdiff(colnames(exprSet),lncRNA_gaincnv_overlap_list[[i]]$sample)] , exprSet[unique(lncRNA_gaincnv_overlap_list[[i]]$X),unique(lncRNA_gaincnv_overlap_list[[i]]$sample)],alternative = "less")$p.value
    }else{
      group_list <- ifelse(colnames(exprSet) %in% lncRNA_gaincnv_overlap_list[[i]]$sample,"case","control")
      group_list <- factor(group_list,levels = c("case","control"))
      temp <- cbind.data.frame(exprSet[unique(lncRNA_gaincnv_overlap_list[[i]]$X),],group_list)
      colnames(temp) <- c("exp","group_list")
      p.value <- FastStatisticsTest(data = temp, variable = "exp", by = "group_list",
                                    test = "wilcox.test", id = NULL, type = "continuous", alternative="greater")
    }
    gr <- c(gr,p.value)
    if(p.value < 0.05){
      gain_sample[[i]] <- unique(lncRNA_gaincnv_overlap_list[[i]]$sample)
    }
  }
  length(lncRNA_gaincnv_overlap_list[gr < 0.05])
  gain_driver_lncRNA <- names(lncRNA_gaincnv_overlap_list[gr < 0.05])
  gain_sample <- gain_sample[-which(sapply(gain_sample, is.null))]
  names(gain_sample) <- gain_driver_lncRNA
  save(gain_driver_lncRNA,gain_sample,file="D:/Rsources/Project/SingleCelllncRNA/Results/4.1DriverLncRNA/gain_driver_lncRNA.rda")
  
  ## 对缺失lncRNA进行差异分析,得到下调的lncRNA  
  lr <- NULL
  loss_sample <- NULL
  for(i in 1:length(lncRNA_losscnv_overlap_list)){
    if(length(unique(lncRNA_losscnv_overlap_list[[i]]$sample)) == 1){
      p.value <- wilcox.test(exprSet[unique(lncRNA_losscnv_overlap_list[[i]]$X),setdiff(colnames(exprSet),lncRNA_losscnv_overlap_list[[i]]$sample)] , exprSet[unique(lncRNA_losscnv_overlap_list[[i]]$X),unique(lncRNA_losscnv_overlap_list[[i]]$sample)],
                             alternative = "greater")$p.value
    }else{
      group_list <- ifelse(colnames(exprSet) %in% lncRNA_losscnv_overlap_list[[i]]$sample,"case","control")
      group_list <- factor(group_list,levels = c("case","control")) 
      temp <- cbind.data.frame(exprSet[unique(lncRNA_losscnv_overlap_list[[1]]$X),],group_list)
      colnames(temp) <- c("exp","group_list")
      #p.value <- wilcox.test(temp[,"exp"]~temp[,"group_list"],alternative = "less")$p.value
      p.value <- FastStatisticsTest(data = temp, variable = "exp", by = "group_list",
                                    test = "wilcox.test", id = NULL, type = "continuous", alternative="less")
    }
    lr <- c(lr,p.value)
    if(p.value < 0.05){
      loss_sample[[i]] <- unique(lncRNA_losscnv_overlap_list[[i]]$sample)
    }
  }
  length(lncRNA_losscnv_overlap_list[lr < 0.05])
  loss_driver_lncRNA <- names(lncRNA_losscnv_overlap_list[lr < 0.05])
  loss_sample <- loss_sample[-which(sapply(loss_sample, is.null))]
  names(loss_sample) <- loss_driver_lncRNA
  save(loss_driver_lncRNA,loss_sample,file="D:/Rsources/Project/SingleCelllncRNA/Results/4.1DriverLncRNA/loss_driver_lncRNA.rda")
}

## Driver lncRNA 差异表达boxplot
{
  load("D:/Rsources/Project/SingleCelllncRNA/Results/4.1DriverLncRNA/loss_driver_lncRNA.rda")
  load("D:/Rsources/Project/SingleCelllncRNA/Results/4.1DriverLncRNA/gain_driver_lncRNA.rda")
  exprSet <- get(load("D:/Rsources/Project/SingleCelllncRNA/Data/TCGA-GBM/Transcriptome Profiling.Gene Expression Quantification.HTSeq - FPKM.rda"))
  cpn <- get(load("D:/Rsources/Project/SingleCelllncRNA/Data/TCGA-GBM/Copy Number Variation.Masked Copy Number Segment.rda"))
  cpn <- as.data.frame(cpn)
  library(SummarizedExperiment)
  exprSet <- assay(exprSet)
  exprSet[1:4,1:4]
  # 这里需要解析TCGA数据库的ID规律，来判断样本归类问题,去除正常样本
  group_list <- ifelse(as.numeric(substr(colnames(exprSet),14,15)) < 10,'tumor','normal')
  exprSet <- exprSet[, group_list == "tumor"] ## 去除5个正常样本
  # 修改表达谱列名和拷贝数变异谱样本名字
  colnames(exprSet) <- substring(colnames(exprSet),1,12)
  cpn$Sample <- substring(cpn$Sample,1,12) 
  
  # cpn和exprSet取样本交集
  common.patients <- intersect(colnames(exprSet),cpn$Sample)
  exprSet <- exprSet[,common.patients]
  gain_sample
  exp <- exprSet["ENSG00000226950",]
  label <- ifelse(colnames(exprSet) %in% c("TCGA-12-1597","TCGA-06-0221","TCGA-14-0817","TCGA-19-1390"),"Amp","w/o Amp")
  gain.data <- cbind.data.frame(label=label,exp=exp)
  exp1 <- exprSet["ENSG00000214783",]
  label1 <- ifelse(colnames(exprSet) %in% c("TCGA-32-1970","TCGA-32-2615"),"Amp","w/o Amp")
  gain.data1 <- cbind.data.frame(label=label1,exp=exp1)
  pdf("D:/Rsources/Project/SingleCelllncRNA/Results/4.1DriverLncRNA/gain.expression.pdf")
  p <- Fastboxplot(data = gain.data,
                    x = "label",
                    y="exp",
                    color="label",
                    facet=NULL,
                    palette="jco",
                    x.lab = "" ,y.lab = "FPKM",title="ENSG00000226950",
                    legend.position="none",
                    x.test.angle = 0,
                    add.jitter=T)
  p1 <- Fastboxplot(data = gain.data1,
                   x = "label",
                   y="exp",
                   color="label",
                   facet=NULL,
                   palette="jco",
                   x.lab = "" ,y.lab = "FPKM",title="ENSG00000214783",
                   legend.position="none",
                   x.test.angle = 0,
                   add.jitter=T)
  print(p)
  print(p1)
  dev.off()
  
  
  label <- ifelse(colnames(exprSet) %in% c("TCGA-14-0871"),"Loss","w/o Loss")
  exp <- exprSet["ENSG00000178722",]
  loss.data <- cbind.data.frame(label=label,exp=exp)
  pdf("D:/Rsources/Project/SingleCelllncRNA/Results/4.1DriverLncRNA/loss.expression.pdf")
  p <- Fastboxplot(data = loss.data,
                   x = "label",
                   y="exp",
                   fill="label",
                   facet=NULL,
                   palette=NULL,
                   x.lab = "" ,y.lab = "FPKM",title="ENSG00000178722",
                   legend.position="none",
                   x.test.angle = 0,
                   add.jitter=T)
  print(p)
  dev.off()
  
}












