#' @description  评估driver lncRNA的预后以及临床诊断效能
#' @author shi jian
#' 
#' 

## 加载TCGA gbm表达谱以及临床数据
load("D:/Rsources/Project/SingleCelllncRNA/Data/TCGA-GBM/Transcriptome Profiling.Gene Expression Quantification.HTSeq - FPKM.rda")
gbm.clinic <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/TCGA-GBM/gbm_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
load("D:/Rsources/Project/SingleCelllncRNA/Data/TCGA-GBM/GBM.subtype.rda")
head(GBM.subtype)
head(gbm.clinic)

## 对临床数据进行清洗和整合
{
  load("D:/Rsources/Project/SingleCelllncRNA/Results/1.singleCellSignature/lnc_sig.rda")
  head(gbm.clinic)
  head(GBM.subtype)
  table(gbm.clinic$Diagnosis.Age,useNA="ifany")
  table(gbm.clinic$Race.Category,useNA="ifany")
  table(gbm.clinic$Subtype,useNA="ifany")
  table(gbm.clinic$Disease.specific.Survival.status,useNA="ifany")
  table(gbm.clinic$Overall.Survival.Status,useNA="ifany")
  table(gbm.clinic$Progression.Free.Status,useNA="ifany")
  table(gbm.clinic$Sex,useNA="ifany")
  table(gbm.clinic$Race.Category,useNA="ifany")
  gbm.clinic$Disease.specific.Survival.status[gbm.clinic$Disease.specific.Survival.status == "0:ALIVE OR DEAD TUMOR FREE"] = 0
  gbm.clinic$Disease.specific.Survival.status[gbm.clinic$Disease.specific.Survival.status == "1:DEAD WITH TUMOR"] = 1
  gbm.clinic$Overall.Survival.Status[gbm.clinic$Overall.Survival.Status == "0:LIVING"] = 0
  gbm.clinic$Overall.Survival.Status[gbm.clinic$Overall.Survival.Status == "1:DECEASED"] = 1
  gbm.clinic$Progression.Free.Status[gbm.clinic$Progression.Free.Status == "0:CENSORED"] = 0
  gbm.clinic$Progression.Free.Status[gbm.clinic$Progression.Free.Status == "1:PROGRESSION"] = 1
  gbm.clinic$Race.Category[gbm.clinic$Race.Category == "Black or African American"] = "Black"
  gbm.subtype <- GBM.subtype[,c("patient","Age..years.at.diagnosis.","Gender","ABSOLUTE.purity","Original.Subtype")]
  gbm.clinic <- merge(gbm.clinic,gbm.subtype,by.x="Patient.ID",by.y="patient",all.x=TRUE)
  save(gbm.clinic, file = "D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/gbm.clinic.rda")

}

## 对signature lncRNA进行单多因素cox分析
{
  ##生成cox分析所需的table
  load("D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/gbm.clinic.rda")
  load("D:/Rsources/Project/SingleCelllncRNA/Results/1.singleCellSignature/lnc_sig.rda")
  library(SummarizedExperiment)
  exprSet <- assay(data)
  colnames(exprSet) <- substring(colnames(exprSet),1,12)
  exprSet[1:4,1:4]
  exprSet <- exprSet[unique(lnc_sig)[unique(lnc_sig) %in% rownames(exprSet)],  ]
  common_patient <- intersect(gbm.clinic$Patient.ID,colnames(exprSet))
  exprSet <- exprSet[,common_patient]
  gbm.clinic <- gbm.clinic[match(common_patient,gbm.clinic$Patient.ID),]
  identical(colnames(exprSet),gbm.clinic$Patient.ID)
  gbm.clinic.exp <- cbind.data.frame(gbm.clinic,t(exprSet))
  gbm.clinic.exp$Overall.Survival.Status <- as.numeric(gbm.clinic.exp$Overall.Survival.Status)
  gbm.clinic.exp$Progression.Free.Status <- as.numeric(gbm.clinic.exp$Progression.Free.Status)
  gbm.clinic.exp$Disease.specific.Survival.status <- as.numeric(gbm.clinic.exp$Disease.specific.Survival.status)
  save(gbm.clinic.exp,file="D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/gbm.clinic.exp.rda")
  
  os_stepwise_cox_result<- getCoxTable(data=gbm.clinic.exp,
                                          time.col = "Overall.Survival..Months.",
                                          status.col = "Overall.Survival.Status",
                                          stepwise=T,direction="both",
                                          cluster = c("Age..years.at.diagnosis.","Gender","Original.Subtype",colnames(t(exprSet))),
                                          control = c(NA,"female","Proneural",rep(NA,length(colnames(t(exprSet))))),
                                          dig=4,
                                          names = "os_stepwise")
  
  os_cox_result<- getCoxTable(data=gbm.clinic.exp,
                               time.col = "Overall.Survival..Months.",
                               status.col = "Overall.Survival.Status",
                               stepwise=F,direction="both",
                               cluster = c("Age..years.at.diagnosis.","Gender","Original.Subtype",colnames(t(exprSet))),
                               control = c(NA,"female","Proneural",rep(NA,length(colnames(t(exprSet))))),
                               dig=4,
                               names = "os")
  
  pfs_stepwise_cox_result<- getCoxTable(data=gbm.clinic.exp,
                                         time.col = "Progress.Free.Survival..Months.",
                                         status.col = "Progression.Free.Status",
                                         stepwise = T,direction="both",
                                         cluster = c("Age..years.at.diagnosis.","Gender","Original.Subtype",colnames(t(exprSet))),
                                         control = c(NA,"female","Proneural",rep(NA,length(colnames(t(exprSet))))),
                                         dig=4,
                                         names = "pfs_stepwise")
  
  
  
  pfs_cox_result<- getCoxTable(data=gbm.clinic.exp,
                                time.col = "Progress.Free.Survival..Months.",
                                status.col = "Progression.Free.Status",
                                stepwise=F,direction="both",
                                cluster = c("Age..years.at.diagnosis.","Gender","Original.Subtype",colnames(t(exprSet))),
                                control = c(NA,"female","Proneural",rep(NA,length(colnames(t(exprSet))))),
                                dig=4,
                                names = "pfs")
  
  dss_stepwise_cox_result <- getCoxTable(data=gbm.clinic.exp,
                                          time.col = "Months.of.disease.specific.survival",
                                          status.col = "Disease.specific.Survival.status",
                                          stepwise=T,direction="both",
                                          cluster = c("Age..years.at.diagnosis.","Gender","Original.Subtype",colnames(t(exprSet))),
                                          control = c(NA,"female","Proneural",rep(NA,length(colnames(t(exprSet))))),
                                          dig=4,
                                          names = "dss_stepwise")
  
  dss_cox_result <- getCoxTable(data=gbm.clinic.exp,
                                 time.col = "Months.of.disease.specific.survival",
                                 status.col = "Disease.specific.Survival.status",
                                 stepwise=F,direction="both",
                                 cluster = c("Age..years.at.diagnosis.","Gender","Original.Subtype",colnames(t(exprSet))),
                                 control = c(NA,"female","Proneural",rep(NA,length(colnames(t(exprSet))))),
                                 dig=4,
                                 names = "dss")
} 

## 对每个driver lncRNA基于表达中位值分组做KM曲线分析
{
  head(gbm.clinic.exp)
  dim(gbm.clinic.exp)
  lab_matrix <- apply(gbm.clinic.exp[,-c(1:15)],2,function(x){
    med <- median(x)
    lab <- ifelse(x > med,"high","low")
    lab <- factor(lab,levels=c("low","high"))
    lab
  })
  gbm.clinic.km <- cbind.data.frame(gbm.clinic.exp[,1:15],lab_matrix)
  
  for(i in 16:63){
    if(length(unique(gbm.clinic.km[,i])) == 2){
      FastSurvKM(clinical.data = gbm.clinic.km,
                 time = "Overall.Survival..Months.",
                 status = "Overall.Survival.Status",
                 marker = colnames(gbm.clinic.km)[i],
                 upper.time=NULL,
                 xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"),#用于转换x轴坐标刻度
                 unit.xlabel = "month",#所用生存时间的单位,
                 surv.median.line = "hv",#是否画出中位生存时间，默认不给出
                 risk.table = c(TRUE, FALSE)[1],#是否显示risk table
                 pval = c(TRUE, FALSE)[1],#是否给出log-rank的p值
                 conf.int = c(FALSE, TRUE)[1],#是否画出置信区间
                 main = colnames(gbm.clinic.km)[i],#主标题名字 
                 survival.event = "Overall Survival",#事件类型
                 inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/TCGA_os")
      
      FastSurvKM(clinical.data = gbm.clinic.km,
                 time = "Progress.Free.Survival..Months.",
                 status = "Progression.Free.Status",
                 marker = colnames(gbm.clinic.km)[i],
                 upper.time=NULL,
                 xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"),#用于转换x轴坐标刻度
                 unit.xlabel = "month",#所用生存时间的单位,
                 surv.median.line = "hv",#是否画出中位生存时间，默认不给出
                 risk.table = c(TRUE, FALSE)[1],#是否显示risk table
                 pval = c(TRUE, FALSE)[1],#是否给出log-rank的p值
                 conf.int = c(FALSE, TRUE)[1],#是否画出置信区间
                 main = colnames(gbm.clinic.km)[i],#主标题名字 
                 survival.event = "Progress Free Survival",#事件类型
                 inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/TCGA_dfs")
      
      FastSurvKM(clinical.data = gbm.clinic.km,
                 time = "Months.of.disease.specific.survival",
                 status = "Disease.specific.Survival.status",
                 marker = colnames(gbm.clinic.km)[i],
                 upper.time=NULL,
                 xscale = c("d_m", "d_y", "m_d", "m_y", "y_d", "y_m"),#用于转换x轴坐标刻度
                 unit.xlabel = "month",#所用生存时间的单位,
                 surv.median.line = "hv",#是否画出中位生存时间，默认不给出
                 risk.table = c(TRUE, FALSE)[1],#是否显示risk table
                 pval = c(TRUE, FALSE)[1],#是否给出log-rank的p值
                 conf.int = c(FALSE, TRUE)[1],#是否画出置信区间
                 main = colnames(gbm.clinic.km)[i],#主标题名字 
                 survival.event = "DSS",#事件类型
                 inputFilePath = "D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/TCGA_dss")
    }
    
    
  }

}

## 评估signature lncRNA的诊断效能（构建线性判别模型，支持向量机模型以及logistic回归模型）
{
  ## 构建模型训练数据 train_data
  load("D:/Rsources/Project/SingleCelllncRNA/Results/1.singleCellSignature/lnc_sig.rda")
  load("D:/Rsources/Project/SingleCelllncRNA/Data/TCGA-GBM/Transcriptome Profiling.Gene Expression Quantification.HTSeq - FPKM.rda")
  ## load("D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/gbm.clinic.rda")
  library(SummarizedExperiment)
  exprSet <- assay(data)
  table(substring(colnames(exprSet),14,15) < 10) # 5个正常样本，169个疾病样本
  exprSet <- exprSet[unique(lnc_sig)[unique(lnc_sig) %in% rownames(exprSet)],]
  label <- ifelse(substring(colnames(exprSet),14,15) < 10,"gbm","normal")
  train_data <- cbind.data.frame(t(exprSet),label)
  dim(train_data) #174 124
  train_data <- train_data[,colSums(train_data[,-49]) > 0]
  
  ## cfs特征选择
  # 预后显著的lncRNA作为诊断marker
  lnc_sig1 <- unique(c("ENSG00000176124","ENSG00000231889","ENSG00000232677","ENSG00000245694","ENSG00000224032","ENSG00000226950","ENSG00000231419","ENSG00000231889","ENSG00000232677","ENSG00000245694","ENSG00000176124","ENSG00000231889","ENSG00000232677","ENSG00000245694"))
  exprSet <- exprSet[lnc_sig1,]
  train_data <- cbind.data.frame(t(exprSet),label)
  dim(train_data) #174 124
  train_data <- train_data[,colSums(train_data[,-8]) > 0]
  library(Biocomb)
  train_data[,ncol(train_data)] <- as.factor(train_data[,ncol(train_data)])
  out <- select.cfs(matrix=train_data)
  train_data <- train_data[,-out$Index]
  out$Biomarker # "ENSG00000176124"
  lnc_sig1 <- setdiff(unique(lnc_sig1)[unique(lnc_sig1) %in% rownames(exprSet)],out$Biomarker)
  
  ## 构建检测数据集GSE7696，GSE108474，GSE16011
  {
    GSE7696_GPL570_exprSet <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE7696/GSE7696-GPL570-expression.MATRIX.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
    GSE7696_GPL570_pData <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE7696/GSE7696-GPL570-pData.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
    GSE7696_GPL570_probeAnnotation <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE7696/GSE7696-GPL570-ProbeAnnotation.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
    
    GSE108474_GPL570_exprSet <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE108474/GSE108474-GPL570-expression.MATRIX.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
    GSE108474_GPL570_pData <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE108474/GSE108474-GPL570-pData.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
    GSE108474_GPL570_probeAnnotation <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE108474/GSE108474-GPL570-ProbeAnnotation.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
    
    GSE16011_GPL8542_exprSet <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE16011/GSE16011-GPL8542-expression.MATRIX.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
    GSE16011_GPL8542_pData <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE16011/GSE16011-GPL8542-pData.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
    GSE16011_GPL8542_probeAnnotation <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GSE16011/GSE16011-GPL8542-ProbeAnnotation.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
  }
  
  ## 表型和表达数据样本顺序一致
  {
    identical(GSE7696_GPL570_pData$geo_accession,colnames(GSE7696_GPL570_exprSet))
    table(GSE7696_GPL570_pData$characteristics_ch1.1)
    GSE7696_GPL570_pData$characteristics_ch1.1[grep("GBM",GSE7696_GPL570_pData$characteristics_ch1.1)] <- "gbm"
    GSE7696_GPL570_pData$characteristics_ch1.1[grep("non-tumoral",GSE7696_GPL570_pData$characteristics_ch1.1)] <- "normal" #四个正常样本
    GSE7696_label <- GSE7696_GPL570_pData$characteristics_ch1.1
    table(GSE7696_label)
    identical(GSE108474_GPL570_pData$geo_accession,colnames(GSE108474_GPL570_exprSet))
    table(GSE108474_GPL570_pData$characteristics_ch1.1)
    GSE108474_GPL570_pData$characteristics_ch1.1[grep("glioblastoma multiforme",GSE108474_GPL570_pData$characteristics_ch1.1)] <- "gbm"
    GSE108474_GPL570_pData$characteristics_ch1.1[grep("normal",GSE108474_GPL570_pData$characteristics_ch1.1)] <- "normal" #四个正常样本
    GSE108474_GPL570_exprSet <- GSE108474_GPL570_exprSet[,which(GSE108474_GPL570_pData$characteristics_ch1.1 == "gbm" | GSE108474_GPL570_pData$characteristics_ch1.1 == "normal")]
    GSE108474_GPL570_exprSet <- GSE108474_GPL570_exprSet[,setdiff(colnames(GSE108474_GPL570_exprSet),c("GSM2899720","GSM2899721","GSM2899722","GSM2899723","GSM2899724","GSM2899725","GSM2899726"))]
    GSE108474_GPL570_pData <- GSE108474_GPL570_pData[which(GSE108474_GPL570_pData$characteristics_ch1.1 == "gbm" | GSE108474_GPL570_pData$characteristics_ch1.1 == "normal"),]
    GSE108474_GPL570_pData <- GSE108474_GPL570_pData[GSE108474_GPL570_pData$geo_accession %in% setdiff(colnames(GSE108474_GPL570_exprSet),c("GSM2899720","GSM2899721","GSM2899722","GSM2899723","GSM2899724","GSM2899725","GSM2899726")),]
    GSE108474_label <- GSE108474_GPL570_pData$characteristics_ch1.1
    table(GSE108474_label)
    identical(GSE16011_GPL8542_pData$geo_accession,colnames(GSE16011_GPL8542_exprSet))
    table(GSE16011_GPL8542_pData$source_name_ch1)
    GSE16011_GPL8542_pData$source_name_ch1[grep("glioma",GSE16011_GPL8542_pData$source_name_ch1)] <- "gbm"
    GSE16011_GPL8542_pData$source_name_ch1[grep("control",GSE16011_GPL8542_pData$source_name_ch1)] <- "normal" #四个正常样本
    GSE16011_label <- GSE16011_GPL8542_pData$source_name_ch1
    head(GSE16011_GPL8542_probeAnnotation)
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
    GSE7696_common_lnc <- intersect(rownames(GSE7696.anno.res$exprSet),lnc_sig1 )
    GSE7696.test.data <- cbind.data.frame(t(GSE7696.anno.res$exprSet[GSE7696_common_lnc,]),GSE7696_label)
    colnames(GSE7696.test.data)[ncol(GSE7696.test.data)] <- "label"
    
    ## GSE108474表达谱探针注释
    GSE108474.anno.res <- probeAnnotation(exprSet = GSE108474_GPL570_exprSet,
                                          probeAnnotation = GSE108474_GPL570_probeAnnotation,
                                          id = "ID",
                                          symbol = "Gene.Symbol",
                                          entrezid = "ENTREZ_GENE_ID",
                                          index = "ENTREZ_GENE_ID",
                                          fromtype = "ENTREZID",
                                          totype = c("entrez","symbol","ensemble")[3])
    GSE108474_common_lnc <- intersect(rownames(GSE108474.anno.res$exprSet),lnc_sig1 )
    GSE108474.test.data <- cbind.data.frame(t(GSE108474.anno.res$exprSet[GSE108474_common_lnc,]),GSE108474_label)
    colnames(GSE108474.test.data)[ncol(GSE108474.test.data)] <- "label"
    ## GSE16011表达谱探针注释
    GSE16011.anno.res <- probeAnnotation(exprSet = GSE16011_GPL8542_exprSet,
                                         probeAnnotation = GSE16011_GPL8542_probeAnnotation,
                                         id = "ID",
                                         symbol = NULL,
                                         entrezid = "ORF",
                                         index = "ORF",
                                         fromtype = "ENTREZID",
                                         totype = c("entrez","symbol","ensemble")[3])
    GSE16011_common_lnc <- intersect(rownames(GSE16011.anno.res$exprSet),lnc_sig1 )
    GSE16011.test.data <- cbind.data.frame(t(GSE16011.anno.res$exprSet[GSE16011_common_lnc,]),GSE16011_label)
    colnames(GSE16011.test.data)[ncol(GSE16011.test.data)] <- "label"
    ## 11个lncRNA
  }

  ## 构建分类器模型以及不同检测数据集中进行ROC曲线评估
  {
    # 线性SVM模型
    GSE7696.svmLinear.res <- FastSVM(train.data = train_data,
                                     test.data = GSE7696.test.data,
                                     response = "label",
                                     cluster = GSE7696_common_lnc,
                                     seed = 123,
                                     method = c("svmLinear","svmRadial","svmPoly","nb","rf")[1],
                                     name = "GSE7696")
    # 朴素贝叶斯模型
    GSE7696.nb.res <- FastSVM(train.data = train_data,
                              test.data = GSE7696.test.data,
                              response = "label",
                              cluster = GSE7696_common_lnc,
                              seed = 123,
                              method = "nb",
                              name = "GSE7696")
    
    # 随机森林模型
    GSE7696.rf.res <- FastSVM(train.data = train_data,
                      test.data = GSE7696.test.data,
                      response = "label",
                      cluster = GSE7696_common_lnc,
                      seed = 123,
                      method = "rf",
                      name = "GSE7696")
    # 支持向量机poly模型
    GSE7696.svmPoly.res <- FastSVM(train.data = train_data,
                           test.data = GSE7696.test.data,
                           response = "label",
                           cluster = GSE7696_common_lnc,
                           seed = 123,
                           method = "svmPoly",
                           name = "GSE7696")
    # 支持向量机Radial模型
    GSE7696.svmRadial.res <- FastSVM(train.data = train_data,
                                     test.data = GSE7696.test.data,
                                     response = "label",
                                     cluster = GSE7696_common_lnc,
                                     seed = 123,
                                     method = "svmRadial",
                                     name = "GSE7696")
    # 线性判别模型
    GSE7696.lda.res <- FastLDA(train.data = train_data,
                       test.data = GSE7696.test.data,
                       response = "label",
                       cluster = GSE7696_common_lnc,
                       method = c("lda","qda","mda","fda","NaiveBayes")[1],
                       name = "GSE7696")
    
    # lasso logistic回归模型
    train_data$label1 <- ifelse(train_data$label == "gbm",1,0)
    GSE7696.test.data$label1 <- ifelse(GSE7696.test.data$label == "gbm",1,0)
    GSE7696.log.res <- FastGlmnet(data = train_data,
                          test.data = GSE7696.test.data,
                          response = "label1",
                          cluster = GSE7696_common_lnc,
                          col = "label",
                          control = "normal",
                          alpha = c(1,0)[1],
                          k = 10,
                          family = c("binomial","cox")[1],
                          optimize.method = c("min","lse")[1],
                          seed = 123,
                          plot.label = T,
                          plot.lwd = 2,
                          line.lty=3,
                          line.lwd=2,
                          line.col="blue",
                          verbose = F,
                          name = "GSE7696")
    save(GSE7696.svmLinear.res,
         GSE7696.nb.res,
         GSE7696.rf.res,
         GSE7696.svmPoly.res,
         GSE7696.svmRadial.res,
         GSE7696.lda.res,
         GSE7696.log.res,
         file="D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/GSE7696.classification.model.rda")
    
    
    # 线性SVM模型
    GSE108474.svmLinear.res <- FastSVM(train.data = train_data,
                                       test.data = GSE108474.test.data,
                                       response = "label",
                                       cluster = GSE108474_common_lnc,
                                       seed = 123,
                                       method = c("svmLinear","svmRadial","svmPoly","nb","rf")[1],
                                       name = "GSE108474")
    # 朴素贝叶斯模型
    GSE108474.nb.res <- FastSVM(train.data = train_data,
                                test.data = GSE108474.test.data,
                                response = "label",
                                cluster = GSE108474_common_lnc,
                                seed = 123,
                                method = "nb",
                                name = "GSE108474")
    
    # 随机森林模型
    GSE108474.rf.res <- FastSVM(train.data = train_data,
                                test.data = GSE108474.test.data,
                                response = "label",
                                cluster = GSE108474_common_lnc,
                                seed = 123,
                                method = "rf",
                                name = "GSE108474")
    # 支持向量机poly模型
    GSE108474.svmPoly.res <- FastSVM(train.data = train_data,
                                     test.data = GSE108474.test.data,
                                     response = "label",
                                     cluster = GSE108474_common_lnc,
                                     seed = 123,
                                     method = "svmPoly",
                                     name = "GSE108474")
    # 支持向量机Radial模型
    GSE108474.svmRadial.res <- FastSVM(train.data = train_data,
                                       test.data = GSE108474.test.data,
                                       response = "label",
                                       cluster = GSE108474_common_lnc,
                                       seed = 123,
                                       method = "svmRadial",
                                       name = "GSE108474")
    # 线性判别模型
    GSE108474.lda.res <- FastLDA(train.data = train_data,
                                 test.data = GSE108474.test.data,
                                 response = "label",
                                 cluster = GSE108474_common_lnc,
                                 method = c("lda","qda","mda","fda","NaiveBayes")[1],
                                 name = "GSE108474")
    
    # lasso logistic回归模型
    train_data$label1 <- ifelse(train_data$label == "gbm",1,0)
    GSE108474.test.data$label1 <- ifelse(GSE108474.test.data$label == "gbm",1,0)
    GSE108474.log.res <- FastGlmnet(data = train_data,
                                    test.data = GSE108474.test.data,
                                    response = "label1",
                                    cluster = GSE108474_common_lnc,
                                    col = "label",
                                    control = "normal",
                                    alpha = c(1,0)[1],
                                    k = 10,
                                    family = c("binomial","cox")[1],
                                    optimize.method = c("min","lse")[1],
                                    seed = 123,
                                    plot.label = T,
                                    plot.lwd = 2,
                                    line.lty=3,
                                    line.lwd=2,
                                    line.col="blue",
                                    verbose = F,
                                    name = "GSE108474")
    
    save(GSE108474.svmLinear.res,
         GSE108474.nb.res,
         GSE108474.rf.res,
         GSE108474.svmPoly.res,
         GSE108474.svmRadial.res,
         GSE108474.lda.res,
         GSE108474.log.res,
         file="D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/GSE108474.classification.model.rda")
    
    # 线性SVM模型
    GSE16011.svmLinear.res <- FastSVM(train.data = train_data,
                                       test.data = GSE16011.test.data,
                                       response = "label",
                                       cluster = GSE16011_common_lnc,
                                       seed = 123,
                                       method = c("svmLinear","svmRadial","svmPoly","nb","rf")[1],
                                       name = "GSE16011")
    # 朴素贝叶斯模型
    GSE16011.nb.res <- FastSVM(train.data = train_data,
                                test.data = GSE16011.test.data,
                                response = "label",
                                cluster = GSE16011_common_lnc,
                                seed = 123,
                                method = "nb",
                                name = "GSE16011")
    
    # 随机森林模型
    GSE16011.rf.res <- FastSVM(train.data = train_data,
                                test.data = GSE16011.test.data,
                                response = "label",
                                cluster = GSE16011_common_lnc,
                                seed = 123,
                                method = "rf",
                                name = "GSE16011")
    # 支持向量机poly模型
    GSE16011.svmPoly.res <- FastSVM(train.data = train_data,
                                     test.data = GSE16011.test.data,
                                     response = "label",
                                     cluster = GSE16011_common_lnc,
                                     seed = 123,
                                     method = "svmPoly",
                                     name = "GSE16011")
    # 支持向量机Radial模型
    GSE16011.svmRadial.res <- FastSVM(train.data = train_data,
                                       test.data = GSE16011.test.data,
                                       response = "label",
                                       cluster = GSE16011_common_lnc,
                                       seed = 123,
                                       method = "svmRadial",
                                       name = "GSE16011")
    # 线性判别模型
    GSE16011.lda.res <- FastLDA(train.data = train_data,
                                 test.data = GSE16011.test.data,
                                 response = "label",
                                 cluster = GSE16011_common_lnc,
                                 method = c("lda","qda","mda","fda","NaiveBayes")[1],
                                 name = "GSE16011")
    
    # lasso logistic回归模型
    train_data$label1 <- ifelse(train_data$label == "gbm",1,0)
    GSE16011.test.data$label1 <- ifelse(GSE16011.test.data$label == "gbm",1,0)
    GSE16011.log.res <- FastGlmnet(data = train_data,
                                    test.data = GSE16011.test.data,
                                    response = "label1",
                                    cluster = GSE16011_common_lnc,
                                    col = "label",
                                    control = "normal",
                                    alpha = c(1,0)[1],
                                    k = 10,
                                    family = c("binomial","cox")[1],
                                    optimize.method = c("min","lse")[1],
                                    seed = 123,
                                    plot.label = T,
                                    plot.lwd = 2,
                                    line.lty=3,
                                    line.lwd=2,
                                    line.col="blue",
                                    verbose = F,
                                    name = "GSE16011")
    
    save(GSE16011.svmLinear.res,
         GSE16011.nb.res,
         GSE16011.rf.res,
         GSE16011.svmPoly.res,
         GSE16011.svmRadial.res,
         GSE16011.lda.res,
         GSE16011.log.res,
         file="D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/GSE16011.classification.model.rda")
    
  }
  ## ROC曲线重绘图
  {
    load("D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/GSE7696.classification.model.rda")
    load("D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/GSE16011.classification.model.rda")
    load("D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/GSE108474.classification.model.rda")
    library(pROC)
    roc.list <- list(LSVM = GSE7696.svmLinear.res$roc.res,
                       LDA = GSE7696.lda.res$roc.res,
                       #LOGISTIC = GSE7696.log.res$Data$roc.res,
                       #RF = GSE7696.rf.res$roc.res,
                       #NB  = GSE7696.nb.res$roc.res,
                       PSVM = GSE7696.svmPoly.res$roc.res,
                       RSVM = GSE7696.svmRadial.res$roc.res
                 )
    roc.list1 <- list(LSVM = GSE16011.svmLinear.res$roc.res,
                      LDA = GSE16011.lda.res$roc.res,
                      #LOGISTIC = GSE16011.log.res$Data$roc.res,
                      #RF = GSE16011.rf.res$roc.res,
                     # NB  = GSE16011.nb.res$roc.res,
                      PSVM = GSE16011.svmPoly.res$roc.res,
                      RSVM = GSE16011.svmRadial.res$roc.res
    )
    roc.list2 <- list(LSVM = GSE108474.svmLinear.res$roc.res,
                      LDA = GSE108474.lda.res$roc.res,
                      #LOGISTIC = GSE108474.log.res$Data$roc.res,
                      #RF = GSE108474.rf.res$roc.res,
                      #NB  = GSE108474.nb.res$roc.res,
                      PSVM = GSE108474.svmPoly.res$roc.res,
                      RSVM = GSE108474.svmRadial.res$roc.res
    )
    p1 <- ggroc(roc.list,legacy.axes = TRUE) + theme_bw() + 
      geom_abline(intercept=0,slope=1,linetype = 2) + 
      ggtitle("GSE7696")
    p2 <- ggroc(roc.list1,legacy.axes = TRUE) + theme_bw() + 
      geom_abline(intercept=0,slope=1,linetype = 2) +
      ggtitle("GSE16011")
    p3 <- ggroc(roc.list2,legacy.axes = TRUE) + theme_bw() + 
      geom_abline(intercept=0,slope=1,linetype = 2) +
      ggtitle("GSE108474")
    pdf("D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/diagnosis.ROC.pdf")
    print(p1)
    print(p2)
    print(p3)
    dev.off()
    
    
  }
  ## AUC曲线下面积绘图
  {
    load("D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/GSE7696.classification.model.rda")
    load("D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/GSE16011.classification.model.rda")
    load("D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/GSE108474.classification.model.rda")
    roc.list <- list(LSVM = GSE7696.svmLinear.res$roc.res,
                     LDA = GSE7696.lda.res$roc.res,
                     #LOGISTIC = GSE7696.log.res$Data$roc.res,
                     #RF = GSE7696.rf.res$roc.res,
                     #NB  = GSE7696.nb.res$roc.res,
                     PSVM = GSE7696.svmPoly.res$roc.res,
                     RSVM = GSE7696.svmRadial.res$roc.res)
    roc.list1 <- list(LSVM = GSE16011.svmLinear.res$roc.res,
                      LDA = GSE16011.lda.res$roc.res,
                      #LOGISTIC = GSE16011.log.res$Data$roc.res,
                      #RF = GSE16011.rf.res$roc.res,
                      # NB  = GSE16011.nb.res$roc.res,
                      PSVM = GSE16011.svmPoly.res$roc.res,
                      RSVM = GSE16011.svmRadial.res$roc.res)
    roc.list2 <- list(LSVM = GSE108474.svmLinear.res$roc.res,
                      LDA = GSE108474.lda.res$roc.res,
                      #LOGISTIC = GSE108474.log.res$Data$roc.res,
                      #RF = GSE108474.rf.res$roc.res,
                      #NB  = GSE108474.nb.res$roc.res,
                      PSVM = GSE108474.svmPoly.res$roc.res,
                      RSVM = GSE108474.svmRadial.res$roc.res)
    library(ggplot2)
    df2 <- data.frame(dataset=c("GSE7696","GSE7696","GSE7696","GSE7696","GSE16011","GSE16011","GSE16011","GSE16011",
                                "GSE108474","GSE108474","GSE108474","GSE108474"),
                      method=c("LSVM","LDA","RSVM","PSVM","LSVM","LDA","RSVM","PSVM","LSVM","LDA","RSVM","PSVM"),
                      auc=c(GSE7696.svmLinear.res$roc.res$auc,GSE7696.lda.res$roc.res$auc,GSE7696.svmPoly.res$roc.res$auc,GSE7696.svmRadial.res$roc.res$auc,
                            GSE16011.svmLinear.res$roc.res$auc,GSE16011.lda.res$roc.res$auc,GSE16011.svmPoly.res$roc.res$auc,GSE16011.svmRadial.res$roc.res$auc,
                            GSE108474.svmLinear.res$roc.res$auc,GSE108474.lda.res$roc.res$auc,GSE108474.svmPoly.res$roc.res$auc,GSE108474.svmRadial.res$roc.res$auc) )
    pdf("D:/Rsources/Project/SingleCelllncRNA/Results/4.2Driver_sinature_clinical/auc.pdf")
    p <- ggplot(data=df2, aes(x=dataset, y=auc, fill=method)) +
      geom_bar(stat="identity", position=position_dodge()) +
      scale_fill_manual(values=mycolor) +
      theme_bw() 
    print(p)
    dev.off()
  }
  
  
  
}

