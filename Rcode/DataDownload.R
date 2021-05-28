#' @description  下载GEO乳腺癌芯片表达数据以及TCGA乳腺癌所有数据
#' @author shi jian

## 读入GEO数据表格,批量下载GEO BRCA芯片数据
{
  GEOid <- read.table(file = "D:/Rsources/Project/SingleCelllncRNA/Data/breast_cancer.txt",sep="\t",header=F,fill=T,quote=NULL,stringsAsFactors=F)
  apply(GEOid,1,function(x){
    FastDownloadGEO(x,"D:/Rsources/Project/SingleCelllncRNA/Data",origin = F)
  })
}

## 下载TCGA-BRCA 表达数据，拷贝数数据以及突变数据 临床数据
{
  FastTCGADownload(project = "TCGA-BRCA",filepath = "D:/Rsources/Rcode/BioinforRCode-master/data/TCGAbiolinks and GDC - Hamornized.csv")
  FastTCGADownload(project = "TCGA-BRCA",filepath = "D:/Rsources/Rcode/BioinforRCode-master/data/TCGAbiolinks and GDC - Copy of Legacy.csv",legacy=T)
  
  ## 下载TCGA-BRCA亚型信息
  brca.subtype <- TCGAquery_subtype(tumor = "brca")
  brca.subtype <- as.data.frame(brca.subtype)
  head(brca.subtype)
  save(brca.subtype,file="D:/Rsources/Project/SingleCelllncRNA/Data/TCGA-BRCA/brca.subtype.rda")
  
  subtypes <- PanCancerAtlas_subtypes()
  subtypes <- as.data.frame(subtypes)
  head(subtypes[subtypes$cancer.type == "BRCA",])
}

## 读入GEO数据表格,批量下载GEO GBM芯片数据
{
  GEOid <- read.table(file = "D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM/GBM.txt",sep="\t",header=F,fill=T,quote=NULL,stringsAsFactors=F)
  apply(GEOid,1,function(x){
    FastDownloadGEO(x,"D:/Rsources/Project/SingleCelllncRNA/Data/GEO-GBM",origin = F)
  })
}

## 下载TCGA-GBM 表达数据，拷贝数数据以及突变数据 临床数据
{
  FastTCGADownload(project = "TCGA-GBM",filepath = "D:/Rsources/Rcode/BioinforRCode-master/data/TCGAbiolinks and GDC - Hamornized.csv")
  FastTCGADownload(project = "TCGA-GBM",filepath = "D:/Rsources/Rcode/BioinforRCode-master/data/TCGAbiolinks and GDC - Copy of Legacy.csv",legacy=T)
  
  query <- GDCquery(project ='TCGA-GBM',data.category = 'Clinical',access = 'open',legacy = FALSE,file.type = "xml")
  clinical <- GDCprepare_clinic(query, clinical.info = "patient")
  clinical.admin <- GDCprepare_clinic(query, clinical.info = "admin")
  clinical.follow_up <- GDCprepare_clinic(query, clinical.info = "follow_up")
  clinical.stage_event <- GDCprepare_clinic(query, clinical.info = "stage_event")
  clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")
  clinical.radiation <- GDCprepare_clinic(query, clinical.info = "radiation")
  clinical.new_tumor_event <- GDCprepare_clinic(query, clinical.info = "new_tumor_event")
  save(clinical,clinical.admin,clinical.follow_up,clinical.stage_event,clinical.drug,clinical.radiation,clinical.new_tumor_event,file=paste0(paste0("./",'TCGA-GBM',"/"),'TCGA-GBM',"-Clinical",".rda")) 
  
  ## 下载TCGA-GBM亚型信息
  GBM.subtype <- TCGAquery_subtype(tumor = "gbm")
  GBM.subtype <- as.data.frame(GBM.subtype)
  head(GBM.subtype)
  save(GBM.subtype,file="D:/Rsources/Project/SingleCelllncRNA/Data/TCGA-GBM/GBM.subtype.rda")
  
  subtypes <- PanCancerAtlas_subtypes()
  subtypes <- as.data.frame(subtypes)
  head(subtypes[subtypes$cancer.type == "GBM",])
  head(subtypes[subtypes$cancer.type == "LGG",])
}

