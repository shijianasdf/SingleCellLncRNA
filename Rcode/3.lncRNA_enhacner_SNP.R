#' @description signature lncRNA,enhancer,SNP共定位刻画
#' @author shi jian

## 读入GBM对应细胞系的增强子数据
gbm_enhancer <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/enhancer/gbm_enhancer.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) 
gbm_enhancer_region <- cbind.data.frame(gbm_enhancer$chromosome,gbm_enhancer$start,gbm_enhancer$end,rep("*",length(gbm_enhancer$chromosome)),gbm_enhancer$is.enhancer.a.super.enhancer)
colnames(gbm_enhancer_region) <- c("chr","start","end","strand","is.superEnhancer")
dim(gbm_enhancer_region) #20304 5

## 提取lncRNA区间信息
load("D:/Rsources/Project/SingleCelllncRNA/Results/1.singleCellSignature/lnc_sig.rda")
gtf.annotation <- getAnnotation(gtf.path="D:/Rsources/Project/SingleCelllncRNA/Data/AnnotationData/gencode.v35.annotation.gtf/gencode.v35.annotation.gtf",
                                select.col=c("seqnames","start","end","width","strand","source","type","gene_id","gene_name","gene_type"),
                                Newcolnames=c("chr","start","end","width","strand","source","type","ENSEMBL","SYMBOL","gene_type"),
                                control.id=NULL) 
save(gtf.annotation,file="D:/Rsources/Project/SingleCelllncRNA/Data/AnnotationData/gtf.annot.rda")
load("D:/Rsources/Project/SingleCelllncRNA/Data/AnnotationData/gtf.annot.rda")
lncRNA <- gtf.annotation[match(unique(lnc_sig),gtf.annotation$ENSEMBL),]
lncRNA_region <- cbind.data.frame(lncRNA$chr,lncRNA$start,lncRNA$end,lncRNA$strand,lncRNA$ENSEMBL)
dim(lncRNA_region) #14059 5
colnames(lncRNA_region) <- c("chr","start","end","strand","ensemble")
lncRNA_region$strand <- rep("*",length(lncRNA_region$chr))
lncRNA_enhancer_overlap <- RegionOverlapping.gr(lncRNA_region,gbm_enhancer_region)
head(lncRNA_enhancer_overlap)
dim(lncRNA_region[unique(lncRNA_enhancer_overlap$Qindex),]) #2778个lncRNA和GBM enhancer有overlap
dim(gbm_enhancer_region[unique(lncRNA_enhancer_overlap$Sindex),]) #3737个enhancer和lncRNA有overlap

## 读入GWAS GBM相关SNP信息
gbm_snp <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/gwas_snp.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
head(gbm_snp)
gbm_snp <- gbm_snp[grep("Glioma|glioma|Glioblastoma|glioblastoma",gbm_snp$DISEASE.TRAIT),]
gbm_snp_region <- unique(data.frame(chr=paste0("chr",gbm_snp$CHR_ID),start=as.numeric(gbm_snp$CHR_POS),end=as.numeric(gbm_snp$CHR_POS)+1,strand=rep("*",length(gbm_snp$CHR_ID)),rs=gbm_snp$SNPS))
dim(gbm_snp_region) #55 5
snp_lncRNA_overlap <- RegionOverlapping.gr(lncRNA_region,gbm_snp_region)
lncRNA_region[unique(snp_lncRNA_overlap$Qindex),]  # 4个lncRNA和GBM 相关SNP有overlap
# chr     start       end strand        ensemble 
# 2328  chr8 128634199 129683770      * ENSG00000229140
# 5164  chr9  21994139  22128103      * ENSG00000240498
# 5585 chr11  81821272  82718299      * ENSG00000245832
# 8509 chr12  75563202  75984015      * ENSG00000258077
gbm_snp[unique(snp_lncRNA_overlap$Sindex),] 
# SNPS           DISEASE.TRAIT   REGION CHR_ID   CHR_POS REPORTED.GENE.S.             MAPPED_GENE
# 244     rs891835                  Glioma  8q24.21      8 129479506           CCDC26                  CCDC26
# 16541 rs59060240            Glioblastoma   7p11.2      7  55080369               NR                    EGFR
# 243    rs4295627                  Glioma  8q24.21      8 129673211           CCDC26                  CCDC26
# 16548 rs12230172 Non-glioblastoma glioma  12q21.2     12  75848895       intergenic              AC078923.1
# 10301  rs6010620     Glioma (high-grade) 20q13.33     20  63678486            RTEL1 "RTEL1-TNFRSF6B, RTEL1"
# 3575   rs2736100                  Glioma  5p15.33      5   1286401             TERT                    TERT
# 13077  rs2736100                  Glioma  5p15.33      5   1286401             TERT                    TERT
# 248    rs4977756                  Glioma   9p21.3      9  22068653 "CDKN2A, CDKN2B"              CDKN2B-AS1
# 16563 rs55705857                  Glioma  8q24.21      8 129633446               NR                  CCDC26
# 16539  rs3851634            Glioblastoma  12q23.3     12 106419124           POLR3B                  POLR3B
# 16620 rs12803321                  Glioma  11q23.3     11 118609400               NR                  PHLDB1


## 扩展0,5000,10000,15000,20000后overlap情况
extends <- c(0,10000,20000,30000,40000)
lncRNA_region_gr_list <- list()
lncRNA_enhancer_ov_list <- list()
lncRNA_snp_ov_list <- list()
for(i in 1:length(extends)){
  lncRNA_region_gr <- extendRegion(lncRNA_region,extend=extends[i])
  lncRNA_enhancer_ov  <- RegionOverlapping.gr(lncRNA_region_gr,gbm_enhancer_region)
  lncRNA_snp_ov  <- RegionOverlapping.gr(lncRNA_region_gr,gbm_snp_region)
  lncRNA_region_gr_list[[i]] <- lncRNA_region_gr
  lncRNA_enhancer_ov_list[[i]] <- lncRNA_enhancer_ov
  lncRNA_snp_ov_list[[i]] <- lncRNA_snp_ov
}
names(lncRNA_region_gr_list) <- extends
names(lncRNA_enhancer_ov_list) <- extends
names(lncRNA_snp_ov_list) <- extends

# 生成可视化数据tt
tt <- NULL; pp<- NULL
for(i in 1:length(extends)){
  tp.i <- nrow(lncRNA_region[unique(lncRNA_enhancer_ov_list[[i]]$Qindex),])/nrow(lncRNA_region)
  tpp.i <- nrow(gbm_snp_region[unique(lncRNA_snp_ov_list[[i]]$Sindex),])/nrow(gbm_snp_region)
  x <- extends[i]
  tt <- rbind.data.frame(tt,data.frame(propotion=tp.i,extend=x))
  pp <- rbind.data.frame(pp,data.frame(propotion=tpp.i,extend=x))
}

# Line plot 
library(ggplot2)
pdf("D:/Rsources/Project/SingleCelllncRNA/Results/3.lncRNAEnhancerSNP/extend.enhancer.propotion.pdf")
lineplotA <- ggplot(data=tt, aes(x=extend, y=propotion))+
              geom_line()+
              geom_point()+
              labs(x="lncRNA extend,bp",y="proportion of lncRNA with enhancer")+
              theme_bw()
print(lineplotA)
dev.off()

library(ggplot2)
pdf("D:/Rsources/Project/SingleCelllncRNA/Results/3.lncRNAEnhancerSNP/extend.snp.propotion.pdf")
lineplotB <- ggplot(data=pp, aes(x=extend, y=propotion))+
              geom_line()+
              geom_point()+
              labs(x="lncRNA extend,bp",y="fraction of GWAS SNP")+
              theme_bw()
print(lineplotB)
dev.off()

## lncRNA和gbm snp关系卡方检验分析
{
  #SNP数据清洗
  snps <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/gwas_snp.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
  snps <- snps[!grepl(";|x",snps$CHR_ID),] 
  snps <- snps[snps$CHR_ID != "",]
  snps_region <- unique(data.frame(chr=paste0("chr",snps$CHR_ID),start=as.numeric(snps$CHR_POS),end=as.numeric(snps$CHR_POS)+1,strand=rep("*",length(snps$CHR_ID)),rs=snps$SNPS))
  #lncRNA上下游extend
  extends <- c(0,10000,20000,30000,40000)
  lncRNA_region_gr_list <- list()
  lncRNA_allsnp_ov_list <- list()
  for(i in 1:length(extends)){
    #lncRNA区间上下游扩展extend
    lncRNA_region_gr <- extendRegion(lncRNA_region,extend=extends[i])
    #lncRNA和snp区间overlap
    lncRNA_allsnp_ov  <- RegionOverlapping.gr(lncRNA_region_gr,snps_region)
    lncRNA_region_gr_list[[i]] <- lncRNA_region_gr
    lncRNA_allsnp_ov_list[[i]] <- lncRNA_allsnp_ov
  }
  names(lncRNA_region_gr_list) <- extends
  names(lncRNA_allsnp_ov_list) <- extends
  
  pos <- grepl("Glioma|glioma|Glioblastoma|glioblastoma",snps$DISEASE.TRAIT)
  is_GBM <- ifelse(pos,"Y","N")
  lncRNA_gbm_snp_p_list <- lapply(lncRNA_allsnp_ov_list,function(x){
    is_ov_lncRNA <- ifelse(snps$SNPS %in% snps_region$rs[unique(x$Sindex)],"Y","N")
    #构建卡方检验数据框
    lncRNA_gbm_snp_chisq <- unique(data.frame(rs=snps$SNPS,is.ov.lncRNA=is_ov_lncRNA,is.gbm=is_GBM))
    #卡方检验分析
    p <- chisq.test(lncRNA_gbm_snp_chisq$is.ov.lncRNA,lncRNA_gbm_snp_chisq$is.gbm)$p.value
    #p <- FastStatisticsTest(data = lncRNA_gbm_snp_chisq, variable = "is.ov.lncRNA", by = "is.gbm",
                            #test = "chisq.test", id = NULL, type = "categorical")
  })
  names(lncRNA_gbm_snp_p_list) <- names(lncRNA_allsnp_ov_list)
  save(lncRNA_gbm_snp_p_list,file="D:/Rsources/Project/SingleCelllncRNA/Results/3.lncRNAEnhancerSNP/lncRNA_gbm_snp_p_list.rda")
}

## lncRNA和SNP关系分析基于随机扰动
{
  library(regioneR)
  library(BSgenome.Hsapiens.UCSC.hg38.masked)
  lsp <- permTest(A=gbm_snp_region[,-4], B=lncRNA_region[,-4], ntimes=1000,genome="hg38",randomize.function=randomizeRegions,
                  evaluate.function=numOverlaps,alternative = "greater")
  summary(lsp)
  pdf("./Results/3.lncRNAEnhancerSNP/lncRNA.snp.association.pdf")
  plot(lsp)
  dev.off()
}
## lncRNA和enhancer关系分析基于随机扰动
{
  library(regioneR)
  library(BSgenome.Hsapiens.UCSC.hg38.masked)
  pt <- permTest(A=gbm_enhancer_region[,-4], B=lncRNA_region[,-4], ntimes=1000,genome="hg38",randomize.function=randomizeRegions,
                 evaluate.function=numOverlaps, alternative = "greater")
  summary(pt)
  pdf("./Results/3.lncRNAEnhancerSNP/lncRNA.enhancer.association.pdf")
  plot(pt)
  dev.off()
}

## 与增强子有overlap的lncRNA区间great富集分析
head(lncRNA_region[unique(lncRNA_enhancer_ov_list[[1]]$Qindex),])
fg <- FastGreat(region=lncRNA_region[unique(lncRNA_enhancer_ov_list[[1]]$Qindex),1:3],adv_upstream = 5,adv_downstream =5)
save(fg,file="./Results/3.lncRNAEnhancerSNP/lncRNA.great.rda")

load("./Results/3.lncRNAEnhancerSNP/lncRNA.great.rda")
