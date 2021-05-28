#' @description signature lncRNA功能刻画
#' @author shi jian

## signature LncRNA共表达分析
load("D:/Rsources/Project/SingleCelllncRNA/Results/1.singleCellSignature/lnc_sig.rda")
exprSet <- read.table(file = "D:/Rsources/Project/SingleCelllncRNA/Data/SingleCell/IDHwtGBM.processed.SS2.logTPM.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
rownames(exprSet) <- exprSet[,1]
exprSet <- as.matrix(exprSet[,-1])
exprSet[1:4,1:4]
dim(exprSet)

## We excluded low-quality cells with less than 200,000 aligned reads or with less than 3000 detected genes
pos <- apply(exprSet,2,function(x){sum(as.numeric(x == 0)) > 3000})
exprSet <- exprSet[,pos]
#abundant expression were selected according to two quality
#filtering criterions: (i) genes expressed in more than three single
#cells for each tumor or detected by the bulk RNA-seq and (ii) genes
#with average expression level more than three [log2(TPM/10+1) > 0.2]
a <- apply(exprSet,1,function(x){mean(x)})
plot(density(a))
fivenum(a)
exprSet <- exprSet[a > 0.15, ]
dim(exprSet)
## 提取lncRNA表达谱
ID2gene <- cbind.data.frame(
  SYMBOL = rownames(exprSet),
  ENSEMBL = convert(rownames(exprSet),fromtype="SYMBOL",totype = "ENSEMBL"),
  ENTREZID = convert(rownames(exprSet),fromtype="SYMBOL",totype = "ENTREZID"),
  type = convert(rownames(exprSet),fromtype="SYMBOL",totype = "gene_type")
)
rownames(exprSet) <- ID2gene$ENSEMBL




lnc_exprSet <- exprSet[which(ID2gene$type == "lncRNA"),  ]
lnc_exprSet[1:4,1:4]
dim(lnc_exprSet)
rowSums(lnc_exprSet)
lnc_exprSet[1:4,1:4]

pro_exprSet <- exprSet[which(ID2gene$type == "protein_coding"),  ]
## 计算lncRNA signature共表达蛋白编码基因(spearman)
lnc.gene.cor.res <- MakeCorrelationMatrix( x = t(pro_exprSet[!is.na(rownames(pro_exprSet)),]),
                                           y = t(lnc_exprSet[unique(lnc_sig),]),
                                           size = 500, 
                                           adjust = "BH",
                                           method = c("pearson", "kendall", "spearman")[3])
dim( lnc.gene.cor.res$r ) #11620   124
## 宽数据转换为长数据
lnc.gene.cor <- outputDF(lnc.gene.cor.res,symmetric= FALSE)
## 筛选共表达显著的蛋白编码基因进行功能注释
dim(lnc.gene.cor[(abs(lnc.gene.cor$r) >= 0.4 & lnc.gene.cor$fdr <= 0.05),])
co_gene <- as.character(lnc.gene.cor[(abs(lnc.gene.cor$r) >= 0.4 & lnc.gene.cor$fdr <= 0.05),]$V1)
geneList <- pro_exprSet[,1]
names(geneList) <- rownames(pro_exprSet)

go.res <- FastGO(genes = co_gene,
                 geneList = geneList,
                 default.universe = F,
                 classlevel = 2:2,
                 OrgDb  = NULL,
                 keyType = NULL,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff  = 0.05,
                 cnet.showCategory = 5,
                 verbose = TRUE,
                 save.path = "GO",
                 names = "love")

tgo.res <- FastGO(genes = co_gene,
                 geneList = geneList,
                 default.universe = T,
                 classlevel = 2:2,
                 OrgDb  = NULL,
                 keyType = NULL,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff  = 0.05,
                 cnet.showCategory = 5,
                 verbose = TRUE,
                 save.path = "GO",
                 names = "llove")
save(go.res,tgo.res,file = "D:/Rsources/Project/SingleCelllncRNA/Results/2.SigLncfunction/go.res.rda")


kegg.res <- FastKEGG(genes = co_gene,
                     geneList = geneList,
                     default.universe = F,
                     organism = 'hsa',
                     pvalueCutoff = list(enrichKEGG = 0.05,
                                         enrichMKEGG = 0.05),
                     qvalueCutoff = 0.05,
                     cnet.showCategory = 10,
                     verbose = T,
                     save.path = "KEGG",
                     names = "love")

lkegg.res <- FastKEGG(genes = co_gene,
                     geneList = geneList,
                     default.universe = T,
                     organism = 'hsa',
                     pvalueCutoff = list(enrichKEGG = 0.05,
                                         enrichMKEGG = 0.05),
                     qvalueCutoff = 0.05,
                     cnet.showCategory = 10,
                     verbose = T,
                     save.path = "KEGG",
                     names = "llove")
save(kegg.res,lkegg.res,file = "D:/Rsources/Project/SingleCelllncRNA/Results/2.SigLncfunction/kegg.res.rda")

## lncRNA signature GSEA富集分析
lnc.gene.cor.list <- split(lnc.gene.cor,lnc.gene.cor$V2)
geneLists <- lapply(lnc.gene.cor.list,function(x){ tt <- x$r;names(tt) <- x$V1;tt })
names(geneLists) <- names(lnc.gene.cor.list)
gsea.res <- NULL
for(i in names(geneLists)[1:19]){
  print(i)
  gsea.res.i <- FastGSEA(geneList = geneLists[[i]],
                         OrgDb = NULL,
                         keyType = NULL,
                         pvalueCutoff = list(gseGO = 0.05,
                                             gseKEGG = 0.05,
                                             gseMKEGG = 0.05),
                         verbose = T,
                         save.path = "GSEA",
                         names = i)
  gsea.res <- c(gsea.res,gsea.res.i)
}

save(gsea.res,file="D:/Rsources/Project/SingleCelllncRNA/Results/2.SigLncfunction/gsea.res.rda")


## 重新绘制enrichplot
{
  load("D:/Rsources/Project/SingleCelllncRNA/Results/2.SigLncfunction/go.res.rda")
  p <- FastEnrichPlot(object = go.res$GOEnrichment$BP,
                      pvalueCutoff = NULL,
                      qvalueCutoff = NULL,
                      id.col = "Description",
                      select = NULL,
                      x.title = "Biological Process terms",
                      y.title = "-log10(FDR)",
                      size = 12,
                      short.cutoff = 46,
                      topshow = 10)
  p1 <- FastEnrichPlot(object = go.res$GOEnrichment$MF,
                       pvalueCutoff = NULL,
                       qvalueCutoff = NULL,
                       id.col = "Description",
                       select = NULL,
                       x.title = "Molecular Functions terms",
                       y.title = "-log10(FDR)",
                       size = 12,
                       short.cutoff = 46,
                       topshow = 10)
  p2 <- FastEnrichPlot(object = go.res$GOEnrichment$CC,
                       pvalueCutoff = NULL,
                       qvalueCutoff = NULL,
                       id.col = "Description",
                       select = NULL,
                       x.title = "Cellular Components terms",
                       y.title = "-log10(FDR)",
                       size = 12,
                       short.cutoff = 46,
                       topshow = 10)
  pdf("D:/Rsources/Project/SingleCelllncRNA/Results/2.SigLncfunction/GO/love/GO/cogenes_enrichment.pdf")
  print(p$Plot)
  print(p1$Plot)
  print(p2$Plot)
  dev.off()
 
}



  
