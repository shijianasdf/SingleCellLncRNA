#' @description 探索single cell表达数据
#' @author shi jian

## 读入单细胞GBM表达数据(TPM),查看数据有多少lncRNA
exprSet <- read.table(file = "D:/Rsources/Project/SingleCelllncRNA/Data/SingleCell/IDHwtGBM.processed.SS2.logTPM.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
table( convert(exprSet$GENE,fromtype="SYMBOL",totype = "gene_type") )
rownames(exprSet) <- exprSet[,1]
exprSet <- as.matrix(exprSet[,-1])
exprSet[1:4,1:4]
dim(exprSet)
## We excluded low-quality cells with less than 200,000 aligned reads or with less than 3000 detected genes
pos <- apply(exprSet,2,function(x){
  sum(as.numeric(x == 0)) > 3000
})
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
table( convert(rownames(exprSet),fromtype="ENSEMBL",totype = "gene_type") )
## 提取lncRNA表达谱
ID2gene <- cbind.data.frame(
                  SYMBOL = rownames(exprSet),
                  ENSEMBL = convert(rownames(exprSet),fromtype="SYMBOL",totype = "ENSEMBL"),
                  ENTREZID = convert(rownames(exprSet),fromtype="SYMBOL",totype = "ENTREZID"),
                  type = convert(rownames(exprSet),fromtype="SYMBOL",totype = "gene_type")
            )
rownames(exprSet) <- ID2gene$ENSEMBL
lnc_exprSet <- exprSet[ which(ID2gene$type == "lncRNA") ,  ]
lnc_exprSet[1:4,1:4]
dim(lnc_exprSet)
rowSums(lnc_exprSet)
lnc_exprSet[1:4,1:4]

## 对lncRNA表达谱进行主成分分析
library(FactoMineR)
library(factoextra)
p.res <- PCA(t(lnc_exprSet),graph = T)
eig.val <- get_eigenvalue(p.res)
pdf("D:/Rsources/Project/SingleCelllncRNA/Results/1.singleCellSignature/cumulativeVarianceAnalysis.pdf")
fe <- fviz_eig(p.res, addlabels = TRUE, ylim = c(0, 10))
fe
dev.off()

## 计算PCA第一主成分，第二主成分，第三主成分与lncRNA表达的spearman相关性
var1 <- get_pca_ind(p.res)
cor.res <- MakeCorrelationMatrix( x = t(lnc_exprSet),
                                  y = var1$coord[,1:3],
                                  size = 500, 
                                  adjust = "BH",
                                  method = c("pearson", "kendall", "spearman")[3] )
log <- abs(cor.res$r) > 0.2 & cor.res$fdr < 0.05
length(rownames(cor.res$p[log[,1],])) ## 第一主成分相关的lncRNA 80个
length(rownames(cor.res$p[log[,2],])) ## 第二主成分相关的lncRNA 52个
lnc_sig <- c( rownames(cor.res$p[log[,1],]), rownames(cor.res$p[log[,2],]) )
save(lnc_sig,file="D:/Rsources/Project/SingleCelllncRNA/Results/1.singleCellSignature/lnc_sig.rda")

## 计算差异表达lncRNA
design <- read.table(file = "D:/Rsources/Project/SingleCelllncRNA/Data/SingleCell/IDHwt.GBM.Metadata.SS2.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
rownames(design) <- gsub("-",".",design$NAME)
design$label <- ifelse(design$CellAssignment == "Macrophage","tumor","normal")
diff.res <- ScreenGenes(object = lnc_exprSet,
                        design = design,
                        contrast.col = "label",
                        contrast.level = c("tumor","normal"),
                        contrast.control = c("normal"),
                        cutoff.p = 0.05,
                        cutoff.q = 0.05,
                        cutoff.logFC = 1,
                        method = "BH",
                        verbose = F)
lnc_sig <- unique(intersect(diff.res[diff.res$lab != "NO",]$X1,lnc_sig))
save(lnc_sig,file = "D:/Rsources/Project/SingleCelllncRNA/Results/1.singleCellSignature/lnc_sig.rda")
save(diff.res,file = "D:/Rsources/Project/SingleCelllncRNA/Results/1.singleCellSignature/diff.res.rda")
load("D:/Rsources/Project/SingleCelllncRNA/Results/1.singleCellSignature/diff.res.rda")

## 差异表达lncRNA火山图展示
{
  this_title <- paste0('The number of up gene is ', nrow(diff.res[diff.res$lab =='UP',]),                      
                      'The number of down gene is ', nrow(diff.res[diff.res$lab =='DOWN',]))
  # Make a basic ggplot2 object with x-y values
  vol <- ggplot(diff.res, aes(x = logFC, y = -log10(qvalue), color = lab))+  
    ggtitle(label = this_title, subtitle = "Colored by fold-change direction") +
    geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
    scale_color_manual(name = "Directionality",
                       values = c(UP = "#008B00", DOWN = "#CD4F39", NO = "darkgray")) +
    theme_bw(base_size = 14) + # change overall theme
    theme(legend.position = "right") + # change the legend
    xlab(expression(log[2]("FoldChange"))) + # Change X-Axis label
    ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
    geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
    scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
  pdf(paste0("lnc_sig_t.test volcano plot.pdf"),8,8);print(vol);dev.off()
}



## 绘制lncRNA在单细胞中异质性的散点图
{
  load("D:/Rsources/Project/SingleCelllncRNA/Results/1.singleCellSignature/lnc_sig.rda")
  meta.data <- read.table("D:/Rsources/Project/SingleCelllncRNA/Data/SingleCell/IDHwt.GBM.Metadata.SS2.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
  seurat_expr <- lnc_exprSet[rownames(lnc_exprSet) %in% lnc_sig,]
  seurat_expr <- 10*(exp(seurat_expr)-1)
  table(meta.data$CellAssignment)
  # identical(meta.data$NAME,colnames(seurat_expr))
  meta.data$NAME <- gsub("-",".",meta.data$NAME,)
  seurat_expr1 <- seurat_expr[,meta.data$NAME]
  # seurat_expr1 <- seurat_expr1[,meta.data$CellAssignment == "Malignant"]

  tsne_data <- Run.Seurat.Cluster(seurat_expr1, 
                      nfeatures= 49, 
                      do.scale= T, 
                      do.center= T, 
                      seed= 1,
                      resolution= 0.2, 
                      perplexity= 30, 
                      n.neighbors= 30,            
                      n.components= 2, 
                      min.dist= 0.5,  # range 0.001 to 0.5
                      negative.sample.rate= 5, 
                      features= lnc_sig, 
                      plot.PC= "shijian.pdf"
  )
  tsne_data <- cbind.data.frame(tsne_data[-dim(tsne_data)[1],],meta.data$CellAssignment)
  colnames(tsne_data)[4] <- "cellType"
  pdf("D:/Rsources/Project/SingleCelllncRNA/Results/1.singleCellSignature/tsne.pdf") 
  p <- ggplot(tsne_data, aes(x= tSNE_1 ,y= tSNE_2)) + 
        geom_point(aes(col= Idents), size = 1.5, alpha= 0.8) + 
    #scale_colour_manual(values=  index_col) + 
        theme_bw()+ 
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),  panel.grid.minor = element_blank() , axis.line = element_line(colour = "black")) +
        xlab("tSNE_1") + 
        ylab("tSNE_2") + stat_ellipse()
  print(p)
  dev.off()
  
  
}
