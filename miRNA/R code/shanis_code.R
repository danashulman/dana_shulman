#help from: https://dibsi-rnaseq.readthedocs.io/en/latest/DE.html
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("tximeta")

library("tximport")
library("DESeq2")
library("EnhancedVolcano")
library("ggplot2")
library("limma")
library("biomaRt")
library("calibrate")
library("edgeR")
library("gridExtra")

setwd("G:/PPMI/RNA_210720/quant/genes_BL")
dir<-"G:/PPMI/RNA_210720"
files_list = list.files()
files <- file.path(dir, "quant/genes_BL",files_list, "")
head(files)
#samples<-data.frame(ID=character(),timePoint=character(),Name=character(),fileName=character())
samples<-as.data.frame(files)
for(i in 1:length(files)){
  c <- NULL
  c <- strsplit(files[i], "[.]")
  samples$ID[i] <- c[[1]][2]
  samples$timePoint[i] <- c[[1]][3]
  samples$names[i] <- paste0(c[[1]][2],"_",c[[1]][3])
}
head(samples)

names(files) <- files_list
files
print(file.exists(files))

# library(tximportData)
# dir_annot <- system.file("extdata", package = "tximportData")
# list.files(dir_annot)
# library(readr)
# tx2gene <- read_csv(file.path(dir_annot, "tx2gene.gencode.v27.csv"))
# head(tx2gene)

# library("biomaRt")
quant_samp <- read.csv("G:/PPMI/RNA_210720/code_try/list_gene.csv") #annotation file (read in!)
length(unique(quant_samp$Name))#34569
length(quant_samp$Name)#34569
for(i in 1:length(quant_samp$Name)){
  c <- NULL
  c <- strsplit(as.character(quant_samp$Name[i]), "[.]")
  quant_samp$X[i] <- c[[1]][1]
}
length(unique(quant_samp$X))#34569

# ensembl <- useMart( "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org" )
# genemap <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
#                   filters = "ensembl_gene_id",
#                   values = quant_samp$X,
#                   mart = ensembl)


tx2gene <- quant_samp[,c(1,2)]
colnames(tx2gene)<-c("tx_name","gene_id")

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi.salmon$counts)
dim(txi.salmon$counts)

setwd("G:/PPMI/RNA_210720/output/BL")
TermPdata <- read.csv("G:/PPMI/RNA_210720/TermPdata_BL.csv")

joined_sampdata <- merge(samples, TermPdata, by.x = "names", by.y = "SampleName", all.x = TRUE, all.y = FALSE)
joined_sampdata <- joined_sampdata[ ,-c(5,6,11,12,17:31)]

table(joined_sampdata$combCurrentDiagnosisTerm)
joined_sampdata$GenderTerm[is.na(joined_sampdata$GenderTerm)]

#joined_sampdata$GenderTerm <- as.factor(joined_sampdata$GenderTerm)
#joined_sampdata$combCurrentDiagnosisTerm <- as.factor(joined_sampdata$combCurrentDiagnosisTerm)

write.csv(joined_sampdata,"joined_sampdata.csv")

## work 1
library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(txi.salmon, colData = joined_sampdata, design = ~GenderTerm + combDiagnosisTerm)
ddsTxi
mean(rowSums(counts(ddsTxi))) #466898.7
mean(rowMeans(counts(ddsTxi))) #427.9548

ddsTxi$combDiagnosisTerm <- relevel(ddsTxi$combDiagnosisTerm, ref = "HC")

#preforming DE analysis
ddsTxi <- DESeq(ddsTxi)

save(ddsTxi, file = "ddsTxi.RData")
countsnorm = counts(ddsTxi, normalize=T)
dim(countsnorm)

plotDispEsts(ddsTxi, main="Dispersion plot")

write.csv(as.data.frame(countsnorm), "countsnorm.csv")

resultsNames(ddsTxi) # check order of the coefficient, sanity check
#create vector with all desired contrasts
contrastList <- c("PDvsHC","SWEDDvsHC","PDvsSWEDD","GR_PDvsGR_UN","GC_PDvsGC_UN")

#prepare for adding annotations
library( "biomaRt" )

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
#ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",version = "Ensembl Genes 100", host = "uswest.ensembl.org")
#listAttributes(ensembl)
#listEnsembl()

#some times doesn't work because the server is full, so switch host (^) and try again.
genemap1 <- getBM(attributes = c("ensembl_gene_id","ensembl_gene_id_version","entrezgene_id","external_gene_name",
                                 "gene_biotype","description","chromosome_name","start_position","end_position"),
                  filters = "ensembl_gene_id",
                  values = rownames(countsnorm),
                  mart = ensembl)

idx_counts <- match(rownames(countsnorm), genemap1$ensembl_gene_id )
countsnorm_DF <- as.data.frame(countsnorm)
countsnorm_DF$ensembl_gene_id <- genemap1$ensembl_gene_id[idx_counts]
countsnorm_DF$ensembl_gene_id_version <- genemap1$ensembl_gene_id_version[idx_counts]
countsnorm_DF$entrez <- genemap1$entrezgene_id[idx_counts]
countsnorm_DF$external_gene_name <- genemap1$external_gene_name[idx_counts]
countsnorm_DF$rnacentral <- genemap1$rnacentral[idx_counts]
countsnorm_DF$gene_biotype <- genemap1$gene_biotype[idx_counts]
countsnorm_DF$description <- genemap1$description[idx_counts]
countsnorm_DF$chromosome_name <- genemap1$chromosome_name[idx_counts]
countsnorm_DF$start_position <- genemap1$start_position[idx_counts]
countsnorm_DF$end_position <- genemap1$end_position[idx_counts]

write.csv(as.data.frame(countsnorm_DF), "countsnormAnnotated.csv")

##create result objects for all contrasts
for (i in 1:length(contrastList)){
  resName <- paste('res', contrastList[i], sep='_')
  ##create res file
  assign(resName, results(ddsTxi, contrast=c("combDiagnosisTerm",strsplit(contrastList[i], "vs")[[1]][1],strsplit(contrastList[i], "vs")[[1]][2]), alpha=0.05))
  summary(get(resName))
  # # # graphs for checking stuff from the beginner guide to DESeq2
  # create bins using the quantile function 
  ##ratio of small p.values to mean normalized counts
  qs <- c( 0, quantile(get(resName)$baseMean[get(resName)$baseMean > 0], 0:7/7 ) )
  # "cut" the genes into the bins
  bins <- cut( get(resName)$baseMean, qs )
  # rename the levels of the bins using the middle point
  levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
  # calculate the ratio of ?p? values less than .01 for each bin
  ratios <- tapply( get(resName)$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
  # 
  # # plot these ratios
  png(paste('qc', contrastList[i],'p.valueVSmeancounts.png', sep='_'), 1000, 1000, pointsize=20)
  barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")
  dev.off()
  
  ##how much filtering is happening automaticly:
  metadata(get(resName))$filterThreshold  #res_DPDvsNDC 0% 0.6536203
  
  png(paste('qc', contrastList[i],'filterThreshold.png', sep='_'), 1000, 1000, pointsize=20)
  plot(metadata(get(resName))$filterNumRej,
       type="b", ylab="number of rejections",
       xlab="quantiles of filter")
  lines(metadata(get(resName))$lo.fit, col="red")
  abline(v=metadata(get(resName))$filterTheta)
  dev.off()
  ##Histogram of p values for all tests. The area shaded in blue indicates the subset of
  #those that pass the filtering, the area in khaki those that do not pass:
  
  use <- get(resName)$baseMean > metadata(get(resName))$filterThreshold
  h1 <- hist(get(resName)$pvalue[!use], breaks=0:50/50, plot=FALSE)
  h2 <- hist(get(resName)$pvalue[use], breaks=0:50/50, plot=FALSE)
  colori <- c(`do not pass`="khaki", `pass`="powderblue")
  
  png(paste('qc', contrastList[i],'p.valHistogram.png', sep='_'), 1000, 1000, pointsize=20)
  barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
          col = colori, space = 0, main = "", ylab="frequency")
  text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
       adj = c(0.5,1.7), xpd=NA)
  legend("topright", fill=rev(colori), legend=rev(names(colori)))
  dev.off()
  
  #create ordered file by p.value
  resNameOrdered <- paste(resName, "Ordered", sep='_')
  assign(resNameOrdered, get(resName)[order(get(resName)$padj),])# order results table by the smallest p value
  get(resNameOrdered)[1:10, ] #present 10 most significant genes
  sum(get(resName)$padj < 0.05, na.rm=TRUE) #How many genes were less than 0.05 
  
  png(paste('qc', contrastList[i],'MAaplot.png', sep='_'), 1000, 1000, pointsize=20)
  DESeq2::plotMA(get(resName), ylim=c(-30,30))
  dev.off()
  
  # #add gene symbol,entrez id, etc. from: https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf
  res_FullAnnot <- NULL
  res_FullAnnot <- as.data.frame(get(resNameOrdered))
  
  idx2 <- match(rownames(res_FullAnnot), genemap1$ensembl_gene_id)
  
  res_FullAnnot$ensembl_gene_id <- genemap1$ensembl_gene_id[idx2]
  res_FullAnnot$ensembl_gene_id_version <- genemap1$ensembl_gene_id_version[idx2]
  res_FullAnnot$entrez <- genemap1$entrezgene_id[idx2]
  res_FullAnnot$external_gene_name <- genemap1$external_gene_name[idx2]
  res_FullAnnot$rnacentral <- genemap1$rnacentral[idx2]
  res_FullAnnot$gene_biotype <- genemap1$gene_biotype[idx2]
  res_FullAnnot$description <- genemap1$description[idx2]
  res_FullAnnot$chromosome_name <- genemap1$chromosome_name[idx2]
  res_FullAnnot$start_position <- genemap1$start_position[idx2]
  res_FullAnnot$end_position <- genemap1$end_position[idx2]
  
  #save res
  write.csv(res_FullAnnot, paste(contrastList[i],'res.csv', sep='_'))
  
  #creat new dataframe of padj < 0.05
  resNamePadj <- paste("Padj", resNameOrdered, sep='_')
  assign(resNamePadj, as.data.frame(res_FullAnnot[which(res_FullAnnot$padj <= 0.05), ]))
  
  write.csv(get(resNamePadj), paste(contrastList[i],'resPadj.csv', sep='_'))
  
  ##plot volcanos
  #volc with of all
  # png(paste('volc',contrastList[i],'With.png', sep='_'), 1000, 1000, pointsize=20)
  # with(get(resName), plot(log2FoldChange, -log10(pvalue), pch=20, main=paste('Blood:',strsplit(contrastList[i], "vs")[[1]][1],'vs.',strsplit(contrastList[i], "vs")[[1]][2], sep=' ')))
  # with(subset(get(resName), padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  # with(subset(get(resName), abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="dark green"))
  # with(subset(get(resName), padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  # with(subset(get(resName), padj<.05 & abs(log2FoldChange)>1), 
  #      legend("topleft", xjust=1, yjust=1, legend=c("NS", 
  #                                                   paste("FDR<",0.05,sep=""), 
  #                                                   paste("|LogFC|>",1,sep=""), 
  #                                                   "both"), 
  #             pch=16, col=c("dark green", "blue", "red", "black")))
  # dev.off()
  
  #volc with of above10
  res_above10 <- res_FullAnnot[which(res_FullAnnot$baseMean >= 10),]
  res_above5 <- res_FullAnnot[which(res_FullAnnot$baseMean >= 5),]
  png(paste('volc_above10',contrastList[i],'With.png', sep='_'), 1000, 1000, pointsize=20)
  with(res_above10, plot(log2FoldChange, -log10(pvalue), pch=20, main=paste('Blood:',strsplit(contrastList[i], "vs")[[1]][1],'vs.',strsplit(contrastList[i], "vs")[[1]][2], sep=' ')))
  with(subset(res_above10, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  with(subset(res_above10, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="dark green"))
  with(subset(res_above10, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  with(subset(res_above10, padj<.05 & abs(log2FoldChange)>1), 
       legend("topleft", xjust=1, yjust=1, legend=c("NS", 
                                                    paste("FDR<",0.05,sep=""), 
                                                    paste("|LogFC|>",1,sep=""), 
                                                    "both"), 
              pch=16, col=c("dark green", "blue", "red", "black")))
  dev.off()
  
  #volc of all ggplot
  ggplot_full<- ggplot(res_FullAnnot, aes(x=log2FoldChange, y=-log10(pvalue))) + 
    geom_point(data = subset(res_FullAnnot, padj > 0.05), color = "black",size = 0.3,shape = 16) +
    geom_point(data = subset(res_FullAnnot, padj <= 0.05), color = "blue",size = 0.3,shape = 16) +
    geom_point(data = subset(res_FullAnnot, abs(log2FoldChange) >= 1), color = "dark green",size = 0.3,shape = 16) +
    geom_point(data = subset(res_FullAnnot, padj <= 0.05 & abs(log2FoldChange) >= 1), color = "red",size = 0.3, shape = 16) +
    xlab("Fold change (log2)") + ylab("-log10(P.Value)") + ggtitle(label = paste(contrastList[i]),subtitle = "ggplot afull volcano") +
    theme(axis.title=element_text(size=5), title = element_text(size=7)) #+ xlim(-50,30) + ylim(0,25) 
  ggsave(paste('volc',contrastList[i],'ggplot_full.png', sep='_'),plot = ggplot_full, width = 7, height = 7,units = c("cm"))
  
  res_above10 <- as.data.frame(res_above10)
  res_above5 <- as.data.frame(res_above5)
  ggplot_above10<- ggplot(res_above10, aes(x=log2FoldChange, y=-log10(pvalue))) +
    geom_point(data = subset(res_above10, padj > 0.05), color = "black",size = 0.3,shape = 16) +
    geom_point(data = subset(res_above10, padj <= 0.05), color = "blue",size = 0.3,shape = 16) +
    geom_point(data = subset(res_above10, abs(log2FoldChange) >= 1), color = "dark green",size = 0.3,shape = 16) +
    geom_point(data = subset(res_above10, padj <= 0.05 & abs(log2FoldChange) >= 1), color = "red",size = 0.3, shape = 16) +
    xlab("Fold change (log2)") + ylab("-log10(P.Value)") + ggtitle(label = paste(contrastList[i]),subtitle = "ggplot above 10 basemean volcano") +
    theme(axis.title=element_text(size=5), title = element_text(size=7)) #+ xlim(-6,6) + ylim(0,10)
  ggsave(paste('volc_above10',contrastList[i],'ggplot.png', sep='_'),plot = ggplot_above10, width = 7, height = 7,units = c("cm"))
  
  ggplot_above10_labs<- ggplot(res_above10, aes(x=log2FoldChange, y=-log10(pvalue))) +
    geom_point(data = subset(res_above10, padj > 0.05), color = "black",size = 0.3,shape = 16) +
    geom_point(data = subset(res_above10, padj <= 0.05), color = "blue",size = 0.3,shape = 16) +
    geom_point(data = subset(res_above10, abs(log2FoldChange) >= 1), color = "dark green",size = 0.3,shape = 16) +
    geom_point(data = subset(res_above10, padj <= 0.05 & abs(log2FoldChange) >= 1), color = "red",size = 0.3, shape = 16) +
    xlab("Fold change (log2)") + ylab("-log10(P.Value)") + ggtitle(label = paste(contrastList[i]),subtitle = "ggplot above 10 basemean volcano") +
    theme(axis.title=element_text(size=5), title = element_text(size=7)) + #xlim(-10,10) + ylim(0,10) +
    geom_text(aes(label=ifelse(padj <= 0.05 & abs(log2FoldChange) >= 1,external_gene_name,"")), color="black", size=1, vjust=0, hjust=-0.1)
  ggsave(paste('volc_above10',contrastList[i],'ggplot_withLabels.png', sep='_'),plot = ggplot_above10_labs, width = 7, height = 7,units = c("cm"))
  
  top10 <- res_above10[1:10,]
  ggplot_above10_labs2<- ggplot(res_above10, aes(x=log2FoldChange, y=-log10(pvalue))) +
    geom_point(data = subset(res_above10, padj > 0.05), color = "black",size = 0.3,shape = 16) +
    geom_point(data = subset(res_above10, padj <= 0.05), color = "blue",size = 0.3,shape = 16) +
    geom_point(data = subset(res_above10, abs(log2FoldChange) >= 1), color = "dark green",size = 0.3,shape = 16) +
    geom_point(data = subset(res_above10, padj <= 0.05 & abs(log2FoldChange) >= 1), color = "red",size = 0.3, shape = 16) +
    xlab("Fold change (log2)") + ylab("-log10(P.Value)") + ggtitle(label = paste(contrastList[i]),subtitle = "ggplot above 10 basemean volcano") +
    theme(axis.title=element_text(size=5), title = element_text(size=7)) + #xlim(-6,6) + ylim(0,10) +
    geom_text(aes(label=ifelse(rownames(res_above10)%in%rownames(top10),external_gene_name,"")), color="black", size=1, vjust=0, hjust=-0.1)
  ggsave(paste('volc_above10',contrastList[i],'ggplot_withLabels2.png', sep='_'),plot = ggplot_above10_labs2, width = 7, height = 7,units = c("cm"))
  
}  

## PCA
vsd <- varianceStabilizingTransformation(ddsTxi, blind=FALSE) #preforme trnsformation 
mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat, batch = vsd$GenderTerm) #correction of BOTH
assay(vsd) <- mat

png("pca_vsd_Diag.png", 500, 500, pointsize=20)
pcaData <- plotPCA(vsd, intgroup = c("combCurrentDiagnosisTerm"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = combCurrentDiagnosisTerm )) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()


png("pca_vsd_Gender.png", 500, 500, pointsize=20)
pcaData <- plotPCA(vsd, intgroup = c("GenderTerm"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = GenderTerm )) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

png("pca_vsd_Ethnicity.png", 500, 500, pointsize=20)
pcaData <- plotPCA(vsd, intgroup = c("EthnicityTerm"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = EthnicityTerm )) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

png("pca_vsd_corr_DiagGender.png", 500, 500, pointsize=20)
pcaData <- plotPCA(vsd, intgroup = c("combCurrentDiagnosisTerm", "GenderTerm"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = combCurrentDiagnosisTerm, shape = GenderTerm )) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

png("pca_vsd_corr_DiagPlusCurrent.png", 500, 500, pointsize=20)
pcaData <- plotPCA(vsd, intgroup = c("combCurrentDiagnosisTerm", "combDiagnosisTerm"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = combCurrentDiagnosisTerm, shape = combDiagnosisTerm )) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()


#plot heatmap
library("RColorBrewer")
library("pheatmap")
library("gplots")
library("genefilter")

#choose top 50 genes with the biggest variance
topVarGenes <- head( order( rowVars( assay(vsd) ), decreasing=TRUE ), 50 )
topVarGenes_100 <- head( order( rowVars( assay(vsd) ), decreasing=TRUE ), 100 )
#identify genes
topVarGenes_list <- assay(vsd)[ topVarGenes, ]
topVarGenes_list <- row.names(topVarGenes_list)
topVarGenes_annot <- countsnorm_DF[intersect(rownames(countsnorm_DF), topVarGenes_list), ]

topVarGenes_list_100 <- assay(vsd)[ topVarGenes_100, ]
topVarGenes_list_100 <- row.names(topVarGenes_list_100)
topVarGenes_annot_100 <- countsnorm_DF[intersect(rownames(countsnorm_DF), topVarGenes_list), ]
write.csv(topVarGenes_annot_100, "topVarGenes_vsd100_annotCounts.csv")

png("heatmap_vsd_corr_geneCluster.png", 1800, 1800, pointsize=26)
heatmap.2( assay(vsd)[ topVarGenes, ], scale="row", labCol =  paste(colnames(vsd),vsd$combCurrentDiagnosisTerm,vsd$GenderTerm,sep="-"), 
           trace="none", dendrogram="column",labRow = paste(rownames(vsd),topVarGenes_annot$external_gene_name, sep="-"),
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           cexRow=0.9, margins=c(7,16))
dev.off()

png("heatmap_vsd_corr_geneCluster_100.png", 1800, 1800, pointsize=26)
heatmap.2( assay(vsd)[ topVarGenes_100, ], scale="row", labCol =  paste(colnames(vsd),vsd$combCurrentDiagnosisTerm,vsd$GenderTerm,sep="-"), 
           trace="none", dendrogram="column",labRow = paste(rownames(vsd),topVarGenes_annot_100$external_gene_name, sep="-"),
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           cexRow=0.9, margins=c(7,16))
dev.off()