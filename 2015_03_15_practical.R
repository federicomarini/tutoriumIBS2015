## ----echo=FALSE----------------------------------------------------------
library("knitr")
opts_chunk$set(fig.width = 6,
               fig.height = 4,
               message = FALSE,
               warning = FALSE,
               size='small',
               comment=NA,
               results='hide',
               fig.show='hide')

## ----eval=FALSE----------------------------------------------------------
## ## historically, htseq-count was used...
## 
## ## now there is featureCounts in the Rsubread package!
## ## http://subread.sourceforge.net/
## setwd("tutorium2015")
## 
## benchmarkData <- list.files(pattern=".*.bam$",full.names=T)
## 
## fcBenchmark <- featureCounts(files=benchmarkData,annot.ext="hg19_genes.gtf",isGTFAnnotationFile=T,nthreads=6)
## # nthreads exploits multicore architectures
## # very quick - and no information loss!
## # lot of other options for dealing with paired-end data, strand-specific data, multi-mapping reads, etc.
## 

## ----, message=FALSE,warning=FALSE,results='hide'------------------------
# source("http://bioconductor.org/biocLite.R")

# biocLite() is the Bioconductor installer function. Run it without any
# arguments to install the core packages or update any installed packages.
# This requires internet connectivity and will take some time!

library("DESeq2")
library("gplots")
library("RColorBrewer")
library("edgeR")
library("topGO")
library("dplyr")
library("RColorBrewer")

## ------------------------------------------------------------------------
## from the .RData object...
load("fcBenchmark.RData")
countMatrix <- fcBenchmark$counts

## or from the text file
countMatrix <- read.table("countMatrix.txt",header=TRUE,sep="\t")

## ------------------------------------------------------------------------
head(countMatrix)
colnames(countMatrix)
class(countMatrix)

## ----,results='markup'---------------------------------------------------
samplesDesign <- read.delim("sampleInfo.tsv")
samplesDesign

## ------------------------------------------------------------------------
## manually renaming the columns...
colnames(countMatrix) <- gsub(x = gsub(x = colnames(countMatrix),
                                       pattern = "_accepted_hits.bam$",
                                       replacement = ""),
                              pattern="^..",
                              replacement = "")

## ------------------------------------------------------------------------
?DESeqDataSetFromMatrix
?SummarizedExperiment

# reminder...
samplesDesign
dds <- DESeqDataSetFromMatrix(countMatrix,
                              colData = data.frame(condition=samplesDesign$condition),
                              design = ~ condition)
dds

colData(dds)
colnames(dds) <- samplesDesign$sampleID
dds$condition

## ----,cache=TRUE---------------------------------------------------------
dds <- DESeq(dds)

## if you have BiocParallel installed...
# library("BiocParallel")
# register(MulticoreParam(4))
# dds <- DESeq(dds,parallel=T,BPPARAM = MulticoreParam(4))

## ----,cache=TRUE---------------------------------------------------------
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
head(dispersions(dds))
dds <- estimateDispersions(dds)
# diagnostic plot: 
# plotDispEsts(dds)
dds <- nbinomWaldTest(dds)

## ----,cache=TRUE---------------------------------------------------------
rld <- rlogTransformation(dds)
# alternatively...
# vsd <- varianceStabilizingTransformation(dds)
par(mfrow = c(1,2))
dds <- estimateSizeFactors(dds)
plot(log2( 1 + counts(dds, normalized=TRUE)[ , c(1,3)] ),
     col=rgb(0,0,0,.2), pch=16, cex=0.3 )
plot(assay(rld)[ , c(1,3)],
     col=rgb(0,0,0,.2), pch=16, cex=0.3 )

## ------------------------------------------------------------------------
normCounts <- counts(dds, normalized=TRUE)

head(countMatrix)
head(normCounts)
head(assay(rld))

## ------------------------------------------------------------------------
# top N genes

N <- 50
select <- order(rowMeans(normCounts),decreasing=TRUE)[1:N]
hmcols <- colorRampPalette((brewer.pal(9, "Oranges")))(255)

heatmap.2(normCounts[select,], col = hmcols,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))
heatmap.2(assay(rld)[select,], col = hmcols,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))

## ------------------------------------------------------------------------
sampleDists <- dist(t(assay(rld)))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc <- hclust(sampleDists)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE )

# principal component plot of the samples
plotPCA(rld, intgroup = c("condition"))

## ------------------------------------------------------------------------
res <- results(dds)

summary(res)

## helper function
deseqresult2tbl <- function(deseqresult) {
  if (class(deseqresult) != "DESeqResults") stop("Not a DESeqResults object.")
  deseqresult <- as.data.frame(deseqresult)
  deseqresult$id <- rownames(deseqresult)
  rownames(deseqresult) <- NULL
  deseqresult <- tbl_df(deseqresult)
  deseqresult <- dplyr::select(deseqresult, id, log2FoldChange:padj)
  deseqresult %>% arrange(padj)
}
deseqresult2tbl(res)

## ------------------------------------------------------------------------
DESeq2::plotMA(res)

FDR <- 0.1
FDR <- 0.05

sigGenes <- subset(res, padj < FDR)

## ------------------------------------------------------------------------
attributes(res)
attr(res,"filterThreshold")
plot(attr(res,"filterNumRej"),ylab ="nr rejections",pch=16)

resNoFilt <- results(dds,independentFiltering = FALSE)
addmargins(table(withIndepFiltering = res$padj < FDR, 
                 withoutFiltering = resNoFilt$padj < FDR))

## ----,eval=FALSE---------------------------------------------------------
## ## simple export
## write.csv(res,file="res_exported.csv",quote=F)
## ## or with the sorted version
## write.csv(deseqresult2tbl(res),file="resSorted_exported.csv",quote=F,row.names=FALSE)
## 
## ## to interactively access the results...
## library("ReportingTools")
## rprt <- HTMLReport(shortName = "analysis",
##                    title = "RNA-Seq analysis",
##                    reportDirectory = "./reports")
## publish(res, rprt,
##         DataSet=dds)
## finish(rprt)
## 

## ------------------------------------------------------------------------
topGOtable <- function(DEgenes,BGgenes, ontology="BP",maxP = 0.001,desc="",annot = annFUN.org,
                       mapping = "org.Mm.eg.db",geneID = "symbol",topTablerows = 200,plotGraph=FALSE, 
                       plotNodes= 10,writeOutput=FALSE, outputFile=""
                       )
{
  DEgenes_input <- factor(as.integer(BGgenes %in% DEgenes))
  names(DEgenes_input) <- BGgenes
  GOdata <- new("topGOdata",ontology = ontology,allGenes = DEgenes_input,
                nodeSize = 10,annot = annot,mapping = mapping,ID = geneID)
  resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
  resultClassic <- runTest(GOdata,algorithm="classic",statistic = "fisher")
  sTab <- GenTable(GOdata,
                   p.value_elim=resultFisher,
                   p.value_classic=resultClassic,
                   orderBy= "p.value_elim",
                   ranksOf= "p.value_classic",
                   topNodes=topTablerows)
  sTabSig <- subset(sTab, as.numeric(p.value_elim) < maxP)
  if(writeOutput) write.table(sTab,file=outputFile,sep="\t",quote=F,col.names=T,row.names=F)
  if(plotGraph) showSigOfNodes(GOdata,topGO::score(resultFisher),firstSigNodes=plotNodes, useInfo="all")
  return(sTab)
}


## ----,cache=TRUE---------------------------------------------------------
geneUniverseExpr <- rownames(normCounts)[rowSums(normCounts) > 0]

# for biological processes
goenrich_DEgenes_BP <- topGOtable(DEgenes = rownames(sigGenes),
                                 BGgenes = geneUniverseExpr, 
                                 ontology="BP",maxP = 0.001, 
                                 geneID = "symbol",mapping = "org.Hs.eg.db",
                                 topTablerows = 100,plotGraph=FALSE)
# for molecular function
goenrich_DEgenes_MF <- topGOtable(DEgenes = rownames(sigGenes),
                               BGgenes = geneUniverseExpr, 
                               ontology="MF",maxP = 0.001, 
                               geneID = "symbol",mapping = "org.Hs.eg.db",
                               topTablerows = 100,plotGraph=FALSE)
goenrich_DEgenes_BP
goenrich_DEgenes_MF

## ----,cache=TRUE---------------------------------------------------------
library("edgeR")
countsEdger <- DGEList(counts=countMatrix,group=samplesDesign$condition)
y <- countsEdger
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)

keep <- rowSums(cpm(y)) >= 1
d <- y[keep,]
et <- exactTest(d)
topTags(et)
resTable <- topTags(et,dim(et)[1],sort.by = "logFC",adjust.method = "BH") # sorting by log fold change
resDF <- as.data.frame(resTable)

de <- decideTestsDGE(et, p=0.05, adjust="BH")
summary(decideTestsDGE(et, p=0.05, adjust="BH"))
detags <- rownames(d)[as.logical(de)]
plotSmear(et, de.tags=detags)
cpm(d)[detags,]
yy <- cpm(d, prior.count=2, log=TRUE)
heatmap.2(as.matrix(dist(t(yy))))

## ----,results='markup'---------------------------------------------------
sessionInfo()

