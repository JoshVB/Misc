library("tidyverse")
library("reshape2")
library("DESeq2")
library("gplots")

# read counts
tmp<-read.table("3col.tsv",header=F)
x<-as.matrix(acast(tmp, V2~V1, value.var="V3", fun.aggregate = sum))
x<-as.data.frame(x)
x$geneid<-sapply((strsplit(rownames(x),"\\|")),"[[",2)
xx<-aggregate(. ~ geneid,x,sum)
rownames(xx)<-xx$geneid
xx$geneid=NULL
xx<-round(xx)

pdf("MDSplot.pdf")
plot(cmdscale(dist(t(xx))), xlab="Coordinate 1", ylab="Coordinate 2", type = "n") ; text(cmdscale(dist(t(xx))), labels=colnames(xx)) 
dev.off()

# curate the samplesheet
samplesheet<-read.table("samplesheet.tsv",header=T,row.names=1)
samplesheet$crtl<-factor(samplesheet$crtl)
samplesheet$hcy<-factor(samplesheet$hcy)

####################################################
# Crtl Vs Hcy
####################################################
y<-xx[,which(colnames(xx) %in% rownames(samplesheet))]
y<-y[which(rowSums(y)/ncol(y)>=(10)),]
dds <- DESeqDataSetFromMatrix(countData = y , colData = samplesheet, design = ~ hcy )
res <- DESeq(dds)
z<- results(res)
vsd <- vst(dds, blind=FALSE)
zz<-cbind(as.data.frame(z),assay(vsd))
dge<-as.data.frame(zz[order(zz$pvalue),])
write.table(dge,file="dge_hcy_deseq.tsv",quote=F,sep="\t")
rnk<-as.data.frame( sign(dge$log2FoldChange) * (-log(dge$pvalue + 1E-307) ))
rownames(rnk)<-row.names(dge)
colnames(rnk)="Score"
write.table(rnk,file="dge_hcy_deseq.rnk",sep='\t',quote=F)

#some plots
pdf("dge_hcy_plots.pdf")
sig<-subset(dge,padj<0.05)
SIG=nrow(sig)
DN=nrow(subset(sig,log2FoldChange<0))
UP=nrow(subset(sig,log2FoldChange>0))
HEADER=paste("ST vs LS:", SIG , "DGEs,", UP ,"upregulated,", DN, "downregulated")
plot(log2(dge$baseMean),dge$log2FoldChange,cex=0.6, xlab="log2 base mean", ylab="log2 fold change" ,pch=19,col="#838383")
points(log2(sig$baseMean),sig$log2FoldChange,cex=0.6,pch=19,col="red")
mtext(HEADER)
top<-head(sig,20)
#text(log2(top$baseMean)+1, top$log2FoldChange, labels = rownames(top),cex=0.7)
#volcano plot
plot(dge$log2FoldChange, -log2(dge$pvalue) ,cex=0.6, xlim=c(-4,6),xlab="log2 fold change", ylab="-log2 p-value" ,pch=19,col="#838383")
points(sig$log2FoldChange, -log2(sig$pvalue),cex=0.6,pch=19,col="red")
#text(top$log2FoldChange+0.5, -log2(top$pvalue), labels = rownames(top),cex=0.7)
mtext(HEADER)
# top N gene heatmap
colfunc <- colorRampPalette(c("blue", "white", "red"))
heatmap.2(  as.matrix(dge[1:50,c(7:ncol(dge))]), col=colfunc(25),scale="row",
 trace="none",margins = c(6,6), cexRow=.4, main="Top 100 genes")
dev.off()

####################################################
# Crtl Vs Hcy
####################################################
