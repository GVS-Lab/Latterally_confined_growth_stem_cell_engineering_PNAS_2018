# load annoation libraries
library(EnsDb.Mmusculus.v79)
edb <- EnsDb.Mmusculus.v79
tx <- transcripts(edb, columns=c("tx_id", "gene_id", "gene_name"))
tx <- transcripts(edb, columns=c("tx_id", "gene_id", "gene_name"
                                 ,"gene_seq_start","gene_seq_end",
                                 "entrezid","seq_name"))
head(tx)
## extract the transcript ids and gene names and ll other informaiton
mapping <- cbind(tx_id=tx$tx_id, name=tx$gene_name, entrez=tx$entrezid,
                 genestart=tx$gene_seq_start,geneend=tx$gene_seq_end,
                 chrname=tx$seq_name,ensgene=tx$gene_id)
mapping<-as.data.frame(mapping)
head(mapping)
#install annotables to get chromosome number
#devtools::install_github("stephenturner/annotables")
library(dplyr)
library(annotables)
temp<-grcm38
temp_anno<-merge(mapping,temp,by="ensgene")
colnames(temp_anno)[2]<-"trnascript_id"

#merge all anotaiton files and aggregate by gene name
rnaseq_mapped_anno <- merge(rnaseq[,1:78],temp_anno,by="trnascript_id")
write.csv(rnaseq_mapped_anno, file="rnaseq_anno_full_tpkm.csv")
rnaseq_mapped_anno[,86]<-as.factor(rnaseq_mapped_anno[,86])


rnaseq_mapped_genes_chr<-aggregate(. ~ chr,data=rnaseq_mapped_genes_new[,c(75:78,86)],mean,na.rm=T)
png(filename="chr_activity_mean.png", units="in", width=3, height=4 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
x<-heatmap.2(as.matrix(rnaseq_mapped_genes_chr[,2:5]), col=greenred(75),
             scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=0.7,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",labRow = rnaseq_mapped_genes_chr[,1],labCol=c("3 Hour","3 Days","6 Days","10 Days"),density.info="none",
             keysize=0.6, key.par=list(mar=c(3.5,1,2.5,0)))
dev.off()
png(filename="Chr_activity_mean_Zscore_summed_barplot.png", units="in", width=1.5, height=1.5 , pointsize=5, res=1200)
ch_zscore<-as.data.frame(x$carpet)
t<-rowSums(ch_zscore,na.rm=T)
t_up<-max(t)+5
t_dw<-min(t)-5
barplot(t,las=1,ylab="Summed zscore",names =c("S1","S2","S3","S4"),ylim=c(t_dw,t_up))
box()
abline(h=0)
dev.off()
png(filename="Chrmean_Expression_summed_barplot.png", units="in", width=1.5, height=1.5 , pointsize=5, res=1200)
par(font.axis=2,font.lab=2)
x<-colSums(rnaseq_mapped_genes_chr[,2:5])
x_up<-max(x)+200
barplot(x,las=1,ylab="Summed Tpkm",names =c("S1","S2","S3","S4"),ylim=c(0,x_up),cex.axis=0.6)
box()
dev.off()
Chr<-c(1:19)
rnaseq_mapped_genes_chr_sub<-subset(rnaseq_mapped_genes_chr,rnaseq_mapped_genes_chr[,1]%in%Chr)
png(filename="chr_activity_mean_sub.png", units="in", width=3, height=4 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
x<-heatmap.2(as.matrix(rnaseq_mapped_genes_chr_sub[,2:5]), col=greenred(75),
             scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=0.7,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",labRow = rnaseq_mapped_genes_chr_sub[,1],labCol=c("3 Hour","3 Days","6 Days","10 Days"),density.info="none",
             keysize=0.6, key.par=list(mar=c(3.5,1,2.5,0)))
dev.off()
distmatrix<-distmatrix_S3-distmatrix_S1
library(corrplot)

distmatrix<-as.data.frame(matrix(nrow=19,ncol=19))
rownames(distmatrix)<-rnaseq_mapped_genes_chr_sub[,1]
colnames(distmatrix)<-rnaseq_mapped_genes_chr_sub[,1]
for(i in c(1:19)){
  distmatrix[,i]<-abs(rnaseq_mapped_genes_chr_sub[,2]-rnaseq_mapped_genes_chr_sub[i,2])/
    ((rnaseq_mapped_genes_chr_sub[,2]+rnaseq_mapped_genes_chr_sub[i,2]))
}
library(corrplot)

png(filename="Mean_Activity_distance_S1.png", units="in", width=3, height=3 , pointsize=6, res=1200)
x<-heatmap.2(as.matrix(distmatrix), col=bluered(75),
             scale="none",key=T, symkey=F, Rowv=F,Colv=FALSE,cexRow=1.2,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",density.info="none",keysize=1, key.par=list(mar=c(3,2,3,0)))
dev.off()
distmatrix->distmatrix_S1

distmatrix<-as.data.frame(matrix(nrow=19,ncol=19))
rownames(distmatrix)<-rnaseq_mapped_genes_chr_sub[,1]
colnames(distmatrix)<-rnaseq_mapped_genes_chr_sub[,1]
for(i in c(1:19)){
  distmatrix[,i]<-abs(rnaseq_mapped_genes_chr_sub[,3]-rnaseq_mapped_genes_chr_sub[i,3])/
    (rnaseq_mapped_genes_chr_sub[,3]+rnaseq_mapped_genes_chr_sub[i,3])
}
library(corrplot)
png(filename="Mean_Activity_distance_S2.png", units="in", width=3, height=3 , pointsize=6, res=1200)
x<-heatmap.2(as.matrix(distmatrix), col=bluered(75),
             scale="none",key=T, symkey=F, Rowv=F,Colv=FALSE,cexRow=1.2,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",density.info="none",keysize=1, key.par=list(mar=c(3,2,3,0)))
dev.off()
distmatrix->distmatrix_S2

distmatrix<-as.data.frame(matrix(nrow=19,ncol=19))
rownames(distmatrix)<-rnaseq_mapped_genes_chr_sub[,1]
colnames(distmatrix)<-rnaseq_mapped_genes_chr_sub[,1]
for(i in c(1:19)){
  distmatrix[,i]<-abs(rnaseq_mapped_genes_chr_sub[,4]-rnaseq_mapped_genes_chr_sub[i,4])/
    (rnaseq_mapped_genes_chr_sub[,4]+rnaseq_mapped_genes_chr_sub[i,4])
}
library(corrplot)
png(filename="Mean_Activity_distance_S3.png", units="in", width=3, height=3 , pointsize=6, res=1200)
x<-heatmap.2(as.matrix(distmatrix), col=bluered(75),
             scale="none",key=T, symkey=F, Rowv=F,Colv=FALSE,cexRow=1.2,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",density.info="none",keysize=1, key.par=list(mar=c(3,2,3,0)))
dev.off()
distmatrix->distmatrix_S3

distmatrix<-as.data.frame(matrix(nrow=19,ncol=19))
rownames(distmatrix)<-rnaseq_mapped_genes_chr_sub[,1]
colnames(distmatrix)<-rnaseq_mapped_genes_chr_sub[,1]
for(i in c(1:19)){
  distmatrix[,i]<-abs(rnaseq_mapped_genes_chr_sub[,5]-rnaseq_mapped_genes_chr_sub[i,5])/
    (rnaseq_mapped_genes_chr_sub[,5]+rnaseq_mapped_genes_chr_sub[i,5])
}
library(corrplot)
png(filename="Mean_Activity_distance_S4.png", units="in", width=3, height=3 , pointsize=6, res=1200)
x<-heatmap.2(as.matrix(distmatrix), col=bluered(75),
             scale="none",key=T, symkey=F, Rowv=F,Colv=FALSE,cexRow=1.2,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",density.info="none",keysize=1, key.par=list(mar=c(3,2,3,0)))
dev.off()
distmatrix->distmatrix_S4





distmatrix<-distmatrix_S3-distmatrix_S1
library(corrplot)
png(filename="Mean_Activity_distance_S3S1.png", units="in", width=3, height=3 , pointsize=6, res=1200)
x<-heatmap.2(as.matrix(distmatrix), col=bluered(75),
             scale="none",key=T, symkey=F, Rowv=F,Colv=FALSE,cexRow=1.2,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",density.info="none",keysize=1, key.par=list(mar=c(3,2,3,0)))
dev.off()
png(filename="hist_S3_S1_chr_dist.png", units="in", width=3, height=3 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
hist(as.vector(as.matrix(distmatrix)), col="gray", breaks=50, xlim=c(-0.5,0.5))
box()
abline(v=-0.2, col="red", lwd=2)
abline(v=0.2, col="red", lwd=2)
dev.off()
temp<-as.data.frame(which(distmatrix< (-0.2), TRUE))
temp[,3]<-paste(rownames(distmatrix)[temp[,1]],colnames(distmatrix)[temp[,2]],sep='-' )
temp[,3]
temp<-as.data.frame(which(distmatrix> (0.2), TRUE))
temp[,3]<-paste(rownames(distmatrix)[temp[,1]],colnames(distmatrix)[temp[,2]],sep='-' )
temp[,3]

distmatrix<-distmatrix_S2-distmatrix_S1
library(corrplot)
png(filename="Mean_Activity_distance_S2S1.png", units="in",  width=3, height=3 , pointsize=6, res=1200)
x<-heatmap.2(as.matrix(distmatrix), col=bluered(75),
             scale="none",key=T, symkey=F, Rowv=F,Colv=FALSE,cexRow=1.2,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",density.info="none",keysize=1, key.par=list(mar=c(3,2,3,0)))
dev.off()
png(filename="hist_S2_S1_chr_dist.png", units="in", width=3, height=3  , pointsize=5, res=1200)
par(font.axis=2,font.lab=2)
hist(as.vector(as.matrix(distmatrix)), col="gray", breaks=50, xlim=c(-0.5,0.5))
box()
abline(v=-0.2, col="red", lwd=2)
abline(v=0.2, col="red", lwd=2)
dev.off()
temp<-as.data.frame(which(distmatrix< (-0.2), TRUE))
temp[,3]<-paste(rownames(distmatrix)[temp[,1]],colnames(distmatrix)[temp[,2]],sep='-' )
temp[,3]
temp<-as.data.frame(which(distmatrix> (0.2), TRUE))
temp[,3]<-paste(rownames(distmatrix)[temp[,1]],colnames(distmatrix)[temp[,2]],sep='-' )
temp[,3]


distmatrix<-distmatrix_S4-distmatrix_S3
library(corrplot)
png(filename="Mean_Activity_distance_S4S3.png", units="in", width=3, height=3 , pointsize=6, res=1200)
x<-heatmap.2(as.matrix(distmatrix), col=bluered(75),
             scale="none",key=T, symkey=F, Rowv=F,Colv=FALSE,cexRow=1.2,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",density.info="none",keysize=1, key.par=list(mar=c(3,2,3,0)))
dev.off()
png(filename="hist_S4_S3_chr_dist.png", units="in", width=3, height=3 , pointsize=5, res=1200)
par(font.axis=2,font.lab=2)
hist(as.vector(as.matrix(distmatrix)), col="gray", breaks=50, xlim=c(-0.5,0.5))
box()
abline(v=-0.2, col="red", lwd=2)
abline(v=0.2, col="red", lwd=2)
dev.off()
temp<-as.data.frame(which(distmatrix< (-0.2), TRUE))
temp[,3]<-paste(rownames(distmatrix)[temp[,1]],colnames(distmatrix)[temp[,2]],sep='-' )
temp[,3]
temp<-as.data.frame(which(distmatrix> (0.2), TRUE))
temp[,3]<-paste(rownames(distmatrix)[temp[,1]],colnames(distmatrix)[temp[,2]],sep='-' )
temp[,3]


distmatrix<-distmatrix_S2-distmatrix_S3
library(corrplot)
png(filename="Mean_Activity_distance_S2S3.png", units="in",  width=3, height=3 , pointsize=6, res=1200)
x<-heatmap.2(as.matrix(distmatrix), col=bluered(75),
             scale="none",key=T, symkey=F, Rowv=F,Colv=FALSE,cexRow=1.2,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",density.info="none",keysize=1, key.par=list(mar=c(3,2,3,0)))
dev.off()
png(filename="hist_S2_S3_chr_dist.png", units="in", width=3, height=3 , pointsize=5, res=1200)
par(font.axis=2,font.lab=2)
hist(as.vector(as.matrix(distmatrix)), col="gray", breaks=50, xlim=c(-0.5,0.5))
box()
abline(v=-0.2, col="red", lwd=2)
abline(v=0.2, col="red", lwd=2)
dev.off()
temp<-as.data.frame(which(distmatrix< (-0.2), TRUE))
temp[,3]<-paste(rownames(distmatrix)[temp[,1]],colnames(distmatrix)[temp[,2]],sep='-' )
temp[,3]
temp<-as.data.frame(which(distmatrix> (0.2), TRUE))
temp[,3]<-paste(rownames(distmatrix)[temp[,1]],colnames(distmatrix)[temp[,2]],sep='-' )
temp[,3]

distmatrix<-distmatrix_S4-distmatrix_S1
library(corrplot)
png(filename="Mean_Activity_distance_S4S1.png", units="in", width=3, height=3 , pointsize=6, res=1200)
x<-heatmap.2(as.matrix(distmatrix), col=bluered(75),
             scale="none",key=T, symkey=F, Rowv=F,Colv=FALSE,cexRow=1.2,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",density.info="none",keysize=1, key.par=list(mar=c(3,2,3,0)))
dev.off()
png(filename="hist_S4_S1_chr_dist.png", units="in", width=3, height=3 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
hist(as.vector(as.matrix(distmatrix)), col="gray", breaks=50, xlim=c(-0.5,0.5))
box()
abline(v=-0.2, col="red", lwd=2)
abline(v=0.2, col="red", lwd=2)
dev.off()
temp<-as.data.frame(which(distmatrix< (-0.2), TRUE))
temp[,3]<-paste(rownames(distmatrix)[temp[,1]],colnames(distmatrix)[temp[,2]],sep='-' )
temp[,3]
temp<-as.data.frame(which(distmatrix> (0.2), TRUE))
temp[,3]<-paste(rownames(distmatrix)[temp[,1]],colnames(distmatrix)[temp[,2]],sep='-' )
temp[,3]


distmatrix<-distmatrix_S2-distmatrix_S4
library(corrplot)
png(filename="Mean_Activity_distance_S2S4.png", units="in", width=3, height=3 , pointsize=6, res=1200)
x<-heatmap.2(as.matrix(distmatrix), col=bluered(75),
             scale="none",key=T, symkey=F, Rowv=F,Colv=FALSE,cexRow=1.2,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",density.info="none",keysize=1, key.par=list(mar=c(3,2,3,0)))
dev.off()
png(filename="hist_S2_S4_chr_dist.png", units="in", width=3, height=3 , pointsize=5, res=1200)
par(font.axis=2,font.lab=2)
hist(as.vector(as.matrix(distmatrix)), col="gray", breaks=50, xlim=c(-0.5,0.5))
box()
abline(v=-0.2, col="red", lwd=2)
abline(v=0.2, col="red", lwd=2)
dev.off()
temp<-as.data.frame(which(distmatrix< (-0.2), TRUE))
temp[,3]<-paste(rownames(distmatrix)[temp[,1]],colnames(distmatrix)[temp[,2]],sep='-' )
temp[,3]
temp<-as.data.frame(which(distmatrix> (0.2), TRUE))
temp[,3]<-paste(rownames(distmatrix)[temp[,1]],colnames(distmatrix)[temp[,2]],sep='-' )
temp[,3]