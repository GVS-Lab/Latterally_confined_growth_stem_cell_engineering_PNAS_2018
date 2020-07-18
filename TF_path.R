library(gskb)
data(mm_TF)
rnaseq_1<-rnaseq_mapped_genes_new[complete.cases(rnaseq_mapped_genes_new), ]
rnaseq_1[,92]<-toupper(rnaseq_1[,85])

tf_sum<-as.data.frame(matrix(nrow=372, ncol=5))
colnames(tf_sum)<-c("TF","S1","S2","S3","S4")

for(i in 1:372){
  temp<-subset(rnaseq_1,rnaseq_1[,92]%in% (mm_TF[[i]][3:length(mm_TF[[i]])]))[,75:78]
  tf_sum[i,1]<-mm_TF[[i]][1]
  tf_sum[i,2:5]<-colSums(temp)
}
library(gplots)
png(filename="TF_activity.png", units="in", width=3, height=7 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
x<-heatmap.2(as.matrix(tf_sum[,2:5]), col=greenred(75),
             scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",labRow = tf_sum[,1],labCol=c("S1","S2","S3","S4"),density.info="none",
             keysize=0.6, key.par=list(mar=c(3.5,1,2.5,0)))
dev.off()

tf_sum[,6]<-tf_sum[,5]/tf_sum[,2]
tf_sum[,7]<-rowMeans(tf_sum[,2:5])
tf_sum_1<-subset(tf_sum, tf_sum[,7]>500)
tf_sum_1<-tf_sum_1[order(tf_sum_1[,6], decreasing=F),]
png(filename="TF_activity_top20_decreasing.png", units="in", width=3, height=5 , pointsize=8, res=1200)
x<-heatmap.2(as.matrix(tf_sum_1[(1:30),2:5]), col=redgreen(75),
             scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",labRow = tf_sum_1[1:30,1],labCol=c("S1","S2","S3","S4"),density.info="none",
             keysize=0.6, key.par=list(mar=c(3.5,1,2.5,0)))
dev.off()
tf_sum_1<-tf_sum_1[order(tf_sum_1[,6], decreasing=T),]
png(filename="TF_activity_top20_ascending.png", units="in", width=3, height=5 , pointsize=8, res=1200)
x<-heatmap.2(as.matrix(tf_sum_1[(1:30),2:5]), col=redgreen(75),
             scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",labRow = tf_sum_1[1:30,1],labCol=c("S1","S2","S3","S4"),density.info="none",
             keysize=0.6, key.par=list(mar=c(3.5,1,2.5,0)))
dev.off()
library(gskb)
data("mm_pathway")
rnaseq_1<-rnaseq_mapped_genes_new[complete.cases(rnaseq_mapped_genes_new), ]
rnaseq_1[,92]<-toupper(rnaseq_1[,85])

path_sum<-as.data.frame(matrix(nrow=892, ncol=5))
colnames(path_sum)<-c("TF","S1","S2","S3","S4")

for(i in 1:892){
  temp<-subset(rnaseq_1,rnaseq_1[,92]%in% (mm_pathway[[i]][3:length(mm_pathway[[i]])]))[,75:78]
  path_sum[i,1]<-mm_pathway[[i]][1]
  path_sum[i,2:5]<-colSums(temp)
}
library(gplots)
png(filename="path_activity.png", units="in", width=3, height=7 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
x<-heatmap.2(as.matrix(path_sum[,2:5]), col=greenred(75),
             scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",labRow = path_sum[,1],labCol=c("S1","S2","S3","S4"),density.info="none",
             keysize=0.6, key.par=list(mar=c(3.5,1,2.5,0)))
dev.off()

path_sum[,6]<-path_sum[,5]/path_sum[,2]
path_sum[,7]<-rowMeans(path_sum[,2:5])
path_sum_1<-subset(path_sum, path_sum[,7]>500)
path_sum_1<-path_sum_1[order(path_sum_1[,6], decreasing=F),]
png(filename="path_activity_decreasing_top20.png", units="in", width=3, height=5 , pointsize=8, res=1200)
x<-heatmap.2(as.matrix(path_sum_1[(1:30),2:5]), col=redgreen(75),
             scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=0.5,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",labRow = path_sum_1[1:30,1],labCol=c("S1","S2","S3","S4"),density.info="none",
             keysize=0.6, key.par=list(mar=c(3.5,1,2.5,0)))
dev.off()


path_sum_1<-path_sum_1[order(path_sum_1[,6], decreasing=T),]
png(filename="path_activity_ascending_top20.png", units="in", width=3, height=5 , pointsize=8, res=1200)
x<-heatmap.2(as.matrix(path_sum_1[(1:30),2:5]), col=redgreen(75),
             scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=0.5,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",labRow = path_sum_1[1:30,1],labCol=c("S1","S2","S3","S4"),density.info="none",
             keysize=0.6, key.par=list(mar=c(3.5,1,2.5,0)))
dev.off()



library(gskb)
data("mm_metabolic")
rnaseq_1<-rnaseq_mapped_genes_new[complete.cases(rnaseq_mapped_genes_new), ]
rnaseq_1[,92]<-toupper(rnaseq_1[,85])

path_sum<-as.data.frame(matrix(nrow=12525, ncol=5))
colnames(path_sum)<-c("GO","S1","S2","S3","S4")

for(i in 1:4603){
  temp<-subset(rnaseq_1,rnaseq_1[,92]%in% (mm_metabolic[[i]][3:length(mm_GO[[i]])]))[,75:78]
  path_sum[i,1]<-mm_metabolic[[i]][1]
  path_sum[i,2:5]<-colSums(temp)
}




path_sum[,6]<-path_sum[,5]/path_sum[,2]
path_sum[,7]<-rowMeans(path_sum[,2:5])
path_sum_1<-subset(path_sum, path_sum[,7]>500)

path_sum_1<-path_sum_1[order(path_sum_1[,6], decreasing=F),]
path_sum_1[1:30,1]
path_sum_1<-path_sum_1[order(path_sum_1[,6], decreasing=T),]
path_sum_1[1:30,1]

