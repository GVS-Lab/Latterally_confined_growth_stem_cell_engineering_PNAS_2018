#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq")
library( "DESeq" )
# prepare count tables and design matrix
counttable<-rnaseq_mapped_genes[,c(50:61)]
rownames(counttable)<-rnaseq_mapped_genes[,1]
head( counttable )
Design = data.frame(row.names = colnames( counttable ),
                    condition = c( "S1", "S2", "S3","S4",
                                   "S1", "S2", "S3","S4",
                                   "S1", "S2", "S3","S4") )
Design
condition = factor( c( "S1", "S2", "S3","S4",
                       "S1", "S2", "S3","S4",
                       "S1", "S2", "S3","S4")  )

cds = newCountDataSet( counttable, condition )
#sizefactors to normalised
cds = estimateSizeFactors( cds )
sizeFactors( cds )
png(filename="DEseq_sizeFactor.png", units="in", width=3, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
barplot(sizeFactors( cds ), ylab="DEseq_sizeFactor", las=2)
box()
dev.off()
norm_counts<-( counts( cds, normalized=TRUE ) )
raw_counts<-( counts( cds, normalized=F ) )
tej<-cor((norm_counts),(raw_counts))

png(filename="DEseq_colsums.png", units="in", width=3, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
barplot(colSums(norm_counts[,1:12]), ylab="summed norm.count in each sample", las=2,cex.axis=0.6)
box()
dev.off()


library("corrplot")
png(filename="norm_vs_raw_DEseq.png", units="in", width=2, height=2 , pointsize=6, res=1200)
corrplot(as.matrix(tej),is.corr=F, tl.col="black",cl.pos="b", cl.ratio=0.3)
dev.off()
tpkm_norm<-as.matrix(rnaseq_mapped_genes[,c(63:74)])
tej<-cor((norm_counts),(tpkm_norm))
library("corrplot")
png(filename="normtpkm_vs_norm_DEseq.png", units="in", width=2, height=2 , pointsize=6, res=1200)
corrplot(as.matrix(tej),is.corr=F, tl.col="black",cl.pos="b", cl.ratio=0.3)
dev.off()

norm_counts<-as.data.frame(norm_counts)
norm_counts[,13]<-rownames(norm_counts)
norm_counts[,14]<-rowMeans(subset(norm_counts, select = c(S1_B1,S1_B2,S1_B3)), na.rm = TRUE)
norm_counts[,15]<-rowMeans(subset(norm_counts, select = c(S2_B2,S2_B2,S3_B3)), na.rm = TRUE)
norm_counts[,16]<-rowMeans(subset(norm_counts, select = c(S3_B3,S3_B2,S3_B3)), na.rm = TRUE)
norm_counts[,17]<-rowMeans(subset(norm_counts, select = c(S4_B3,S4_B2,S4_B3)), na.rm = TRUE)
colnames(norm_counts)[14:17]<-c("S1","S2","S3","S4")

#The package uses the relationship between the data’s variance (or dispersion) and its mean. 
cds = estimateDispersions( cds )
str( fitInfo(cds) )
png(filename="dispersion_deseq.png", units="in", width=2, height=2 , pointsize=5, res=1200)
plotDispEsts( cds )
dev.off()

#calculate the significantly differing genes. Plot the pvalue distribution, MAplot,hist fo foldchange, 
#heatmap of these genes. write the sigificant genes to a csv file. Store the number and names of genes

#S1S2
S1S2 = nbinomTest( cds, "S1","S2" )
png(filename="S1S2_deseq_maplot.png", units="in", width=2, height=2 , pointsize=5, res=1200)
plotMA(S1S2)
dev.off()
png(filename="S1S2_deseq_pvalue.png", units="in", width=2, height=2 , pointsize=5, res=1200)
hist(S1S2$pval, breaks=25,  main="", col="gray")
box()
dev.off()
S1S2Sig = S1S2[ S1S2$padj < 0.1, ]
S1S2Sig<-S1S2Sig[complete.cases(S1S2Sig), ]
S1S2_num<-nrow(S1S2Sig)
S1S2_gene<-(S1S2Sig[,1])
write.csv(S1S2Sig, file="S1S2sig_deseq_norm.csv")
png(filename="S1S2_deseq_foldchange.png", units="in", width=2, height=2 , pointsize=5, res=1200)
hist((S1S2Sig[,6]), col="gray", las=1, cex.axis=0.8,xlab="log2Foldchange",breaks=50,
     main="S1S2DESeqnorm")
box()
dev.off()
temp<-subset(norm_counts,norm_counts[,13] %in% S1S2Sig[,1])


#S1S3
S1S3 = nbinomTest( cds, "S1","S3" )
png(filename="S1S3_deseq_maplot.png", units="in", width=2, height=2 , pointsize=5, res=1200)
plotMA(S1S3)
dev.off()
png(filename="S1S3_deseq_pvalue.png", units="in", width=2, height=2 , pointsize=5, res=1200)
hist(S1S3$pval, breaks=25,  main="", col="gray")
box()
dev.off()
S1S3Sig = S1S3[ S1S3$padj < 0.1, ]
S1S3Sig<-S1S3Sig[complete.cases(S1S3Sig), ]
S1S3_num<-nrow(S1S3Sig)
S1S3_gene<-(S1S3Sig[,1])
write.csv(S1S3Sig, file="S1S3sig_deseq_norm.csv")
png(filename="S1S3_deseq_foldchange.png", units="in", width=2, height=2 , pointsize=5, res=1200)
hist((S1S3Sig[,6]), col="gray", las=1, cex.axis=0.8,xlab="log2Foldchange",breaks=50,
     main="S1S3DESeqnorm")
box()
dev.off()
temp<-subset(norm_counts,norm_counts[,13] %in% S1S3Sig[,1])

#S1S4
S1S4 = nbinomTest( cds, "S1","S4" )
png(filename="S1S4_deseq_maplot.png", units="in", width=2, height=2 , pointsize=5, res=1200)
plotMA(S1S4)
dev.off()
png(filename="S1S4_deseq_pvalue.png", units="in", width=2, height=2 , pointsize=5, res=1200)
hist(S1S4$pval, breaks=25,  main="", col="gray")
box()
dev.off()
S1S4Sig = S1S4[ S1S4$padj < 0.1, ]
S1S4Sig<-S1S4Sig[complete.cases(S1S4Sig), ]
S1S4_num<-nrow(S1S4Sig)
S1S4_gene<-(S1S4Sig[,1])
write.csv(S1S4Sig, file="S1S4sig_deseq_norm.csv")
png(filename="S1S4_deseq_foldchange.png", units="in", width=2, height=2 , pointsize=5, res=1200)
hist((S1S4Sig[,6]), col="gray", las=1, cex.axis=0.8,xlab="log2Foldchange",breaks=50,
     main="S1S4DESeqnorm")
box()
dev.off()
temp<-subset(norm_counts,norm_counts[,13] %in% S1S4Sig[,1])

#S2S3
S2S3 = nbinomTest( cds, "S2","S3" )
png(filename="S2S3_deseq_maplot.png", units="in", width=2, height=2 , pointsize=5, res=1200)
plotMA(S2S3)
dev.off()
png(filename="S2S3_deseq_pvalue.png", units="in", width=2, height=2 , pointsize=5, res=1200)
hist(S2S3$pval, breaks=25,  main="", col="gray")
box()
dev.off()
S2S3Sig = S2S3[ S2S3$padj < 0.1, ]
S2S3Sig<-S2S3Sig[complete.cases(S2S3Sig), ]
S2S3_num<-nrow(S2S3Sig)
S2S3_gene<-(S2S3Sig[,1])
write.csv(S2S3Sig, file="S2S3sig_deseq_norm.csv")
png(filename="S2S3_deseq_foldchange.png", units="in", width=2, height=2 , pointsize=5, res=1200)
hist((S2S3Sig[,6]), col="gray", las=1, cex.axis=0.8,xlab="log2Foldchange",breaks=50,
     main="S2S3DESeqnorm")
box()
dev.off()
temp<-subset(norm_counts,norm_counts[,13] %in% S2S3Sig[,1])

#S2S4
S2S4 = nbinomTest( cds, "S2","S4" )
png(filename="S2S4_deseq_maplot.png", units="in", width=2, height=2 , pointsize=5, res=1200)
plotMA(S2S4)
dev.off()
png(filename="S2S4_deseq_pvalue.png", units="in", width=2, height=2 , pointsize=5, res=1200)
hist(S2S4$pval, breaks=25,  main="", col="gray")
box()
dev.off()
S2S4Sig = S2S4[ S2S4$padj < 0.1, ]
S2S4Sig<-S2S4Sig[complete.cases(S2S4Sig), ]
S2S4_num<-nrow(S2S4Sig)
S2S4_gene<-(S2S4Sig[,1])
write.csv(S2S4Sig, file="S2S4sig_deseq_norm.csv")
png(filename="S2S4_deseq_foldchange.png", units="in", width=2, height=2 , pointsize=5, res=1200)
hist((S2S4Sig[,6]), col="gray", las=1, cex.axis=0.8,xlab="log2Foldchange",breaks=50,
     main="S2S4DESeqnorm")
box()
dev.off()
temp<-subset(norm_counts,norm_counts[,13] %in% S2S4Sig[,1])

#S3S4
S3S4 = nbinomTest( cds, "S3","S4" )
png(filename="S3S4_deseq_maplot.png", units="in", width=2, height=2 , pointsize=5, res=1200)
plotMA(S3S4)
dev.off()
png(filename="S3S4_deseq_pvalue.png", units="in", width=2, height=2 , pointsize=5, res=1200)
hist(S3S4$pval, breaks=25,  main="", col="gray")
box()
dev.off()
S3S4Sig = S3S4[ S3S4$padj < 0.1, ]
S3S4Sig<-S3S4Sig[complete.cases(S3S4Sig), ]
S3S4_num<-nrow(S3S4Sig)
S3S4_gene<-(S3S4Sig[,1])
write.csv(S3S4Sig, file="S3S4sig_deseq_norm.csv")
png(filename="S3S4_deseq_foldchange.png", units="in", width=2, height=2 , pointsize=5, res=1200)
hist((S3S4Sig[,6]), col="gray", las=1, cex.axis=0.8,xlab="log2Foldchange",breaks=50,
     main="S3S4DESeqnorm")
box()
dev.off()
temp<-subset(norm_counts,norm_counts[,13] %in% S3S4Sig[,1])


DE_norm_number = matrix( c(0,S1S2_num,S1S3_num,S1S4_num,
        S1S2_num,0,S2S3_num,S2S4_num,
        S1S3_num,S2S3_num,0,S3S4_num,
        S1S4_num,S2S4_num,S3S4_num,0), 
     nrow=4, ncol=4) 
DE_norm_number<-as.data.frame(DE_norm_number)
colnames(DE_norm_number)<-c("S1","S2","S3","S4")
rownames(DE_norm_number)<-c("S1","S2","S3","S4")
library(corrplot)
png(filename="Number_DE_genes_DESEQ_norm.png", units="in", width=2, height=2 , pointsize=5, res=1200)
corrplot(as.matrix(DE_norm_number),is.corr=F, tl.col="black",cl.pos="r", cl.ratio=0.3, method="color")
dev.off()

#ME
a1<-length(subset(S1S2_gene,S1S2_gene %in% candidate_genes[,1]))
a2<-length(subset(S1S3_gene,S1S3_gene %in% candidate_genes[,1]))
a3<-length(subset(S1S4_gene,S1S4_gene %in% candidate_genes[,1]))
b1<-length(subset(S2S3_gene,S2S3_gene %in% candidate_genes[,1]))
b2<-length(subset(S2S4_gene,S2S4_gene %in% candidate_genes[,1]))
c1<-length(subset(S3S4_gene,S2S4_gene %in% candidate_genes[,1]))
DE_norm_number_me = matrix( c(0,a1,a2,a3,
                              a1,0,b1,b2,
                              a2,b1,0,c1,
                              a3,b2,c1,0), 
                         nrow=4, ncol=4) 
DE_norm_number_me<-as.data.frame(DE_norm_number_me)
colnames(DE_norm_number_me)<-c("S1","S2","S3","S4")
rownames(DE_norm_number_me)<-c("S1","S2","S3","S4")
library(corrplot)
png(filename="Number_DE_genes_DESEQ_norm_me.png", units="in", width=2, height=2 , pointsize=5, res=1200)
corrplot(as.matrix(DE_norm_number_me),is.corr=F, tl.col="black",cl.pos="r", cl.ratio=0.3, method="color")
dev.off()

#ES
a1<-length(subset(S1S2_gene,S1S2_gene %in% candidate_genes[,2]))
a2<-length(subset(S1S3_gene,S1S3_gene %in% candidate_genes[,2]))
a3<-length(subset(S1S4_gene,S1S4_gene %in% candidate_genes[,2]))
b1<-length(subset(S2S3_gene,S2S3_gene %in% candidate_genes[,2]))
b2<-length(subset(S2S4_gene,S2S4_gene %in% candidate_genes[,2]))
c1<-length(subset(S3S4_gene,S2S4_gene %in% candidate_genes[,2]))
DE_norm_number_es = matrix( c(0,a1,a2,a3,
                              a1,0,b1,b2,
                              a2,b1,0,c1,
                              a3,b2,c1,0), 
                            nrow=4, ncol=4) 
DE_norm_number_es<-as.data.frame(DE_norm_number_es)
colnames(DE_norm_number_es)<-c("S1","S2","S3","S4")
rownames(DE_norm_number_es)<-c("S1","S2","S3","S4")
library(corrplot)
png(filename="Number_DE_genes_DESEQ_norm_es.png", units="in", width=2, height=2 , pointsize=5, res=1200)
corrplot(as.matrix(DE_norm_number_es),is.corr=F, tl.col="black",cl.pos="r", cl.ratio=0.3, method="color")
dev.off()

#iPSC
a1<-length(subset(S1S2_gene,S1S2_gene %in% candidate_genes[,3]))
a2<-length(subset(S1S3_gene,S1S3_gene %in% candidate_genes[,3]))
a3<-length(subset(S1S4_gene,S1S4_gene %in% candidate_genes[,3]))
b1<-length(subset(S2S3_gene,S2S3_gene %in% candidate_genes[,3]))
b2<-length(subset(S2S4_gene,S2S4_gene %in% candidate_genes[,3]))
c1<-length(subset(S3S4_gene,S2S4_gene %in% candidate_genes[,3]))
DE_norm_number_ip = matrix( c(0,a1,a2,a3,
                              a1,0,b1,b2,
                              a2,b1,0,c1,
                              a3,b2,c1,0), 
                            nrow=4, ncol=4) 
DE_norm_number_ip<-as.data.frame(DE_norm_number_ip)
colnames(DE_norm_number_ip)<-c("S1","S2","S3","S4")
rownames(DE_norm_number_ip)<-c("S1","S2","S3","S4")
library(corrplot)
png(filename="Number_DE_genes_DESEQ_norm_ip.png", units="in", width=2, height=2 , pointsize=5, res=1200)
corrplot(as.matrix(DE_norm_number_ip),is.corr=F, tl.col="black",cl.pos="r", cl.ratio=0.3, method="color")
dev.off()



# log transform unnormalised
exon <- log2(norm_counts[, 1:12])
is.na(exon) <- do.call(cbind,lapply(exon, is.infinite)) 
exon.pca <- prcomp(na.omit(exon),scale= F) 
t<-as.data.frame(exon.pca$rotation)
png(filename="pca_DESeq_norm_log_scaled.png", units="in", width=2, height=2 , pointsize=5, res=1200)
plot(t[,2]~t[,1],pch=18,col=c("purple","blue","red","green4"), cex.axis=1,las=1, cex=1.2, 
     ylab="PC2",xlab="PC1")
legend(-0.2805761,0.08247002,c("S1","S2","S3","S4"), lty=c(1,1), lwd=c(2.5,2.5),col=c("purple","blue","red","green4")) 
dev.off()

library(corrplot)
png(filename="Sample_distances_log_DESeqnorm.png", units="in", width=4, height=4 , pointsize=5, res=1200)
par(mfrow=c(2,2))
t<-as.matrix(dist(t(exon), upper=T),method = "euclidean")
corrplot(t, is.corr=F, tl.col="black",method="color")
plot(hclust(dist(t(exon), method="euclidean")))
t<-as.matrix(dist(t(exon), upper=F,method = "maximum"))
corrplot(t, is.corr=F, tl.col="black",method="color")
plot(hclust(dist(t(exon), method="maximum")))
dev.off()

png(filename="DESeq_norm_histogram.png", units="in", width=2, height=2, pointsize=4, res=1200)
par(font.axis=2,font.lab=2)
par(mfrow=c(2,2))
hist(log(norm_counts[,14], base=2), col="gray", las=1, cex.axis=0.8,xlab="DESeq_norm",main="S1",xlim=c(-3,20))
box()
hist(log(norm_counts[,15], base=2), col="gray", las=1, cex.axis=0.8,xlab="DESeq_norm",main="S2",xlim=c(-3,20))
box()
hist(log(norm_counts[,16], base=2), col="gray", las=1, cex.axis=0.8,xlab="DESeq_norm",main="S3",xlim=c(-3,20))
box()
hist(log(norm_counts[,17], base=2), col="gray", las=1, cex.axis=0.8,xlab="DESeq_norm",main="S4",xlim=c(-3,20))
box()
dev.off()

png(filename="DESeq_norm_density.png", units="in", width=2, height=2, pointsize=4, res=1200)
par(font.axis=2,font.lab=2)
d <- density(log(norm_counts[,17], base=2)) # returns the density data 
plot(d, col="green4", las=1, lwd=1, main="")
d <- density(log(norm_counts[,16], base=2)) # returns the density data 
lines(d,col="red",lwd=1)
d <- density(log(norm_counts[,15], base=2)) # returns the density data 
lines(d,col="blue",lwd=1)
d <- density(log(norm_counts[,14], base=2)) # returns the density data 
lines(d,col="purple",lwd=1)
legend(12,0.04176209,c("S1","S2","S3","S4"), lty=c(1,1), lwd=c(2.5,2.5),col=c("purple","blue","red","green4")) 
dev.off()





#################################################################################################

Me<-subset(norm_counts,norm_counts[,13] %in% candidate_genes[,1])
ES<-subset(norm_counts,norm_counts[,13] %in% candidate_genes[,2])
iPSC<-subset(norm_counts,norm_counts[,13] %in% candidate_genes[,3])

library("gplots")
png(filename="Me_heatmap_DESeq_count.png", units="in", width=2, height=4 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
x<-heatmap.2(as.matrix(Me[,14:17]), col=greenred(75),
             scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=0.9,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",labRow = Me[,13],labCol=c("3 Hour","3 Days","6 Days","10 Days"),density.info="none",
             keysize=0.6, key.par=list(mar=c(3.5,1,2.5,0)))
dev.off()
png(filename="Me_Zscore_summed_barplot_DESeq.png", units="in", width=1.5, height=1.5 , pointsize=5, res=1200)
par(font.axis=2,font.lab=2)
me_zscore<-as.data.frame(x$carpet)
t<-rowSums(me_zscore)
t_up<-max(t)+5
t_dw<-min(t)-5
barplot(t,las=1,ylab="Summed zscore",names =c("S1","S2","S3","S4"),ylim=c(t_dw,t_up))
box()
abline(h=0)
dev.off()
png(filename="Me_Expression_summed_barplot_DESeq.png", units="in", width=1.5, height=1.5 , pointsize=5, res=1200)
par(font.axis=2,font.lab=2)
x<-colSums(Me[,14:17])
x_up<-max(x)+200
barplot(x,las=1,ylab="Summed Tpkm",names =c("S1","S2","S3","S4"),ylim=c(0,x_up),cex.axis=0.7)
box()
dev.off()


png(filename="ES_heatmap_DESeq_count.png", units="in", width=2, height=6 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
x<-heatmap.2(as.matrix(ES[,14:17]), col=greenred(75),
             scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=0.9,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",labRow = ES[,13],labCol=c("3 Hour","3 Days","6 Days","10 Days"),density.info="none",
             keysize=0.6, key.par=list(mar=c(3.5,1,6,0)))
dev.off()
png(filename="ES_Zscore_summed_barplot_DESeq.png", units="in", width=1.5, height=1.5 , pointsize=5, res=1200)
par(font.axis=2,font.lab=2)
es_zscore<-as.data.frame(x$carpet)
t<-rowSums(es_zscore)
t_up<-max(t)+5
t_dw<-min(t)-5
barplot(t,las=1,ylab="Summed zscore",names =c("S1","S2","S3","S4"),ylim=c(t_dw,t_up))
box()
abline(h=0)
dev.off()
png(filename="ES_Expression_summed_barplot_DESeq.png", units="in", width=1.5, height=1.5 , pointsize=5, res=1200)
par(font.axis=2,font.lab=2)
x<-colSums(ES[,14:17])
x_up<-max(x)+200
barplot(x,las=1,ylab="Summed Tpkm",names =c("S1","S2","S3","S4"),ylim=c(0,x_up))
box()
dev.off()


png(filename="iPSC_heatmap_DESeq_count.png", units="in", width=2, height=6 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
x<-heatmap.2(as.matrix(iPSC[,14:17]), col=greenred(75),
             scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=0.9,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",labRow = iPSC[,13],labCol=c("3 Hour","3 Days","6 Days","10 Days"),density.info="none",
             keysize=0.6, key.par=list(mar=c(3.5,1,6,0)))
dev.off()
png(filename="iPSC_Zscore_summed_barplot_DESeq.png", units="in", width=1.5, height=1.5 , pointsize=5, res=1200)
iPSC_zscore<-as.data.frame(x$carpet)
t<-rowSums(iPSC_zscore, na.rm = T)
t_up<-max(t)+5
t_dw<-min(t)-5
barplot(t,las=1,ylab="Summed zscore",names =c("S1","S2","S3","S4"),ylim=c(t_dw,t_up))
box()
abline(h=0)
dev.off()
png(filename="iPSC_Expression_summed_barplot_DESeq.png", units="in", width=1.5, height=1.5 , pointsize=5, res=1200)
par(font.axis=2,font.lab=2)
x<-colSums(iPSC[,14:17])
x_up<-max(x)+200
barplot(x,las=1,ylab="Summed Tpkm",names =c("S1","S2","S3","S4"),ylim=c(0,x_up))
box()
dev.off()
# top genes
n<-round(nrow(norm_counts)*0.1,0)
S1_topgenes<-norm_counts[order(norm_counts[,14]),]
S1_topgenes<-tail(S1_topgenes,n)
S2_topgenes<-norm_counts[order(norm_counts[,15]),]
S2_topgenes<-tail(S2_topgenes,n)
S3_topgenes<-norm_counts[order(norm_counts[,16]),]
S3_topgenes<-tail(S3_topgenes,n)
S4_topgenes<-norm_counts[order(norm_counts[,17]),]
S4_topgenes<-tail(S4_topgenes,n)
m<-nrow(subset(S1_topgenes,S1_topgenes[,13] %in% candidate_genes[,1]))
e<-nrow(subset(S1_topgenes,S1_topgenes[,13] %in% candidate_genes[,2]))
i<-nrow(subset(S1_topgenes,S1_topgenes[,13] %in% candidate_genes[,3]))
S1_fractions<-c(m,e,i)
m<-nrow(subset(S2_topgenes,S2_topgenes[,13] %in% candidate_genes[,1]))
e<-nrow(subset(S2_topgenes,S2_topgenes[,13] %in% candidate_genes[,2]))
i<-nrow(subset(S2_topgenes,S2_topgenes[,13] %in% candidate_genes[,3]))
S2_fractions<-c(m,e,i)
m<-nrow(subset(S3_topgenes,S3_topgenes[,13] %in% candidate_genes[,1]))
e<-nrow(subset(S3_topgenes,S3_topgenes[,13] %in% candidate_genes[,2]))
i<-nrow(subset(S3_topgenes,S3_topgenes[,13] %in% candidate_genes[,3]))
S3_fractions<-c(m,e,i)
m<-nrow(subset(S4_topgenes,S4_topgenes[,13] %in% candidate_genes[,1]))
e<-nrow(subset(S4_topgenes,S4_topgenes[,13] %in% candidate_genes[,2]))
i<-nrow(subset(S4_topgenes,S4_topgenes[,13] %in% candidate_genes[,3]))
S4_fractions<-c(m,e,i)
top_10_percent<-as.data.frame(S1_fractions)
top_10_percent[,2]<-S2_fractions
top_10_percent[,3]<-S3_fractions
top_10_percent[,4]<-S4_fractions
row.names(top_10_percent)<-c("Mesenchymal","ES","iPSC")
colnames(top_10_percent)<-c("S1","S2","S3","S4")
library("corrplot")
png(filename="top_genes_DESEQ.png", units="in", width=2, height=2 , pointsize=6, res=1200)
corrplot(as.matrix(top_10_percent),is.corr=F, tl.col="black",cl.pos="b", cl.ratio=0.3)
dev.off()
temp<-top_10_percent
colnames(temp)<-c("3 Hour","3 Days","6 Days","10 Days")
png(filename="top_genes_final_DESEQ.png", units="in", width=3, height=3 , pointsize=8, res=1200)
par(font.axis=2,font.lab=2)
corrplot(as.matrix(temp),is.corr=F, tl.col="black",cl.pos="b", cl.ratio=0.3,method="color",addgrid.col="black")
dev.off()

#epitect
rna_epi<-subset(norm_counts,norm_counts[,13] %in% epitect[,1])
library("gplots")
png(filename="epitect_rna_DESEQnorm_heatmap.png", units="in", width=2, height=5 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
x<-heatmap.2(as.matrix(rna_epi[,14:17]), col=greenred(75),
             scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=0.9,
             trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="",cex.axis=0.7,
             ylab = "",labRow = rna_epi[,13],labCol=c("3 Hour","3 Days","6 Days","10 Days"),density.info="none",
             keysize=0.6, key.par=list(mar=c(3.5,1,2.5,0)))
dev.off()

te<-subset(epitect,epitect[,6]=="Epithelial")[,1]
te_rna<-subset(norm_counts,norm_counts[,13] %in% te)[,c(1,14:17)]
te_rna[,6]<-"Epithelial"
tm<-subset(epitect,epitect[,6]=="Mesenchymal")[,1]
tm_rna<-subset(norm_counts,norm_counts[,13] %in% tm)[,c(1,14:17)]
tm_rna[,6]<-"Mesenchymal"
tp<-subset(epitect,epitect[,6]=="Proliferation")[,1]
tp_rna<-subset(norm_counts,norm_counts[,13] %in% tp)[,c(1,14:17)]
tp_rna[,6]<-"Proliferation"
tr<-subset(epitect,epitect[,6]=="Reprogramming")[,1]
tr_rna<-subset(norm_counts,norm_counts[,13] %in% tr)[,c(1,14:17)]
tr_rna[,6]<-"Reprogramming"
rna_epitect_groups<-rbind(te_rna,tm_rna,tp_rna,tr_rna)

temp_rna<-aggregate(. ~ V6, data = rna_epitect_groups[,-1],sum)
temp_epitect<-aggregate(. ~ X.1, data = epitect[,-1],sum)
png(filename="Epi_RNA_genes_sum_DESEQ.png", units="in", width=3, height=6 , pointsize=6, res=1200)
par(mfrow=c(4,2),font.axis=2,font.lab=2)
k<-as.matrix(temp_rna[1,-1])
barplot(k, las=1, horiz = T, main="Epithelial_RNAseq")
box()
k<-as.matrix(temp_epitect[1,-1])
barplot(k, las=1,horiz = T,main="Epithelial_Epitect")
box()
k<-as.matrix(temp_rna[2,-1])
barplot(k, las=1, horiz = T, main="Mesenchymal_RNAseq")
box()
k<-as.matrix(temp_epitect[2,-1])
barplot(k, las=1,horiz = T,main="Mesenchymal_Epitect")
box()
k<-as.matrix(temp_rna[3,-1])
barplot(k, las=1, horiz = T, main="Proliferation_RNAseq")
box()
k<-as.matrix(temp_epitect[3,-1])
barplot(k, las=1,horiz = T,main="Proliferation_Epitect")
box()
k<-as.matrix(temp_rna[4,-1])
barplot(k, las=1, horiz = T, main="Reprogramming_RNAseq")
box()
k<-as.matrix(temp_epitect[4,-1])
barplot(k, las=1,horiz = T,main="Reprogramming_Epitect")
box()
dev.off()

temp_rna<-aggregate(. ~ V6, data = rna_epitect_groups[,-1],mean)
temp_epitect<-aggregate(. ~ X.1, data = epitect[,-1],mean)
png(filename="Epi_RNA_genes_mean_DESEQ.png", units="in", width=3, height=6 , pointsize=6, res=1200)
par(mfrow=c(4,2),font.axis=2,font.lab=2)
k<-as.matrix(temp_rna[1,-1])
barplot(k, las=1, horiz = T, main="Epithelial_RNAseq")
box()
k<-as.matrix(temp_epitect[1,-1])
barplot(k, las=1,horiz = T,main="Epithelial_Epitect")
box()
k<-as.matrix(temp_rna[2,-1])
barplot(k, las=1, horiz = T, main="Mesenchymal_RNAseq")
box()
k<-as.matrix(temp_epitect[2,-1])
barplot(k, las=1,horiz = T,main="Mesenchymal_Epitect")
box()
k<-as.matrix(temp_rna[3,-1])
barplot(k, las=1, horiz = T, main="Proliferation_RNAseq")
box()
k<-as.matrix(temp_epitect[3,-1])
barplot(k, las=1,horiz = T,main="Proliferation_Epitect")
box()
k<-as.matrix(temp_rna[4,-1])
barplot(k, las=1, horiz = T, main="Reprogramming_RNAseq")
box()
k<-as.matrix(temp_epitect[4,-1])
barplot(k, las=1,horiz = T,main="Reprogramming_Epitect")
box()
dev.off()



