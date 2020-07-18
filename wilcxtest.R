d<-as.matrix(rnaseq_mapped_genes[,c(1,63:74)])
row.names(d)<-d[,1]

wil_S1S2<-apply(d, 1, function(row) unlist(wilcox.test(as.numeric(row[c(2,6,10)]),as.numeric(row[c(3,7,11)]))[c("p.value","statistic")]))
wil_S1S2<-as.data.frame(t(wil_S1S2))
wil_S1S2[,3]<-as.character(rownames(wil_S1S2))
S1S2_num<-nrow(subset(wil_S1S2,wil_S1S2[,1]<0.1))
S1S2_gene<-subset(wil_S1S2,wil_S1S2[,1]<0.1)[,3]

wil_S1S3<-apply(d, 1, function(row) unlist(wilcox.test(as.numeric(row[c(2,6,10)]),as.numeric(row[c(4,8,12)]))[c("p.value","statistic")]))
wil_S1S3<-as.data.frame(t(wil_S1S3))
wil_S1S3[,3]<-as.character(rownames(wil_S1S3))
S1S3_num<-nrow(subset(wil_S1S3,wil_S1S3[,1]<0.1))
S1S3_gene<-subset(wil_S1S3,wil_S1S3[,1]<0.1)[,3]

wil_S1S4<-apply(d, 1, function(row) unlist(wilcox.test(as.numeric(row[c(2,6,10)]),as.numeric(row[c(5,9,13)]))[c("p.value","statistic")]))
wil_S1S4<-as.data.frame(t(wil_S1S4))
wil_S1S4[,3]<-as.character(rownames(wil_S1S4))
S1S4_num<-nrow(subset(wil_S1S4,wil_S1S4[,1]<0.1))
S1S4_gene<-subset(wil_S1S4,wil_S1S4[,1]<0.1)[,3]

wil_S2S3<-apply(d, 1, function(row) unlist(wilcox.test(as.numeric(row[c(3,7,11)]),as.numeric(row[c(4,8,12)]))[c("p.value","statistic")]))
wil_S2S3<-as.data.frame(t(wil_S2S3))
wil_S2S3[,3]<-as.character(rownames(wil_S2S3))
S2S3_num<-nrow(subset(wil_S2S3,wil_S2S3[,1]<0.1))
S2S3_gene<-subset(wil_S2S3,wil_S2S3[,1]<0.1)[,3]

wil_S2S4<-apply(d, 1, function(row) unlist(wilcox.test(as.numeric(row[c(3,7,11)]),as.numeric(row[c(5,9,13)]))[c("p.value","statistic")]))
wil_S2S4<-as.data.frame(t(wil_S2S4))
wil_S2S4[,3]<-as.character(rownames(wil_S2S4))
S2S4_num<-nrow(subset(wil_S2S4,wil_S2S4[,1]<0.1))
S2S4_gene<-subset(wil_S2S4,wil_S2S4[,1]<0.1)[,3]

wil_S3S4<-apply(d, 1, function(row) unlist(wilcox.test(as.numeric(row[c(4,8,12)]),as.numeric(row[c(5,9,13)]))[c("p.value","statistic")]))
wil_S3S4<-as.data.frame(t(wil_S3S4))
wil_S3S4[,3]<-as.character(rownames(wil_S3S4))
S3S4_num<-nrow(subset(wil_S3S4,wil_S3S4[,1]<0.1))
S3S4_gene<-subset(wil_S3S4,wil_S3S4[,1]<0.1)[,3]

Wilcox_number = matrix( c(0,S1S2_num,S1S3_num,S1S4_num,
                           S1S2_num,0,S2S3_num,S2S4_num,
                           S1S3_num,S2S3_num,0,S3S4_num,
                           S1S4_num,S2S4_num,S3S4_num,0), 
                         nrow=4, ncol=4) 
Wilcox_number<-as.data.frame(Wilcox_number)
colnames(Wilcox_number)<-c("S1","S2","S3","S4")
rownames(Wilcox_number)<-c("S1","S2","S3","S4")
library(corrplot)
png(filename="Number_DE_genes_Wilcox_number.png", units="in", width=2, height=2 , pointsize=5, res=1200)
corrplot(as.matrix(Wilcox_number),is.corr=F, tl.col="black",cl.pos="r", cl.ratio=0.3, method="color")
dev.off()


a1<-length(subset(S1S2_gene,S1S2_gene %in% candidate_genes[,1]))
a2<-length(subset(S1S3_gene,S1S3_gene %in% candidate_genes[,1]))
a3<-length(subset(S1S4_gene,S1S4_gene %in% candidate_genes[,1]))
b1<-length(subset(S2S3_gene,S2S3_gene %in% candidate_genes[,1]))
b2<-length(subset(S2S4_gene,S2S4_gene %in% candidate_genes[,1]))
c1<-length(subset(S3S4_gene,S2S4_gene %in% candidate_genes[,1]))
Wilcox_number_me = matrix( c(0,a1,a2,a3,
                              a1,0,b1,b2,
                              a2,b1,0,c1,
                              a3,b2,c1,0), 
                            nrow=4, ncol=4) 
Wilcox_number_me<-as.data.frame(Wilcox_number_me)
colnames(Wilcox_number_me)<-c("S1","S2","S3","S4")
rownames(Wilcox_number_me)<-c("S1","S2","S3","S4")
library(corrplot)
png(filename="Number_DE_genes_Wilcox_number_me.png", units="in", width=2, height=2 , pointsize=5, res=1200)
corrplot(as.matrix(Wilcox_number_me),is.corr=F, tl.col="black",cl.pos="r", cl.ratio=0.3, method="color")
dev.off()

#ES
a1<-length(subset(S1S2_gene,S1S2_gene %in% candidate_genes[,2]))
a2<-length(subset(S1S3_gene,S1S3_gene %in% candidate_genes[,2]))
a3<-length(subset(S1S4_gene,S1S4_gene %in% candidate_genes[,2]))
b1<-length(subset(S2S3_gene,S2S3_gene %in% candidate_genes[,2]))
b2<-length(subset(S2S4_gene,S2S4_gene %in% candidate_genes[,2]))
c1<-length(subset(S3S4_gene,S2S4_gene %in% candidate_genes[,2]))
Wilcox_number_es = matrix( c(0,a1,a2,a3,
                              a1,0,b1,b2,
                              a2,b1,0,c1,
                              a3,b2,c1,0), 
                            nrow=4, ncol=4) 
Wilcox_number_es<-as.data.frame(Wilcox_number_es)
colnames(Wilcox_number_es)<-c("S1","S2","S3","S4")
rownames(Wilcox_number_es)<-c("S1","S2","S3","S4")
library(corrplot)
png(filename="Number_DE_genes_Wilcox_number_es.png", units="in", width=2, height=2 , pointsize=5, res=1200)
corrplot(as.matrix(Wilcox_number_es),is.corr=F, tl.col="black",cl.pos="r", cl.ratio=0.3, method="color")
dev.off()

#iPSC
a1<-length(subset(S1S2_gene,S1S2_gene %in% candidate_genes[,3]))
a2<-length(subset(S1S3_gene,S1S3_gene %in% candidate_genes[,3]))
a3<-length(subset(S1S4_gene,S1S4_gene %in% candidate_genes[,3]))
b1<-length(subset(S2S3_gene,S2S3_gene %in% candidate_genes[,3]))
b2<-length(subset(S2S4_gene,S2S4_gene %in% candidate_genes[,3]))
c1<-length(subset(S3S4_gene,S2S4_gene %in% candidate_genes[,3]))
Wilcox_number_ip = matrix( c(0,a1,a2,a3,
                              a1,0,b1,b2,
                              a2,b1,0,c1,
                              a3,b2,c1,0), 
                            nrow=4, ncol=4) 
Wilcox_number_ip<-as.data.frame(Wilcox_number_ip)
colnames(Wilcox_number_ip)<-c("S1","S2","S3","S4")
rownames(Wilcox_number_ip)<-c("S1","S2","S3","S4")
library(corrplot)
png(filename="Number_DE_genes_Wilcox_number_ip.png", units="in", width=2, height=2 , pointsize=5, res=1200)
corrplot(as.matrix(Wilcox_number_ip),is.corr=F, tl.col="black",cl.pos="r", cl.ratio=0.3, method="color")
dev.off()
