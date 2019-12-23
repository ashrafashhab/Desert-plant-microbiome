#!/bin/bash
OLD_IFS=$IFS; 
IFS=$'\n';
# PLEASE READ BEFORE USE:

# for using this script make sure you download the follwing
# 1- bowtie for removing Phi-X sequences 
# 2- PEAR for merging both R1 and R2 reads
# 3- Mothur for additional quality control step
# 4- Qiime for making OTU classified sequences table
# 5- R for data analysis and making graph (if you wish, the script will run withouit)
# 6- rdp_classifier-2.2.jar for OTUS classification

# also make sure you have a map (text) file (MOST_core.files) in the folder where the row_files located formated as: short_file_name, R1 and R2, space seperated with no header 
# also make sure you have a map (text) file (MOST_core.files) in the files where the mothurQC will be created formated as: #sampleID(S001), InputFileName (mothurQC_S001.fasta (as will be created)), Description. tab seperated with header 

#Filting OTU table
rich_dense_biom  = system.file("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/rdp_F649_R889.bac.filtered.sample.biom", "rich_dense_otu_table.biom",  package="phyloseq")
rich_sparse_biom = system.file("extdata", "rich_sparse_otu_table.biom", package="phyloseq")


otufile = system.file("extdata", "/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/core_all.txt.gz", package="phyloseq")
mapfile = system.file("extdata", "/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/F649_R889_map.files", package="phyloseq")
trefile = system.file("extdata", "/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/rep_set.tre", package="phyloseq")
rs_file = system.file("extdata", "/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/rep_set/seqs_rep_set.fasta", package="phyloseq")
qiimedata = import_qiime(otufile, mapfile, trefile, rs_file)


source activate qiime1 
filter_samples_from_otu_table.py \
-i /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/rdp_F649_R889.bac.filtered.sample.biom \
-m /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/F649_R889_map.files \
-s INC:core_all \
-o /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/core_all.biom



biom convert \
-i /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/core_all.biom \
-o /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/core_all.txt \
--to-json \
--table-type='OTU Table'






sed 1d /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/core_all.txt>L6_0.txt;
sed 's/ /_/g' L6_0.txt>L6_1.txt;
sed s/'#OTU_ID'/'OTU_ID'/ L6_1.txt>L6_2.txt;
sed 's/;/\t/g' L6_2.txt>L6_3.txt;
sed s/'ConsensusLineage'/'Kingdome\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies'/ L6_3.txt>L6_4.txt;
sed 's/-/_/g' L6_4.txt>L6_5.txt;
sed 's/\[//g' L6_5.txt>L6_6.txt;
sed 's/\]//g' L6_6.txt>L6_7.txt;
sed 's/_s__//g' L6_7.txt>L6_8.txt;
sed 's/_g__//g' L6_8.txt>L6_9.txt;
sed 's/_f__//g' L6_9.txt>L6_10.txt;
sed 's/_o__//g' L6_10.txt>L6_11.txt;
sed 's/_c__//g' L6_11.txt>L6_12.txt;
sed 's/_//g' L6_12.txt>L6_13.txt;
sed 's/k__//g' L6_13.txt>/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/core_all_R.txt;
rm L6_0.txt;
rm L6_1.txt;
rm L6_2.txt;
rm L6_3.txt;
rm L6_4.txt;
rm L6_5.txt;
rm L6_6.txt;
rm L6_7.txt;
rm L6_8.txt;
rm L6_9.txt;
rm L6_10.txt;
rm L6_11.txt;
rm L6_12.txt;
rm L6_13.txt;
rm /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/core_all.txt;
source deactivate qiime1 
#################################################################endo vs exo 
R --save <<RSCRIPT 
library(MASS)
library(vegan)
library(ggplot2)
library(plotrix)
library("plot3D")

otu<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/core_all_R.txt",header=T,sep="\t",fill=T,row.names=1)
otu2<-otu[,1:(dim(otu)[2]-7)]
otus<-otu2
t_otus <- as.data.frame(t(otus))
t_otus<-decostand(t_otus,method="total")

map<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/F649_R889_map_r2.files",row.names=1,sep="\t",header=T)
map2<-map[map$INC=="core_all",]
x<-subset(t_otus[rownames(map2),])
p<-cbind(map2,x)
write.table(p,"/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",sep="\t",quote=FALSE)
p<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",header=T,row.names=1,sep="\t")

p<-p[p$Face=="South",]
write.table(p,"/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",sep="\t",quote=FALSE)
p<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",header=T,row.names=1,sep="\t")
pp<-p[,((dim(map2)[2])+1):dim(p)[2]]
ppp<-pp[,!!colSums(pp)]
NMDS_bray<-metaMDS(ppp,distance="bray",k=2,trymax=10000,autotransform=F)
l<-anosim(dat=(as.matrix((vegdist(ppp,distance="bray")))),grouping=p$Endo_exo)
zz<-scores(NMDS_bray,display="site")

p_endo<-p[p$Endo_exo=="Endosphere",]
write.table(p_endo,"/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",sep="\t",quote=FALSE)
p_endo<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",header=T,row.names=1,sep="\t")
pp_endo<-p_endo[,((dim(map2)[2])+1):dim(p_endo)[2]]
ppp_endo<-pp_endo[,!!colSums(pp_endo)]
l_endo<-anosim(dat=(as.matrix((vegdist(ppp_endo,distance="bray")))),grouping=p_endo$Species)


colfunc<-colorRampPalette(c("blue", "red"))
pdf("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/CCA_endo.pdf")
ab<-cca(ppp_endo~VPD+rain+hum+temp+GI+RI,p_endo)
plot(ab,type="n",scaling=1)
with(ppp_endo,points(ab,display="species",col="gray",pch=3,lwd=0.1,scaling=1))
with(ppp_endo,points(ab,col=colfunc(2)[p_endo$Species],display="sites",scaling=1,pch=unique(as.numeric(p_endo$Species))))
fit<-envfit(ab~VPD+rain+hum+temp+GI+RI,p_endo,strata=p_endo$Species)
plot(fit, cex=1.2, axis=TRUE,p.max=0.05)
#text(ab, display="bp",scaling=1)
legend(x=-6.4,y=4.2,legend=unique(p_endo$Species),col=colfunc(2)[unique(p_exo$Species)])
#legend(x=-6.4,y=4.2,legend=unique(p_endo$Species),pch=unique(as.numeric(p_endo$Species)))

dev.off()







p_exo<-p[p$Endo_exo=="Exosphere",]
write.table(p_exo,"/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",sep="\t",quote=FALSE)
p_exo<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",header=T,row.names=1,sep="\t")
pp_exo<-p_exo[,((dim(map2)[2])+1):dim(p_exo)[2]]
ppp_exo<-pp_exo[,!!colSums(pp_exo)]
l_exo<-anosim(dat=(as.matrix((vegdist(ppp_exo,distance="bray")))),grouping=p_exo$Species)

colfunc<-colorRampPalette(c("blue", "red"))
#pdf("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/Endo_Exo.pdf")
plot(NMDS_bray,dis="site",type="n")
points(NMDS_bray,col=colfunc(2)[p$Endo_exo],pch=as.numeric(p$Species),cex=1.5,lwd=1.5)
legend(x=1.4,y=1.5,legend=unique(p$Endo_exo),fill=colfunc(2)[unique(p$Endo_exo)])
legend(x=1.45,y=1.1,legend=unique(p$Species),pch=unique(as.numeric(p$Species)))
text(y=-1.5,x=2,label=paste("Stress=",round(NMDS_bray$stress,digits=4)))
text(y=1.3,x=-0.2,label=paste("Endo vs. Exo p-value=",l$signif))
text(y=-1.3,x=-0.6,label=paste("p-value=",l_endo$signif),col="blue")
text(y=-0.9,x=0.8,label=paste("p-value=",l_exo$signif),col="red")
title("Endosphere vs. Exosphere bacterial communities (South canopy face)",font.main=1,cex.main=1,adj=0)
dev.off()
#################################################################different canopy

otu<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/core_all_R.txt",header=T,sep="\t",fill=T,row.names=1)
otu2<-otu[,1:(dim(otu)[2]-7)]
otus<-otu2
t_otus <- as.data.frame(t(otus))
t_otus<-decostand(t_otus,method="total")

map<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/F649_R889_map_r2.files",row.names=1,sep="\t",header=T)
map2<-map[map$INC=="core_all",]
x<-subset(t_otus[rownames(map2),])
p<-cbind(map2,x)
write.table(p,"/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",sep="\t",quote=FALSE)
p<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",header=T,row.names=1,sep="\t")

p_exo<-p[p$Endo_exo=="Exosphere",]
write.table(p_exo,"/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",sep="\t",quote=FALSE)
p_exo<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",header=T,row.names=1,sep="\t")

pp_exo<-p_exo[,((dim(map2)[2])+1):dim(p_exo)[2]]
ppp_exo<-pp_exo[,!!colSums(pp_exo)]


NMDS_bray<-metaMDS(ppp_exo,distance="bray",k=2,trymax=10000,autotransform=F)
zz<-scores(NMDS_bray,display="site")
l_canopy<-anosim(dat=(as.matrix((vegdist(ppp_exo,distance="bray")))),grouping=p_exo$Face)
l_species<-anosim(dat=(as.matrix((vegdist(ppp_exo,distance="bray")))),grouping=p_exo$Species)

colfunc<-colorRampPalette(c("blue", "red","green"))
pdf("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/canopy.pdf")
plot(NMDS_bray,dis="site",type="n")
points(NMDS_bray,col=colfunc(3)[p_exo$Face],pch=as.numeric(p_exo$Species),cex=1.5,lwd=1.5)
legend(x=1.55,y=1.76,legend=unique(p_exo$Face),fill=colfunc(3)[unique(p_exo$Face)])
legend(x=1.37,y=1.1,legend=unique(p_exo$Species),pch=unique(as.numeric(p_exo$Species)))
text(y=-1.5,x=1.8,label=paste("Stress=",round(NMDS_bray$stress,digits=4)))
text(y=1.3,x=-0.9,label=paste("Canopy p-value=",l_canopy$signif))
text(y=1.15,x=-0.9,label=paste("Species p-value=",l_species$signif))
title("Exosphere different canopy bacterial communities",font.main=1,cex.main=1,adj=0)
dev.off()


#################################################################abiotic









library(MASS)
library(vegan)
library(ggplot2)
library(plotrix)
library(plot3D)
library (Heatplus)
library(reshape)


otu<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/core_all_R.txt",header=T,sep="\t",fill=T,row.names=1)
rownames(otu)<-paste(otu$Family,row.names(otu),sep="_")
#otu<-aggregate(otu[,1:(dim(otu)[2]-7)],by=list(otu$Family),FUN=sum)
#rownames(otu)<-otu[,1]
#otu<-otu[,-1]
#otu<-otu[-1,]
otu2<-otu[,1:(dim(otu)[2]-7)]
otus<-otu2
t_otus <- as.data.frame(t(otus))
t_otus<-decostand(t_otus,method="total")
map<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/F649_R889_map_r2.files",row.names=1,sep="\t",header=T)
#map2<-map[map$INC=="core_all",]
#map2<-map[map$Endo_exo=="Exosphere",]
map2<-map[map$INC=="core_all" & map$Endo_exo=="Exosphere" & map$Species=="tortilis" & map$Month=="1" & map$Tree=="23",]
#map2<-map[map$INC=="core_all" & map$Face=="South",]
x<-subset(t_otus[rownames(map2),])
p<-cbind(map2,x)
write.table(p,"/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",sep="\t",quote=FALSE)
p<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",header=T,row.names=1,sep="\t")
rownames(p)<-p$Description
p_exoS<-p
write.table(p_exoS,"/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",sep="\t",quote=FALSE)
p_exoS<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",header=T,row.names=1,sep="\t")
p_exoS<-p_exoS[order(row.names(p_exoS)),]
#cbind(row.names(p_exoS),p_exoS$VPD)
pp_exo<-p_exoS[,((dim(map2)[2])+1):dim(p_exoS)[2]]
ppp_exo<-pp_exo[,!!colSums(pp_exo)]
ppp_exo<-decostand(ppp_exo,method="total")
write.table(ppp_exo,"/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",sep="\t",quote=FALSE)
ppp_exo<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",header=T,row.names=1,sep="\t")
maxab <- apply(ppp_exo, 2, max)
n1 <- names(which(maxab < 0.05))
data.prop.1 <- ppp_exo[, -which(names(ppp_exo) %in% n1)]
data.dist <- vegdist(ppp_exo, method = "bray")
row.clus <- hclust(data.dist, "complete")
data.dist.g <- vegdist(t(data.prop.1), method = "bray")
col.clus <- hclust(data.dist.g, "complete")
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)
#pdf("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/exo_raddiana.pdf")
plot(annHeatmap2(as.matrix(data.prop.1),col = colorRampPalette(c("lightyellow", "red"), space = "rgb")(72), breaks = 71,legend=3,labels=list(Col=list(nrow=18))))
#dev.off()



###############################
otu<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/core_all_R.txt",header=T,sep="\t",fill=T,row.names=1)
rownames(otu)<-paste(otu$Family,row.names(otu),sep="_")
otu2<-otu[,1:(dim(otu)[2]-7)]
otus<-otu2
t_otus <- as.data.frame(t(otus))
t_otus<-decostand(t_otus,method="total")
map<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/F649_R889_map_r2.files",row.names=1,sep="\t",header=T)
map1<-map[map$INC=="core_all" & map$Endo_exo=="Exosphere" & map$Species=="tortilis" & map$Month=="1" & map$Tree=="23",]
x1<-subset(t_otus[rownames(map1),])
p<-cbind(map1,x1)
write.table(p,"/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",sep="\t",quote=FALSE)
p<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",header=T,row.names=1,sep="\t")
rownames(p)<-p1$Description
p_exoS<-p
write.table(p_exoS,"/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",sep="\t",quote=FALSE)
p_exoS<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",header=T,row.names=1,sep="\t")
p_exoS<-p_exoS[order(row.names(p_exoS)),]
#cbind(row.names(p_exoS),p_exoS$VPD)
pp_exo<-p_exoS[,((dim(map2)[2])+1):dim(p_exoS)[2]]
ppp_exo<-pp_exo[,!!colSums(pp_exo)]
ppp_exo<-decostand(ppp_exo,method="total")
write.table(ppp_exo,"/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",sep="\t",quote=FALSE)
ppp_exo<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",header=T,row.names=1,sep="\t")
maxab <- apply(ppp_exo, 2, max)
n1 <- names(which(maxab < 0.05))
data.prop.1 <- ppp_exo[, -which(names(ppp_exo) %in% n1)]
data.dist <- vegdist(ppp_exo, method = "bray")
row.clus <- hclust(data.dist, "complete")
data.dist.g <- vegdist(t(data.prop.1), method = "bray")
col.clus <- hclust(data.dist.g, "complete")
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)
#pdf("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/tortilis_map$month=="1".pdf")
plot(annHeatmap2(as.matrix(data.prop.1),col = colorRampPalette(c("lightyellow", "red"), space = "rgb")(72), breaks = 71,legend=3,labels=list(Col=list(nrow=18))))
dev.off()













otu<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/core_all_R.txt",header=T,sep="\t",fill=T,row.names=1)
#otu2<-otu[,1:(dim(otu)[2]-7)]
map<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/F649_R889_map_r2.files",row.names=1,sep="\t",header=T)
map2<-map[map$INC=="core_all" & map$Face=="South",]



t_otu<-t(otu2)
x<-subset(t_otus[rownames(map2),])


#p_exoS<-p[p$Endo_exo=="Exosphere",]
#heatmap(as.matrix(data.prop.1))

pdf("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/heat_exo.pdf")
plot(annHeatmap2(as.matrix(data.prop.1),col = colorRampPalette(c("lightyellow", "red"), space = "rgb")(72), breaks = 61,legend=3,labels=list(Col=list(nrow=18))))
dev.off()



plot(annHeatmap2(as.matrix(data.prop.1),col = colorRampPalette(c("lightyellow", "red"), space = "rgb")(60), breaks = 61, dendrogram = list(Row = list(dendro = as.dendrogram(row.clus)),Col=list(dendro=as.dendrogram(col.clus))),legend=2,labels=list(Col=list(nrow=18))))
,ann=list(Row=list(data=ann.dat))))


otumat=as.matrix(otu2)
taxmat=as.matrix(otu2<-otu[,(dim(otu)[2]-6):(dim(otu)[2])])
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
physeq = phyloseq(OTU, TAX)


mydataAbundence <- apply(mydata,1,sum)
mydataAbundence <- sort(mydataAbundence, decreasing = T)
top20 <- mydataAbundence[1:20]
mydataTop20 <- mydata[row.names(mydata) %in% names(top20),]


##############################################3
library(MASS)
library(vegan)
library(ggplot2)
library(plotrix)
library(plot3D)
library (Heatplus)
otu<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/core_all_R.txt",header=T,sep="\t",fill=T,row.names=1)

rownames(otu)<-paste(row.names(otu),otu$Family)

otu2<-otu[,1:(dim(otu)[2]-7)]
otus<-otu2
t_otus <- as.data.frame(t(otus))
#t_otus<-decostand(t_otus,method="total")

map<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/F649_R889_map_r2.files",row.names=1,sep="\t",header=T)

map2<-map[map$INC=="core_all",]
#map2<-map[map$INC=="core_all" & map$Face=="South",]
x<-subset(t_otus[rownames(map2),])
p<-cbind(map2,x)
write.table(p,"/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",sep="\t",quote=FALSE)
p<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",header=T,row.names=1,sep="\t")
rownames(p)<-p$Description
p_exoS<-p[p$Endo_exo=="Exosphere",]
#p_exoS<-p
write.table(p_exoS,"/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",sep="\t",quote=FALSE)
p_exoS<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",header=T,row.names=1,sep="\t")
p_exoS<-p_exoS[order(row.names(p_exoS)),]

#cbind(row.names(p_exoS),p_exoS$VPD)

p_exo<-p_exoS
pp_exo<-p_exoS[,((dim(map2)[2])+1):dim(p_exoS)[2]]
ppp_exo<-pp_exo[,!!colSums(pp_exo)]
ppp_exo<-decostand(ppp_exo,method="total")


write.table(ppp_exo,"/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",sep="\t",quote=FALSE)
ppp_exo<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",header=T,row.names=1,sep="\t")


x<-rda(ppp_exo)
xx<-summary(x)
xxx<-as.data.frame(xx$sites)
scatter3d(xxx[,1],xxx[,2],xxx[,3],clab=c("pc1"),xlab="NMDS1",ylab="NMDS2",zlab="pc1",pch=as.numeric(p_exo$Face),groups=p_exo$season,grid=F,surface=F)



library("car")

NMDS_bray<-metaMDS(ppp_exo,distance="bray",k=3,trymax=10000,autotransform=F)
zz<-scores(NMDS_bray,display="site")
l_canopy<-anosim(dat=(as.matrix((vegdist(ppp_exo,distance="bray")))),grouping=p_exo$Face)
l_species<-anosim(dat=(as.matrix((vegdist(ppp_exo,distance="bray")))),grouping=p_exo$Species)

colfunc<-colorRampPalette(c("blue", "red","green"))
pdf("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/canopy.pdf")
plot(NMDS_bray,dis="site",type="n")
points(NMDS_bray,col=colfunc(3)[p_exo$Face],pch=as.numeric(p_exo$Species),cex=1.5,lwd=1.5)
legend(x=1.55,y=1.76,legend=unique(p_exo$Face),fill=colfunc(3)[unique(p_exo$Face)])
legend(x=1.37,y=1.1,legend=unique(p_exo$Species),pch=unique(as.numeric(p_exo$Species)))
text(y=-1.5,x=1.8,label=paste("Stress=",round(NMDS_bray$stress,digits=4)))
text(y=1.3,x=-0.9,label=paste("Canopy p-value=",l_canopy$signif))
text(y=1.15,x=-0.9,label=paste("Species p-value=",l_species$signif))
title("Exosphere different canopy bacterial communities",font.main=1,cex.main=1,adj=0)
dev.off()


l_minhum<-anosim(dat=(as.matrix((vegdist(ppp_exo,distance="bray")))),grouping=p_exo$minhum)
l_maxhum<-anosim(dat=(as.matrix((vegdist(ppp_exo,distance="bray")))),grouping=p_exo$maxhum)
l_mintemp<-anosim(dat=(as.matrix((vegdist(ppp_exo,distance="bray")))),grouping=p_exo$min_temp)
l_maxtemp<-anosim(dat=(as.matrix((vegdist(ppp_exo,distance="bray")))),grouping=p_exo$max_temp)
l_month<-anosim(dat=(as.matrix((vegdist(ppp_exo,distance="bray")))),grouping=p_exo$Month)

pdf("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/canopy_month.pdf")
scatter3D(zz[,1],zz[,2],p_exo$Month,clab=c("Month"),xlab="NMDS1",ylab="NMDS2",zlab="Month",pch=as.numeric(p_exo$Face))
legend("topleft",legend=unique(p_exo$Face),pch=unique(as.numeric(p_exo$Face)))
text(y=-0.45,x=-0.32,label=paste("Canopy p-value=",l_canopy$signif))
text(y=-0.49,x=-0.32,label=paste("Species p-value=",l_species$signif))
text(y=-0.53,x=-0.31,label=paste("Month p-value=",l_month$signif))
dev.off()

pdf("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/canopy_min_tem.pdf")
scatter3D(zz[,1],zz[,2],p_exo$min_temp,clab=c("min_temp","celsius"),xlab="NMDS1",ylab="NMDS2",zlab="min_temp",pch=as.numeric(p_exo$Face))
legend("topleft",legend=unique(p_exo$Face),pch=unique(as.numeric(p_exo$Face)))
text(y=-0.45,x=-0.32,label=paste("Canopy p-value=",l_canopy$signif))
text(y=-0.49,x=-0.32,label=paste("Species p-value=",l_species$signif))
text(y=-0.53,x=-0.31,label=paste("min_temp p-value=",l_mintemp$signif))
dev.off()

pdf("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/canopy_max_tem.pdf")
scatter3D(zz[,1],zz[,2],p_exo$max_temp,clab=c("max_temp","celsius"),xlab="NMDS1",ylab="NMDS2",zlab="min_temp",pch=as.numeric(p_exo$Face))
legend("topleft",legend=unique(p_exo$Face),pch=unique(as.numeric(p_exo$Face)))
text(y=-0.45,x=-0.32,label=paste("Canopy p-value=",l_canopy$signif))
text(y=-0.49,x=-0.32,label=paste("Species p-value=",l_species$signif))
text(y=-0.53,x=-0.31,label=paste("max_temp p-value=",l_maxtemp$signif))
dev.off()


pdf("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/canopy_maxhum.pdf")
scatter3D(zz[,1],zz[,2],p_exo$maxhum,clab=c("max_hum","%"),xlab="NMDS1",ylab="NMDS2",zlab="max_hum",pch=as.numeric(p_exo$Face))
legend("topleft",legend=unique(p_exo$Face),pch=unique(as.numeric(p_exo$Face)))
text(y=-0.45,x=-0.32,label=paste("Canopy p-value=",l_canopy$signif))
text(y=-0.49,x=-0.32,label=paste("Species p-value=",l_species$signif))
text(y=-0.53,x=-0.31,label=paste("max_hum p-value=",l_maxhum$signif))
dev.off()

pdf("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/canopy_mimhum.pdf")
scatter3D(zz[,1],zz[,2],p_exo$minhum,clab=c("min_hum","%"),xlab="NMDS1",ylab="NMDS2",zlab="min_hum",pch=as.numeric(p_exo$Face))
legend("topleft",legend=unique(p_exo$Face),pch=unique(as.numeric(p_exo$Face)))
text(y=-0.45,x=-0.32,label=paste("Canopy p-value=",l_canopy$signif))
text(y=-0.49,x=-0.32,label=paste("Species p-value=",l_species$signif))
text(y=-0.53,x=-0.31,label=paste("min_hum p-value=",l_minhum$signif))
dev.off()

################################################################################################3

#################################################################abiotic_tortilis

otu<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/core_all_R.txt",header=T,sep="\t",fill=T,row.names=1)
otu2<-otu[,1:(dim(otu)[2]-7)]
otus<-otu2
t_otus <- as.data.frame(t(otus))
t_otus<-decostand(t_otus,method="total")

map<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/F649_R889_map_r.files",row.names=1,sep="\t",header=T)
map2<-map[map$INC=="core_all",]
x<-subset(t_otus[rownames(map2),])
p<-cbind(map2,x)
write.table(p,"/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",sep="\t",quote=FALSE)
p<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",header=T,row.names=1,sep="\t")

p_exo<-p[p$Endo_exo=="Exosphere",]
write.table(p_exo,"/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",sep="\t",quote=FALSE)
p_exo<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",header=T,row.names=1,sep="\t")


p_exo<-p_exo[p_exo$Species=="A. tortilis",]
write.table(p_exo,"/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",sep="\t",quote=FALSE)
p_exo<-read.table("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/p.txt",header=T,row.names=1,sep="\t")



pp_exo<-p_exo[,((dim(map2)[2])+1):dim(p_exo)[2]]
ppp_exo<-pp_exo[,!!colSums(pp_exo)]


NMDS_bray<-metaMDS(ppp_exo,distance="bray",k=2,trymax=10000,autotransform=F)
zz<-scores(NMDS_bray,display="site")
l_canopy<-anosim(dat=(as.matrix((vegdist(ppp_exo,distance="bray")))),grouping=p_exo$Face)


colfunc<-colorRampPalette(c("blue", "red","green"))
pdf("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/canopy_tortilis.pdf")
plot(NMDS_bray,dis="site",type="n")
points(NMDS_bray,col=colfunc(3)[p_exo$Face],pch=as.numeric(p_exo$Face),cex=1.5,lwd=1.5)
legend(x=1.55,y=1.76,legend=unique(p_exo$Face),fill=colfunc(3)[unique(p_exo$Face)])
#legend(x=1.37,y=1.1,legend=unique(p_exo$Face),pch=unique(as.numeric(p_exo$Face)))
text(y=-1.5,x=1.8,label=paste("Stress=",round(NMDS_bray$stress,digits=4)))
text(y=1.3,x=-0.9,label=paste("Canopy p-value=",l_canopy$signif))
title("A. tortilis exosphere different canopy bacterial communities",font.main=1,cex.main=1,adj=0)
dev.off()


pdf("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/tortilis_canopy_min_tem.pdf")
scatter3D(zz[,1],zz[,2],p_exo$min_temp,clab=c("min_temp","celsius"),xlab="NMDS1",ylab="NMDS2",zlab="min_temp",pch=as.numeric(p_exo$Face))
legend("topleft",legend=unique(p_exo$Face),pch=unique(as.numeric(p_exo$Face)))
dev.off()

pdf("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/tortilis_canopy_max_tem.pdf")
scatter3D(zz[,1],zz[,2],p_exo$max_temp,clab=c("max_temp","celsius"),xlab="NMDS1",ylab="NMDS2",zlab="min_temp",pch=as.numeric(p_exo$Face))
legend("topleft",legend=unique(p_exo$Face),pch=unique(as.numeric(p_exo$Face)))
dev.off()


pdf("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/tortilis_canopy_maxhum.pdf")
scatter3D(zz[,1],zz[,2],p_exo$maxhum,clab=c("max_hum","%"),xlab="NMDS1",ylab="NMDS2",zlab="max_hum",pch=as.numeric(p_exo$Face))
legend("topleft",legend=unique(p_exo$Face),pch=unique(as.numeric(p_exo$Face)))
dev.off()

pdf("/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/tortilis_canopy_mimhum.pdf")
scatter3D(zz[,1],zz[,2],p_exo$minhum,clab=c("min_hum","%"),xlab="NMDS1",ylab="NMDS2",zlab="min_hum",pch=as.numeric(p_exo$Face))
legend("topleft",legend=unique(p_exo$Face),pch=unique(as.numeric(p_exo$Face)))
dev.off()
########################3





plot(p_exo$max_temp,zz[,1],type="n",main="Dor_Site_month metaMDS_13")
points(x=p_exo$max_temp,y=zz[,1], col=colfunc(3)[p_exo$Face],pch=as.numeric(p_exo$Species))

plot(p_exo$min_temp,zz[,1],type="n",main="Dor_Site_month metaMDS_13")
points(x=p_exo$min_temp,y=zz[,1], col=colfunc(3)[p_exo$Face],pch=as.numeric(p_exo$Species))

plot(p_exo$maxhum,zz[,1],type="n",main="Dor_Site_month metaMDS_13")
points(x=p_exo$minhum,y=zz[,1], col=colfunc(3)[p_exo$Face],pch=as.numeric(p_exo$Species))

plot(p_exo$minhum,zz[,1],type="n",main="Dor_Site_month metaMDS_13")
points(x=p_exo$minhum,y=zz[,1], col=colfunc(3)[p_exo$Face],pch=as.numeric(p_exo$Species))



#NMDS_bray<-metaMDS(p[,14:dim(p)[2]],distance="bray",k=2,trymax=10000,autotransform=F)
#NMDS_euclid<-metaMDS(p[,14:dim(p)[2]],distance="euclid",k=2,trymax=10000,autotransform=F)
#NMDS_manhattan<-metaMDS(p[,14:dim(p)[2]],distance="manhattan",k=2,trymax=10000,autotransform=F)
#NMDS_horn<-metaMDS(p[,14:dim(p)[2]],distance="horn",k=2,trymax=10000,autotransform=F)








plot(NMDS_euclid,dis="site",type="n")
points(NMDS_euclid,col=colfunc(3)[pp[,1]],pch=as.numeric(pp[,1]))

plot(NMDS_manhattan,dis="site",type="n")
points(NMDS_manhattan,col=colfunc(3)[pp[,1]],pch=as.numeric(pp[,1]))

plot(NMDS_horn,dis="site",type="n")
points(NMDS_horn,col=colfunc(3)[pp[,1]],pch=as.numeric(pp[,1]))




plot(p$Month,zz[,1],type='n')
points(x=p$Month,y=zz[,1], col=colfunc(2)[p$Species],pch=as.numeric(p$Species)





sort_otu_table.py \
-i /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/Tortilis_exosphere_300_North_YY.biom \
-m /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/F649_R889_map.files \
-s Month \
-o /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/Tortilis_exosphere_300_North_YY_M.biom

summarize_taxa.py -i \
/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/Tortilis_exosphere_300_North_YY_M.biom \
-o /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/

plot_taxa_summary.py \
-i /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/Tortilis_exosphere_300_North_YY_M_L2.txt \
-o /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/Tortilis_exosphere_300_North_YY_M_Figures_L2/



filter_samples_from_otu_table.py \
-i /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/rdp_F649_R889.bac.filtered.sample.biom \
-m /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/F649_R889_map.files \
-s Treatment1:exosphere \
-o /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/exosphere.biom



sort_otu_table.py \
-i /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/Tortilis_exosphere_300_Center_YY.biom \
-m /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/F649_R889_map.files \
-s Month \
-o /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/Tortilis_exosphere_300_Center_YY_M.biom

summarize_taxa.py -i \
/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/Tortilis_exosphere_300_Center_YY_M.biom \
-o /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/
 
plot_taxa_summary.py \
-i /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/Tortilis_exosphere_300_Center_YY_M_L2.txt \
-o /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/Tortilis_exosphere_300_Center_YY_M_Figures_L2/

filter_samples_from_otu_table.py \
-i /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/rdp_F649_R889.bac.filtered.sample.biom \
-m /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/F649_R889_map.files \
-s Sample_type:Tortilis_exosphere_300_South_YY \
-o /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/Tortilis_exosphere_300_South_YY.biom

sort_otu_table.py \
-i /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/Tortilis_exosphere_300_South_YY.biom \
-m /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/F649_R889_map.files \
-s Month \
-o /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/Tortilis_exosphere_300_South_YY_M.biom

summarize_taxa.py -i \
/home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/Tortilis_exosphere_300_South_YY_M.biom \
-o /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/


plot_taxa_summary.py \
-i /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/Tortilis_exosphere_300_South_YY_M_L2.txt \
-o /home/ashraf/Documents/Acacia/pear/bowtie2/mothurQC/splited_16S/F649_R889/denovo/Tortilis_exosphere_300_South_YY_M_Figures_L2/





done;
OLD_IFS=$IFS
