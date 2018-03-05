####################################
#having removed treatment samples###
####################################


#  PCA figure
setwd("~/Documents/Barnacles/notreat/mafs/")
library(stringr)

#bring in MAF files estimated for each pool of barnacles

warnings()
file_list<-list.files(pattern='*bam.remsecondary.removeddup.bam.mafs')
file_list
table<-c()

for (file in file_list){
  sampleID<-str_match(file,"\\[([A-Z]*\\d_\\d)\\]")[,2]
  print(sampleID)
  df<-read.table(file,header=T,sep="\t",row.names=NULL)
  head(df)
  df2<-data.frame(paste(df$chromo,df$position,sep="_"),df$phat)
  colnames(df2)[2] <- paste(sampleID,"_maf",sep="")
  colnames(df2)[1] <- "chromo_position"
  if (length(table) ==0){table<-df2
  } else {
    table<-merge(table,df2,by="chromo_position",all=TRUE) }
}  

warnings()  
colnames(table)

table2<-t(table[-1])
table2<-table2[ , colSums(is.na(table2)) == 0]

#annotation file to get groupings for colours etc
annot <- read.table('Barnacles.mindepth10ind.minind30.siteq30.nomapqual.majmin2.maf005.sitesNotDE.bam.remsecondary.removeddup.bam.annot', sep="\t", header=T, row.names=NULL); # note that plink cluster files are usually tab-separated instead
ansPop = annot$Group2
pCols = c("red","blue")
popCol = pCols[as.numeric(as.factor(ansPop))]
names(popCol) = ansPop
table(popCol,ansPop)

?prcomp
PCA.prcomp<-prcomp(table2,center=T) #center data here

#standard PCA
PCA.prcomp
screeplot(PCA.prcomp)
summary(PCA.prcomp)
#porportion explain
prop<-PCA.prcomp$sdev^2 / sum(PCA.prcomp$sdev^2)
plot(PCA.prcomp$x)
title <- paste("PC1"," (",signif(prop[1], digits=3)*100,"%)"," / PC2"," (",signif(prop[2], digits=3)*100,"%)",sep="",collapse="")
title

#data to plot
X<-as.data.frame(PCA.prcomp$x)
rownames(X)
X$pop<-annot$Group2

library(ggrepel)

cols<-c("red","blue")

X$names<-as.character(annot$Group1)
plottets<-ggplot(x=X$PC1, y=X$PC2) + 
  geom_point(aes(x=X$PC1, y=X$PC2,color=X$pop)) + 
  scale_colour_manual(values=cols) +
  ggtitle(title) + 
  theme_bw ()+ 
  geom_text_repel(aes(X$PC1,X$PC2,label=X$names),size=6,box.padding = 0.5,point.padding = 1) + #seperate package called ggrepel that has function geom_text_repel that can be useful for moving labels about nicely
  labs(x="PC1",y="PC2") +
  theme (legend.key.size=unit (0.05,'cm'))

plot(plottets)


#pictures of the average FST on sliding window
fst_slide <- read.table('../popoolation/BCT_WST_w500.maxcov50000min50.fst.clean', sep="\t", header=T, row.names=NULL); # note that plink cluster files are usually tab-separated instead

library(stringr)
fst_slide$fst_values_BCT.WST_split<-str_split_fixed(fst_slide$fst_values_BCT.WST, "=", 2)[,2]
plot(fst_slide$scaffold,fst_slide$fst_values_BCT.WST_split)


#install.packages("qqman")
#install.packages("manhattanly")
library(qqman)
library(manhattanly)
library(ggplot2)
?manhattan

fst_slide$scaffold
l=unique(fst_slide$scaffold)

fst_slide$CHR<-as.numeric(factor(fst_slide$scaffold, levels=l))
fst_slide$BP<-fst_slide$position_middle_window
fst_slide$P<-as.numeric(fst_slide$fst_values_BCT.WST_split)
fst_slide$SNP<-fst_slide$position_middle_window
manhattan(fst_slide,logp=FALSE)
manhattan(fst_slide)
#doesnt work on the larger sets too many 'chromosomes' so need to use ggplot

ggplot() + 
  geom_point(aes(x=fst_slide$CHR,y=fst_slide$P)) + 
  ggtitle("Fst Slide") + 
  theme_bw ()+ 
  labs(x="scaffold",y="Fst") #+
#  scale_y_reverse( )


fet <- read.table('../popoolation/BCT_WST.maxcov50000min50.fet', sep="\t", header=T, row.names=NULL) # note that plink cluster files are usually tab-separated instead
fet$scaffold
l=unique(fet$scaffold)

fet$CHR<-as.numeric(factor(fet$scaffold, levels=l))
fet$BP<-fet$position
fet$P<-as.numeric(str_split_fixed(fet$fisher_values_.log10.p.value._BCT.WST, "=", 2)[,2]) #fisher
#fet$P<-as.numeric(str_split_fixed(fet$fst_values_BCT.WST, "=", 2)[,2]) #fst
fet$SNP<-fet$position
#manhattan(fet,logp=FALSE)

plot(fet$CHR,fet$P,ylim=c(320, 0))

#fisher
ggplot() + 
  geom_point(aes(x=fet$CHR,y=fet$P)) + 
  ggtitle("Fisher Test -log10(pvalue)") + 
  theme_bw ()+ 
  labs(x="scaffold",y="-log10(p)") #+
  #scale_y_reverse( )

#fst
ggplot() + 
  geom_point(aes(x=fet$CHR,y=fet$P)) + 
  ggtitle("Fst Each Site") + 
  theme_bw ()+ 
  labs(x="scaffold",y="Fst") 

setwd("~/Documents/Barnacles/notreat/mafs/")

#plot of maf?
maf <- read.table('../popoolation/BCT_WST.max50000in10.freq_rc', sep="\t", header=T, row.names=NULL,stringsAsFactors = F) # note that plink cluster files are usually tab-separated instead
maf$MAF_BCT<-sapply(maf$maa_BCT, function(x) eval(parse(text=x)))
maf$MAF_WST<-sapply(maf$maa_WST, function(x) eval(parse(text=x)))

test<-as.matrix(cbind(maf$MAF_BCT, maf$MAF_WST))
library(RColorBrewer)

pcolors<-brewer.pal(10,"BrBG")

#my_palette <- colorRampPalette(c("darkturquoise", "darkgoldenrod2"))(n = 1000)
library(pheatmap)
#file is too big to run

#pheatmap(test,col = pcolors,treeheight_row=0,treeheight_col=0)
#pheatmap(test,col = pcolors,treeheight_row=0,treeheight_col=0,labels_col=c("BCT_MAF", "WST_MAF"))
#pheatmap(test,col = my_palette,treeheight_row=0,treeheight_col=0, cluster_cols=FALSE, cluster_rows=FALSE)




#####################################
##########exploring 2018#############
#####################################


#fst calculations using the stampp package
library(StAMPP)
setwd("~/Documents/Barnacles/notreat/mafs/")

#just want to compare between BCT and WST so bring in two column FST table
maf.file<-read.table('Barnacles.BCTWST.notreat.mindepth10wihtin.minind2.siteq30.nomapqual.majmin3.mafs.both.clean',row.names=NULL,sep='\t',header=T)

maf.file<-t(maf.file)
#write.table(maf.file,'Barnacles.BCTWST.notreat.mindepth10wihtin.minind2.siteq30.nomapqual.majmin3.mafs.both.clean.trans',row.names=F,sep='\t',quote=F,header=NULL)

colnames(maf.file) =maf.file[1, ]
maf.file = maf.file[-1, ] #chromo
maf.file = maf.file[-1, ] #position
maf.file = maf.file[-1, ] #allele 1
maf.file = maf.file[-1, ] #allele 2
dets<-data.frame(c('Sample','BCT','WST'),c('Pop','BCT','WST'),c('pop.num',1,2),c('ploidy',6,6),c('format','BiA','BiA'))
colnames(dets) <- dets[1,]

names(dets)[1] <- "Sample"
names(dets)[2] <- "Pop"
names(dets)[3] <- "pop.num"
names(dets)[4] <- "ploidy"
names(dets)[5] <- "format"
dets = dets[-1, ] #chromo
#write.table(maf.file,'Barnacles.BCTWST.notreat.mindepth10wihtin.minind2.siteq30.nomapqual.majmin3.mafs.both.clean.trans',row.names=F,sep='\t',quote=F,header=NULL)



maf.file2 <- data.frame(dets,  maf.file)
write.table(maf.file2,'Barnacles.BCTWST.notreat.mindepth10wihtin.minind2.siteq30.nomapqual.majmin3.mafs.both.clean.trans',row.names=F,sep='\t',quote=F)
maf.file3<-read.table('Barnacles.BCTWST.notreat.mindepth10wihtin.minind2.siteq30.nomapqual.majmin3.mafs.both.clean.trans',row.names=NULL,sep='\t',header=T)



maf.file.D.pop <- stamppNeisD(maf.file3, TRUE)
maf.file.D.fst <- stamppFst(maf.file3)
maf.file.D.fst


data(potato.mini, package="StAMPP")
potato.freq <- stamppConvert(potato.mini, "r")
potato.D.pop <- stamppNeisD(potato.freq, TRUE)

maf.file.upright<-read.table('Barnacles.BCTWST.notreat.mindepth10wihtin.minind2.siteq30.nomapqual.majmin3.mafs.both.clean',row.names=NULL,sep='\t',header=T)
head(maf.file.upright)
maf.file.upright$abs<-abs(maf.file.upright$BCT_MAF-maf.file.upright$WST_MAF)
maf.file.upright.nofixed<-maf.file.upright[!(maf.file.upright$abs> 0.95),]
maf.file.upright.nofixed$abs<-NULL
maf.file.upright.nofixed<-t(maf.file.upright.nofixed)
colnames(maf.file.upright.nofixed) =maf.file.upright.nofixed[1, ]
maf.file.upright.nofixed= maf.file.upright.nofixed[-1, ] #chromo
maf.file.upright.nofixed = maf.file.upright.nofixed[-1, ] #position
maf.file.upright.nofixed = maf.file.upright.nofixed[-1, ] #allele 1
maf.file.upright.nofixed = maf.file.upright.nofixed[-1, ] #allele 2
maf.file.upright.nofixed2 <- data.frame(dets,  maf.file.upright.nofixed)
write.table(maf.file.upright.nofixed2,'Barnacles.BCTWST.notreat.mindepth10wihtin.minind2.siteq30.nomapqual.majmin3.mafs.both.clean.trans.notfixed',row.names=F,sep='\t',quote=F)
maf.file.upright.nofixed3<-read.table('Barnacles.BCTWST.notreat.mindepth10wihtin.minind2.siteq30.nomapqual.majmin3.mafs.both.clean.trans.notfixed',row.names=NULL,sep='\t',header=T)
maf.file.upright.nofixed3.D.pop <- stamppNeisD(maf.file.upright.nofixed3, TRUE)

maf.file.upright.nofixed3.D.fst <- stamppFst(maf.file.upright.nofixed3)

maf.file.upright.nofixed3.D.fst



################################
######exploring data 2017#######
################################
#this was before we decided to throw out the treatment groups that act weird

setwd("~/Documents/Barnacles/")

#lazy PCA from angsd
covar <- read.table('Barnacles.mindepth20ind.minind32.siteq30.nomapqual.majmin2.maf005.sitesNotDE.geno.PCA', stringsAsFact=F);
covar <- read.table('Barnacles.nobcttreat.mindepth10across.minind10.siteq20.nomapqual.majmin2.maf005.geno.PCA', stringsAsFact=F);

#read table identifying groups and populations
annot <- read.table('Barnacles.mindepth10across.minind20.siteq30.nomapqual.majmin3.maf005.sites.PCA.geno.annot', sep="\t", header=T, row.names=NULL); # note that plink cluster files are usually tab-separated instead
annot <- read.table('Barnacles.nobcttreat.mindepth10across.minind10.siteq20.nomapqual.majmin2.maf005.geno.annot', sep="\t", header=T, row.names=NULL); # note that plink cluster files are usually tab-separated instead

a<-c("1","2")
comp<-as.numeric(a)

#turn eigen into % accounted
eig <- eigen(covar, symm=TRUE);
eig$val <- eig$val/sum(eig$val);
cat(signif(eig$val, digits=3)*100,"\n");

#associate groups with PCA values
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))
PC$Pop <- factor(annot$CLUSTER)
PC$Group1 <- factor(annot$Group1)
PC$Group2 <- factor(annot$Group2)

title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)",sep="",collapse="")

#axis labels
x_axis = paste("PC",comp[1],sep="")
y_axis = paste("PC",comp[2],sep="")

#use grid to more easily modify the legend size etc
library(grid)
library(ggplot2)
#plot
cols<-c("red","blue")
shap<-c(15,2,4,16,7)
plottets<-ggplot() + geom_point(data=PC, aes_string(x=x_axis, y=y_axis)) + ggtitle(title) + theme_bw ()+ theme (legend.key.size=unit (0.05,'cm'))
plot(plottets)
#lets add colour by population and shape by clusters identifyed in the ngsADMIX plot
plottets<-ggplot() + 
  geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Group2",shape="Pop")) + 
  scale_colour_manual(values=cols) +
  scale_shape_manual(values=shap) +
  ggtitle(title) + 
  theme_bw ()+ 
  theme (legend.key.size=unit (0.05,'cm'))

plot(plottets)

library(ggrepel)


PC$names<-as.character(annot$Group1)
plottets<-ggplot() + 
  geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Group2",shape="Pop")) + 
  scale_colour_manual(values=cols) +
  scale_shape_manual(values=shap) +
  ggtitle(title) + 
  theme_bw ()+ 
  geom_text_repel(aes(PC$PC1,PC$PC2,label=PC$Group1),size=4) + #seperate package called ggrepel that has function geom_text_repel that can be useful for moving labels about nicely
  theme (legend.key.size=unit (0.05,'cm'))

plot(plottets)

#historgrams:

table <- read.table('Barnacles.mindepth10ind.minind30.siteq30.nomapqual.majmin3.maf005.sites.counts2.counts', stringsAsFact=F,header=T);
hist(table$WST2_5,xlim=c(10,500),breaks=1000)

hist(table$BCT1_2,xlim=c(0,500),breaks=1000)

library(dplyr)
table.wst<-select(table,contains("WST"))
table.bct<-select(table,contains("BCT"))

table.wst$mean<-rowMeans(table.wst, na.rm = FALSE, dims = 1)
table.bct$mean<-rowMeans(table.bct, na.rm = FALSE, dims = 1)
hist(table.wst$mean,xlim=c(0,500),breaks=1000)
hist(table.bct$mean,xlim=c(0,500),breaks=1000)

table.oddballs<-data.frame(table$BCT1_3,table$BCT2_1,table$BCT2_2,table$BCT2_3)
table.oddballs$mean<-rowMeans(table.oddballs, na.rm = FALSE, dims = 1)
hist(table.oddballs$mean,xlim=c(0,500),breaks=1000)

table.notoddballs<-select(table,contains("BCT"))

table.notoddballs$BCT1_3<-NULL
table.notoddballs$BCT2_1<-NULL
table.notoddballs$BCT2_2<-NULL
table.notoddballs$BCT2_3<-NULL
table.notoddballs$mean<-rowMeans(table.notoddballs, na.rm = FALSE, dims = 1)
hist(table.notoddballs$mean,xlim=c(0,500),breaks=1000)

#########################
##########PCA############
#########################

#PCA from mafs
library(stringr)
setwd("~/Documents/Barnacles/mafs_nonDE/")

#######create table of counts#######

file_list<-list.files(pattern='*mafs')
file_list
table<-c()

for (file in file_list){
  sampleID<-str_match(file,"(^.*?-\\[(.*?)\\])")[,3]
  print(sampleID)
  df<-read.table(file,header=T,sep="\t",row.names=NULL)
  df2<-df[,c(1,2)]
  colnames(df2)[colnames(df2)=="phat"] <- paste(sampleID,"_maf",sep="")
  if (length(table) ==0){table<-df2
  } else {
    table<-merge(table,df2,by="chromo_position_major_minor") }
}  
warnings()


table2<-t(table[-1])

pca = pcaGenotypes(table2)
annot <- read.table('../Barnacles.mindepth10across.minind20.siteq30.nomapqual.majmin3.maf005.sites.PCA.geno.annot', sep="\t", header=T, row.names=NULL); # note that plink cluster files are usually tab-separated instead
ansPop = annot$Group2
pCols = c("red","blue")
popCol = pCols[as.numeric(as.factor(ansPop))]
names(popCol) = ansPop
table(popCol,ansPop)
plotPCA(pca, popCol)



PCA.prcomp<-prcomp(table2)

PCA.prcomp
screeplot(PCA.prcomp)
summary(PCA.prcomp)
prop<-PCA.prcomp$sdev^2 / sum(PCA.prcomp$sdev^2)
plot(PCA.prcomp$x)
title <- paste("PC1"," (",signif(prop[1], digits=3)*100,"%)"," / PC2"," (",signif(prop[2], digits=3)*100,"%)",sep="",collapse="")
title

X<-as.data.frame(PCA.prcomp$x)
rownames(X)
X$pop<-annot$Group2

library(ggrepel)

cols<-c("red","blue")

X$names<-as.character(annot$Group1)
plottets<-ggplot(x=X$PC1, y=X$PC2) + 
  geom_point(aes(x=X$PC1, y=X$PC2,color=X$pop)) + 
  scale_colour_manual(values=cols) +
  ggtitle(title) + 
  theme_bw ()+ 
  geom_text_repel(aes(X$PC1,X$PC2,label=X$names),size=4) + #seperate package called ggrepel that has function geom_text_repel that can be useful for moving labels about nicely
  theme (legend.key.size=unit (0.05,'cm'))

plot(plottets)

#
