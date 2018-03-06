#how to look at the comparisons between expression and snp datasets
#expression

#markdown etc is available on github:
#https://github.com/eddydowle/barnacles

library(dplyr)

#6th March 2018
####testing for evidence of selection on DE transcipts###

#######################################################
##############across all DE transcripts################
#######################################################

expressiondat<-read.csv("~/Documents/github/barnacles/annot.results.csv",header=T)


dcorrected_BCT<-read.table("~/Documents/Barnacles/notreat/popoolation/barnacles.BCT.notreat.sort.gtf.min4max50000.downsample.dcorrected.clean",header=T)
dcorrected_WST<-read.table("~/Documents/Barnacles/notreat/popoolation/barnacles.WST.notreat.sort.gtf.min4max50000.downsample.dcorrected.clean",header=T)

#some of these got annotated as having two genes on the same scaffold Im just going to average by scafold ~not sure if that is the 'right' way to go about this
dcorrected_BCT<-dcorrected_BCT %>% group_by(.,SeqName) %>% summarise(.,dcorrected_ave = mean(dcorrected, na.rm = T),snp_count=sum(number_snps))
dcorrected_WST<-dcorrected_WST %>% group_by(.,SeqName) %>% summarise(.,dcorrected_ave = mean(dcorrected, na.rm = T),snp_count=sum(number_snps))

dcorrected_BCT %>% filter(., SeqName %in% expressiondat$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-1.188611 

dcorrected_BCT %>% filter(., !SeqName %in% expressiondat$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-0.6417929 

dcorrected_WST %>% filter(., SeqName %in% expressiondat$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-1.422189 

dcorrected_WST %>% filter(., !SeqName %in% expressiondat$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-0.8071944


a<- dcorrected_BCT %>% filter(., SeqName %in% expressiondat$SeqName) %>% ungroup() %>% select(.,dcorrected_ave)
hist(a$dcorrected_ave)

b<-dcorrected_BCT %>% filter(., !SeqName %in% expressiondat$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) 
hist(b$dcorrected_ave)

t.test(a$dcorrected_ave,b$dcorrected_ave)

#Welch Two Sample t-test

#data:  a$dcorrected_ave and b$dcorrected_ave
#t = -20.982, df = 2021, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.5979277 -0.4957094
#sample estimates:
#  mean of x  mean of y 
#-1.1886114 -0.6417929 

a<- dcorrected_WST %>% filter(., SeqName %in% expressiondat$SeqName) %>% ungroup() %>% select(.,dcorrected_ave)
hist(a$dcorrected_ave)

b<-dcorrected_WST %>% filter(., !SeqName %in% expressiondat$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) 
hist(b$dcorrected_ave)

t.test(a$dcorrected_ave,b$dcorrected_ave)

#Welch Two Sample t-test

#data:  a$dcorrected_ave and b$dcorrected_ave
#t = -15.42, df = 955.63, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.6932628 -0.5367260
#sample estimates:
#  mean of x  mean of y 
#-1.4221888 -0.8071944 


########7th March 2018########

#testing on just transcipts DE in WST population#


##########################################
################question 1################
##########################################


#1) transcripts that were only differentially regulated in the White Sea (that's the more freeze tolerant) population

expressionWST<-read.csv("~/Documents/github/barnacles/allwsdiffreg.csv",header=T)


dcorrected_BCT<-read.table("~/Documents/Barnacles/notreat/popoolation/barnacles.BCT.notreat.sort.gtf.min4max50000.downsample.dcorrected.clean",header=T)
dcorrected_WST<-read.table("~/Documents/Barnacles/notreat/popoolation/barnacles.WST.notreat.sort.gtf.min4max50000.downsample.dcorrected.clean",header=T)

#some of these got annotated as having two genes on the same scaffold Im just going to average by scafold ~not sure if that is the 'right' way to go about this
dcorrected_BCT<-dcorrected_BCT %>% group_by(.,SeqName) %>% summarise(.,dcorrected_ave = mean(dcorrected, na.rm = T),snp_count=sum(number_snps))
dcorrected_WST<-dcorrected_WST %>% group_by(.,SeqName) %>% summarise(.,dcorrected_ave = mean(dcorrected, na.rm = T),snp_count=sum(number_snps))

dcorrected_BCT %>% filter(., SeqName %in% expressionWST$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-1.189827 

dcorrected_BCT %>% filter(., !SeqName %in% expressionWST$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-0.7021807

dcorrected_WST %>% filter(., SeqName %in% expressionWST$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-1.445604 

dcorrected_WST %>% filter(., !SeqName %in% expressionWST$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-0.8319091

a<- dcorrected_BCT %>% filter(., SeqName %in% expressionWST$SeqName) %>% ungroup() %>% select(.,dcorrected_ave)
hist(a$dcorrected_ave)

b<-dcorrected_BCT %>% filter(., !SeqName %in% expressionWST$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) 
hist(b$dcorrected_ave)

t.test(a$dcorrected_ave,b$dcorrected_ave)

#Welch Two Sample t-test

#data:  a$dcorrected_ave and b$dcorrected_ave
#t = -2.9455, df = 41.176, p-value = 0.005284
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.8219491 -0.1533430
#sample estimates:
#  mean of x  mean of y 
#-1.1898267 -0.7021807 

a<- dcorrected_WST %>% filter(., SeqName %in% expressionWST$SeqName) %>% ungroup() %>% select(.,dcorrected_ave)
hist(a$dcorrected_ave)

b<-dcorrected_WST %>% filter(., !SeqName %in% expressionWST$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) 
hist(b$dcorrected_ave)

t.test(a$dcorrected_ave,b$dcorrected_ave)

#Welch Two Sample t-test

#data:  a$dcorrected_ave and b$dcorrected_ave
#t = -8.3487, df = 252.35, p-value = 4.594e-15
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.7584620 -0.4689284
#sample estimates:
#  mean of x  mean of y 
#-1.4456042 -0.8319091 

##########################################
################question 2################
##########################################


#2)Anything labelled with "transcription" as a function.  I think these are likely transcription factors, and in particular, XLOC_006827 and XLOC_01766 are upregulated in both the BC and WS populations in the 4 hours following freezing.

colnames(expressionWST)

expressionWST_transciption<-expressionWST %>% filter(.,grepl('transcription',InterPro.GO.Names) | grepl('transcription',GO.Names)) 

#there is only one transcipt in this comparison though for BCT coverage isnt high enough for these transcripts for a Tajima D calculation
dcorrected_BCT %>% filter(., SeqName %in% expressionWST_transciption$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-1.157307

dcorrected_BCT %>% filter(., !SeqName %in% expressionWST_transciption$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-0.7035338

dcorrected_WST %>% filter(., SeqName %in% expressionWST_transciption$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-2.479258 

dcorrected_WST %>% filter(., !SeqName %in% expressionWST_transciption$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-0.8406581

#there isnt enough of these transcription transcripts covered at high enough coverage in BCT so we cant do a t-test
a<- dcorrected_BCT %>% filter(., SeqName %in% expressionWST_transciption$SeqName) %>% ungroup() %>% select(.,dcorrected_ave)
hist(a$dcorrected_ave)
a
b<-dcorrected_BCT %>% filter(., !SeqName %in% expressionWST_transciption$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) 
hist(b$dcorrected_ave)
b
t.test(a$dcorrected_ave,b$dcorrected_ave)


#WST
a<- dcorrected_WST %>% filter(., SeqName %in% expressionWST_transciption$SeqName) %>% ungroup() %>% select(.,dcorrected_ave)
hist(a$dcorrected_ave)
a
b<-dcorrected_WST %>% filter(., !SeqName %in% expressionWST_transciption$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) 
hist(b$dcorrected_ave)
b
t.test(a$dcorrected_ave,b$dcorrected_ave)

#Welch Two Sample t-test

#data:  a$dcorrected_ave and b$dcorrected_ave
#t = -5.6627, df = 8.013, p-value = 0.0004716
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -2.3056946 -0.9715055
#sample estimates:
#  mean of x  mean of y 
#-2.4792582 -0.8406581 

##########################################
################question 3################
##########################################

#3) Anything labelled as a "macrophage mannose receptor".  I think these are potentially novel antifreeze proteins--they often blast to fish Type 2 antifreeze proteins.  

colnames(expressionWST)

expressionWST %>% filter(.,grepl('macrophage mannose receptor',Description)) 

expressionWST_MMR <- expressionWST %>% filter(.,grepl('macrophage mannose receptor',Description)) 


#there is only one transcipt in this comparison though for BCT coverage isnt high enough for these transcripts for a Tajima D calculation
dcorrected_BCT %>% filter(., SeqName %in% expressionWST_MMR$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-1.108875

dcorrected_BCT %>% filter(., !SeqName %in% expressionWST_MMR$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-0.7035097

dcorrected_WST %>% filter(., SeqName %in% expressionWST_MMR$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-1.320399

dcorrected_WST %>% filter(., !SeqName %in% expressionWST_MMR$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-0.8414786

#there isnt enough of these transcription transcripts covered at high enough coverage in BCT so we cant do a t-test
a<- dcorrected_BCT %>% filter(., SeqName %in% expressionWST_MMR$SeqName) %>% ungroup() %>% select(.,dcorrected_ave)
hist(a$dcorrected_ave)
a
b<-dcorrected_BCT %>% filter(., !SeqName %in% expressionWST_MMR$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) 
hist(b$dcorrected_ave)
b
#t.test(a$dcorrected_ave,b$dcorrected_ave)


#WST also only has four observations for MMR so this isnt going to work either.
a<- dcorrected_WST %>% filter(., SeqName %in% expressionWST_MMR$SeqName) %>% ungroup() %>% select(.,dcorrected_ave)
hist(a$dcorrected_ave)
a
b<-dcorrected_WST %>% filter(., !SeqName %in% expressionWST_MMR$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) 
hist(b$dcorrected_ave)
b
t.test(a$dcorrected_ave,b$dcorrected_ave)

#Welch Two Sample t-test

#data:  a$dcorrected_ave and b$dcorrected_ave
#t = -0.97935, df = 3.0017, p-value = 0.3996
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -2.034701  1.076860
#sample estimates:
#  mean of x  mean of y 
#-1.3203988 -0.8414786 


####################################################
####################data exploring##################
####################################################

####################FST and FET###################

expressiondat<-read.csv("~/Documents/github/barnacles/annot.results.csv",header=T)
head(expressiondat)
#snps

snps<-read.table("~/Documents/Barnacles/notreat/popoolation/barnacles.BCT.notreat.sort.gtf.min4max50000.downsample.dcorrected.clean",header=T)
#can be strict or not
#fst<-read.table("~/Documents/Barnacles/notreat/popoolation/BCT_WST.max50000min10.fst",header=T)
fst<-read.table("~/Documents/Barnacles/notreat/popoolation/BCT_WST.max50000min50.fst",header=T)

#Just looking at the value of the top snp in a scaffold here just a rough look
fst_top <- fst %>% group_by(.,SeqName) %>% slice(which.max(fst_values_BCT.WST))

fst_top %>% filter(., SeqName %in% expressiondat$SeqName) %>% ungroup() %>% select(.,fst_values_BCT.WST) %>% colMeans()
#fst_values_BCT.WST 
#0.9777981 

#strict
#$fst_values_BCT.WST 
#0.960366 


fst_top %>% filter(., !SeqName %in% expressiondat$SeqName) %>% ungroup() %>% select(.,fst_values_BCT.WST) %>% colMeans()
#fst_values_BCT.WST 
#0.9808638 
#actually slightly higher in those not DE

#strict
#fst_values_BCT.WST 
#0.9480188 

dcorrected_BCT<-read.table("~/Documents/Barnacles/notreat/popoolation/barnacles.BCT.notreat.sort.gtf.min4max50000.downsample.dcorrected.clean",header=T)
dcorrected_WST<-read.table("~/Documents/Barnacles/notreat/popoolation/barnacles.WST.notreat.sort.gtf.min4max50000.downsample.dcorrected.clean",header=T)

#some of these got annotated as having two genes on the same scaffold Im just going to average by scafold ~not sure if that is the 'right' way to go about this
dcorrected_BCT<-dcorrected_BCT %>% group_by(.,SeqName) %>% summarise(.,dcorrected_ave = mean(dcorrected, na.rm = T),snp_count=sum(number_snps))
dcorrected_WST<-dcorrected_WST %>% group_by(.,SeqName) %>% summarise(.,dcorrected_ave = mean(dcorrected, na.rm = T),snp_count=sum(number_snps))

dcorrected_BCT %>% filter(., SeqName %in% expressiondat$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-1.188611 

dcorrected_BCT %>% filter(., !SeqName %in% expressiondat$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-0.6417929 

dcorrected_WST %>% filter(., SeqName %in% expressiondat$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-1.422189 

dcorrected_WST %>% filter(., !SeqName %in% expressiondat$SeqName) %>% ungroup() %>% select(.,dcorrected_ave) %>% colMeans()
#dcorrected_ave 
#-0.8071944




fet<-read.table("~/Documents/Barnacles/notreat/popoolation/BCT_WST.maxcov50000min10.fet",header=T)
#fet<-read.table("~/Documents/Barnacles/notreat/popoolation/BCT_WST.maxcov50000min50.fet",header=T)

#have to drop values that came in as Inf-less than ideal
#I think the tajimas D values are going to be most interesting as the differentiation is so high I think it just sweeps any signal out
fet_top <- fet %>% group_by(.,SeqName) %>% filter(., !fisher_values_.log10.p.value._BCT.WST==Inf) %>%  slice(which.max(fisher_values_.log10.p.value._BCT.WST))

fet_top %>% filter(., SeqName %in% expressiondat$SeqName) %>% ungroup() %>% select(.,fisher_values_.log10.p.value._BCT.WST) %>% colMeans()
#fisher_values_.log10.p.value._BCT.WST 
#63.74633 

#strict
#fisher_values_.log10.p.value._BCT.WST 
#122.4163 

fet_top %>% filter(., !SeqName %in% expressiondat$SeqName) %>% ungroup() %>% select(.,fisher_values_.log10.p.value._BCT.WST) %>% colMeans()
#fisher_values_.log10.p.value._BCT.WST 
#50.62169 

#strict
#fisher_values_.log10.p.value._BCT.WST 
#123.2578 



