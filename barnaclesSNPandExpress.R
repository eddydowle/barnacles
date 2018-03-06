#how to look at the comparisons between expression and snp datasets
#expression

library(dplyr)

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














  