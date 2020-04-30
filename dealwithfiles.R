########Loading packages##################################
library(rtracklayer)
library(GenomicRanges)
library(clusterProfiler)
library(org.Dm.eg.db)
library(tidyverse)

########Loading the RNAseq file and update ID##########################
setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/raw data/rnaseq")
fbrnaseq<-read.csv(paste(getwd(),"/","rnaseq_switch_fb_foxo.csv",sep=''))
rnaid<-read.csv(paste(getwd(),"/","switchfb_rnaseq_619to632_ID.csv",sep='')) #change ID from 6.19 to 6.32 on flybase Id converter ##14555 object to 14556 object

colnames(fbrnaseq)[1]<-"previousID"
colnames(rnaid)[1]<-"previousID"

rnaiddup<-as.character(rnaid$updatedID[duplicated(rnaid$updatedID)]) #duplicated "FBgn0286778" "FBgn0286809" found
rnaidspl<-as.character(rnaid$updatedID[grep("::", rnaid$updatedID)]) #none split found
a<-as.character(rnaid$previousID[duplicated(rnaid$previousID)]) #FBgn0036125
a1<-subset(rnaid,rnaid$previousID %in% a)
a2<-as.character(a1$converted_id) 
rnaidremove<-c(rnaiddup,rnaidspl,a2)

rnaidclearid<-as.character(setdiff(rnaid$updatedID,rnaidremove))
rnaidclear<-subset(rnaid,rnaid$updatedID %in% rnaidclearid) #only left 14552 object

rnafinal<-left_join(rnaidclear,fbrnaseq,by="previousID") 

#convert the flybase ID to NCBI refseq for analyze. (the BETA only recognize [a-z][a-z]/_//d* type, which is NW_123456 something) #fake
#switchcg <- bitr(rnafinal1$FLYBASE, fromType="FLYBASE", toType="REFSEQ", OrgDb="org.Dm.eg.db") #1.81% of input gene IDs are fail to map ##what more, 14552 flybase bitr to 59576 NCBI Refseq
#switchcg1 <- bitr(rnafinal1$FLYBASE, fromType="FLYBASE", toType="SYMBOL", OrgDb="org.Dm.eg.db") #This give gene symbol ##0.14% unmapped ###same as GENENAME
#colnames(rnafinal1)[1]<-"FLYBASE"
#afterid<-left_join(switchcg,rnafinal1,by="FLYBASE")
#afterid1<-afterid[,-1]

rnaextract<-rnafinal[,c(3,6,10)] #only output flybase ID, Log2fd, padj
rnaextrac1<-rnaextract[complete.cases(rnaextract),] #only 13822 object left.
colnames(rnaextrac1)[1:3]<-c("ID","status","value")

write.csv(rnaextrac1, paste(getwd(),"/","fb_switch_de.csv",sep=''),row.names = FALSE)

#change the chr name of the fasta file manually
#change the atac-peak file chr name manually
#change the DE file to tab-deliminated txt version, add the #at the beginning.
#change the chr name of the gtf file

#############change the switch gut de data#########################################
setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch experiment/raw data/rnaseq")
gutrnaseq<-read.csv(paste(getwd(),"/","switch_gut_foxo_rnaseq.csv",sep=''))
rnaid<-read.csv(paste(getwd(),"/","gut_switch_6.19to6.32_flybaseid.csv",sep='')) #change ID from 6.19 to 6.32 on flybase Id converter ##14254 object to 14253 object

colnames(gutrnaseq)[1]<-"previousID"
colnames(rnaid)[1]<-"previousID"

rnaiddup<-as.character(rnaid$converted_id[duplicated(rnaid$converted_id)]) #"FBgn0286778" "FBgn0286809"
rnaidspl<-as.character(rnaid$converted_id[grep("::", rnaid$converted_id)]) #empty
a<-as.character(rnaid$previousID[duplicated(rnaid$previousID)]) #FBgn0036125
a1<-subset(rnaid,rnaid$previousID %in% a)
a2<-as.character(a1$converted_id) 
rnaidremove<-c(rnaiddup,rnaidspl,a2)

rnaidclearid<-as.character(setdiff(rnaid$converted_id,rnaidremove))
rnaidclear<-subset(rnaid,rnaid$converted_id %in% rnaidclearid) #only left 14249 object

rnafinal<-left_join(rnaidclear,gutrnaseq,by="previousID") 

rnaextract<-rnafinal[,c(3,6,10)] #only output flybase ID, Log2fd, padj
rnaextrac1<-rnaextract[complete.cases(rnaextract),] #only 13551 object left.
colnames(rnaextrac1)[1:3]<-c("ID","status","value")

write.csv(rnaextrac1, paste(getwd(),"/","gut_switch_de.csv",sep=''),row.names = FALSE)

#############change the annotation data##############################
setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch experiment/BETA")
anno<-import(paste(getwd(),"/","dmel-all-r6.32.gtf",sep=''))
seqlevelsStyle(anno)<-"UCSC"
anno@seqnames


#deal with the annotation file by extracting the flybase ID out of the metadata


export.gff3(anno,paste(getwd(),"/","dmel-all-r6.32.gff3",sep='')) #the gff3 will include all the metadata in the annotation.
#delete the header of gff3 annotation mannually (left#r6.32 only)

#deal with the annotation file by extracting the flybase ID out of the metadata
#export the annotation like colnames(anno2)[1:6]<-c("name","chrom","strand","txStart","txEnd","name2"), delete"", add#

#still meeting index out of range problem, (think is due to the rna-seq and annotation file> so import these two and setdiff)

de<-read.table("fb_switch_dev1.txt")
an<-read.table("r6.32minimumv3.txt")
k<-setdiff(de$V1,an$V1) #really found two FBgn0267117 FBgn0262818 (lncRNA:CR45557,CG43189), how? for now delete both 
k1<-setdiff(an$V1,de$V1) #only 4049 object, is that normal?

anu<-unique(an$V1) #17869 object
angene<-subset(an,an$V6 == "gene") #also 17869 object, then i should use this 
angene1<-angene[,c(6,2,3,4,5,1)]
colnames(angene1)[1:6]<-c("name","chrom","strand","txStart","txEnd","name2")

#chromosome<-c("chr2L","chr2R","chr3L","chr3R","chr4","chrX","chrY")
#angene2<-subset(angene1,angene1$chrom %in% chromosome) #17792 object left

gname<-angene1$name2
angene2<-cbind(gname,angene1[,-1])

write.csv(angene2,paste(getwd(),"/","r6.32minimumv3.csv",sep=''),row.names = FALSE)

#now check what's wrong the atac-seq peak.
#atac<-read.table("foxoatacv6.bed") 
#2L 48113-23513653
#2R 1203-25284572
#3L 28036-28109545
#3R 10716-31841689
#4 270045-1295229
#X 4866-23538562
#Y 22683-3667352
# the problem is...should i only use the significant peak? no, just use all #no still i think only use those significant (or even just those significant upregulated one), since the score itself make nosense.
## one possibility is that i need to use raw count rather than log2fc as score, indicated in the template file.
setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch experiment/raw data")

atacraw<-read.csv("ATACFBallResults.csv") #3447 peak
ata1<-atacraw[,c(9,10,11,2,3,6,7)]
chr<-c("2L","2R","3L","3R","4","X","Y")
ata2<-subset(ata1,ata1$Chr %in% chr) #2074 peak
ata2$Chr<-paste("chr",ata2$Chr, sep='')
strandn<-1:2074
strand<-paste("undefined_",strandn, sep='')
ata3<-cbind(ata2,strand)

#collect only significant peak
atasig<-subset(ata3,ata3$padj<0.1)#79 left
atasigp<-subset(ata3,ata3$pvalue<0.05) #244 left

#only up 
atasigpup<-subset(atasigp,atasigp$log2FoldChange>0)#186 keft

atafinal<-atasigpup[,c(1,2,3,8,4)]
write.csv(atafinal,paste(getwd(),"/","atacsigpup.csv",sep=''),row.names = FALSE)
#Then delete header, convert to tab-deliminated txt, then change to .bed 

#deal with the gut peak
setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch experiment/raw data/atacpeak")

atacraw<-read.csv("GuallResults.csv") #2906 peak
ata1<-atacraw[,c(9,10,11,2,3,6,7)]
chr<-c("2L","2R","3L","3R","4","X","Y")
ata2<-subset(ata1,ata1$Chr %in% chr) #2316 peak
ata2$Chr<-paste("chr",ata2$Chr, sep='')
strandn<-1:2316
strand<-paste("undefined_",strandn, sep='')
ata3<-cbind(ata2,strand)

#collect only significant peak
atasig<-subset(ata3,ata3$padj<0.2)#0 left
atasigp<-subset(ata3,ata3$pvalue<0.05) #125 left

#only up 
atasigpup<-subset(atasigp,atasigp$log2FoldChange>0)#31 keft
#only down
atasigpdown<-subset(atasigp,atasigp$log2FoldChange<0)#94 keft

atafinal<-atasigpdown[,c(1,2,3,8,4)]
write.csv(atafinal,paste(getwd(),"/","Gutsigpdown.csv",sep=''),row.names = FALSE)
#Then delete header, convert to tab-deliminated txt, then change to .bed 


########################################Try the foxo CHIP peak with the foxo RNA-seq#########################
#adam-2014-chipseq-r5 #Check that pervious conversion is fine. #location of all the peaks detected on ChIP-chip of GFP-dFOXO from induced S1106>GFP-dfoxo flies (genome release 5)
#nazif-2011-chipchip-r4 #flybase r4>r6 #foxo binding in will type, whole flies, genome release 4
r4<-read.csv("nazif 2011 peak.csv")
position<-paste(r4[,1],":",r4$Start,"..",r4$End,se)
colnames(r4)[1]<-"Chromosome"
r41<-cbind(r4,position)
write.csv(r41,paste(getwd(),"/","nazif2011.csv", sep=''))
r6<-read.csv("nazif2011.csv")
Chromosome<-str_extract(r6$R6,"^[\\s\\S]*?:")
start<-str_extract(r6$R6,":[\\s\\S]*[.]")
end<-str_extract(r6$R6,"[.][\\s\\S]*")
start<-sub("[.]",'',start)
start<-sub("[.]",'',start)
start<-sub(":",'',start)
end<-sub("[.]",'',end)
end<-sub("[.]",'',end)
Chromosome<-sub(":",'',Chromosome)
Chromosome<-paste("chr",Chromosome,sep='')
r61<-cbind(Chromosome,start,end)
write.csv(r61,paste(getwd(),"/","nazif2011chippeakr6.csv",sep=''),row.names = FALSE)

#############adam-2014-rnaseq#############################################
#fb
setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch experiment/raw data/rnaseq/unknown")
fbrnaseq<-read.csv(paste(getwd(),"/","naziffbrename.csv",sep=''))
rnaid<-read.csv(paste(getwd(),"/","2014_fb_FlyBase_IDs.csv",sep='')) 

colnames(fbrnaseq)[1]<-"previousID"
colnames(rnaid)[1]<-"previousID"

rnaiddup<-as.character(rnaid$converted_id[duplicated(rnaid$converted_id)]) #empty
rnaidspl<-as.character(rnaid$converted_id[grep("::", rnaid$converted_id)]) #empty
a<-as.character(rnaid$previousID[duplicated(rnaid$previousID)]) #FBgn0036125 is that mean the new way of split?
a1<-subset(rnaid,rnaid$previousID %in% a)
a2<-as.character(a1$converted_id) #"FBgn0286979" "FBgn0286980"

rnaidremove<-c(rnaiddup,rnaidspl,a2)

rnaidclearid<-as.character(setdiff(rnaid$converted_id,rnaidremove))
rnaidclear<-subset(rnaid,rnaid$converted_id %in% rnaidclearid) #7585 objects

rnafinal<-left_join(rnaidclear,fbrnaseq,by="previousID") 

rnaextract<-rnafinal[,c(2,6,11)] #only output flybase ID, Log2fd, padj
rnaextrac1<-rnaextract[complete.cases(rnaextract),] #7585 objects
colnames(rnaextrac1)[1:3]<-c("ID","status","value")

write.csv(rnaextrac1, paste(getwd(),"/","2014_fb.csv",sep=''),row.names = FALSE)

#############nazif-2014-rnaseq#############################################
#gut
fbrnaseq<-read.csv(paste(getwd(),"/","nazifgutrename.csv",sep=''))
rnaid<-read.csv(paste(getwd(),"/","2014_gut_FlyBase_IDs.csv",sep='')) 

colnames(fbrnaseq)[1]<-"previousID"
colnames(rnaid)[1]<-"previousID"

rnaiddup<-as.character(rnaid$converted_id[duplicated(rnaid$converted_id)]) #empty
rnaidspl<-as.character(rnaid$converted_id[grep("::", rnaid$converted_id)]) #empty
a<-as.character(rnaid$previousID[duplicated(rnaid$previousID)]) #empty
a1<-subset(rnaid,rnaid$previousID %in% a)
a2<-as.character(a1$converted_id) 

rnaidremove<-c(rnaiddup,rnaidspl,a2)

rnaidclearid<-as.character(setdiff(rnaid$converted_id,rnaidremove))
rnaidclear<-subset(rnaid,rnaid$converted_id %in% rnaidclearid) #7388 objects

rnafinal<-left_join(rnaidclear,fbrnaseq,by="previousID") 

rnaextract<-rnafinal[,c(3,6,11)] #only output flybase ID, Log2fd, padj
rnaextrac1<-rnaextract[complete.cases(rnaextract),] #7388 objects
colnames(rnaextrac1)[1:3]<-c("ID","status","value")

write.csv(rnaextrac1, paste(getwd(),"/","2014_gut.csv",sep=''),row.names = FALSE)


#############adam-2019-rnaseq#############################################
#fb
fbrnaseq<-read.csv(paste(getwd(),"/","2019_foxo_fb.csv",sep=''))
rnaid<-read.csv(paste(getwd(),"/","2019_fb_FlyBase_IDs.csv",sep='')) 

colnames(fbrnaseq)[1]<-"previousID"
colnames(rnaid)[1]<-"previousID"

rnaiddup<-as.character(rnaid$converted_id[duplicated(rnaid$converted_id)]) #"FBgn0286778" "FBgn0286809"
rnaidspl<-as.character(rnaid$converted_id[grep("::", rnaid$converted_id)]) #empty
a<-as.character(rnaid$previousID[duplicated(rnaid$previousID)]) #FBgn0036125
a1<-subset(rnaid,rnaid$previousID %in% a)
a2<-as.character(a1$converted_id) # "FBgn0286979" "FBgn0286980"

rnaidremove<-c(rnaiddup,rnaidspl,a2)

rnaidclearid<-as.character(setdiff(rnaid$converted_id,rnaidremove))
rnaidclear<-subset(rnaid,rnaid$converted_id %in% rnaidclearid) #11064 objects

rnafinal<-left_join(rnaidclear,fbrnaseq,by="previousID") 

rnaextract<-rnafinal[,c(3,6,9)] #only output flybase ID, Log2fd, padj
rnaextrac1<-rnaextract[complete.cases(rnaextract),] #11012 objects
colnames(rnaextrac1)[1:3]<-c("ID","status","value")

write.csv(rnaextrac1, paste(getwd(),"/","2019_fb.csv",sep=''),row.names = FALSE)

#############adam-2019-rnaseq#############################################
#gut
fbrnaseq<-read.csv(paste(getwd(),"/","2019_foxo_gut.csv",sep=''))
rnaid<-read.csv(paste(getwd(),"/","2019_gut_FlyBase_IDs.csv",sep='')) 

colnames(fbrnaseq)[1]<-"previousID"
colnames(rnaid)[1]<-"previousID"

rnaiddup<-as.character(rnaid$converted_id[duplicated(rnaid$converted_id)]) #"FBgn0286778" "FBgn0286809"
rnaidspl<-as.character(rnaid$converted_id[grep("::", rnaid$converted_id)]) #empty
a<-as.character(rnaid$previousID[duplicated(rnaid$previousID)]) #empty
a1<-subset(rnaid,rnaid$previousID %in% a)
a2<-as.character(a1$converted_id) #empty

rnaidremove<-c(rnaiddup,rnaidspl,a2)

rnaidclearid<-as.character(setdiff(rnaid$converted_id,rnaidremove))
rnaidclear<-subset(rnaid,rnaid$converted_id %in% rnaidclearid) #10362 objects

rnafinal<-left_join(rnaidclear,fbrnaseq,by="previousID") 

rnaextract<-rnafinal[,c(3,6,9)] #only output flybase ID, Log2fd, padj
rnaextrac1<-rnaextract[complete.cases(rnaextract),] #10346 objects
colnames(rnaextrac1)[1:3]<-c("ID","status","value")

write.csv(rnaextrac1, paste(getwd(),"/","2019_gut.csv",sep=''),row.names = FALSE)
