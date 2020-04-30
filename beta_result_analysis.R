library(tidyverse)
library(clusterProfiler)
library(topGO)
library(org.Dm.eg.db)
library(veccompare)
library(KEGGprofile)
library(DOSE)
library(pathview)
library(pheatmap)

#load to BETA results and flybase BATCH result --gut --chip2014 --microarray 2014
setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/BETA/result/foxo/2014gut-2014chip/test1")

#gutbeta<-read.csv("2014_chip_microarray_gut_beta.csv")
#gutbatch<-read.csv("2014_chip_microarray_gut_beta_batch.csv")

#colnames(gutbeta)<-c("score","id")
#colnames(gutbatch)[1]<-"id"

#gutbeta1<-left_join(gutbatch,gutbeta, by="id") #384 gene reported 

#write.csv(gutbeta1, paste(getwd(),"/","2014_chip_microarray_gut_beta_final.csv", sep=''), row.names = FALSE)
gutbeta1<-read.csv("2014_chip_microarray_gut_beta_final.csv") #here if us both direction the number of genes reported is:726
#load to BETA results and flybase BATCH result --fb --chip2014 --microarray 2014
setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/BETA/result/foxo/2014fb-2014chip/test1/")

fbbeta<-read.csv("2014_chip_microarray_fb_beta.csv")
fbbatch<-read.csv("2014_chip_microarray_fb_beta_batch.csv")

fbbeta<-fbbeta[,c(5,6,7)]
colnames(fbbeta)<-c("score","direction","id")
colnames(fbbatch)[1]<-"id"

fbbeta1<-left_join(fbbatch,fbbeta, by="id") #32 gene reported\

rm(fbbeta,fbbatch)

#write.csv(fbbeta1, paste(getwd(),"/","2014_chip_microarray_fb_beta_final.csv", sep=''), row.names = FALSE)



#load the Original analyze done by Nazif
###################################################################################
# ORGANISE ALL DATASETS FROM BOTH EXPERIMENTS
###################################################################################

setwd("F:/Datacollection/Nazif switch/Programming/post-alignment_analysis/RNAseq/gene_lists_sources/chip_chips/")
chip2011 <-read.csv("chipseq_genelist_alic2011.csv", fileEncoding="UTF-8-BOM", header=F)# WHOLEFLY genes bound to FOXO from 2011 chip-chip alic
chip2014 <- read.csv("chipseq_genelist_alic2014.csv", fileEncoding="UTF-8-BOM") # WHOLEFLY genes bound to FOXO from 2014 chip-chip alic

names(chip2011) <- c("id")

#common_chips <- intersect (chip2011$id, chip2014$id) # 411 GENES IN COMMON BETWEEN BOTH CHIP-CHIPS 2011-2014

setwd("F:/Datacollection/Nazif switch/Programming/post-alignment_analysis/RNAseq/gene_lists_sources/rnaseq_microarray_foxo/")
a_gut <- read_csv("microarray_foxoinduction_gut.csv") # gut genes from microarray of Nazif, UAS-foxo, 1weekold
a_fb <- read_csv("microarray_foxoinduction_fb.csv") # fb genes from microarray of Nazif, UAS-foxo, 1weekold
setwd("F:/Datacollection/Nazif switch/Programming/post-alignment_analysis/RNAseq/gene_lists_sources/final_lists_july2019/adams_RNAseq/")
r_gut <- read_csv("s106uasfoxo_gut.csv") # gut genes from RNAseq of Adam, UAS-foxo, 1weekold
r_fb <- read_csv("s106uasfoxo_fb.csv") # fb genes from RNAseq of Adam, UAS-foxo, 1weekold

###################################################################################
#TRANSFORM THE NAMES AND MAKE THEM ALL THE SAME FROM ALL DATASETS
###################################################################################
names(a_gut) <- tolower(names(a_gut))
names(r_gut) <- tolower(names(r_gut))
names(a_fb) <- tolower(names(a_fb))
names(r_fb) <- tolower(names(r_fb))

###################################################################################
#FILTER THE SIGNIFICANT ONES FROM RNASEQ AND MICROARRAY AND THEN MERGE BOTH DATASETS TO COMPARE THEM TO CHIPSEQ
###################################################################################
a_gut <- dplyr::select(a_gut, id, logfc, padj) %>%
  filter(padj<0.1) #change by mengjia
r_gut <- dplyr::select(r_gut, id, log2fc, padj) %>%
  filter(padj<0.1) #change by mengjia

names(a_gut) <- c("id", "logfc", "padj")
names(r_gut) <- c("id", "logfc", "padj")
#a_r_gut <- rbind(a_gut, r_gut)
#WE ALSO HAVE TO DELETE THE ONES THAT ARE DUPLICATED
#a_r_gut <- distinct(a_r_gut,id, .keep_all= TRUE)

a_fb <- dplyr::select(a_fb, id, logfc, padj) %>%
  filter(padj<0.1) #mengjia change to 0,2
r_fb <- dplyr::select(r_fb, id, log2fc, padj) %>%
  filter(padj<0.1) #mengjia change to 0.2

names(a_fb) <- c("id", "logfc", "padj")
names(r_fb) <- c("id", "logfc", "padj")

#a_r_fb <- rbind(a_fb, r_fb)
#WE ALSO HAVE TO DELETE THE ONES THAT ARE DUPLICATED
#a_r_fb <- distinct(a_r_fb,id, .keep_all= TRUE)

#remove useless object

###################################################################################
#WE NOW GET THE LIST OF GENES THAT ARE BOTH IN RNA-MICROARRAY AND CHIPSEQ
###################################################################################

common_genes_gut_2011 <- intersect(a_r_gut$id, chip2011$id) # 155 GENES COMPARING WITH CHIP-2011
common_genes_gut_2014 <- intersect(a_r_gut$id, chip2014$id) # 230 GENES COMPARING WITH CHIP-2014
common_common_gut <- intersect(common_genes_gut_2011, common_genes_gut_2014) # 57 COMMON GENES IN ALL 3 DATASETS

common_genes_fb_2011 <- intersect(a_r_fb$id, chip2011$id) # 194 GENES COMPARING WITH CHIP-2011
common_genes_fb_2014 <- intersect(a_r_fb$id, chip2014$id) # 264 GENES COMPARING WITH CHIP-2014
common_common_fb <- intersect(common_genes_fb_2011, common_genes_fb_2014) # 49 COMMON GENES IN ALL 3 DATASETS

common_genes_gut_14chip_14rna <-intersect(a_gut$id, chip2014$id)#221
common_genes_fb_14chip_14rna <-intersect(a_fb$id, chip2014$id)#9

##########################################################
#Then compare with the BETA gut upregulated
##########################################################
betagutcommon<-intersect(gutbeta1$id,common_genes_gut_14chip_14rna) # 141 OVERLAPPING in both list 

# IS THIS OVERLAP SIGNIFICANT?

a <- 912 # Number of gene with padj<0.2
b <- (length(gutbeta1$id)) # 384 GENES from BETA result
c <- length(common_genes_gut_14chip_14rna) # 221 from pervious overlap result
d <- length(betagutcommon) # 141 OVERLAPPING in both list

dhyper(d, b, (a-b), c, log = FALSE) # IS THIS OVERLAP SIGNIFICANT??  4.65369e-14

##########################################################
#Then compare with the BETA gut both direction
##########################################################
betagutcommon<-intersect(gutbeta1$id,common_genes_gut_14chip_14rna) # 205 OVERLAPPING in both list 

# IS THIS OVERLAP SIGNIFICANT?

a <- 912 # Number of gene with padj<0.2
b <- (length(gutbeta1$id)) # 726 GENES from BETA result
c <- length(common_genes_gut_14chip_14rna) # 221 from pervious overlap result
d <- length(betagutcommon) # 204 OVERLAPPING in both list

dhyper(d, b, (a-b), c, log = FALSE) # IS THIS OVERLAP SIGNIFICANT??  1.159007e-09

########################################################
#Then compare with the BETA fb
########################################################
betafbcommon<-intersect(fbbeta1$id,common_genes_fb_14chip_14rna) # 6 OVERLAPPING (with gut) in both list(mistaked by gutbeta1$id with common_genes_fb_14chip_14rna), 8 OVERLAPPING with fb
#Try with common_genes_fb_2014: betafbcommon1<-intersect(fbbeta1$id,common_genes_fb_2014) #still the 8 Overlapping 
#If try with common_genes_fb_2014 only 1, try with common_common_fb still 1


# IS THIS OVERLAP SIGNIFICANT?

a <- 912 # Number of gene with padj<0.2
b <- (length(fbbeta1$id)) # 32 GENES from BETA result
c <- length(common_genes_fb_14chip_14rna) # 9 from pervious overlap result
d <- length(betafbcommon) # 8 OVERLAPPING in both list

dhyper(d, b, (a-b), c, log = FALSE) # IS THIS OVERLAP SIGNIFICANT?? 8.006347e-12

########################################################
#The mistake of fb/gut i made inspire to look at fb/gut_beta
########################################################
beta_common<-intersect(fbbeta1$id,gutbeta1$id) #23 reported

# IS THIS OVERLAP SIGNIFICANT?

a <- 912 # Number of gene with padj<0.2
b <- (length(fbbeta1$id)) # 32 GENES from fbBETA result
c <- length(gutbeta1$id) # 726 GENES from gutBETA result
d <- length(beta_common) # 23 OVERLAPPING in both list

dhyper(d, b, (a-b), c, log = FALSE) # IS THIS OVERLAP SIGNIFICANT?? 0.09034418

beta_common_table<-subset(fbbeta1,fbbeta1$id %in% beta_common)
setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/BETA/result")
write.csv(beta_common_table,"fb_gut_beta_common.csv",row.names = FALSE)

#################################################################################
#NOW DEAL WITH THE SWITCH RESULT
#################################################################################
#--fb --atac --switchrnaseq --only upregulated, as downregulated only have 5 genes
#setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/BETA/result/switch/fb/test1/raw")
#switch_fb_beta<-read.csv("fb_atac_switch_beta_uptarget.csv")
#switchfbbatch<-read.csv("fb_atac_switch_beta_uptarget_batch.csv")

#switch_fb_beta<-switch_fb_beta[,c(5,7)]

#colnames(switch_fb_beta)<-c("score","id")
#colnames(switchfbbatch)[1]<-"id"

#switch_fb_beta1<-left_join(switchfbbatch,switch_fb_beta, by="id") #211 gene reported with upregulated, 216 with both direction
#write.csv(switch_fb_beta1, paste(getwd(),"/","fb_atac_switch_beta_uptarget_final.csv", sep=''), row.names = FALSE)

setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/BETA/result/switch/fb/test1")
switch_fb_beta1<-read.csv("fb_atac_switch_beta_final.csv") #216 gene reported, meanwhile in the microarray results, padj<0.2 have 912

#--gut --atac --switchrnaseq --both upregulated/downregulated
setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/BETA/result/switch/gut/sigp")
switch_gut_beta<-read.csv("gut_atac_sigp.csv")
switchgutbatch<-read.csv("gut_atac_sigp_batch.csv")

switch_gut_beta<-switch_gut_beta[,c(5,6,7)]
colnames(switch_gut_beta)<-c("score","direction","id")
colnames(switchgutbatch)[1]<-"id"

switch_gut_beta1<-left_join(switchgutbatch,switch_gut_beta, by="id") #17 gene reported

rm(switch_gut_beta,switchgutbatch)
#write.csv(switch_gut_beta1, paste(getwd(),"/","gut_atac_switch_beta_final.csv", sep=''), row.names = FALSE)


#According to NAZIF suggestion change the padj=0.1 on microarray DE
#setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/BETA/result/switch/fb/padj01/raw")

#switch_fb_beta_01<-read.csv("fb_atac_switch_beta_padj01.csv")
#switchfbbatch01<-read.csv("fb_atac_switch_beta_padj01_batch.csv")

#switch_fb_beta_01<-switch_fb_beta_01[,c(5,6,7)]

#colnames(switch_fb_beta_01)<-c("score","direction","id")
#colnames(switchfbbatch01)[1]<-"id"

#switch_fb_beta_01_1<-left_join(switchfbbatch01,switch_fb_beta_01, by="id") #155 genes both direction, |#461 gene padj<0.1 in the microarray
#write.csv(switch_fb_beta_01_1, paste(getwd(),"/","fb_atac_switch_beta_01_final.csv", sep=''), row.names = FALSE)

setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/BETA/result/switch/fb/padj01")
switch_fb_beta_01_1<-read.csv("fb_atac_switch_beta_padj01_final.csv")
#################################################################################
#COMBINATION STUDY OF SWITCH AND FOXO LONG-TERM
################################################################################# 

#Gut_2014chip_2014microarray_atac_switchrnaseq
gutcommon<-intersect(switch_gut_beta1$id, gutbeta1$id) # 2 OVERLAPPING in both list£º "FBgn0063493=Glutathione S transferase E7" "FBgn0033913=CG8468" #5 overlappping after both direction include
#Three additional (-) are "FBgn0040349=CG3699(perdicted NADPH activity"  "FBgn0024957=cytosolic aconitase" "FBgn0032775=CG17544 (Acyl-CoA oxidase)"
#So all three are involved in the tricarboxylic acid (TCA)?
#Are these three both downregulated in these two list? FBgn0040349 upregulated in switch, FBgn0024957 and FBgn0032775 both downregulated

#fb_2014chip_2014microarray_atac_switchrnaseq
fbcommon<-intersect(switch_fb_beta1$id, fbbeta1$id) # 1 OVERLAPPING in both list£º "FBgn0261985=Protein tyrosine phosphatase Meg"

ataccommon<-intersect(switch_gut_beta1$id, switch_fb_beta1) #nothing

###############################################################################
#ALSO LOOK AT THE PEAK
###############################################################################

#setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/BETA/result/switch/fb/padj01/raw")

#switch_fb_beta_01_peak<-read.csv("fb_atac_switch_beta_padj01_peak.csv")
#length(unique(switch_fb_beta_01_peak$start))#188 peaks load, 104 unique peak

#switchfbbatch01_peak<-read.csv("fb_atac_switch_beta_padj01_peak_batch.csv")

#switch_fb_beta_01_peak<-switch_fb_beta_01_peak[,c(2,3,5,6,7)]

#colnames(switch_fb_beta_01_peak)<-c("start","end","id","distance","score")
#colnames(switchfbbatch01_peak)[1]<-"id"

#switch_fb_beta_01_peak_1<-left_join(switchfbbatch01_peak,switch_fb_beta_01_peak, by="id") 

#switch_fb_beta_01_peak_1_unique<-switch_fb_beta_01_peak_1[!duplicated(switch_fb_beta_01_peak_1$distance),] #FBgn0013949 was deleted due to same distance
#deleted<-subset(switch_fb_beta_01_peak_1, switch_fb_beta_01_peak_1$id == "FBgn0013949")
#deleted<-as.vector(deleted[1,])
#switch_fb_beta_01_peak_1_unique<-rbind(switch_fb_beta_01_peak_1_unique,deleted)#still 188 peaks, good
#write.csv(switch_fb_beta_01_peak_1_unique, paste(getwd(),"/","fb_atac_switch_beta_padj01_peak_final.csv", sep=''), row.names = FALSE)

setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/BETA/result/switch/fb/padj01")
switch_fb_beta_01_peak_1_unique<-read.csv("fb_atac_switch_beta_padj01_peak_final.csv")
###############################################################################
#Focusd on several group of genes i found interesting
##############################################################################
#purine<-c("FBgn0000052", "FBgn0003204", "FBgn0039241", "FBgn0000053", "FBgn0032781", "FBgn0004888","FBgn0027493") #Purine metabolism related gene in fb_switch_atac results, padj0.2
purine<-c("FBgn0000052",  "FBgn0039241", "FBgn0000053", "FBgn0032781", "FBgn0004888","FBgn0027493") #Purine metabolism related gene in fb_switch_atac results, padj0.1
ribosome<-c("FBgn0004867", "FBgn0003279", "FBgn0015288", "FBgn0014026", "FBgn0029095", "FBgn0000100", "FBgn0016978", "FBgn0034497", "FBgn0028646")#Purine metabolism related gene in fb_switch_atac results,

purine1<- bitr(purine, fromType="FLYBASE", toType="UNIPROT", OrgDb="org.Dm.eg.db")
ribosome1<- bitr(ribosome, fromType="FLYBASE", toType="FLYBASECG", OrgDb="org.Dm.eg.db")

#####################################################################
#The overlap between these with previous RNA-seq list
#####################################################################
#load the gene 
#Go back to load the foxo microarray/foxo rnaseq sig in both fat body and gut
setwd("F:/Datacollection/Nazif switch/Programming/sugar")
sugar_before<-read.csv("DE_8xsugar_before.csv") #8x sugar before switch(1week), whole flies, 6254DE<0.1padj,r6 , here abondon those duplicated ID
sugar_after<-read.csv("DE_8xsugar_after.csv") #8x sugar after switch(1week+1week), whole flies, 64DE<0.1padj,r6
colnames(sugar_before)[1]<-"id"
colnames(sugar_after)[1]<-"id"


#purine intersect with previous lists 
a_gut_purine<-intersect(purine,a_gut$id) #"FBgn0000052" "FBgn0027493"
a_fb_purine<-intersect(purine,a_fb$id) #null
r_gut_purine<-intersect(purine,r_gut$id) #null
r_fb_purine<-intersect(purine,r_fb$id) #null
sugar_before_purine<-intersect(purine,sugar_before$id)#"FBgn0032781"
sugar_after_purine<-intersect(purine,sugar_after$id)#null

#intersect of all fb-beta padj0.1 perdicted target with previous list.
a_gut_beta<-intersect(switch_fb_beta_01_1$id,a_gut$id)#14 ("FBgn0086450" "FBgn0053144" "FBgn0053080" "FBgn0027493" "FBgn0000052" "FBgn0031069" "FBgn0014455" "FBgn0003507" "FBgn0037720" "FBgn0034501" "FBgn0037057" "FBgn0010100" "FBgn0020545"
#"FBgn0039507)
a_fb_beta<-intersect(switch_fb_beta_01_1$id,a_fb$id) #null
r_gut_beta<-intersect(switch_fb_beta_01_1$id,r_gut$id) #null
r_fb_beta<-intersect(switch_fb_beta_01_1$id,r_fb$id) #9 "FBgn0037146" "FBgn0026428" "FBgn0003279" "FBgn0027560" "FBgn0004885" "FBgn0035812" "FBgn0035811" "FBgn0004867" "FBgn0013949"
sugar_before_beta<-intersect(switch_fb_beta_01_1$id,sugar_before$id)#82 (of course, since 6254 gene) 
sugar_after_beta<-intersect(switch_fb_beta_01_1$id,sugar_after$id)#null
#Seem like not very meaningful?



############################################################
#PLOT THE BETA-PERDICTED PEAK AND BETA-PERDICTED GENE
############################################################





overlap<-list()
ID2<-vector()
object2<-vector()

for (i in 1:length(whichcomponent)){           #choose either whichgrange or whichcomponent, and change three position
  num<-whichcomponent[i]
  overlapped<-findOverlaps(foxogRange, filetableGrange[[num]])
  overlap[[i]]<-as.vector(as.character(overlapped@from))
  ID2[i]<-ID[num]
  object2[i]<-object1[num]
}

overlapdata<-data.frame(matrix(NA, nrow=length(whichcomponent), ncol=length(foxogRange)))
a<-as.vector(rep(0, times=length(foxogRange)))


for (i in 1:length(overlap)){
  overlapdata[i,]<-replace(a, as.numeric(overlap[[i]]), 1)
}

colnames(overlapdata)<-c(1:length(foxogRange))
overlapdata<-cbind(ID2,object2,overlapdata) 

overlapdata1<-overlapdata[order(pmatch(overlapdata$object2,complex, duplicates.ok = TRUE)),, drop=FALSE]

forheatmap<-as.matrix(overlapdata1[,-c(1,2)])
rownames(forheatmap)<-overlapdata1$object2
colnames(forheatmap)<-foxopeakid


pheatmap(forheatmap, color=c( "gray70","violetred2"), fontsize=9, fontsize_row=6, cellwidth = 7, cellheight = 5, cluster_cols=TRUE, cluster_rows=FALSE)






















#############################################################################
#Copy the guillermo's  GO analysis
############################################################################

#switch_FB_padj0.1
setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/BETA/result/switch/fb/padj01")
switch_fb_beta_01_1<-read.csv("fb_atac_switch_beta_01_final.csv")

switchfbcg <- bitr(switch_fb_beta_01_1$id, fromType="FLYBASE", toType="FLYBASECG", OrgDb="org.Dm.eg.db")
switchfbcg2 <- switch_fb_beta_01_1 %>%left_join(switchfbcg, c("id" = "FLYBASE"))
switchfbcg3 <- as.vector(t(switchfbcg2[,4]))
names(switchfbcg3)<-switchfbcg2$FLYBASECG
switchfbcg3<-switchfbcg3[complete.cases(switchfbcg3)]


topDiffGenes <- function (allScore) 
{
  return(allScore < 0.1)
}


###################################
#build the topGO objects
#calculate the enriched GO terms for their BP, MF or CC, and their significance using Fisher's exact test
###################################

BP_switch <- new("topGOdata",
                 description = "Simple session", ontology = "BP",
                 allGenes = switchfbcg3, geneSel = topDiffGenes,
                 nodeSize = 10,
                 annot=annFUN.org, mapping="org.Dm.eg.db", ID = "alias")

BP_switch

Fisher_BP_switch <- runTest(BP_switch, algorithm = "classic", statistic = "fisher")

topBP_switch <- GenTable(BP_switch, classicFisher = Fisher_BP_switch, topNodes = 100)  #if genes number too low, may not have 20
topBP_switch

MF_switch <- new("topGOdata",
                 description = "Simple session", ontology = "MF",
                 allGenes = switchfbcg3, geneSel = topDiffGenes,
                 nodeSize = 10,
                 annot=annFUN.org, mapping="org.Dm.eg.db", ID = "alias")

MF_switch

Fisher_MF_switch <- runTest(MF_switch, algorithm = "classic", statistic = "fisher")

topMF_switch <- GenTable(MF_switch, classicFisher = Fisher_MF_switch, topNodes = 20)
topMF_switch

CC_switch <- new("topGOdata",
                 description = "Simple session", ontology = "CC",
                 allGenes = switchfbcg3, geneSel = topDiffGenes,
                 nodeSize = 10,
                 annot=annFUN.org, mapping="org.Dm.eg.db", ID = "alias")

CC_switch

Fisher_CC_switch <- runTest(CC_switch, algorithm = "classic", statistic = "fisher")

topCC_switch <- GenTable(CC_switch, classicFisher = Fisher_CC_switch, topNodes = 20)
topCC_switch

write.csv(topBP_switch, file = "BP_fb.csv")
write.csv(topMF_switch, file = "MF_fb.csv")
write.csv(topCC_switch, file = "CC_fb.csv")

###################################
# to identify which of your significant genes correspond to each GO category
#this is an example that says that from the general list of genes
#you have 10 genes in the first GO of each category
###################################
allGO_BP <- genesInTerm(BP_switch)
allGO_BP ["GO:0019827"]

allGO_MF <- genesInTerm(MF_switch)
allGO_MF ["GO:0019827"]

allGO_CC <- genesInTerm(CC_switch)
allGO_CC ["GO:0005828"]
