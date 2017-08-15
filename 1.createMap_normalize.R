#d rename the headers
setwd('/home/isai/Documents/results/primary_metabolites_phosphate_context/primary_metabolites_phosphate_context')
###Read data and assign new shared names to each column
# rep1<-read.csv("replicate1.csv",header=T,row.names=1)
# # 
# colnames(rep1)<-c("Bio1_Exudate+/-","Bio1_Exudate+/-","Bio1_Control-","Bio1_Exudate-/+","Bio1_Control+","Bio1_Control-","Bio1_Control-"
#                   ,"Bio1_Control+","Bio1_Exudate+/-","Bio1_Control-","Bio1_Exudate-/+","Bio1_Control+","Bio1_Control+"
#                   ,"Bio1_Exudate-/+","Bio1_Exudate-/+","Bio1_Exudate+/-")
rep2<-read.csv("replicate2.csv",header=T,row.names=1)
colnames(rep2)<-c("B1_T1_Control_+","B1_T2_Control_+","B1_T1_Control_-","B1_T2_Control_-",
                  "B2_T1_Control_+","B2_T2_Control_+","B2_T1_Control_-","B2_T2_Control_-",
                  "B1_T1_Exudate_+/+","B1_T2_Exudate_+/+","B1_T1_Exudate_-/-","B1_T2_Exudate_-/-",
                  "B1_T1_Exudate_+/-","B1_T2_Exudate_+/-","B1_T1_Exudate_-/+","B1_T2_Exudate_-/+",
                  "B2_T1_Exudate_+/+","B2_T2_Exudate_+/+","B2_T1_Exudate_-/-","B2_T2_Exudate_-/-",
                  "B2_T1_Exudate_+/-","B2_T2_Exudate_+/-","B2_T1_Exudate_-/+","B2_T2_Exudate_-/+")

#Create the map file

Map<-as.data.frame(matrix(data = unlist(strsplit(x = colnames(rep2),split = "_")),ncol = 4,byrow = T))
Map<-cbind(colnames(rep2),Map)
colnames(Map)<-c("Colname","BioRep","TechRep","Type","Specific")
#Scale the variables according to columns so transforming to relative abundance in each sample 
#scaled_rep1 <- scale(x = rep1,center=F,scale = colSums(rep1))
scaled_rep2 <- scale(x = rep2,center=F,scale = colSums(rep2))

write.table(x = scaled_rep2,file = "raw_rep2.tsv",append = F,quote = F,sep = "\t",col.names = T,row.names = T)
write.table(x = Map,file = "metadata_phosphate_metabolites_oct27_2016.tsv",append = F,quote = F,sep = "\t",col.names = T,row.names = F)


##Create the normalized versions
controlminus<-rowMeans(scaled_rep2[,which(colnames(scaled_rep2)%in%as.vector(subset(Map,Type=="Control" & Specific=="-")$Colname))])
controlplus <-rowMeans(scaled_rep2[,which(colnames(scaled_rep2)%in%as.vector(subset(Map,Type=="Control" & Specific=="+")$Colname))])

##Subset the columns of the different conditions that have a plant
plants.plus.plus<-scaled_rep2[,which(colnames(scaled_rep2)%in%as.vector(subset(Map,Type!="Control" & Specific=="+/+")$Colname))]
plants.plus.minus<-scaled_rep2[,which(colnames(scaled_rep2)%in%as.vector(subset(Map,Type!="Control" & Specific=="+/-")$Colname))]
plants.minus.plus<-scaled_rep2[,which(colnames(scaled_rep2)%in%as.vector(subset(Map,Type!="Control" & Specific=="-/+")$Colname))]
plants.minus.minus<-scaled_rep2[,which(colnames(scaled_rep2)%in%as.vector(subset(Map,Type!="Control" & Specific=="-/-")$Colname))]


normalized.plants.plus.plus<-log2(plants.plus.plus)-log2(controlplus)
normalized.plants.plus.minus<-log2(plants.plus.minus)-log2(controlminus)
normalized.plants.minus.minus<-log2(plants.minus.minus)-log2(controlminus)
normalized.plants.minus.plus<-log2(plants.minus.plus)-log2(controlplus)




rep2.normalized<-cbind(normalized.plants.minus.minus,normalized.plants.plus.plus,normalized.plants.minus.plus,normalized.plants.plus.minus)
Map.normalized<-Map[match(colnames(rep2.normalized),Map$Colname),]
rownames(Map.normalized)<-Map.normalized$Colname

write.table(x = rep2.normalized,file = "normalized_rep2.tsv",append = F,quote = F,sep = "\t",col.names = T,row.names = T)
write.table(x = Map.normalized,file = "normalized_metadata_phosphate_metabolites_oct27_2016.tsv",quote=F,sep="\t",col.names = T,row.names = F)

