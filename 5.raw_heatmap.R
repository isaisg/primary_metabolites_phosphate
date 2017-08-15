library(gplots)
library(RColorBrewer)
library(vegan)
setwd('/home/isai/Documents/results/primary_metabolites_phosphate_context/primary_metabolites_phosphate_context/')
figures <- "figures/"

filename_Tab <- "raw_rep2.tsv"
filename_Map <- "metadata_phosphate_metabolites_oct27_2016.tsv"
Tab <- as.matrix(read.table(file = filename_Tab,header = T,sep = "\t",row.names = 1,check.names = F))
Map <- read.table(file = filename_Map,header=T,sep = "\t",check.names = F)

Map$Specific<-factor(Map$Specific,levels=c("-","+","-/-","-/+","+/-","+/+"))



Tab <-Tab[,match(Map$Colname,colnames(Tab))]


# #Define the colors to use 
#Define colors
#Define colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
paleta = gg_color_hue(n)
#Colsidecolors
paleta<-c("gray","gray",paleta[c(2,1,3,4)])
####ColSideColors
colsidecolors<- paleta[as.numeric(Map$Specific)]


#Heatmap.2 
paletahm<-(brewer.pal(n = 9,name = "PiYG"))

# #Corrplot
Tab.normalized <-Tab
Tab.normalized.cor<-1-(cor(t(Tab.normalized)))
Tab.normalized.cor.hclust<-hclust(as.dist(Tab.normalized.cor),method = "ward.D")
compounds.order<-Tab.normalized.cor.hclust$labels[Tab.normalized.cor.hclust$order]
# #Now determine the order of the columns 
Tab.normalized.distsamples<-hclust(vegdist(t(Tab.normalized),method = "bray"),method = "ward.D")

samples.order<-Tab.normalized.distsamples$labels[Tab.normalized.distsamples$order]

Tab.normalized.ordered<-Tab.normalized[match(compounds.order,rownames(Tab.normalized)),
                                       match(samples.order,colnames(Tab.normalized))]

filename <- paste(figures,"heatmap_raw_plants_allmetabolites_metabolitesbycor_samplesbyeuclidean_ward.svg",sep = "")
svg(filename =filename,width = 12,height = 12,pointsize = 12 )
heatmap.2(x = Tab.normalized,trace = "none",
          col =paletahm,Rowv = as.dendrogram(Tab.normalized.cor.hclust),
          Colv = rev(as.dendrogram(Tab.normalized.distsamples)),scale = "row",dendrogram = "both"
          ,cexRow = 0.2,density.info = "none",ColSideColors = colsidecolors,cexCol = 1,margins = c(10,10))
dev.off()

########Only known compounds
Tab<-Tab[1:76,]
Tab.normalized <-Tab
Tab.normalized.cor<-1-(cor(t(Tab.normalized)))
Tab.normalized.cor.hclust<-hclust(as.dist(Tab.normalized.cor),method = "ward.D")
compounds.order<-Tab.normalized.cor.hclust$labels[Tab.normalized.cor.hclust$order]
# #Now determine the order of the columns 
Tab.normalized.distsamples<-hclust(vegdist(t(Tab.normalized),method = "bray"),method = "ward.D")

samples.order<-Tab.normalized.distsamples$labels[Tab.normalized.distsamples$order]

Tab.normalized.ordered<-Tab.normalized[match(compounds.order,rownames(Tab.normalized)),
                                       match(samples.order,colnames(Tab.normalized))]

filename <- paste(figures,"heatmap_raw_exudates_knownmetabolites_metabolitesbycor_samplesbyeuclidean_ward.svg",sep = "")
svg(filename =filename,width = 12,height = 12,pointsize = 12 )
heatmap.2(x = Tab.normalized,trace = "none",
          col =paletahm,Rowv = as.dendrogram(Tab.normalized.cor.hclust),
          Colv = rev(as.dendrogram(Tab.normalized.distsamples)),scale = "row",dendrogram = "both"
          ,cexRow = 0.7,density.info = "none",ColSideColors = colsidecolors,cexCol = 1,margins = c(10,10))

dev.off()

