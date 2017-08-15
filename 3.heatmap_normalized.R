library(corrplot)
library(RColorBrewer)
library(gplots)
#Read data
setwd('/home/isai/Documents/results/primary_metabolites_phosphate_context/primary_metabolites_phosphate_context/')
figures <- "figures/"

Tab.normalized <- as.matrix(read.table("normalized_rep2.tsv",header = T,
                                       sep = "\t",row.names = 1,check.names = F))
Map.normalized <- read.table("normalized_metadata_phosphate_metabolites_oct27_2016.tsv",header=T,
                             sep="\t",check.names = F)

#Match
Tab.normalized <-Tab.normalized[,match(Map.normalized$Colname,colnames(Tab.normalized))]

##Colors
#Define colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
paleta = gg_color_hue(n)
#Colsidecolors
paleta<-paleta[c(2,1,3,4)]
colsidecolors<-paleta[as.numeric(Map.normalized$Specific)]

# #Corrplot
Tab.normalized.cor<-1-(cor(t(Tab.normalized)))
Tab.normalized.cor.hclust<-hclust(as.dist(Tab.normalized.cor),method = "ward.D")
#plot(Tab.normalized.cor.hclust)
compounds.order<-Tab.normalized.cor.hclust$labels[Tab.normalized.cor.hclust$order]
# #Now determine the order of the columns 
Tab.normalized.distsamples<-hclust(dist(t(Tab.normalized)),method = "ward.D")

samples.order<-Tab.normalized.distsamples$labels[Tab.normalized.distsamples$order]
#Heatmap.2 test

paletahm<-(brewer.pal(n = 9,name = "PiYG"))


filename <- paste(figures,"heatmap_normalized_plants_allmetabolites_metabolitesbycor_samplesbyeuclidean_ward.svg",sep = "")
svg(filename =filename,width = 12,height = 12,pointsize = 12 )
heatmap.2(x = Tab.normalized,trace = "none",
          col =paletahm,Rowv = as.dendrogram(Tab.normalized.cor.hclust),
          Colv = as.dendrogram(Tab.normalized.distsamples),scale = "row",dendrogram = "both"
          ,cexRow = 0.2,density.info = "none",ColSideColors = colsidecolors,cexCol = 1,margins = c(10,10))

dev.off()

filename <- paste(figures,"heatmap_normalized_plants_allmetabolites_metabolitesbycor_samplesbyeuclidean_ward.pdf",sep = "")
pdf(file=filename,width = 12,height = 12,pointsize = 12 )
heatmap.2(x = Tab.normalized,trace = "none",
          col =paletahm,Rowv = as.dendrogram(Tab.normalized.cor.hclust),
          Colv = as.dendrogram(Tab.normalized.distsamples),scale = "row",dendrogram = "both"
          ,cexRow = 0.2,density.info = "none",ColSideColors = colsidecolors,cexCol = 1,margins = c(10,10))

dev.off()




######Only characterizable compounds
Tab.normalized <- as.matrix(read.table("normalized_rep2.tsv",header = T,
                                       sep = "\t",row.names = 1,check.names = F))
Map.normalized <- read.table("normalized_metadata_phosphate_metabolites_oct27_2016.tsv",header=T,
                             sep="\t",check.names = F)

#Match
Tab.normalized <-Tab.normalized[,match(Map.normalized$Colname,colnames(Tab.normalized))]

Tab.normalized<-Tab.normalized[1:76,]
colsidecolors<-paleta[as.numeric(Map.normalized$Specific)]

# #Corrplot
Tab.normalized.cor<-1-(cor(t(Tab.normalized)))

Tab.normalized.cor.hclust<-hclust(as.dist(Tab.normalized.cor),method = "ward.D")
compounds.order<-Tab.normalized.cor.hclust$labels[Tab.normalized.cor.hclust$order]
# #Now determine the order of the columns 
Tab.normalized.distsamples<-hclust(dist(t(Tab.normalized)),method = "ward.D")
samples.order<-Tab.normalized.distsamples$labels[Tab.normalized.distsamples$order]
#Heatmap.2 test
paletahm<-(brewer.pal(n = 9,name = "PiYG"))


filename <- paste(figures,"heatmap_normalized_exudates_knownmetabolites_metabolitesbycor_samplesbyeuclidean_ward.svg",sep = "")
svg(filename =filename,width = 12,height = 12,pointsize = 12 )
heatmap.2(x = Tab.normalized,trace = "none",
          col =paletahm,Rowv = as.dendrogram(Tab.normalized.cor.hclust),
          Colv = as.dendrogram(Tab.normalized.distsamples),scale = "row",dendrogram = "both"
          ,cexRow = 0.7,density.info = "none",ColSideColors = colsidecolors,cexCol = 1,margins = c(10,10))

dev.off()

filename <- paste(figures,"heatmap_normalized_exudates_knownmetabolites_metabolitesbycor_samplesbyeuclidean_ward.pdf",sep = "")
pdf(file =filename,width = 12,height = 12,pointsize = 12 )
heatmap.2(x = Tab.normalized,trace = "none",
          col =paletahm,Rowv = as.dendrogram(Tab.normalized.cor.hclust),
          Colv = as.dendrogram(Tab.normalized.distsamples),scale = "row",dendrogram = "both"
          ,cexRow = 0.7,density.info = "none",ColSideColors = colsidecolors,cexCol = 1,margins = c(10,10))

dev.off()


