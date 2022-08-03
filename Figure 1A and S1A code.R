library(Seurat)
library(ComplexHeatmap)
library(RColorBrewer)
library(xlsx)
library(pheatmap)

library(ggplot2)
library(ggplotify)
library(pheatmap)
library(patchwork)


psor <- c("165", "173", "194", "199", "211", "222", "234", "235")
atopic <- c("170", "198", "230", "231", "232", "236")
normal <- c("150", "154", "155", "169", "195", "204", "207")

All <- c("165", "173", "194", "199", "211", "222", "234", "235",
         "170", "198", "230", "231", "232", "236",
         "150", "154", "155", "169", "195", "204", "207")

PsorNorm <- c("165", "173", "194", "199", "211", "222", "234", "235",
              "150", "154", "155", "169", "195", "204", "207")


#Subset Donors you want
Idents(our) <- our$donor
our <- subset(our, idents = PsorNorm)




#### Subset only 6 clusters
Idents(our) <- our$ID

our2 <- subset(our, idents = c(1:6))
our12 <- subset(our, idents = c(1:8, 11, 12, 13, 14))




########################## Figure 1A
########calculate the ADT or RNA expression

PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}


######RNA data
calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}







#######ADT data
calc_helper <- function(object,genes){
  counts = object[['ADT']]@scale.data
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>1)/ncells
  }else{return(NA)}
}


genes <- c("CD3", "CD45RA", "CD4", "CD8", "CD69")
genes <- c("ITGAE", "CXCR6", "FOXP3", "TIGIT", "GZMA", "CCL5", "IL7R", "TNFRSF18", "SELL", "CCR7")

Adt <- PrctCellExpringGene(our2, genes=genes, group.by="ID")
Rna <- PrctCellExpringGene(our2, genes=genes, group.by="ID")

#CD45RA, CD4, CD8, CD69 [ADT] 
#CD103, CXCR6, FOXP3, TIGIT, GMZA, GMZB, CCL4, CCL5, IL7R, TNFRSF18, AREG [transcript]





## Write and reconfigure the dataframe (lazy way) use ADT and RNAs
write.xlsx(Adt, "/Users/chris/Desktop/New Paper/New/Figures/Figure 1/AdtPre.xlsx")
write.xlsx(Rna, "/Users/chris/Desktop/New Paper/New/Figures/Figure 1/RnaPre.xlsx")

Adt <- read.xlsx("/Users/chris/Desktop/New Paper/New/Figures/Figure 1/AdtPost.xlsx", sheetIndex = 1)
Rna <- read.xlsx("/Users/chris/Desktop/New Paper/New/Figures/Figure 1/RnaPost.xlsx", sheetIndex = 1)

rownames(Adt) <- Adt[,1]# change first column to rownames
Adt[,1] <- NULL# remove the first column
colnames(Adt) <- c("Tcm", "Trm1","Treg","Trm2",
                   "CTLex","CTLac")
               
rownames(Rna) <- Rna[,1]# change first column to rownames
Rna[,1] <- NULL# remove the first column
colnames(Rna) <- c("Tcm", "Trm1","Treg","Trm2",
                   "CTLex","CTLac")

colnames(Rna)#double check colnames




breaksList = as.numeric(seq(0, 1, by = .01))
col = colorRampPalette(brewer.pal(n = 4, name = "Reds"))(length(breaksList))

#Test  -- ended up using this one
col = colorRampPalette(brewer.pal(n = 9, name = "Reds"))(length(breaksList))


## Take out SELL and CCR7
Rna <- Rna[1:8, ]

### Swap dataframe columns
Adt <- Adt[,c(1, 2, 4,3,5,6)]
Rna <- Rna[,c(1, 2, 4,3,5,6)]


R <- as.ggplot(pheatmap(as.matrix(t(Rna)),
         kmeans_k = NA,
         color = col,
          breaks = breaksList,
         #cellheight = 10,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         #scale = "column",
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "RNA",
         legend = TRUE,
         cellwidth = 40,
         cellheight = 40
         
))

A <- as.ggplot(pheatmap(as.matrix(t(Adt)),
         kmeans_k = NA,
         color = col,
         breaks = breaksList,
         #cellheight = 10,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         #scale = "column",
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = FALSE,
         show_colnames = TRUE,
         main = "ADT",
         legend = FALSE,
         cellwidth = 40,
         cellheight = 40
         
))

A+R





########1B stacked bar plot

## Calculate frequencies per sample per cluster
Freq <- table(our2$ID, our2$donor)
Freq <- Freq[1:6, PsorNorm]

## Get table into right format
write.xlsx(Freq, "/Users/chris/Desktop/New Paper/New/Figures/Figure 1/Freq1b.xlsx")




###### FINAL

AO <- read.xlsx("/Users/chris/Desktop/New Paper/New/Figures/Figure 1/Freq1b.xlsx", sheetIndex = 1)
AO <- AO[, c(1,2,4)]
AO <- AO[1:35, ]
AO$Cluster <- as.factor(AO$Cluster)
#AO$Group <- factor(AO$Group, levels = levelsX)


#Turn your 'treatment' column into a character vector
AO$Donor <- as.character(AO$Donor)
#Then turn it back into a factor with the levels in the correct order
AO$Donor <- factor(AO$Donor, levels=unique(AO$Donor))



#aaas <- get_palette(palette = "Blues", 11)
#aaas[11] <- "#FFFFFF"

ggplot(AO, aes(fill= factor(Cluster, levels = c("Tcm", "Trm1", "Trm2", "Treg", "CTLex", "CTLac")), y=Freq2, 
               x=Donor)) + 
  scale_y_continuous()+
  geom_bar(position=position_stack(reverse = FALSE), stat="identity", color = "black")+
  scale_fill_manual(values = rev(my_colors))+ylab("Frequency")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.title.x = element_blank(), 
                     axis.text.x = element_text(face="bold", color = "black", size=14, angle = 45),
                     axis.title.y = element_text(face="bold", color = "black", size=17),
                     axis.text.y = element_text(face="bold", color = "black", size=13),
                     legend.text=element_text(face="bold", size = 9),
                     legend.title = element_blank())

#levelsX = c("Psoriasis 1", "Psoriasis 2", "Psoriasis 3", "Psoriasis 4",
#            "Psoriasis 5", "Psoriasis 6", "Psoriasis 7", "Psoriasis 8",
#           "Normal 1", "Normal 2", "Normal 3", "Normal 4",
#           "Normal 5", "Normal 6", "Normal 7")





######!!!!!!!!!!!
my_colors <- RColorBrewer::brewer.pal(11, "Accent")
my_colors <- RColorBrewer::brewer.pal(11, "Accent")
######!!!!!!!!!!!

my_colors <- colorRampPalette(brewer.pal(8, "Blues"))(11)
my_colors <- colorRampPalette(brewer.pal(9, "Blues"))(11)

myplot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


test <- get_palette(palette = "Oranges", 10)

nb.cols <- 22
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

NewMako <- c( "#0B0405FF", "#180D16FF", "#241628FF", "#2E1E3CFF", "#372850FF", "#3D3266FF", "#413D7BFF",
              "#40498EFF", "#3B5799FF", "#37659EFF", "#3573A1FF", "#3481A4FF", "#348FA7FF", "#359CAAFF",
              "#38AAACFF", "#40B7ADFF", "#4EC4ADFF", "#66D0ADFF", "#8AD9B1FF", "#AAE2BEFF", "#C6EBD1FF",
              "#FFFFFF")

levels2 <- c("0/ND", "2", "1", "-1", "-2", "-3", "-4-7", "-8-11", "-12-15", "-16-19",
             "-20-25")

Y <- AO[AO$Indel == 0, ]











