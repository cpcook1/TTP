library(Seurat)
library(dplyr)
library(xlsx)
library(iCellR)
library(cowplot)

Immune <- readRDS("G:/Immune.rds")
DimPlot(Immune)
unique(Immune$ident)
#[1] DETC  cDC2  d¦ÃdT  Mac   ILC2  Thet  M/MdM mDC   LC    M/B   cDC1  NK    Neu  
#####Levels: Mac DETC M/MdM d¦ÃdT cDC2 ILC2 Thet M/B Neu cDC1 LC mDC NK

gdT <- subset(Immune, idents = "d¦ÃdT")

VEH <- subset(gdT, stim == "VEH")
IMQ <- subset(gdT, stim == "IMQ")

objectName=VEH
#objectName=IMQ

my.data <- as.data.frame(as.matrix(objectName@assays$RNA@data))
myUMAP <- as.data.frame(objectName@reductions$umap@cell.embeddings)
myTSNE <- as.data.frame(objectName@reductions$tsne@cell.embeddings)
myPCA <- as.data.frame(objectName@reductions$pca@cell.embeddings)
myCond <- as.data.frame(objectName@meta.data)

dim(my.data)
dim(myTSNE)
dim(myUMAP)

head(my.data)[1:5]
head(myUMAP)

MyCols <- row.names(myUMAP)
MyCols <- gsub("-",".",MyCols)
MyCols <- gsub("_",".",MyCols)
MyCols <- paste0(objectName@meta.data$tech,"_",MyCols,sep="")

colnames(my.data) <- MyCols
rownames(myPCA) <- MyCols
rownames(myUMAP) <- MyCols

my.obj <- make.obj(my.data)
my.obj@main.data <- my.obj@raw.data
my.obj
my.obj@umap.data <- myUMAP
my.obj@tsne.data <- myTSNE
my.obj@pca.data <- myPCA
my.obj@metadata <- myCond

my.obj <- run.pc.tsne(my.obj, dims = 1:20)

### get the clusters 
MyClust <- objectName@meta.data$new_ID
MyClust <- as.numeric(MyClust)
MyClust <- as.data.frame(MyClust)
colnames(MyClust) <- "clusters"
rownames(MyClust) <- rownames(myPCA)

my.obj@best.clust <- MyClust
my.obj

my.obj <- run.impute(my.obj, dims = 1:10, nn = 10, data.type = "pca")

###########CORRELATION ANALYSIS FOR TWO FEATURES
#PsoI <-my.obj
NMLI <- my.obj

library(ggpubr)
Pso <- as.data.frame(t(PsoI@imputed.data))
Nml <- as.data.frame(t(NMLI@imputed.data))

Pso$Status <- "IMQ"
Nml$Status <- "VEH"

library(plyr)
library(dplyr)

my_data <- bind_rows(as.data.frame(Pso), as.data.frame(Nml))
my_data <- as.data.frame(my_data)


ggplot(my_data, aes(x=Cxcl13, y=Zfp36l2, color=Status)) +
  geom_point() +
  scale_color_manual(values=c("brown4", "midnightblue", "coral4"))+
  #scale_color_brewer(palette="Dark2") + 
  theme_classic()+
  stat_cor(aes(x=Cxcl13, y=Zfp36l2),inherit.aes = FALSE, label.y = 3)+
  ylim(0, 3)

ggplot(my_data, aes(x=Il22, y=Zfp36l2, color=Status)) +
  geom_point() +
  scale_color_manual(values=c("brown4", "midnightblue", "coral4"))+
  #scale_color_brewer(palette="Dark2") + 
  geom_smooth(method=lm,se=FALSE, col='black', size=1) +
  theme_classic()+
  stat_cor(aes(x=Il22, y=Zfp36l2),inherit.aes = FALSE, label.y = 3.3)+
  ylim(0, 3.5)

ggplot(my_data, aes(x=Il26, y=Zfp36l2, color=Status)) +
  geom_point() +
  scale_color_manual(values=c("brown4", "midnightblue", "coral4"))+
  #scale_color_brewer(palette="Dark2") + 
  theme_classic()+
  stat_cor(aes(x=Il26, y=Zfp36l2),inherit.aes = FALSE, label.y = 0.8)+
  ylim(0, 0.8)

ggplot(my_data, aes(x=Ifng, y=Zfp36l2, color=Status)) +
  geom_point() +
  scale_color_manual(values=c("brown4", "midnightblue", "coral4"))+
  #scale_color_brewer(palette="Dark2") + 
  theme_classic()+
  stat_cor(aes(x=Ifng, y=Zfp36l2),inherit.aes = FALSE, label.y = 3.3)+
  ylim(0, 3.5)

ggplot(my_data, aes(x=Tnf, y=Zfp36l2, color=Status)) +
  geom_point() +
  scale_color_manual(values=c("brown4", "midnightblue", "coral4"))+
  #scale_color_brewer(palette="Dark2") + 
  theme_classic()+
  stat_cor(aes(x=Tnf, y=Zfp36l2),inherit.aes = FALSE, label.y = 3.3)+
  ylim(0, 3.5)
