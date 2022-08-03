library(Seurat)
library(SeuratObject)
library(openxlsx)
library(iCellR)


Tmix <- readRDS("G:/Chris_ZFP36/Tmix_newID.rds")

Tmix <- subset(Tmix, percent.mt<20)
VlnPlot(Tmix, features = "percent.mt")

DimPlot(Tmix)
Idents(Tmix) <- Tmix$donor
Tpso <- subset(Tmix, idents=c("150","154","155","169","195","204","207","165","173","194","199","211","222","234","235"))
Tpso<- RenameIdents(object = Tpso, `150` = "HC1", `154` = "HC2", `155` = "HC3", `169` = "HC4", `195` = "HC5", `204` = "HC6", `207` = "HC7",
                    `165` = "Pso1",`173` = "Pso2", `194` = "Pso3", `199` = "Pso4", `211` = "Pso5",`222` = "Pso6", `234` = "Pso7",`235` = "Pso8")
DimPlot(Tpso, reduction = "umap", label = T, pt.size = .1)+NoLegend()
Tpso[["Sample"]] <- Idents(object = Tpso)
Idents(Tpso) <- Tpso$ID################22
library("psych")

Pso1 <- subset(Tpso, dis=="Pso")
NML <- subset(Tpso, dis=="NML")

Pso1 <- subset(Pso1, idents="2")
NML <- subset(NML, idents="2")

#######################Pso
objectName=Pso1
my.data <- as.data.frame(as.matrix(objectName@assays[["RNA"]]@data))
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

remove<- c("ABCF2.1", "ATXN7.1", "CCDC39.1", "COG8.1", "CYB561D2.1", "DIABLO.1", "EMG1.1", "H2BFS.1", "HSPA14.1",
           "IGF2.1", "LINC01238.1", "LINC01505.1", "LINC02203.1", "MATR3.1", "PDE11A.1", "PINX1.1", "POLR2J3.1", 
           "PRSS50.1", "RGS5.1", "SCO2.1", "SOD2.1", "TBCE.1", "TMSB15B.1", "TXNRD3NB.1")

my.data = my.data[!row.names(my.data)%in%remove, ]

my.obj <- make.obj(my.data)
my.obj@main.data <- my.obj@raw.data
my.obj
my.obj@umap.data <- myUMAP
my.obj@tsne.data <- myTSNE
my.obj@pca.data <- myPCA
my.obj@metadata <- myCond

my.obj <- run.pc.tsne(my.obj, dims = 1:20)

### get the clusters 
MyClust <- objectName@meta.data$Ident2
MyClust <- as.numeric(MyClust)
MyClust <- as.data.frame(MyClust)
colnames(MyClust) <- "clusters"
rownames(MyClust) <- rownames(myPCA)

my.obj@best.clust <- MyClust
my.obj

my.obj <- run.impute(my.obj, dims = 1:10, nn = 10, data.type = "pca")

PsoI <-my.obj


#######################NML
objectName=NML
my.data <- as.data.frame(as.matrix(objectName@assays[["RNA"]]@data))
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

remove<- c("ABCF2.1", "ATXN7.1", "CCDC39.1", "COG8.1", "CYB561D2.1", "DIABLO.1", "EMG1.1", "H2BFS.1", "HSPA14.1",
           "IGF2.1", "LINC01238.1", "LINC01505.1", "LINC02203.1", "MATR3.1", "PDE11A.1", "PINX1.1", "POLR2J3.1", 
           "PRSS50.1", "RGS5.1", "SCO2.1", "SOD2.1", "TBCE.1", "TMSB15B.1", "TXNRD3NB.1")

my.data = my.data[!row.names(my.data)%in%remove, ]

my.obj <- make.obj(my.data)
my.obj@main.data <- my.obj@raw.data
my.obj
my.obj@umap.data <- myUMAP
my.obj@tsne.data <- myTSNE
my.obj@pca.data <- myPCA
my.obj@metadata <- myCond

my.obj <- run.pc.tsne(my.obj, dims = 1:20)

### get the clusters 
MyClust <- objectName@meta.data$Ident2
MyClust <- as.numeric(MyClust)
MyClust <- as.data.frame(MyClust)
colnames(MyClust) <- "clusters"
rownames(MyClust) <- rownames(myPCA)

my.obj@best.clust <- MyClust
my.obj

my.obj <- run.impute(my.obj, dims = 1:10, nn = 10, data.type = "pca")

NMLI <- my.obj

###########################COMINED AND ASSIGNED
library(ggpubr)
Pso <- as.data.frame(t(PsoI@imputed.data))
Nml <- as.data.frame(t(NMLI@imputed.data))

Nml$Status <- "NML"
Pso$Status <- "Pso"

###############################COMINED
library(plyr)
library(dplyr)
my_data <- bind_rows(as.data.frame(Pso), as.data.frame(Nml))
my_data <- as.data.frame(my_data)

ggplot(my_data, aes(x=CXCL13, y=ZFP36L2, color=Status)) +
  geom_point() +
  scale_color_manual(values=c("brown4", "midnightblue", "coral4"))+
  #scale_color_brewer(palette="Dark2") + 
  theme_classic()+
  stat_cor(aes(x=CXCL13, y=ZFP36L2),inherit.aes = FALSE, label.y = 4)+
  ylim(0, 4)

ggplot(my_data, aes(x=IL22, y=ZFP36L2, color=Status)) +
  geom_point() +
  scale_color_manual(values=c("brown4", "midnightblue", "coral4"))+
  #scale_color_brewer(palette="Dark2") + 
  theme_classic()+
  stat_cor(aes(x=IL22, y=ZFP36L2),inherit.aes = FALSE, label.y = 4)+
  ylim(0, 4)

ggplot(my_data, aes(x=IL26, y=ZFP36L2, color=Status)) +
  geom_point() +
  scale_color_manual(values=c("brown4", "midnightblue", "coral4"))+
  #scale_color_brewer(palette="Dark2") + 
  theme_classic()+
  stat_cor(aes(x=IL26, y=ZFP36L2),inherit.aes = FALSE, label.y =4)+
  ylim(0, 4)

ggplot(my_data, aes(x=IFNG, y=ZFP36L2, color=Status)) +
  geom_point() +
  scale_color_manual(values=c("brown4", "midnightblue", "coral4"))+
  #scale_color_brewer(palette="Dark2") + 
  theme_classic()+
  stat_cor(aes(x=IFNG, y=ZFP36L2),inherit.aes = FALSE, label.y = 4)+
  ylim(0, 4)

ggplot(my_data, aes(x=TNF, y=ZFP36L2, color=Status)) +
  geom_point() +
  scale_color_manual(values=c("brown4", "midnightblue", "coral4"))+
  #scale_color_brewer(palette="Dark2") + 
  theme_classic()+
  stat_cor(aes(x=TNF, y=ZFP36L2),inherit.aes = FALSE, label.y = 4)+
  ylim(0, 4)
