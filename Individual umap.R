library(Seurat)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(openxlsx)

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
Idents(Tpso) <- Tpso$ID
DimPlot(Tpso, reduction = "umap", label = T, pt.size = .1, split.by = "Sample", ncol = 4)+NoLegend()+
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(17))+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family= "italic",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        plot.title=element_text(color="black", size=14, face="bold.italic"),
        legend.text=element_text(family= "italic", size=12),
        legend.key=element_blank())