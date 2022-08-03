#########New Fig 1 Volcanos
library(ggplot2)
library(ggrepel)



Chris <- readRDS("/Users/Chris/Downloads/Chris.rds")

###Volcano plots for new paper*****

################Trm1
Skin165 <- read.xlsx("/Users/chris/Desktop/DEG/AllPSORvsNorm2.xlsx", sheetIndex = 1)
Skin165 <- Skin165[, 1:6]

Skin165$threshold <- ""
Skin165$threshold <- ifelse(-log10(Skin165$p_val_adj)>0.05, TRUE, FALSE)

Skin165$genelabels <- ""
Skin165$genelabels <- ifelse(Skin165$NA.=="CCL22"
                             |Skin165$NA.=="ZFP36L2"
                             |Skin165$NA.=="ZFP36"
                             |Skin165$NA.=="LGALS1"
                             |Skin165$NA.=="SOCS1"
                             |Skin165$NA.=="SOCS3"
                             |Skin165$NA.=="CD69"
                             |Skin165$NA.=="IL17F"
                             |Skin165$NA.=="IL26"
                             |Skin165$NA.=="IFNG"
                             |Skin165$NA.=="CD3G"
                             |Skin165$NA.=="CD2"
                             |Skin165$NA.=="LAYN"
                             |Skin165$NA.=="CTLA4"
                             |Skin165$NA.=="CD82"
                             |Skin165$NA.=="CD96"
                             |Skin165$NA.=="CCR5"
                             |Skin165$NA.=="SOX4"
                             |Skin165$NA.=="ZEB2"
                             |Skin165$NA.=="CXCL13",
                             TRUE, FALSE)

Skin165$genecolour[Skin165$avg_log2FC > 0] <- '#F8766D'
Skin165$genecolour[Skin165$avg_log2FC < 0] <- '#619CFF'
Skin165$genecolour[Skin165$p_val_adj > 0.05] <- 'gray'

Skin165$order <- ""
Skin165$order <- ifelse(Skin165$NA.=="CCL22"
                        |Skin165$NA.=="ZFP36L2"
                        |Skin165$NA.=="ZFP36"
                        |Skin165$NA.=="LGALS1"
                        |Skin165$NA.=="SOCS1"
                        |Skin165$NA.=="SOCS3"
                        |Skin165$NA.=="CD69"
                        |Skin165$NA.=="IL17F"
                        |Skin165$NA.=="IL26"
                        |Skin165$NA.=="IFNG"
                        |Skin165$NA.=="CD3G"
                        |Skin165$NA.=="CD2"
                        |Skin165$NA.=="LAYN"
                        |Skin165$NA.=="CTLA4"
                        |Skin165$NA.=="CD82"
                        |Skin165$NA.=="CD96"
                        |Skin165$NA.=="CCR5"
                        |Skin165$NA.=="SOX4"
                        |Skin165$NA.=="ZEB2"
                        |Skin165$NA.=="CXCL13",
                        TRUE, FALSE)

Skin165$order <-as.numeric(factor(Skin165$order))

df_layer_1 <- Skin165[Skin165$order=="1",]
df_layer_2 <- Skin165[Skin165$order=="2",]

ggplot() +
  geom_point(data=df_layer_1, aes(avg_log2FC, -log10(p_val_adj), alpha = .9 ), colour=df_layer_1$genecolour) +
  geom_point(data=df_layer_2, aes(avg_log2FC, -log10(p_val_adj)), 
             shape = 21, colour= "black", fill = df_layer_2$genecolour, size=4, 
             stroke = .35, show.legend = FALSE)+
  
  geom_text_repel(aes(df_layer_2$avg_log2FC, -log10(df_layer_2$p_val_adj)), label= ifelse(df_layer_2$genelabels, as.character(df_layer_2$NA.), ""),
                  box.padding = unit(.42, "lines"), hjust = .9, vjust=.6, max.overlaps = 10000000)+
  theme(legend.title = element_blank(), legend.position = "none", text=element_text(size=20))+
  theme_bw()+
  geom_vline(xintercept=0, col="gray", linetype = "longdash")+
  geom_hline(yintercept=1.3, col="gray", linetype = "longdash")+
  ylim(0, 350)+
  xlim(-6, 7)+
  theme(legend.position = "none")+
  
  ggtitle("Trm Psoriasis Vs. Normal")

### Cant get GZMA or CCR7






################Tcm1
Skin165 <- read.xlsx("/Users/chris/Desktop/DEG/TcmAllPsor.xlsx", sheetIndex = 1)
Skin165 <- Skin165[, 1:6]

Skin165$threshold <- ""
Skin165$threshold <- ifelse(-log10(Skin165$p_val_adj)>0.05, TRUE, FALSE)

Skin165$genelabels <- ""
Skin165$genelabels <- ifelse(Skin165$NA.=="KLF2"
                             |Skin165$NA.=="ABLIM1"
                             |Skin165$NA.=="KLRB1"
                             |Skin165$NA.=="LTB"
                             |Skin165$NA.=="JUND"
                             |Skin165$NA.=="CD3G"
                             |Skin165$NA.=="CCR7"
                             |Skin165$NA.=="GZMA"
                             |Skin165$NA.=="CXCL8"
                             |Skin165$NA.=="ZFP36"
                             |Skin165$NA.=="ZFP36L2"
                             |Skin165$NA.=="CCL22",
                               TRUE, FALSE)

Skin165$genecolour[Skin165$avg_log2FC > 0] <- '#F8766D'
Skin165$genecolour[Skin165$avg_log2FC < 0] <- '#619CFF'
Skin165$genecolour[Skin165$p_val_adj > 0.05] <- 'gray'

Skin165$order <- ""
Skin165$order <- ifelse(Skin165$NA.=="KLF2"
                        |Skin165$NA.=="ABLIM1"
                        |Skin165$NA.=="KLRB1"
                        |Skin165$NA.=="LTB"
                        |Skin165$NA.=="JUND"
                        |Skin165$NA.=="CD3G"
                        |Skin165$NA.=="CCR7"
                        |Skin165$NA.=="GZMA"
                        |Skin165$NA.=="CXCL8"
                        |Skin165$NA.=="ZFP36"
                        |Skin165$NA.=="ZFP36L2"
                        |Skin165$NA.=="CCL22",
                        TRUE, FALSE)

Skin165$order <-as.numeric(factor(Skin165$order))

df_layer_1 <- Skin165[Skin165$order=="1",]
df_layer_2 <- Skin165[Skin165$order=="2",]

ggplot() +
  geom_point(data=df_layer_1, aes(avg_log2FC, -log10(p_val_adj), alpha = .9 ), colour=df_layer_1$genecolour) +
  geom_point(data=df_layer_2, aes(avg_log2FC, -log10(p_val_adj)), colour=df_layer_2$genecolour, size=4, show.legend = FALSE)+
  
  geom_text_repel(aes(df_layer_2$avg_log2FC, -log10(df_layer_2$p_val_adj)), label= ifelse(df_layer_2$genelabels, as.character(df_layer_2$NA.), ""),
                  box.padding = unit(.42, "lines"), hjust = .9, vjust=.8, max.overlaps = 10000000)+
  theme(legend.title = element_blank(), legend.position = "none", text=element_text(size=20))+
  theme_bw()+
  geom_vline(xintercept=0, col="gray", linetype = "longdash")+
  geom_hline(yintercept=1.3, col="gray", linetype = "longdash")+
  ylim(0, 350)+
  xlim(-6, 7)+
  theme(legend.position = "none")+
  
  ggtitle("Trm Psoriasis Vs. Normal")






################eTreg1
Skin165 <- read.xlsx("/Users/chris/Desktop/DEG/AllPSORvsNorm2.xlsx", sheetIndex = 1)
Skin165 <- Skin165[, 1:6]

Skin165$threshold <- ""
Skin165$threshold <- ifelse(-log10(Skin165$p_val_adj)>0.05, TRUE, FALSE)

Skin165$genelabels <- ""
Skin165$genelabels <- ifelse(Skin165$NA.=="KLF2"
                             |Skin165$NA.=="ABLIM1"
                             |Skin165$NA.=="KLRB1"
                             |Skin165$NA.=="LTB"
                             |Skin165$NA.=="JUND"
                             |Skin165$NA.=="CD3G"
                             |Skin165$NA.=="CCR7"
                             #|Skin165$NA.=="GZMA"
                             |Skin165$NA.=="CXCL8"
                             |Skin165$NA.=="ZFP36"
                             |Skin165$NA.=="ZFP36L2"
                             |Skin165$NA.=="CCL22"
                             |Skin165$NA.=="BTG1"
                             |Skin165$NA.=="PRDM1"
                             |Skin165$NA.=="TSC22D3"
                             |Skin165$NA.=="TXNIP"
                             #|Skin165$NA.=="CNOT6L"
                             |Skin165$NA.=="SOCS1"
                             #|Skin165$NA.=="FOSL2"
                             #|Skin165$NA.=="IL10RA"
                             |Skin165$NA.=="NR3C1",
                             TRUE, FALSE)

Skin165$genecolour[Skin165$avg_log2FC > 0] <- '#F8766D'
Skin165$genecolour[Skin165$avg_log2FC < 0] <- '#619CFF'
Skin165$genecolour[Skin165$p_val_adj > 0.05] <- 'gray'

Skin165$order <- ""
Skin165$order <- ifelse(Skin165$NA.=="KLF2"
                        |Skin165$NA.=="ABLIM1"
                        |Skin165$NA.=="KLRB1"
                        |Skin165$NA.=="LTB"
                        |Skin165$NA.=="JUND"
                        |Skin165$NA.=="CD3G"
                        |Skin165$NA.=="CCR7"
                        #|Skin165$NA.=="GZMA"
                        |Skin165$NA.=="CXCL8"
                        |Skin165$NA.=="ZFP36"
                        |Skin165$NA.=="ZFP36L2"
                        |Skin165$NA.=="CCL22"
                        |Skin165$NA.=="BTG1"
                        |Skin165$NA.=="PRDM1"
                        |Skin165$NA.=="TSC22D3"
                        |Skin165$NA.=="TXNIP"
                        #|Skin165$NA.=="CNOT6L"
                        |Skin165$NA.=="SOCS1"
                        #|Skin165$NA.=="FOSL2"
                        #|Skin165$NA.=="IL10RA"
                        |Skin165$NA.=="NR3C1",
                        TRUE, FALSE)

Skin165$order <-as.numeric(factor(Skin165$order))

df_layer_1 <- Skin165[Skin165$order=="1",]
df_layer_2 <- Skin165[Skin165$order=="2",]

ggplot() +
  geom_point(data=df_layer_1, aes(avg_log2FC, -log10(p_val_adj), alpha = .9 ), colour=df_layer_1$genecolour) +
  geom_point(data=df_layer_2, aes(avg_log2FC, -log10(p_val_adj)), colour=df_layer_2$genecolour, size=4, show.legend = FALSE)+
  
  geom_text_repel(aes(df_layer_2$avg_log2FC, -log10(df_layer_2$p_val_adj)), label= ifelse(df_layer_2$genelabels, as.character(df_layer_2$NA.), ""),
                  box.padding = unit(.7, "lines"), hjust = .9, vjust=.3, max.overlaps = 10000000)+
  theme(legend.title = element_blank(), legend.position = "none", text=element_text(size=20))+
  theme_bw()+
  geom_vline(xintercept=0, col="gray", linetype = "longdash")+
  geom_hline(yintercept=1.3, col="gray", linetype = "longdash")+
  ylim(0, 350)+
  xlim(-5, 5)+
  theme(legend.position = "none")+
  
  ggtitle("Tcm Psoriasis Vs. Normal")






