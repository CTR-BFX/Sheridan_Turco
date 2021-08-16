#!/usr/local/bin/Rscript
# CTR_myt25_0006
# Human trophoblast models RNASeq data
# TSC-3D, TSC-2D, TSC-EVT, TSC-SCT, TO and TO-EVT
#
# Analysis Performed by Russell S. Hamilton & Xiaohui Zhao
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Russell S. Hamilton (rsh46@cam.ac.uk) & Xiaohui Zhao (xz289@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library("tidyr")
  library("dplyr")
  library("methods")
  library("utils")
  library("ggplot2")
  library("ggrepel")
  library("cowplot")
  library("Matrix")
  library("matrixStats")
  library("useful")
  library("reshape")
  library("reshape2")
  library("DESeq2")
  library("biomaRt")
  library("ggforce")
  library("pheatmap")
  library('RColorBrewer')
  library("scales")
  library("ggbeeswarm")
  library("BiocParallel")
  library("ggalt")
  library("ComplexHeatmap")
  library("biomaRt")
})

NUMCORES      <- 3
register(MulticoreParam(NUMCORES))

Project <- "CTR_myt25_0006"
baseDir <- "/storage/CTR-Projects/CTR_myt25/CTR_myt25_0006"
setwd(baseDir)
print(baseDir)

TOPNUM <- 2000
l2fc <- 1
significance <- 0.05
elementTextSize <- 6

message("+-------------------------------------------------------------------------------+")
message("+ Use ensEMBL Annotations+")
message("+-------------------------------------------------------------------------------+")
# attributePages(mart)
# if(interactive()){
#   mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
#  attributeSummary(mart)}
#BiocManager::install('grimbough/biomaRt')
ensembl    =  useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host = 'ensembl.org')
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description',
                                 'chromosome_name'), mart = ensembl, useCache = FALSE) 

head(ensEMBL2id)
nrow(ensEMBL2id)

message("+-------------------------------------------------------------------------------+")
message("+       RNA-Seq                                     ----------------------------+")
message("+-------------------------------------------------------------------------------+")

allLanes_counts           <- read.table(paste0(baseDir, "/reverse_stranded_Analysis/", 
                                               "featurecounts.merged.counts.tsv"), 
                                        header=TRUE, row.names=1)


allLanes_counts           <- allLanes_counts[ ,2:ncol(allLanes_counts)]
corner(allLanes_counts)
colnames(allLanes_counts)

meta.tbl.ori              <- read.csv("SampleTable_4Lanes_2Run.csv", header=T)

meta.tbl            <- as.data.frame(colnames(allLanes_counts), drop=F)
colnames(meta.tbl)  <- c("sampleID")
meta.tbl$sampleName <- meta.tbl$sampleID
meta.tbl$sampleName <- gsub(".L[0-9].*", "", meta.tbl$sampleName)
meta.tbl$origin     <- meta.tbl$sampleName
meta.tbl$origin     <- gsub("^[A-Z]*[0-9]*.", "", meta.tbl$origin)
meta.tbl$individual <- as.character(unlist(lapply(meta.tbl$sampleName, function(x) strsplit(x, split="_")[[1]][1])))
rownames(meta.tbl)  <- meta.tbl$sampleID
meta.tbl$lane       <- as.character(rep(c(1,2), length=272))
meta.tbl            <- data.frame(meta.tbl[ ,2:ncol(meta.tbl)])
meta.tbl$Group      <- ifelse(meta.tbl$origin=="EVT", "TSC-EVT", meta.tbl$origin)
meta.tbl$Group      <- ifelse(meta.tbl$origin=="Okae_TOM", "TSC-3D", meta.tbl$Group)
meta.tbl$Group      <- ifelse(meta.tbl$origin=="STB", "TSC-STC", meta.tbl$Group)
meta.tbl$Group      <- ifelse(meta.tbl$origin=="TSC", "TSC-2D", meta.tbl$Group)
meta.tbl$Group      <- ifelse(meta.tbl$origin=="org_EVT", "TO-EVT", meta.tbl$Group)
meta.tbl$Group      <- ifelse(meta.tbl$origin=="org_TOM", "TO", meta.tbl$Group)

head(meta.tbl)

table( paste0(meta.tbl$origin) )
table( paste0(meta.tbl$individual) )
table( paste0(meta.tbl$individual, "_", meta.tbl$origin) )

write.csv(meta.tbl, file = paste0(baseDir, "/reverse_stranded_Analysis/SampleTable_4Lanes_2Run.csv"))

message("+-------------------------------------------------------------------------------+")
message("+                           DESeq2 Processing                                   +")
message("+-------------------------------------------------------------------------------+")

dds <- DESeqDataSetFromMatrix(countData=allLanes_counts, colData=meta.tbl, design=~Group)
dds <- DESeq(dds, parallel=TRUE)

# colapseReplicates

vsd <- vst(dds,     blind=F)
colData(vsd)

message("+-------------------------------------------------------------------------------+")
message("+--------                   Custom PCA                                 ---------+")
message("+-------------------------------------------------------------------------------+")

pcaData    <- plotPCA(vsd, ntop=TOPNUM, intgroup=c("Group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(paste0(baseDir, "/reverse_stranded_Analysis/", Project, "_Fig.PCA.T2000.pdf"),width=10,height=10)
par(bg=NA)
ggplot(pcaData, aes(PC1, PC2, color=Group)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(Group), group=paste0(Group), 
                        label=Group), alpha=0.1, label.fontsize = 6, label.buffer = unit(1, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("PCA (Top ", TOPNUM, " Most Variable Genes)"),
       subtitle = "2 Sequencing Runs, 4 lanes each (8 points per sample)",
       caption = "CTR_myt25_0006: Margherita Turco" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='none', aspect.ratio=1) 
dev.off()


resultsNames(vsd)
## [1] "Intercept"                     "origin_3D.fluffy_vs_2D.fluffy"
## [3] "origin_EVT_vs_2D.fluffy"       "origin_Okae.TOM_vs_2D.fluffy" 
## [5] "origin_org.EVT_vs_2D.fluffy"   "origin_org.TOM_vs_2D.fluffy"  
## [7] "origin_STB_vs_2D.fluffy"       "origin_TSC_vs_2D.fluffy" 

message("+-------------------------------------------------------------------------------+")
message("+--- Remove fluffy samples, 48                                         ---------+")
message("+-------------------------------------------------------------------------------+")

meta.tbl.rm <- meta.tbl[-grep("*.fluffy", meta.tbl$origin), ]

write.csv(meta.tbl.rm, file = paste0(baseDir, "/reverse_stranded_Analysis/SampleTable_4Lanes_2Run_rmfluffy.csv"), row.names=F)

allLanes_counts.rm <- allLanes_counts[,-grep("*.fluffy", colnames(allLanes_counts))]
dds.rm <- DESeqDataSetFromMatrix(countData=allLanes_counts.rm, colData=meta.tbl.rm, design=~Group)
dds.rm <- DESeq(dds.rm, parallel=TRUE)

# colapseReplicates

vsd.rm <- vst(dds.rm,     blind=F)
colData(vsd.rm)

save(dds, dds.rm, dds.rm.collapse, file = paste0(baseDir, "/reverse_stranded_Analysis/",Project, "-OriginalData_DESeq_dds.RData"))

message("+-------------------------------------------------------------------------------+")
message("+--------                   Custom PCA  for Rm fluffy                  ---------+")
message("+-------------------------------------------------------------------------------+")

pcaData.rm    <- plotPCA(vsd.rm, ntop=TOPNUM, intgroup=c("Group"), returnData=TRUE)
percentVar.rm <- round(100 * attr(pcaData.rm, "percentVar"))

pdf(paste0(baseDir, "/reverse_stranded_Analysis/", Project, "_Fig.PCA.rmFluffy.T2000.pdf"),width=10,height=10)
par(bg=NA)
ggplot(pcaData.rm, aes(PC1, PC2, color=Group)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(Group), group=paste0(Group), 
                        label=Group), alpha=0.1, label.fontsize = 6, label.buffer = unit(1, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  xlab(paste0("PC1: ",percentVar.rm[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.rm[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("PCA (Top ", TOPNUM, " Most Variable Genes)"),
       subtitle = "2 Sequencing Runs, 4 lanes each (8 points per sample)",
       caption = "CTR_myt25_0006: Margherita Turco" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='none', aspect.ratio=1) 
dev.off()



pdf(paste0(baseDir, "/reverse_stranded_Analysis/", Project, "_Fig.PCA.rmFluffy.T2000.rmLab.pdf"),width=10,height=10)
par(bg=NA)
ggplot(pcaData.rm, aes(PC1, PC2, color=Group)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(Group), group=paste0(Group)), alpha=0.1) +
  geom_point(size=2, alpha=0.75) +
  xlab(paste0("PC1: ",percentVar.rm[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.rm[2],"% variance")) + 
  coord_fixed() +
  theme_bw()  +
  labs(title = paste0("PCA (Top ", TOPNUM, " Most Variable Genes)"),
       subtitle = "2 Sequencing Runs, 4 lanes each (8 points per sample)")+
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='none', aspect.ratio=1,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) 
dev.off()


resultsNames(vsd.rm)
## [1] "Intercept"                     "origin_3D.fluffy_vs_2D.fluffy"
## [3] "origin_EVT_vs_2D.fluffy"       "origin_Okae.TOM_vs_2D.fluffy" 
## [5] "origin_org.EVT_vs_2D.fluffy"   "origin_org.TOM_vs_2D.fluffy"  
## [7] "origin_STB_vs_2D.fluffy"       "origin_TSC_vs_2D.fluffy" 

dds.rm.collapse <- collapseReplicates(dds.rm, dds.rm$sampleName, renameCols = TRUE) 
meta.tbl.rm.collapse <- colData(dds.rm.collapse)
vsd.rm.collapse <- vst(dds.rm.collapse,     blind=F)
colData(vsd.rm.collapse)

pcaData.rm.collapse    <- plotPCA(vsd.rm.collapse, ntop=TOPNUM, intgroup=c("Group"), returnData=TRUE)
percentVar.rm.collapse <- round(100 * attr(pcaData.rm.collapse, "percentVar"))


pcaData.rm.collapse$newname <- unlist(lapply(as.character(pcaData.rm.collapse$name), function(x) strsplit(x, split="[_]")[[1]][1]))
pcaData.rm.collapse$origin  <- meta.tbl.rm.collapse$origin

pdf(paste0(baseDir, "/reverse_stranded_Analysis/", Project, "-Fig2A.PCA.rmFluffy.T2000.collapsed.labels.pdf"),width=10,height=10)
par(bg=NA)
ggplot(pcaData.rm.collapse, aes(PC1, PC2, color=Group, label=newname)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(Group), group=paste0(Group), 
                        label=Group), alpha=0.1, label.fontsize = 6, label.buffer = unit(1, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  geom_text_repel(data = pcaData.rm.collapse, aes(label = newname)) +
  xlab(paste0("PC1: ",percentVar.rm[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.rm[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("PCA (Top ", TOPNUM, " Most Variable Genes)"),
       subtitle = "2 Sequencing Runs, 4 lanes each (collapsed)",
       caption = "CTR_myt25_0006: Margherita Turco" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='none', aspect.ratio=1) 
dev.off()

message("+-------------------------------------------------------------------------------+")
message("+--------                  Fig.2A PCA plot                             ---------+")
message("+-------------------------------------------------------------------------------+")


pdf(paste0(baseDir, "/reverse_stranded_Analysis/Fig.2A_PCA.pdf"),width=10,height=10)
par(bg=NA)
ggplot(pcaData.rm.collapse, aes(PC1, PC2, color=origin)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(origin), group=paste0(origin)), alpha=0.1) +
  geom_point(size=2, alpha=0.75) +
  xlab(paste0("PC1: ",percentVar.rm[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.rm[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='none', aspect.ratio=1) 
dev.off()


message("+-------------------------------------------------------------------------------+")
message("+---       Pairwise DESeq Analysis, revised 05/05/21, add TSC2DvsTroph.org -----+")
message("+-------------------------------------------------------------------------------+")

control.names  <- c("*_TSC_*", "*_org_TOM_*", "*_org_TOM_*", "TSC", "TSC", "*_org_EVT_*", "*_org_TOM_*", "*_org_TOM_*")
compare.names  <- c("*_Okae_TOM.*", "*_Okae_TOM_*", "*_org_EVT_*", "EVT", "STB", "EVT", "STB", "*_TSC_*")

control.Tnames <- c("TSC", "org_TOM","org_TOM", "TSC", "TSC", "org_EVT", "org_TOM", "org_TOM")
compare.Tnames <- c("Okae_TOM", "Okae_TOM","org_EVT", "EVT", "STB", "EVT", "STB", "TSC")

control.TnamesN <- c("TSC-2D", "TO","TO", "TSC-2D", "TSC-2D", "TO-EVT", "TO", "TO")
compare.TnamesN <- c("TSC-3D", "TSC-3D","TO-EVT", "TSC-EVT", "TSC-SCT", "TSC-EVT", "TSC-SCT", "TSC-2D")


control.colors <- c("darkgoldenrod4", "cyan", "cyan", "darkgoldenrod4", "darkgoldenrod4", "chartreuse3", "cyan", "cyan")
compare.colors <- c("magenta", "magenta", "chartreuse3", "firebrick3", "deepskyblue1","firebrick3","deepskyblue1","darkgoldenrod4")

models  <- paste0(compare.Tnames, "vs", control.Tnames)
modelsN <- paste0(compare.TnamesN, "vs", control.TnamesN)


message("+---       DESeq analysis and data outputs                                  ---------+")


for(i in 1:length(models)){ 
  ## subset the counts matrix
  colIndex <- c(grep(paste0("",control.names[i],""), colnames(allLanes_counts)),
                grep(paste0("",compare.names[i],""), colnames(allLanes_counts)))
  if(i==4){
    colIndex1<- grep(paste0("",compare.names[3],""), colnames(allLanes_counts))
    colIndex <- colIndex[colIndex%in%colIndex1==F]
  }
  if(i==6){
    colIndex <- grep(paste0("",compare.names[4],""), colnames(allLanes_counts))
  }
  geneCounts <- allLanes_counts[, colIndex]
  geneCounts <- geneCounts[, order(colnames(geneCounts))]
  
  ## subset the meta table
  meta.tbl.sub <- subset(meta.tbl, origin==control.Tnames[i] | origin==compare.Tnames[i])
  meta.tbl.sub <- meta.tbl.sub[order(rownames(meta.tbl.sub)),]
  meta.tbl.sub$origin <- as.factor(meta.tbl.sub$origin)
  meta.tbl.sub$origin <- relevel(meta.tbl.sub$origin, control.Tnames[i])
  
  ## DESeq analysis
  dds.ori <- DESeqDataSetFromMatrix(countData=geneCounts, colData=meta.tbl.sub, design=~origin)
  
  dds.ori <- DESeq(dds.ori, parallel=TRUE)
  vsd.ori <- vst(dds.ori,     blind=F)
  colData(vsd.ori)
  dds.collapse     <- collapseReplicates(dds.ori, dds.ori$sampleName, renameCols = TRUE)
  dds.collapse     <- DESeq(dds.collapse,  parallel=TRUE)
  resultsNames(dds.collapse)
  design(dds.collapse)
  colData(dds.collapse)
  vsd.collapse       <- vst(dds.collapse,  blind=F) 
  
  ## DEGs list 
  
  res     <- lfcShrink(dds=dds.collapse, coef=2, type="apeglm", parallel=TRUE)
  res.sig <- as.data.frame(subset(res, abs(log2FoldChange) >= l2fc & padj < significance))
  res.sig$external_gene_name <- rownames(res.sig)
  res.sig.ann <-  merge(res.sig, ensEMBL2id, by = "external_gene_name")
  numsig      <- dim(res.sig)[1]
  numsigE     <- length(unique(res.sig.ann$external_gene_name))
  write.csv(res.sig, file = paste0(baseDir, "/reverse_stranded_Analysis/",Project, "-", models[i], "_l2fc1_p0.05_N", numsig, ".csv"),
            row.names=F, quote=T)
  
  ## save RData for recalling
  save(dds.ori, vsd.ori, dds.collapse, vsd.collapse, res, file = paste0(baseDir, "/reverse_stranded_Analysis/",Project, "-", models[i], "_DESeq_Object.RData"))
  
}

message("+---       PCA/Volcano/Heatmap plot customised function                          --------+")
customPCA_fn_group   <- function(vsd.ori, vsd.collapse, dds.ori, dds.collapse, TOPNUM, ensEMBL2id, control, 
                                 controlN, treatN,ctrlColor,trtColor){
  ## origin vst
  
  rv.ori     <- rowVars(assay(vsd.ori))
  select.ori <- order(rv.ori, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv.ori)))]
  pca.ori    <- prcomp(t(assay(vsd.ori)[select.ori, ]))
  
  pc1var.ori <- round(summary(pca.ori)$importance[2,1]*100, digits=2)
  pc2var.ori <- round(summary(pca.ori)$importance[2,2]*100, digits=2)
  pc1lab.ori <- paste0("PC1 (",as.character(pc1var.ori),"%)")
  pc2lab.ori <- paste0("PC2 (",as.character(pc2var.ori),"%)")
  
  sampleTBL.ori <- as.data.frame(colData(dds.ori))
  Group         <- ifelse(sampleTBL.ori$origin==control, controlN, treatN)
  scores.ori    <- data.frame(sampleName=rownames(sampleTBL.ori), pca.ori$x, 
                              individual=sampleTBL.ori$individual, 
                              origin=sampleTBL.ori$origin,
                              lane=sampleTBL.ori$lane,
                              Group=Group)
  
  scores.ori$origin  <- relevel(scores.ori$origin, paste0("", control,""))
 
  ## collapse vst
  rv.collapse     <- rowVars(assay(vsd.collapse))
  select.collapse <- order(rv.collapse, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv.collapse)))]
  pca.collapse    <- prcomp(t(assay(vsd.collapse)[select.collapse, ]))
  
  pc1var.collapse <- round(summary(pca.collapse)$importance[2,1]*100, digits=2)
  pc2var.collapse <- round(summary(pca.collapse)$importance[2,2]*100, digits=2)
  pc1lab.collapse <- paste0("PC1 (",as.character(pc1var.collapse),"%)")
  pc2lab.collapse <- paste0("PC2 (",as.character(pc2var.collapse),"%)")
  
  sampleTBL.collapse <- as.data.frame(colData(dds.collapse))
  Group              <- ifelse(sampleTBL.collapse$origin==control, controlN, treatN)
  scores.collapse    <- data.frame(sampleName=rownames(sampleTBL.collapse), pca.collapse$x, 
                                   individual=sampleTBL.collapse$individual, 
                                   origin=as.factor(sampleTBL.collapse$origin),
                                   Group=Group)
  scores.collapse$origin <- relevel(scores.collapse$origin, paste0("", control,""))

  

 
  plt.ori <- ggplot(scores.ori, aes(PC1, PC2, color=Group, fill=Group, label=individual, shape=lane)) +
    geom_mark_ellipse(aes(fill = NULL, color=paste0(Group), group=paste0(Group), label=Group),alpha=0.1,
                      label.fontsize=10,label.buffer = unit(1, 'mm')) +
    geom_point(size=2, alpha=0.75) +
    scale_shape_manual(name="Lane", values = c(1, 2)) + 
    scale_color_manual(values=c(ctrlColor, trtColor)) + 
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    theme_bw() +
    labs(title = paste0("Origin PCA (Top ", TOPNUM, " Most Variable Genes)"),
         subtitle = "2 Sequencing Runs, 4 lanes each") +
    theme(text = element_text(size=12), 
          plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
          aspect.ratio=1, legend.position="none", panel.border = element_blank(),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
           
  
  
  
  plt.collapse <- ggplot(scores.collapse, aes(PC1, PC2, color=Group, fill=Group, label=individual)) +
    geom_mark_ellipse(aes(fill = NULL, color=paste0(Group), group=paste0(Group), label=Group),alpha=0.1,
                      label.fontsize=10,label.buffer = unit(1, 'mm')) +
    geom_point(size=2, alpha=0.75) +
    scale_color_manual(values=c(ctrlColor, trtColor)) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    theme_bw() +
    labs(title = paste0("Collapse PCA (Top ", TOPNUM, " Most Variable Genes)"),
         subtitle = "2 Sequencing Runs, 4 lanes each") +
    theme(text = element_text(size=12), 
          plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
          legend.title = element_blank(),
          legend.key = element_rect(fill = "white",colour = "white"),
          legend.position='none', 
          aspect.ratio=1,
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  plt.collapse.rmL <- ggplot(scores.collapse, aes(PC1, PC2, color=Group, fill=Group)) +
    geom_mark_ellipse(aes(fill = NULL, color=paste0(Group), group=paste0(Group), label=Group),alpha=0.1,
                      label.fontsize=10,label.buffer = unit(1, 'mm')) +
    geom_point(size=2, alpha=0.75) +
    scale_color_manual(values=c(ctrlColor, trtColor)) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    labs(title = paste0("Collapse PCA (Top ", TOPNUM, " Most Variable Genes)")) +
    coord_fixed() +
    theme_bw() +
    theme(text = element_text(size=12), 
          plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
          legend.title = element_blank(),
          legend.key = element_rect(fill = "white",colour = "white"),
          legend.position='right', 
          aspect.ratio=1,
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  
  loadings                 <- as.data.frame(pca.collapse$rotation)
  loadings$external_gene_name <- rownames(loadings)
  loadings                 <- merge(loadings, ensEMBL2id, by="external_gene_name")
  
  pca.1         <- loadings[ order(loadings$PC1,decreasing=TRUE), ]
  pca.1.25      <- pca.1[c(1:25),]
  pca.1.25.plot <- ggplot(data=pca.1.25, aes(x=factor(external_gene_name,levels=unique(external_gene_name)), y=PC1)) + 
    geom_point(size = 3 ) + xlab("") + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  
  pca.2         <- loadings[ order(loadings$PC2,decreasing=TRUE), ]
  pca.2.25      <- pca.2[c(1:25),]
  pca.2.25.plot <- ggplot(data=pca.2.25, aes(x=factor(external_gene_name,levels=unique(external_gene_name)), y=PC2)) + 
    geom_point(size = 3 ) + xlab("") + 
    theme()+ 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  
  
  return(list(plt.ori, plt.collapse, plt.collapse.rmL, pca.1.25.plot, pca.2.25.plot) )
  
}

functionPlotDEVolcano <- function(results, sig_cut, logfc_cut, title,  xrange, yrange, topN, xlabel, ylabel) {
  
  results       <- as.data.frame(results)
  results       <- results[order(-results$log2FoldChange),]
  results$genes <- rownames(results)
  results.up    <- subset(results, log2FoldChange > logfc_cut & padj<= sig_cut)
  results.down  <- subset(results, log2FoldChange < logfc_cut & padj<= sig_cut)
  
  volc.plt <- ggplot(data=results, aes(x=log2FoldChange, y=-log10(padj), label=genes)) +
    geom_vline(xintercept = logfc_cut,     colour="black", linetype = "dashed", alpha=0.5) +
    geom_vline(xintercept = -(logfc_cut),  colour="black", linetype = "dashed", alpha=0.5) +
    geom_hline(yintercept = -log10(sig_cut), colour="black", linetype = "dashed", alpha=0.5) +
    
    geom_point(data=subset(results, abs(log2FoldChange) < logfc_cut | padj > sig_cut), alpha=0.75, size=0.3, colour="grey") +
    geom_point(data=subset(results, padj<=sig_cut & log2FoldChange >= logfc_cut),      alpha=0.75, size=0.8, colour="red") +
    geom_point(data=subset(results, padj<=sig_cut & log2FoldChange <= -(logfc_cut)),   alpha=0.75, size=0.8, colour="blue") +
    
    geom_text_repel( data= results.up[order(results.up$padj, -results.up$log2FoldChange),][1:topN,],
                     show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=3 ) +
    geom_text_repel( data= results.down[order(results.down$padj, results.down$log2FoldChange),][1:topN,],
                     show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=3 ) +
    xlab(xlabel) + ylab(ylabel) +
    scale_x_continuous(limits=c(xrange[1],xrange[2]), breaks=seq(xrange[1],xrange[2],xrange[3])) +
    scale_y_continuous(limits=c(yrange[1],yrange[2]), breaks=seq(yrange[1],yrange[2],yrange[3])) +
    theme(aspect.ratio=1) +
    ggtitle(title) +
    theme_update(plot.title = element_text(size=16, face="bold", hjust=0.5),
                 axis.title.x = element_text(size=12, face= "bold"),
                 axis.text.x = element_text(size=12, face="bold"),
                 axis.title.y.left = element_text(size=12, face= "bold"),
                 axis.text.y = element_text(size=12, face="bold")) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "white", colour = NA)) 
  
  
  return(volc.plt)
  
}
message("+---       PCA/Volcano/Heatmap plot                                              --------+")

for(i in 1:8){
  load(paste0(baseDir, "/reverse_stranded_Analysis/",Project, "-", models[i], "_DESeq_Object.RData"))
  
  ## PCA for origin and replicates collpase, PC1+PC2  
  
  pcaData.ori    <- plotPCA(vsd.ori, ntop=TOPNUM, intgroup=c("origin"), returnData=TRUE)
  percentVar     <- round(100 * attr(pcaData.ori, "percentVar"))
  
  pcaData.collapse    <- plotPCA(vsd.collapse, ntop=TOPNUM, intgroup=c("origin"), returnData=TRUE)
  percentVar.collapse <- round(100 * attr(pcaData.collapse, "percentVar"))
  
  pcas <- customPCA_fn(vsd.ori, vsd.collapse, dds.ori, dds.collapse, TOPNUM, models[i], 
                       ensEMBL2id, control.Tnames[i])
  
  pdf(paste0(baseDir, "/reverse_stranded_Analysis/", Project, "_", models[i], "_ori_collapse_Fig.PCA.pdf"),width=10,height=10)
  par(bg=NA)
  plot_grid(pcas[[1]], pcas[[2]], pcas[[3]], pcas[[4]], nrow=2, align="hv")
  dev.off()
  
  ## Volcano Plot
  
  title <- paste0(compare.Tnames[i], " vs ", control.Tnames[i])
  xrange <- c(-15,15,3)
  if(i!=3){
    yrange<- c(0, 70, 10)
  }else{
    yrange <-c(0,160,40)
  }
  topN <- 20
  
  xlabel <- paste0("log2FC", "(", models[i],")")
  ylabel <- bquote("-log"[10]~"(adj.p.value)")
  
  plt_volt <- functionPlotDEVolcano(res, significance, l2fc, title,  xrange, yrange, 20, xlabel, ylabel)
  
  pdf(paste0(baseDir, "/reverse_stranded_Analysis/", Project, "-", models[i], "_Fig.Volcano.pdf"))
  print(plt_volt)
  dev.off()
  
}

message("+--- PCA plot for TSC-3D vs TO,PCA for origin and replicates collpase, PC1+PC2  --------+")

i=2
load(paste0(baseDir, "/reverse_stranded_Analysis/",Project, "-", models[i], "_DESeq_Object.RData"))

pcas <- customPCA_fn_group(vsd.ori, vsd.collapse, dds.ori, dds.collapse, TOPNUM, 
                     ensEMBL2id, control.Tnames[i], control.TnamesN[i], compare.TnamesN[i], 
                     control.colors[i], compare.colors[i])


pdf(paste0(baseDir, "/reverse_stranded_Analysis/SFig3B_", modelsN[i], "_ori_collapse_Fig.PCA.pdf"),width=12,height=10)
par(bg=NA)
plot_grid(pcas[[1]], pcas[[3]], pcas[[2]], pcas[[4]], nrow=2, align="hv")
dev.off()

pdf(paste0(baseDir, "/reverse_stranded_Analysis/SFig3B_", modelsN[i], "_collapse_Fig.PCA.pdf"),width=6,height=10)
par(bg=NA)
plot_grid(pcas[[3]], pcas[[4]], nrow=2)
dev.off()



message("+---------------------------------------------------------------------------------+")
message("+-------  Identify specific markers for each model against the others          ---+")
message("+---------------------------------------------------------------------------------+")

metatbl.rm.new     <- meta.tbl.rm
cellNames          <- unique(meta.tbl.rm$origin)
newmodels          <- paste0(cellNames, "vsOthers")

for(i in 1:length(cellNames)){
  metatbl.rm.new$neworigin <- ifelse(metatbl.rm.new$origin==cellNames[i], cellNames[i], "Control")
  dds.rm.new               <- DESeqDataSetFromMatrix(countData=allLanes_counts.rm, colData=metatbl.rm.new, design=~neworigin)
  dds.rm.new               <- DESeq(dds.rm.new, parallel=TRUE)
  vsd.rm.new               <- vst(dds.rm.new,     blind=F)
  colData(vsd.rm.new)
  # colapseReplicates
  dds.collapse.rm.new      <- collapseReplicates(dds.rm.new, dds.rm.new$sampleName, renameCols = TRUE)
  dds.collapse.rm.new      <- DESeq(dds.collapse.rm.new,  parallel=TRUE)
  resultsNames(dds.collapse.rm.new)
  design(dds.collapse.rm.new)
  colData(dds.collapse.rm.new)
  vsd.collapse.rm.new      <- vst(dds.collapse.rm.new,  blind=F) 
  
  ## DEGs list 
  
  res <- lfcShrink(dds=dds.collapse.rm.new, coef=2, type="apeglm", parallel=TRUE)
  res.sig <- as.data.frame(subset(res, abs(log2FoldChange) >= l2fc & padj < significance))
  res.sig$external_gene_name <- rownames(res.sig)
  res.sig.ann <-  merge(res.sig, ensEMBL2id, by = "external_gene_name")
  numsig <- dim(res.sig)[1]
  numsigE <- length(unique(res.sig.ann$external_gene_name))
  print(numsigE)
  write.csv(res.sig, file = paste0(baseDir, "/reverse_stranded_Analysis/",Project, "-", newmodels[i], "_l2fc1_p0.05_N", numsig, ".csv"),
            row.names=F, quote=T)
  
  ## save RData for recalling
  save(dds.rm.new, vsd.rm.new, dds.collapse.rm.new, vsd.collapse.rm.new, res, 
       file = paste0(baseDir, "/reverse_stranded_Analysis/",Project, "-", newmodels[i], "_DESeq_Object.RData"))
  
}

message("+---- Generate a heatmap with selected markers and also write csv file for all of the counts ------+")

MononuclearT  <- c("LRP5", "GATA2", "GATA3", "KRT7", "TFAP2A", "TFAP2C", "EPCAM", "TEAD4", "TP63",
                   "ITGA6", "CDH1", "EGFR", "GCM1", "HLA-A", "HLA-B")
SyncytioT     <- c("CGA", "CGB", "PSG1", "CSH1", "SDC1", "ERVW-1", "ERVFRD-1", "CD46", "TFRC",
                   "PGF", "PAPPA2", "GDF15")
ColumnNiche   <- c("ITGA2", "PECAM1", "RARRES3", "SOX13", "EMP3", "ITGB6", "NOTCH1")
ExtravillousT <- c("ITGA5", "HLA-G", "MMP2", "CD9", "ITGA1", "LAIR2", "PLAC8", "MCAM", "NCAM1")

library(ComplexHeatmap)
load(paste0(baseDir, "/reverse_stranded_Analysis/",Project, "-OriginalData_DESeq_dds.RData"))
rld.collapse.rm <- rlogTransformation(dds.rm.collapse, blind=F)
vsd.collapse.rm <- vst(dds.rm.collapse, blind=F)

DEGs.rldmat <- as.data.frame(assay(rld.collapse.rm))
DEGs.vsdmat <- as.data.frame(assay(vsd.collapse.rm))
DEGs.rawdat <- counts(dds.rm.collapse, normalized = FALSE)
write.csv(DEGs.rldmat, file = paste0(baseDir, "/reverse_stranded_Analysis/",Project, "-RLDcounts_N33121_rmFluffy.csv"), row.names=T, quote=T)
write.csv(DEGs.rldmat, file = paste0(baseDir, "/reverse_stranded_Analysis/",Project, "-Rawcounts_N33121_rmFluffy.csv"), row.names=T, quote=T)

## prepare the RNASeq heatmap matrix
Allmarkers <- c(MononuclearT, SyncytioT, ColumnNiche, ExtravillousT )
DEGs.rldmat.sel <- DEGs.rldmat[rownames(DEGs.rldmat)%in%Allmarkers==T, ]

## PCAM1
TSC.mean         <- rowMeans(DEGs.rldmat.sel[,grep("*_TSC", colnames(DEGs.rldmat.sel))])
TSC3D.mean       <- rowMeans(DEGs.rldmat.sel[,grep("*_Okae_TOM", colnames(DEGs.rldmat.sel))])
TSCEVT.mean      <- rowMeans(DEGs.rldmat.sel[,grep("*_EVT", colnames(DEGs.rldmat.sel))[1:5]])
TSCSTB.mean      <- rowMeans(DEGs.rldmat.sel[,grep("*_STB", colnames(DEGs.rldmat.sel))])
TrophOrg.mean    <- rowMeans(DEGs.rldmat.sel[,grep("*_org_TOM", colnames(DEGs.rldmat.sel))])
TrophOrgEVT.mean <- rowMeans(DEGs.rldmat.sel[,grep("*_org_EVT", colnames(DEGs.rldmat.sel))])

mplt.mat <- cbind(TSC.mean, TSC3D.mean, TSCEVT.mean, TSCSTB.mean, TrophOrg.mean, TrophOrgEVT.mean)

breaksList.deg = seq(0, 21, by = 1)
ScaleCols.deg <- colorRampPalette(colors = c("purple4","white","darkgreen"))(length(breaksList.deg))

cnames <- Allmarkers[Allmarkers %in% rownames(mplt.mat)]
mplt.mat.ord <- mplt.mat[cnames,]
mplt.mat.length <- unlist(lapply(list(MononuclearT, SyncytioT, ColumnNiche, ExtravillousT), function(x) length(x)))
colnames(mplt.mat.ord) <- c("TSC-2D", "TSC-3D", "TSC-EVT", "TSC-STB", "Troph.org", "Troph.org-EVT")

## Fig2B remove HLA-A and HLA-B
mplt.mat.sub <- mplt.mat.ord[-c(14:15),]
mplt.mat.sublen <- c(13,12,7,9)
pht_list.sub = Heatmap(mplt.mat.sub, 
                       col = ScaleCols.deg, 
                       name = "DEGs", 
                       show_row_names=T,
                       show_column_names = T, 
                       width = unit(6, "cm"),
                       heatmap_legend_param = list(title = "log2 Normalised \nExpression",
                                                   title_position = "topcenter",
                                                   title_gp = gpar(fontsize = 8),
                                                   labels_gp = gpar(fontsize = 6),
                                                   direction = "horizontal",
                                                   legend_width = unit(4, "cm")),
                       cluster_rows = T,
                       show_row_dend = T,
                       column_title="",
                       row_split=factor(rep(c("Mononuclear \nTrophoblast", "Syncytio\nTrophoblast", "Column \nNiche", "Extravillous \nTrophoblast"), 
                                            mplt.mat.sublen),level=c("Mononuclear \nTrophoblast", "Syncytio\nTrophoblast", "Column \nNiche", "Extravillous \nTrophoblast")),
                       row_order = c(1:dim(mplt.mat.sub)[1]),
                       column_order=factor(c("TSC-2D", "TSC-3D", "TSC-EVT", "TSC-STB", "Troph.org", "Troph.org-EVT"),
                                           levels=c("TSC-2D", "TSC-3D", "TSC-EVT", "TSC-STB", "Troph.org", "Troph.org-EVT")),
                       column_names_side = "top", 
                       cluster_row_slices = T,
                       row_title_rot = 90,
                       row_names_gp = gpar(fontsize = 8),
                       row_title_gp = gpar(col = "black", fontsize = 8),
                       column_names_gp = gpar(fontsize = 8),
                       cluster_columns = T, 
                       column_dend_height = unit(1, "cm")
                       
                       
)


pdf(paste0(baseDir, "/reverse_stranded_Analysis/Fig2B_SelMarkers_Heatmap.pdf"))
draw(pht_list.sub, heatmap_legend_side = "bottom")
dev.off()



message("+-------   Heatmap plot for each cell type selected top up-regulated 20 genes  ------------------+")

selNum            <- 25
selHeat.sigG      <- list()

for(i in 1:length(cellNames)){
  #load(file = paste0(baseDir, "/reverse_stranded_Analysis/",Project, "-", newmodels[i], "_DESeq_Object.RData"))
  dat              <- read.csv(file = paste0(baseDir, "/reverse_stranded_Analysis/",Project, "-", newmodels[i], "_l2fc1_p0.05_N", numsigs[i], ".csv"))
  dat              <- dat[order(-dat$log2FoldChange, -dat$padj), ]
  datsub           <- dat[1:selNum,]
  selHeat.sigG[[i]]<- datsub$external_gene_name
  selHeat.sigG
}


load(paste0(baseDir, "/reverse_stranded_Analysis/",Project, "-OriginalData_DESeq_dds.RData"))
rld.collapse.rm <- rlogTransformation(dds.rm.collapse, blind=F)
vsd.collapse.rm <- vst(dds.rm.collapse, blind=F)

DEGs.rldmat <- as.data.frame(assay(rld.collapse.rm))
DEGs.vsdmat <- as.data.frame(assay(vsd.collapse.rm))
groupNames  <- c("TSC-EVT", "TSC-3D", "TSC-SCT", "TSC-2D", "TSC-EVT",
                 "TSC-3D", "TSC-SCT", "TSC-2D", "TSC-EVT", "TSC-3D",
                 "TSC-SCT", "TSC-2D", "TSC-EVT", "TSC-3D", "TSC-SCT",
                 "TSC-2D", "TSC-EVT", "TSC-3D", "TSC-SCT", "TSC-2D",
                 "TO-EVT", "TO", "TO-EVT", "TO", "TO-EVT", "TO", "TO-EVT", "TO")
sampleSN    <- unlist(lapply(colnames(DEGs.vsdmat), function(x) strsplit(x, split="_")[[1]][1]))      

selHeatmat  <- list(); selHeatmatMean  <- list()
for(i in 1:6){
  selHeatmat[[i]]  <- DEGs.vsdmat[rownames(DEGs.vsdmat)%in%selHeat.sigG[[i]], ]
  colnames(selHeatmat[[i]]) <- paste0(groupNames, "_", sampleSN)
  TSC_EVT <- rowMeans(selHeatmat[[i]][,grep("TSC-EVT_", colnames(selHeatmat[[i]]))]) 
  TO_EVT  <- rowMeans(selHeatmat[[i]][,grep("TO-EVT_", colnames(selHeatmat[[i]]))]) 
  TSC_SCT <- rowMeans(selHeatmat[[i]][,grep("TSC-SCT_", colnames(selHeatmat[[i]]))]) 
  TSC_2D  <- rowMeans(selHeatmat[[i]][,grep("TSC-2D_", colnames(selHeatmat[[i]]))]) 
  TSC_3D  <- rowMeans(selHeatmat[[i]][,grep("TSC-3D_", colnames(selHeatmat[[i]]))]) 
  TO      <- rowMeans(selHeatmat[[i]][,grep("TO_", colnames(selHeatmat[[i]]))]) 
  selHeatmatMean[[i]] <- cbind(TSC_EVT, TO_EVT, TSC_SCT, TSC_2D, TSC_3D, TO) 
  colnames(selHeatmatMean[[i]]) <- c("TSC-EVT", "TO-EVT", "TSC-SCT", "TSC-2D", "TSC-3D", "TO")
}

breaksList.deg = seq(4, 18, by = 1)
ScaleCols.deg <- colorRampPalette(colors = c("purple4","white","darkgreen"))(length(breaksList.deg))

pullheatmat <- rbind(selHeatmat[[6]], selHeatmat[[3]], selHeatmat[[2]], selHeatmat[[4]], selHeatmat[[1]], selHeatmat[[5]])
column.ord <- c(grep("TO_", colnames(pullheatmat)), grep("TSC-SCT", colnames(pullheatmat)), grep("TSC-3D", colnames(pullheatmat)),
                grep("TSC-2D", colnames(pullheatmat)), grep("TSC-EVT", colnames(pullheatmat)), grep("TO-EVT", colnames(pullheatmat)))
pullheatmat <- pullheatmat[,column.ord]

ham.names  <-  unlist(lapply(colnames(pullheatmat), function(x) strsplit(x, "_")[[1]][1]))
ham        <- HeatmapAnnotation(Group = ham.names, 
                col = list(Group = c("TO" = "cyan", "TSC-SCT" = "deepskyblue1", "TSC-3D" = "magenta",                                                         
                 "TSC-2D" = "darkgoldenrod4", "TSC-EVT"= "firebrick3", "TO-EVT"= "chartreuse3")), 
                annotation_name_side = "left")

allpht_list = Heatmap(as.matrix(pullheatmat), 
                      col = ScaleCols.deg, 
                      name = "Norm.Expr", 
                      show_row_names=T,
                      show_column_names = T, 
                      width = unit(10, "cm"),
                      top_annotation=ham,
                      cluster_rows = T,
                      show_row_dend = F,
                      column_title=NULL,
                      row_split=factor(rep(c("TO", "TSC-SCT", "TSC-3D", "TSC-2D", "TSC-3EVT", "TO-EVT"), each=25),
                                       levels=c("TO", "TSC-SCT", "TSC-3D", "TSC-2D", "TSC-3EVT", "TO-EVT")),
                      row_order = c(1:dim(pullheatmat)[1]),
                      column_names_side = "bottom", 
                      cluster_row_slices = T,
                      row_title_rot = 90,
                      row_names_gp = gpar(fontsize = 6),
                      row_title_gp = gpar(col = "black", fontsize = 8),
                      column_names_gp = gpar(fontsize =6, col= "black"),
                      cluster_columns = T,
                      column_split = rep(c("TO", "TSC-SCT", "TSC-3D", "TSC-2D", "TSC-3EVT", "TO-EVT"), c(4,5,5,5,5,4)),
                      #column_title_gp = gpar(fontsize =5.5, fill= rainbow(6)),
                      cluster_column_slices=T
)
pdf(paste0(baseDir, "/reverse_stranded_Analysis/SFig3A_Heatmap_sel25_upreg_new.pdf"), width=8, height=12)
draw(allpht_list, heatmap_legend_side = "right", merge_legend=T)
dev.off()

message("+-------------------------------------------------------------------------------+")
message("+---                         GeneOntology Analysis                      --------+")
message("+-------------------------------------------------------------------------------+")

suppressPackageStartupMessages({
  library("clusterProfiler")
  library("DOSE")
  library("GSEABase")
  library("AnnotationHub")
  library("org.Hs.eg.db")
  library("BiocParallel")
  library("gage")
  library("gageData")
  library("enrichplot")
  library("reshape2")
  library("enrichplot")
  library("ggraph")
  library("ggforce")
  library("reactome.db")
  library("ReactomePA")
  library("Biostrings")
})

hub <- AnnotationHub()
query(hub, "Homo_sapiens")
register(MulticoreParam(2))
ensEMBL2id.GO <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description', "entrezgene_id"), mart = ensembl) 

message("+--- Data sorting for clusterProfiler input format --------------+")

sel.column.names <-c("entrezgene_id", "ensembl_gene_id", "external_gene_name", "log2FoldChange")

resdf.ckegg <- list(); resdf.cBP <- list(); resdf.cMF <- list(); resdf.cCC <- list(); resdf.cPathway <- list();
resdf.ego <- list(); resdf.ego.UP <- list(); resdf.ego.DOWN <- list(); 
resdf.KEGG <- list(); resdf.KEGG.UP <- list(); resdf.KEGG.DOWN <- list()
resdf.gse <- list(); resdf.gse.UP <- list(); resdf.gse.DOWN <- list(); 
resdf.KEGG.gse <- list(); resdf.KEGG.gse.UP <- list(); resdf.KEGG.gse.DOWN <- list(); 

for(i in 1:8){
  load(paste0(baseDir, "/reverse_stranded_Analysis/",Project, "-", models[i], "_DESeq_Object.RData"))
  resdf.dat <- as.data.frame(subset(res, abs(log2FoldChange)>=l2fc & padj < significance))
  resdf.dat$external_gene_name <- rownames(resdf.dat)
  resdf.ann <- merge(resdf.dat, ensEMBL2id.GO, by = "external_gene_name")
  resdf     <- resdf.ann[,colnames(resdf.ann)%in%sel.column.names==T]
  resdf     <- resdf[order(-resdf$log2FoldChange),c(4,3,1,2)]
  colnames(resdf) <- c("ENTREZID", "ENSEMBL", "SYMBOL", "L2FC") 
  resdf     <- resdf[-which(duplicated(resdf$SYMBOL)), ]
  print(dim(resdf)) 
  print(length(unique(resdf$SYMBOL))) 
  print(length(unique(resdf$ENTREZID))) 
  
  resdf.rmNA <- subset(resdf[!is.na(resdf$ENTREZID),])
  resdf.rmNA <- resdf.rmNA %>%
    mutate(rank = rank(L2FC,  ties.method = "random")) 
  resdf.list <- resdf.rmNA$L2FC; names(resdf.list) <- resdf.rmNA$ENTREZID
  resdf.list <- resdf.list[!duplicated(names(resdf.list))]
  ## split up/down regulated genes
  resdf.up   <- subset(resdf, L2FC > 0)  
  resdf.down <- subset(resdf, L2FC < 0) 
  
  resdf.up.rmNA   <- subset(resdf.up[!is.na(resdf.up$ENTREZID),])
  resdf.down.rmNA <- subset(resdf.down[!is.na(resdf.down$ENTREZID),])
  resdf.up.rmNA <- resdf.up.rmNA %>%
    mutate(rank = rank(L2FC,  ties.method = "random"))
  resdf.down.rmNA <- resdf.down.rmNA %>%
    mutate(rank = rank(L2FC,  ties.method = "random"))
  
  resdf.uplist <- resdf.up.rmNA$L2FC; names(resdf.uplist) <- resdf.up.rmNA$ENTREZID
  resdf.uplist <- resdf.uplist[!duplicated(names(resdf.uplist))]
  resdf.downlist <- resdf.down.rmNA$L2FC; names(resdf.downlist) <- resdf.down.rmNA$ENTREZID
  resdf.downlist <- resdf.downlist[!duplicated(names(resdf.downlist))]
  
  
  ## use CompareCluster to generate the up/down GO/KEGG pathways
  
  resdf.gcsample <- list(X1=unique(resdf.up$ENTREZID[!is.na(resdf.up$ENTREZID)]), 
                         X2=unique(resdf.down$ENTREZID[!is.na(resdf.down$ENTREZID)]))
  
  save(resdf.gcsample, resdf, file = paste0(baseDir, "/reverse_stranded_Analysis/",Project,"-", models[i],"_GOData_list.RData"))
  
  ## Compare Up and Down regulated pathways for biological theme comparison
  resdf.ckegg[[i]] <- compareCluster(geneCluster = resdf.gcsample, fun = "enrichKEGG",
                                     organism="hsa",pvalueCutoff=0.05)
  
  resdf.cBP[[i]]   <- compareCluster(geneCluster = resdf.gcsample,  fun = "enrichGO",
                                   OrgDb=org.Hs.eg.db, pvalueCutoff=0.05,
                                   pAdjustMethod="BH", ont="BP", readable=T)
  resdf.cMF[[i]]   <- compareCluster(geneCluster = resdf.gcsample, fun = "enrichGO",
                                   OrgDb=org.Hs.eg.db,pvalueCutoff=0.05,
                                   pAdjustMethod="BH", ont="MF", readable=T)
  resdf.cCC[[i]]   <- compareCluster(geneCluster = resdf.gcsample,  fun = "enrichGO",
                                   OrgDb=org.Hs.eg.db, pvalueCutoff=0.05,
                                   pAdjustMethod="BH",ont="CC",readable=T)
  #
  resdf.cPathway[[i]] <- compareCluster(geneCluster = resdf.gcsample,  fun = "enrichPathway",
                                        pvalueCutoff=0.05, pAdjustMethod="BH", readable=T)
  
  ## enrich overrepresented pathways
  
  resdf.ego[[i]]        <- enrichGO(gene = resdf$SYMBOL, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',
                               ont = "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05 )
  resdf.ego.UP[[i]]     <- enrichGO(gene = resdf.up$SYMBOL, OrgDb = org.Hs.eg.db,keyType= 'SYMBOL',
                               ont = "ALL",pAdjustMethod = "BH", pvalueCutoff  = 0.05 )
  resdf.ego.DOWN[[i]]   <- enrichGO(gene = resdf.down$SYMBOL, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',
                                   ont = "ALL", pAdjustMethod = "BH", pvalueCutoff  = 0.05 )

  resdf.KEGG[[i]]       <- enrichKEGG(resdf$ENTREZID, organism = "hsa") 
  resdf.KEGG[[i]]       <- setReadable(resdf.KEGG[[i]], OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  resdf.KEGG.UP[[i]]    <- enrichKEGG(resdf.up$ENTREZID,  organism = "hsa")  
  resdf.KEGG.UP[[i]]    <- setReadable(resdf.KEGG.UP[[i]], OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  resdf.KEGG.DOWN[[i]]  <- enrichKEGG(resdf.down$ENTREZID, organism = "hsa")  
  resdf.KEGG.DOWN[[i]]  <- setReadable(resdf.KEGG.DOWN[[i]], OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  ## geneset enrichments
  resdf.gse[[i]]        <- gseGO(geneList = resdf.list, OrgDb  = org.Hs.eg.db, ont = "ALL", nPermSimple = 10000,
                                minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
  resdf.gse.UP[[i]]     <- gseGO(geneList = resdf.uplist, OrgDb  = org.Hs.eg.db, ont = "ALL", nPermSimple = 10000,
                                minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
  resdf.gse.DOWN[[i]]   <- gseGO(geneList = resdf.downlist, OrgDb  = org.Hs.eg.db, ont = "ALL", nPermSimple = 10000,
                                minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
  
  resdf.KEGG.gse[[i]]      <- gseKEGG(geneList = resdf.list, organism = "hsa", keyType = "kegg",nPermSimple = 10000,
                                      minGSSize = 10, maxGSSize = 800, pvalueCutoff = 0.05, pAdjustMethod = "BH")
  resdf.KEGG.gse.UP[[i]]   <- gseKEGG(geneList = resdf.uplist, organism = "hsa", keyType = "kegg",nPermSimple = 10000,
                                      minGSSize = 10, maxGSSize = 800, pvalueCutoff = 0.05, pAdjustMethod = "BH")
  resdf.KEGG.gse.DOWN[[i]] <- gseKEGG(geneList = resdf.downlist, organism = "hsa", keyType = "kegg", nPermSimple = 10000,
                                      minGSSize = 10, maxGSSize = 800, pvalueCutoff = 0.05, pAdjustMethod = "BH")
}

save(resdf.ckegg, resdf.cBP, resdf.cMF, resdf.cCC, resdf.cPathway, resdf.ego, resdf.ego.UP, resdf.ego.DOWN, 
     resdf.KEGG, resdf.KEGG.UP, resdf.KEGG.DOWN, resdf.gse, resdf.gse.UP, resdf.gse.DOWN, resdf.KEGG.gse, 
     resdf.KEGG.gse.UP, resdf.KEGG.gse.DOWN,
     file=paste0(baseDir, "/reverse_stranded_Analysis/",Project, "-GO_output_PathwayList_New_July.RData"))


message("+------------       TSC-3D vs TO GeneOntology semantic similarity Analysis        -----------------------+")
message("+--   semantic similarity among GO terms  rrvgo package (1.2.0), esp BP       ---------------------------+")
message("+Reduce and visualize lists of Gene Ontology terms by identifying redudance based on semantic similarity.+")

library("pathview")
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("org.Hs.eg.db")
library("rrvgo")
library("GOSemSim")

load(file=paste0(baseDir, "/reverse_stranded_Analysis/",Project, "-GO_output_PathwayList_New_July.RData"))
GOSeSimdata.BP   <- godata(OrgDb = "org.Hs.eg.db", keytype = "ENTREZID", ont="BP", computeIC = TRUE)

GOSemSim_Data_function <- function(goterm, index, gofun, model, goName, threshcut, Project, max.overlap){
  ego.all       <- goterm[[index]]
  ego.dat       <- as.data.frame(ego.all)
  ego.dat       <- subset(ego.dat, ONTOLOGY==gofun) 
  ego.dat.dim   <- dim(ego.dat)[1]
  
  simMatrix     <- calculateSimMatrix(ego.dat$ID, orgdb="org.Hs.eg.db", ont=gofun, semdata=GOSeSimdata.BP, method="Rel")
  scores        <- setNames(-log10(ego.dat$qvalue), ego.dat$ID)
  reducedTerms  <- reduceSimMatrix(simMatrix, scores, threshold=threshcut, orgdb="org.Hs.eg.db")
  maxclust      <- max(reducedTerms$cluster)
  
  message("+---------               Supplementary data for Fig. 4A                                  --------------+")
  
  write.csv(reducedTerms, file = paste0(baseDir, "/reverse_stranded_Analysis/",Project, "-", model, "_GOSeSim_", goName, "_", gofun, "_N", ego.dat.dim, "_thresh_", threshcut, "_cluster", maxclust, "_summaryTable.csv"), row.names=F)
  save(reducedTerms, simMatrix, file = paste0(baseDir, "/reverse_stranded_Analysis/",Project,"-", model,"_GOSeSim_", goName, "_", gofun, ".RData"))
  
  options(ggrepel.max.overlaps = max.overlap)
  
  message("+---------               Supplementary Fig for Fig. 4A                                  --------------+")
  
  pdf(paste0(baseDir, "/reverse_stranded_Analysis/",Project,"-", model, "_GOSeSim_", goName, "_", gofun, "_", ego.dat.dim, "_thresh_", threshcut, "_cluster", maxclust, "_scatterPlot.pdf"))
  scatterPlot(simMatrix, reducedTerms, size="score", labelSize = 2)
  dev.off()
}

allBP  <- GOSemSim_Data_function(resdf.ego, 2, "BP", models[2], "allBP", 0.9, Project, 160)
BP.UP  <- GOSemSim_Data_function(resdf.ego.UP, 2, "BP", models[2], "BPUP", 0.9, Project, 160)
BP.DW  <- GOSemSim_Data_function(resdf.ego.DOWN, 2, "BP", models[2], "BPDW", 0.9, Project, 50)


message("+----------- Regenerate Fig.3A with big font size (biological pathway, reactome) -----------------------+")

cPathway    <- resdf.cPathway[[2]]
cPathway@compareClusterResult$Cluster <- ifelse(cPathway@compareClusterResult$Cluster=="X1", "Up", "Down")
bold12.text <- element_text(face = "bold", size = 12)

pdf(paste0(baseDir, "/reverse_stranded_Analysis/Fig3A_Reactome_BiologicalPathway_UP36_DW15_DotPlot.pdf"), width= 10 , height=15)
dotplot(cPathway,  showCategory = 40, includeAll = T, title ="TSC-3D vs TO", font.size=12) + 
  theme(axis.text.x = bold12.text, axis.text.y = bold12.text, title =element_text(size=12, face='bold'),plot.title = element_text(hjust=0.5)) 
dev.off()

message("+---------               Supplementary data for Fig. 3A                                  --------------+")

write.csv(as.data.frame(cPathway), file = paste0(baseDir, "/reverse_stranded_Analysis/Supplementary_Data1_Fig3A_Reactome_BiologicalPathways_UP36_DW15_summary.csv"), row.names=F)

message("+---------     Fig. 4A  cNetworking for selected upregulated Biological process          --------------+")
message("+-------         TSC 3D vs Troph organoid TF list                                     -----------------+")

sigDEGs.files  <- paste0(baseDir, "/reverse_stranded_Analysis/CTR_myt25_0006-Okae_TOMvsorg_TOM_l2fc1_p0.05_N4037.csv")
dat            <- read.csv(sigDEGs.files)
updat          <- subset(dat, log2FoldChange >0)
colnames(updat)[6] <- "HGNC symbol"
TFs.newlist    <- read_excel("Human_TFList.xlsx", sheet=1)
newTF.sub      <- as.data.frame(subset(TFs.newlist, TFs.newlist[,4]=="Yes"|TFs.newlist[,19]=="Yes"|TFs.newlist[,21]=="Yes"|TFs.newlist[,25]=="Yes"))
updatTF.mer    <- merge(updat, TFs.newlist, by ="HGNC symbol")

write.csv(updatTF.mer, file = paste0(baseDir, "/reverse_stranded_Analysis/TSC3D_Trophorg_Up2257_overlap_TF_N328.csv"), row.names=F, quote=T)
##
load(paste0(baseDir, "/reverse_stranded_Analysis/",Project,"-", models[2], "_GOData_list.RData"))
cat1           <- c("extracellular matrix organization", "cell-substrate adhesion","regulation of cell-cell adhesion", "tissue migration")

BP.plt.up      <- cnetplot(nresdf.BP.UP, categorySize="pvalue", foldChange=geneList, 
                         layout="kk", colorEdge=T, node_label = "gene",
                         showCategory = cat1, 
                         cex_label_gene=0.8, cex_gene=1, font.size=1, lwd=2,
                         cex_category = 1)
BP.plt.nup    <- BP.plt.up + scale_colour_gradient(low ="darkorchid3",
                                                  high = "darkolivegreen3",
                                                  space = "Lab",
                                                  na.value = "grey50",
                                                  guide = "colourbar",
                                                  aesthetics = "colour")
pdf(paste0(baseDir, "/reverse_stranded_Analysis/Fig4A_TSC3DvsTrophorg_sel_BP_Uponly_cell_shape_organisation_cnetplot_colorscheme.pdf"), width=14, height=12)
print(BP.plt.nup)
dev.off()

message("+------------------   SFig.xx  individual genes boxplot                            -----------------+")
message("+- TFs  relating to HLA-structrue, similar as Fig 2C. order  TSC-2D, TSC-3D and TO------------------+")

## install.packages("ggpubr")
library(ggpubr)
library(plotly)

load(paste0(baseDir, "/reverse_stranded_Analysis/",Project, "-OriginalData_DESeq_dds.RData"))
DDS         <- dds.rm.collapse[,colData(dds.rm.collapse)$Group%in%c("TSC-2D", "TSC-3D", "TO")]
genes.HLAs  <- c("NLRC5", "RFX5", "RFXAP", "RFXANK", "TAPBP", "PDIA3", "CALR", "TAP1", "TAP2", "CIITA")


makeGeneCountPlot <- function(DDS, ensEMBL2id, CONDITION, gene2plot) {
  #
  # Plot the log2 normalised read counts for a specified gene
  #
  
  t2            <- plotCounts(DDS, gene=gene2plot, intgroup=c(CONDITION), normalized=TRUE, returnData=TRUE)
  colnames(t2)  <- c("count", "condition")
  t2$count      <- log2(t2$count)
  t2$condition  <- factor(t2$condition, levels=c("TSC-2D", "TSC-3D", "TO"))

  my_comparisons <- list( c("TSC-3D", "TSC-2D"), c("TO", "TSC-3D"), c("TO", "TSC-2D") )
  plt.cont      <-  ggplot(t2, aes(x=condition, y=count, fill=condition)) + 
    geom_boxplot(width = 0.5, color=c("darkgoldenrod4", "magenta", "cyan"),  alpha=0.5, outlier.shape=NA) + 
    geom_point(aes(fill=condition),position=position_jitterdodge(),size=1.2, alpha=0.5) +
    scale_fill_manual(name="Condition", values = c("darkgoldenrod4", "magenta", "cyan")) +
    scale_color_manual(name="Condition", values = c("darkgoldenrod4", "magenta", "cyan")) +
    theme(text = element_text(size=elementTextSize), legend.position="none") +
    ggtitle(gene2plot) + 
    xlab("") + ylab("log2(Normalised count)") +
    #ylim(-2, 16) +
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5, face="bold"),legend.position = "none",
          axis.text.x = element_text(size=elementTextSize1,angle = 45, hjust = 1, face="bold")) +
    stat_compare_means(comparisons = my_comparisons)
  print(paste("Created plot for", gene2plot), sep=" ")
  
  plt.cont
  
}
elementTextSize1 <- 10
plt.list <- list()
for(i in 1:length(genes.HLAs)){
  plt.list[[i]] <- makeGeneCountPlot(DDS, ensEMBL2id, "Group", genes.HLAs[i])
  plt.list
}

pdf(paste0(baseDir, "/reverse_stranded_Analysis/TFs_HLAs_2D_3D_Org.pdf"), height=20, width=6)

plot_grid(plt.list[[1]], plt.list[[2]],plt.list[[3]],plt.list[[4]],plt.list[[5]],
          plt.list[[6]], plt.list[[7]],plt.list[[8]],plt.list[[9]],plt.list[[10]],
          nrow=4, ncol=3)
dev.off()



message("+-------------------------------------------------------------------------------+")
message("+                              END OF SCRIPT                                    +")
message("+-------------------------------------------------------------------------------+")
