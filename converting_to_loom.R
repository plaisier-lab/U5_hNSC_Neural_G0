##########################################################
## ccAF:  converting_to_loom.R                          ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @Author:  Chris Plaisier, Samantha O'Connor          ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

# docker run -it -v '<where you have U5_hNSC_Neural_G0>:/files' cplaisier/scrna_seq_velocity
# R 3.6.2 Linux (Docker hub instance:  cplaisier/scrna_seq_velocity)
# Seurat 3.1.2

library(Seurat)
library(plyr)
library(dplyr)

library(org.Hs.eg.db)

# Set working directory
#setwd('/files/ccAF')

## WT
WT.data = Read10X(data.dir = "data/U5_hNSC/WT/filtered_gene_bc_matrices/hg19/")
WT = CreateSeuratObject(counts = WT.data, project = "WT", min.cells = 3, min.features = 200)
WT[["percent.mt"]] = PercentageFeatureSet(WT, pattern = "^MT-")
pdf('WT_mt_vs_Count.pdf')
FeatureScatter(WT, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()
WT = subset(WT, subset = nCount_RNA < 18000 & percent.mt < 6)
WT = NormalizeData(WT)
WT = FindVariableFeatures(WT, selection.method = "vst", nfeatures = 3000)
all.genes = rownames(WT)
WT = ScaleData(WT, features = all.genes)
WT[['clusts_named']] = clusts_named[tmp5[paste(colnames(WT),'_1',sep='')],'clusts_named']
WT[['clusts_num']] = clusts_num[tmp5[paste(colnames(WT),'_1',sep='')],'clusts_named']
Idents(WT) = 'clusts_named'
WT.markers = FindAllMarkers(WT, only.pos=F, min.pct=0.25, logfc.threshold = 0.2)
WT_marker_genes_0_3 = (WT.markers %>% group_by(cluster) %>% filter(avg_logFC>0.3))$gene
WT_marker_genes_0_5 = (WT.markers %>% group_by(cluster) %>% filter(avg_logFC>0.5))$gene
WT_marker_gene_top50 = (WT.markers %>% group_by(cluster) %>% top_n(50,avg_logFC))$gene
#as.loom(WT, assay='RNA', filename='WT_6_1_2020.loom')

## sgTAOK1
taok1.data = Read10X(data.dir = "data/U5_hNSC/sgTAOK1/filtered_gene_bc_matrices/hg19/")
taok1 = CreateSeuratObject(counts = taok1.data, project = "sgTAOK1", min.cells = 3, min.features = 200)
taok1[["percent.mt"]] = PercentageFeatureSet(taok1, pattern = "^MT-")
pdf('taok1_mt_vs_Count.pdf')
FeatureScatter(taok1, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()
taok1 = subset(taok1, subset = nCount_RNA < 18000 & percent.mt < 10)
taok1 = NormalizeData(taok1)
taok1 = FindVariableFeatures(taok1, selection.method = "vst", nfeatures = 3000)
all.genes = rownames(taok1)
taok1 = ScaleData(taok1, features = all.genes)
taok1[['clusts_named']] = clusts_named[tmp5[paste(colnames(taok1),'_2',sep='')],'clusts_named']
taok1[['clusts_num']] = clusts_num[tmp5[paste(colnames(taok1),'_2',sep='')],'clusts_named']
Idents(taok1) = 'clusts_named'
taok1.markers = FindAllMarkers(taok1, only.pos=F, min.pct=0.25, logfc.threshold = 0.2)
taok1_marker_genes_0_3 = (taok1.markers %>% group_by(cluster) %>% filter(avg_logFC>0.3))$gene
taok1_marker_genes_0_5 = (taok1.markers %>% group_by(cluster) %>% filter(avg_logFC>0.5))$gene
taok1_marker_gene_top50 = (taok1.markers %>% group_by(cluster) %>% top_n(50,avg_logFC))$gene
#as.loom(taok1, assay='RNA', filename='sgTAOK1_6_1_2020.loom')

## Highly varaible gene list WT int sgTAOK1 (n = 1584), uses vst method
hvg1 = intersect(WT@assays$RNA@var.features, taok1@assays$RNA@var.features)
write.csv(hvg1,'highlyVarGenes_WT_sgTAOK1_1584.csv')

## Integrate the datasets
int_anchors = FindIntegrationAnchors(object.list = list(WT = WT, taok1 = taok1), dims=1:30, anchor.features = 10000)
cellcycle_int_2 = IntegrateData(anchorset = int_anchors, dims = 1:30, features.to.integrate=unique(c(hvg1, mg1, mg2, mg3, hvg2, mg4, mg5, mg6, int_anchors@anchor.features)))
DefaultAssay(cellcycle_int_2) = "integrated"
as.loom(cellcycle_int_2, assay='integrated', filename='cellcycle_int_integrated_V3_6_18_2020.loom')

