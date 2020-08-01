##########################################################
## ccAF:  makeDatasetsForClassification.R               ##
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

### Docker command to start-up Seurat capable analysis container
# docker run -it -v '/home/cplaisier/Dropbox (ASU):/files' cplaisier/seurat_v3.0

### Packages required to run analyses
library(dplyr)
library(plyr)
library(Seurat)
library(ranger)

### Wang and Verhaak et al. GBM subtype classifier
# R CMD INSTALL ssgsea.GBM.classification_1.0.tar.gz
library(ssgsea.GBM.classification)
source('ssgsea.GBM.classification/R/runSsGSEAwithPermutationR3.R')

### Modified functions to classify cells from Seurat V2
ClassifyCells2 <- function (object, classifier, training.genes = NULL, training.classes = NULL, 
    new.data = NULL, ...) 
{
    features <- attr(classifier$terms,'term.labels')
    genes.to.add <- setdiff(x = features, y = rownames(x = new.data))
    data.to.add <- matrix(data = 0, nrow = length(x = genes.to.add), 
        ncol = ncol(x = new.data))
    rownames(x = data.to.add) <- genes.to.add
    colnames(data.to.add) <- colnames(new.data)
    new.data <- rbind(new.data, data.to.add)
    new.data <- new.data[features, ]
    new.data <- as.matrix(x = t(x = new.data))
    message("Running Classifier ...")
    prediction <- predict(classifier, new.data, type='response')
    new.classes <- prediction$predictions
    return(new.classes)
}

ClassifyCells3 <- function (classifier, new.data = NULL, ...) 
{
    features <- classifier$forest$independent.variable.names
    genes.to.add <- setdiff(x = features, y = rownames(x = new.data))
    data.to.add <- matrix(data = 0, nrow = length(x = genes.to.add), 
        ncol = ncol(x = new.data))
    rownames(x = data.to.add) <- genes.to.add
    colnames(data.to.add) <- colnames(new.data)
    new.data <- rbind(new.data, data.to.add)
    new.data <- new.data[features, ]
    new.data <- as.matrix(x = t(x = new.data))
    message("Running Classifier ...")
    prediction <- predict(classifier, new.data)
    new.classes <- prediction$predictions
    return(new.classes)
}


##########################
### GSE70630 GBM IDH+/- ##
##########################
# Load up data
setwd('/files/scGlioma/GSE70630')
d1 = read.table(gzfile('GSE70630_OG_processed_data_v2.txt.gz'),sep='\t',row.names=1,header=T)
d1 = d1[which(apply(d1,1,sum)!=0),]
dim(d1)
n1 = sapply(colnames(d1), function(x) { strsplit(x,'_')[[1]][1] } )
n1[n1=='X93'] = 'MGH93'
n1[n1=='X97'] = 'MGH97'
GSE70630 = CreateSeuratObject(counts=d1)
GSE70630[['patient']] = as.factor(sapply(sapply(colnames(d1), function(x) { strsplit(x,"_")[[1]][1] }), function(x) { gsub('X','MGH',toupper(strsplit(x,"\\.")[[1]][1])) }))

# Find variable features, scale data, and do dimensionality reduction
#GSE70630 = NormalizeData(object = GSE70630)
GSE70630 = FindVariableFeatures(object = GSE70630)
GSE70630 = ScaleData(object = GSE70630)
GSE70630 = RunPCA(object = GSE70630, features = VariableFeatures(GSE70630))

# Determined that 14 PCs looked reasonable
pdf('GSE70630_PCA.pdf')
DimHeatmap(GSE70630, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(GSE70630, dims = 16:30, cells = 500, balanced = TRUE)
dev.off()

# Project into two dimensions and cluster the data
GSE70630 = RunUMAP(GSE70630, dims=1:14)
GSE70630 = FindNeighbors(GSE70630, dims=1:14)
GSE70630 = FindClusters(GSE70630, resolution = 0.5)

# Predict our cell_cycle
clusts.WT.GSE70630 = ClassifyCells3(GSE70630,classifier=ccAF,new.data=as.matrix(GSE70630@assays$RNA@data))
clusts.WT.GSE70630 = revalue(clusts.WT.GSE70630, c('G1/Differentiation'='Neural G0'))
names(clusts.WT.GSE70630) = colnames(GSE70630@assays$RNA@data)
GSE70630[['cell_cycle']] = clusts.WT.GSE70630
table(GSE70630[['cell_cycle']])

# Classify using Wang et al., 2017 (Verhaak) classifier
#res1 = runSsGSEAwithPermutation('GSE70630.gct',1000)
res1 = read.table('data/ssGSEA.GBM.classification/p_result_GSE70630.gct.txt', sep='\t', header=T, row.names=1)
min1 = apply(res1[,4:6],1,min)
calls1 = c('PN','CL','MS')[sapply(rownames(res1), function(x){ which(res1[x,4:6]==min1[x])[[1]] })]
names(calls1) = rownames(res1)
GSE70630[['subtype']] = calls1[colnames(GSE70630)]
table(GSE70630[[c('subtype','cell_cycle')]])
table(GSE70630[[c('subtype','patient')]])

# Run the Seurat cell_cycle classifeir against the data
GSE70630 = CellCycleScoring(GSE70630, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

# Plot all the cell data
pdf('GSE70630_allCells.pdf')
DimPlot(GSE70630, reduction='umap', label=T)
DimPlot(GSE70630, group.by='patient', reduction='umap', label=T)
DimPlot(GSE70630, group.by='cell_cycle', reduction='umap', label=T)
FeaturePlot(GSE70630, features = c('CD3','PTPRC','CD14','CX3CR1','AIF1','MBP','PLP1','GFAP','PRC1'))
dev.off()

# Removed:
# Immune cells = 8
# Oligodendrocytes = 12
GSE70630_clean = subset(GSE70630, cells = rownames(GSE70630@meta.data)[which(!GSE70630[['RNA_snn_res.0.5']][,1] %in% c(8,12))])
GSE70630_clean = RunUMAP(GSE70630_clean, dims=1:14)
dim(GSE70630) # => 4347 cells total
dim(GSE70630_clean) # => 4047 after our filtering, and there were 4044 malignant cells in the publication
as.loom(GSE70630_clean, assay='RNA', filename='GSE70630.loom')

# Final plot
pdf('GSE70630_clean.pdf')
DimPlot(GSE70630_clean, reduction='umap', label=T)
DimPlot(GSE70630_clean, group.by='patient', reduction='umap', label=T)
DimPlot(GSE70630_clean, group.by='cell_cycle', reduction='umap', label=T)
DimPlot(GSE70630_clean, group.by='subtype', reduction='umap', label=T)
FeaturePlot(GSE70630_clean, features = c('CD3','PTPRC','CD14','CX3CR1','AIF1','MBP','PLP1','GFAP','PRC1'))
dev.off()


######################################
### GSE89567 GBM IDH+/-            ###
### Only include grade III tumors: ###
### MGH103, MGH42-45, MGH56, MGH57 ###
######################################
setwd('/files/scGlioma/GSE89567')
d1 = read.table(gzfile('GSE89567_IDH_A_processed_data.txt.gz'),sep='\t',row.names=1,header=T)
dim(d1)
d1 = d1[which(apply(d1,1,sum)!=0),]
dim(d1)
n1 = sapply(colnames(d1), function(x) { strsplit(x,'_')[[1]][1] } )
n1 = sapply(n1, function(x) { strsplit(x,'\\.')[[1]][1] } )
n1[n1=='X57'] = 'MGH57'
n1[n1=='mgh103'] = 'MGH103'
GSE89567 = CreateSeuratObject(counts=d1)
GSE89567[['patient']] = as.factor(sapply(sapply(colnames(d1), function(x) { strsplit(x,"_")[[1]][1] }), function(x) { gsub('X','MGH',toupper(strsplit(x,"\\.")[[1]][1])) }))

#GSE89567 = NormalizeData(object = GSE89567)
GSE89567 = FindVariableFeatures(object = GSE89567)
GSE89567 = ScaleData(object = GSE89567)
GSE89567 = RunPCA(object = GSE89567, features = VariableFeatures(GSE89567))
pdf('GSE89567_PCA.pdf')
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 16:30, cells = 500, balanced = TRUE)
dev.off()

GSE89567 = RunUMAP(GSE89567, dims=1:16)
GSE89567 = FindNeighbors(GSE89567, dims=1:16)
GSE89567 = FindClusters(GSE89567, resolution = 0.5)

# Predict our cell_cycle
clusts.WT.GSE89567 = ClassifyCells3(GSE89567,classifier=ccAF,new.data=as.matrix(GSE89567@assays$RNA@data))
clusts.WT.GSE89567 = revalue(clusts.WT.GSE89567, c('G1/Differentiation'='Neural G0'))
names(clusts.WT.GSE89567) = colnames(GSE89567@assays$RNA@data)
GSE89567[['cell_cycle']] = clusts.WT.GSE89567
table(GSE89567[['cell_cycle']])

# Classify using Wang et al., 2017 (Verhaak) classifier
#res1 = runSsGSEAwithPermutation('GSE89567.gct',1000)
res1 = read.table('data/ssGSEA.GBM.classification/p_result_GSE89567.gct.txt', sep='\t', header=T, row.names=1)
min1 = apply(res1[,4:6],1,min)
calls1 = c('PN','CL','MS')[sapply(rownames(res1), function(x){ which(res1[x,4:6]==min1[x])[[1]] })]
names(calls1) = rownames(res1)
GSE89567[['subtype']] = calls1[colnames(GSE89567)]
table(GSE89567[[c('subtype','cell_cycle')]])
table(GSE89567[[c('subtype','patient')]])

# Run the Seurat cell_cycle classifeir against the data
GSE89567 = CellCycleScoring(GSE89567, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

pdf('GSE89567_allCells.pdf')
DimPlot(GSE89567, reduction='umap', label=T)
DimPlot(GSE89567, group.by='patient', reduction='umap', label=T)
DimPlot(GSE89567, group.by='cell_cycle', reduction='umap', label=T)
FeaturePlot(GSE89567, features = c('CD3','PTPRC','CD14','CX3CR1','AIF1','MBP','PLP1','GFAP','PRC1'))
dev.off()

# Removed:
# Immune cells = 4, 7 & 11
# Oligodendrocytes = 13
GSE89567_clean = subset(GSE89567, cells = intersect(rownames(GSE89567@meta.data)[which(!GSE89567[['patient']][,1] %in% c('MGH107NEG','MGH107POS','MGH61','MGH64'))], rownames(GSE89567@meta.data)[which(!GSE89567[['RNA_snn_res.0.5']][,1] %in% c(4,7,11,13))]))
GSE89567_clean = RunUMAP(GSE89567_clean, dims=1:16)
as.loom(GSE89567_clean, assay='RNA', filename='GSE89567.loom')

pdf('GSE89567_clean.pdf')
DimPlot(GSE89567_clean, reduction='umap', label=T)
DimPlot(GSE89567_clean, group.by='patient', reduction='umap', label=T)
DimPlot(GSE89567_clean, group.by='cell_cycle', reduction='umap', label=T)
DimPlot(GSE89567_clean, group.by='Phase', reduction='umap', label=T)
DimPlot(GSE89567_clean, group.by='subtype', reduction='umap', label=T)
FeaturePlot(GSE89567_clean, features = c('CD3','PTPRC','CD14','CX3CR1','AIF1','MBP','PLP1','GFAP','PRC1'))
dev.off()


#####################
### GSE131928 GBM ###
#####################
library(ggplot2)
library(cowplot)

# Load phenotype information from paper
pheno1 = read.csv('patientData.csv',header=T,row.names=1)

# Load 10X data
d1 = read.table(gzfile('GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv'),sep='\t',row.names=1,header=T)
GSE131928_10X = CreateSeuratObject(counts = d1)
GSE131928_10X = FindVariableFeatures(GSE131928_10X, selection.method='vst', nfeatures=2000)
GSE131928_10X@assays$RNA@var.features = c(GSE131928_10X@assays$RNA@var.features,'PRC1')

# Get patient names for 10X
names1 = colnames(d1)
cur1 = c('X105_B1','X105_B2','X105_C1','X105_C2','X105_D1','X105_D2')
final1 = c('X105B','X105B','X105C','X105C','X105D','X105D')
for(i in 1:length(cur1)) {
    for(j in 1:length(names1)) {
        names1[j] = gsub(cur1[i],final1[i], names1[j])
    }
}
n1 = sapply(names1, function(x) { paste('MGH',strsplit(gsub('X','',x),'_')[[1]][1],sep='') } )
names(n1) = colnames(d1)
GSE131928_10X[['patient']] = n1

# Remove MGH114 and MGH118 because we have no clinical information as per the paper,
# and MGH105 because it is an outlier on 10X (still some cells in Smartseq2)
noPhenos = c("MGH114", "MGH118","MGH105A","MGH105B","MGH105C","MGH105D")
GSE131928_10X = subset(GSE131928_10X, cells = rownames(GSE131928_10X@meta.data)[which(!GSE131928_10X[['patient']][,1] %in% noPhenos)])

d1 = read.csv('/files/ccAF/classifyPrimaryCells/results/ccAF_results_GSE131928_10X.csv',header=T,row.names=1)

GSE131928_10X[['ccAF']] = d1[,'Predictions',drop=F]
GSE131928_10X[['subtype']] = d1[,'subtype',drop=F]

##################
## All Patients ##
##################
patients = c("MGH102", "MGH143", "MGH115", "MGH124", "MGH125")
pcs = list()
pcs$MGH102 = 12
pcs$MGH115 = 6 # 8
pcs$MGH124 = 14
pcs$MGH125 = 11
pcs$MGH143 = 10
res = list()
res$MGH102 = 0.3
res$MGH115 = 0.2
res$MGH124 = 0.5
res$MGH125 = 0.4
res$MGH143 = 0.35

for(pat1 in c('MGH143')) {
#for(pat1 in patients) {

        patSeurat = subset(GSE131928_10X, cells = rownames(GSE131928_10X@meta.data)[which(GSE131928_10X[['patient']][,1] == pat1)])
        patSeurat = subset(patSeurat, cells = rownames(patSeurat@meta.data)[which(!is.na(patSeurat[['ccAF']][,1]))])
        
        patSeurat = ScaleData(patSeurat)
        patSeurat = RunPCA(patSeurat, features = VariableFeatures(object = patSeurat))
        patSeurat@meta.data[,'ccAF'] = factor(patSeurat@meta.data[,'ccAF'], levels=c('G1/other','Neural G0','G1','Late G1','S','S/G2','G2/M','M/Early G1'))

        pdf(paste(pat1,'_PCA.pdf',sep=''),width = 8.5, height=11)
        DimHeatmap(patSeurat, dims = 1:15, cells = 500, balanced = TRUE)
        dev.off()
        
        patSeurat = FindNeighbors(patSeurat, dims = 1:pcs[[pat1]])
        patSeurat = FindClusters(patSeurat, resolution = 0.1)
        patSeurat = FindClusters(patSeurat, resolution = 0.2)
        patSeurat = FindClusters(patSeurat, resolution = 0.3)
        patSeurat = FindClusters(patSeurat, resolution = 0.35)
        patSeurat = FindClusters(patSeurat, resolution = 0.4)
        patSeurat = FindClusters(patSeurat, resolution = 0.5)
        patSeurat = FindClusters(patSeurat, resolution = 0.6)
        patSeurat = FindClusters(patSeurat, resolution = 0.7)
        patSeurat = FindClusters(patSeurat, resolution = 0.8)
        patSeurat = RunUMAP(patSeurat, dims = 1:pcs[[pat1]])
        
        pdf(paste(pat1,'_UMAP_res_test.pdf',sep=''),width = 11, height=8)
        a = DimPlot(patSeurat, reduction = "umap", group.by=c('RNA_snn_res.0.1','RNA_snn_res.0.2','RNA_snn_res.0.3','RNA_snn_res.0.4','RNA_snn_res.0.5','RNA_snn_res.0.6','RNA_snn_res.0.7','RNA_snn_res.0.8','ccAF','subtype'))
        print(a)
        dev.off()
        
        # Get marker genes
        patSeurat = FindClusters(patSeurat, resolution = res[[pat1]])
        png(paste(pat1,'_UMAP_res_Final.png',sep=''),width =6, height=2.5, units='in',res=300)
        a = DimPlot(patSeurat, reduction = "umap", group.by=c('ccAF',paste('RNA_snn_res.',res[[pat1]],sep='')))
        print(a)
        dev.off()
        png(paste(pat1,'_featurePlot_res_Final.png',sep=''),width =6, height=2.5, units='in',res=300)
        a = FeaturePlot(patSeurat, features = c('MBP','ETNPPL'), reduction = "umap")
        print(a)
        dev.off()
        pdf(paste(pat1,'_featurePlot_Stem_res_Final.pdf',sep=''),width =8, height=3)
        a = VlnPlot(patSeurat, features = c('SOX2','NANOG','POU5F1','PROM1'), group.by='ccAF', pt.size=0.25, ncol=4)
        print(a)
        dev.off()
        pdf(paste(pat1,'_featureScatter_res_Final.pdf',sep=''),width = 8, height=3)
       
        a = FeatureScatter(patSeurat, feature1 = 'SOX2', feature2 = 'NANOG', group.by='ccAF')+theme(legend.position='none')
        b = FeatureScatter(patSeurat, feature1 = 'SOX2', feature2 = 'POU5F1', group.by='ccAF')+theme(legend.position='none')
        c = FeatureScatter(patSeurat, feature1 = 'SOX2', feature2 = 'PROM1', group.by='ccAF')+theme(legend.position='none')
        d = CombinePlots(list(a,b,c), ncol=3)
        print(d)
        dev.off()
}

GSE131928_10X@meta.data[,'ccAF'] = factor(GSE131928_10X@meta.data[,'ccAF'], levels=c('G1/other','Neural G0','G1','Late G1','S','S/G2','G2/M','M/Early G1'))
GSE131928_10X_clean = subset(GSE131928_10X, cells = rownames(GSE131928_10X@meta.data)[which(GSE131928_10X[['patient']][,1] %in% patients)])
GSE131928_10X_clean = subset(GSE131928_10X_clean, cells = rownames(GSE131928_10X_clean@meta.data)[which(!is.na(GSE131928_10X_clean[['ccAF']][,1]))])
GSE131928_10X_clean = subset(GSE131928_10X_clean, cells = rownames(GSE131928_10X_clean@meta.data)[which(GSE131928_10X_clean[['ccAF']][,1]!='G1/other')])
pdf('GSE131928_10X_featurePlot_Stem_res_Final.pdf',width =8, height=3)
a = VlnPlot(GSE131928_10X_clean, features = c('SOX2','NANOG','POU5F1','PROM1'), group.by='ccAF', pt.size=0.25, ncol=4)
print(a)
a = VlnPlot(GSE131928_10X_clean, features = c('SOX2','NANOG','POU5F1','PROM1'), group.by='patient', pt.size=0.25, ncol=4)
print(a)
a = VlnPlot(GSE131928_10X_clean, features = 'PROM1', group.by='ccAF', pt.size=0.25,split.by='patient')
print(a)
dev.off()

pdf('GSE131928_10X_featureScatter_res_Final.pdf',width = 8, height=3)
a = FeatureScatter(GSE131928_10X_clean, feature1 = 'SOX2', feature2 = 'NANOG', group.by='ccAF')+theme(legend.position='none')+xlim(0,7000)
b = FeatureScatter(GSE131928_10X_clean, feature1 = 'SOX2', feature2 = 'POU5F1', group.by='ccAF')+theme(legend.position='none')+xlim(0,7000)
c = FeatureScatter(GSE131928_10X_clean, feature1 = 'SOX2', feature2 = 'PROM1', group.by='ccAF')+theme(legend.position='none')+xlim(0,7000)
d = CombinePlots(list(a,b,c), ncol=3)
print(d)
dev.off()

GSE131928_10X_clean = ScaleData(GSE131928_10X_clean, features=rownames(GSE131928_10X))

pdf('GSE131928_10X_DotPlot_res_Final.pdf',width = 8, height=3)
a = DotPlot(GSE131928_10X_clean, features = c('NANOG','POU5F1','PROM1'), group.by='ccAF')
print(a)
a = DotPlot(GSE131928_10X_clean, features = c('NANOG','POU5F1','PROM1'), group.by='patient')
print(a)
a = DoHeatmap(GSE131928_10X_clean, features = c('SOX2','NANOG','POU5F1','PROM1'), group.by='ccAF')
print(a)
a = DoHeatmap(GSE131928_10X_clean, features = c('SOX2','NANOG','POU5F1','PROM1'), group.by='patient')
print(a)
dev.off()

## Load up Smartseq2 data
d2 = as.sparse(read.table(gzfile('GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv'),sep='\t',row.names=1,header=T))
GSE131928_Smartseq2 = CreateSeuratObject(counts = d2)
GSE131928_Smartseq2 = NormalizeData(GSE131928_Smartseq2)
GSE131928_Smartseq2 = FindVariableFeatures(GSE131928_Smartseq2, slection.method='vst', nfeatures=2000)
n2 = sapply(colnames(d2), function(x) { strsplit(x,'\\.')[[1]][1] } )
GSE131928_Smartseq2[['patient']] = n2

# Remove pediatric patients from Smartseq2 dataset
pediatric = c("BT749", "BT771", "BT830", "MGH85", "BT1160", "BT1187", "BT786", "BT920", "MGH106CD3posP1")
GSE131928_Smartseq2 = subset(GSE131928_Smartseq2, cells = rownames(GSE131928_Smartseq2@meta.data)[which(!GSE131928_Smartseq2[['patient']][,1] %in% pediatric)])

d2 = read.csv('/files/ccAF/classifyPrimaryCells/results/ccAF_results_GSE131928_Smartseq2.csv',header=T,row.names=1)

GSE131928_Smartseq2[['ccAF']] = d2[,'Predictions',drop=F]
GSE131928_Smartseq2[['subtype']] = d2[,'subtype',drop=F]

## Merging the datasets using anchors
reference.list = list(GSE131928_10X, GSE131928_Smartseq2)
GSE131928_anchors = FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
GSE131928_integrated = IntegrateData(anchorset = GSE131928_anchors, dims = 1:30)

# Make the integrated data the default assay, scale, and reduce dimensionality
DefaultAssay(GSE131928_integrated) = "integrated"
GSE131928_integrated = ScaleData(GSE131928_integrated)
GSE131928_integrated = RunPCA(GSE131928_integrated, npcs=30)
GSE131928_integrated = RunUMAP(GSE131928_integrated, reduction='pca', dims=1:30)

# Add key phenotypic information
GSE131928_integrated[['patient']] = c(n1,n2)
GSE131928_integrated[['platform']] = c(rep('10X',ncol(GSE131928_10X)), rep('Smartseq2',ncol(GSE131928_Smartseq2)))
GSE131928_integrated[['Name']] =  pheno1[GSE131928_integrated@meta.data$patient,'Name']
GSE131928_integrated[['Type']] =  pheno1[GSE131928_integrated@meta.data$patient,'Type']
GSE131928_integrated[['Sex']] =  as.factor(pheno1[GSE131928_integrated@meta.data$patient,'Sex'])
GSE131928_integrated[['Recurrent']] =  pheno1[GSE131928_integrated@meta.data$patient,'Primary.or.Recurrence']
GSE131928_integrated[['IDH']] =  pheno1[GSE131928_integrated@meta.data$patient,'IDH']
GSE131928_integrated[['EGFR']] =  pheno1[GSE131928_integrated@meta.data$patient,'EGFR']
GSE131928_integrated[['MET']] =  pheno1[GSE131928_integrated@meta.data$patient,'MET']
GSE131928_integrated[['MGMT']] =  pheno1[GSE131928_integrated@meta.data$patient,'MGMT']

# Run the Seurat cell_cycle classifier against the data
GSE131928_integrated = CellCycleScoring(GSE131928_integrated, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
GSE131928_integrated[['Phase']] = factor(GSE131928_integrated[['Phase']][,1], levels=levels(GSE131928_integrated[['Phase']][,1])[c(1,3,2)])

# Load up the subtype computations from previous analyses
res_all = read.csv('data/ssGSEA.GBM.classification/res_All_GSE131928.csv',header=T,row.names=1)
min1 = apply(res_all[,4:6],1,min)
calls1 = c('PN','CL','MS')[sapply(rownames(res_all), function(x){ which(res_all[x,4:6]==min1[x])[[1]] })]
names(calls1) = rownames(res_all)
GSE131928_integrated[['subtype']] = calls1

# Compute clustering based on data
GSE131928_integrated = FindNeighbors(GSE131928_integrated, dims=1:30)
GSE131928_integrated = FindClusters(GSE131928_integrated, resolution=0.5)

# Plot 
pdf('GSE131928_clusters_UMAP.pdf',width = 8.5, height=11)
DimPlot(GSE131928_integrated, reduction = 'umap', label=T)+theme(legend.position='bottom', legend.text=element_text(size=8))+ggtitle(label = 'Clustering')
FeaturePlot(GSE131928_integrated, features =c('CD3','PTPRC','CD14','CX3CR1','AIF1','MBP','PLP1','GFAP','PRC1','CD31','ETNPPL'))
dev.off()

# Remove non-tumor cells
# Immune = 6,9,10,12,16,20,21
# Oligodendrocytes = 14
clustsToRemove = c(2,8,11,12,16,17,18,19)
GSE131928_integrated_tumor = subset(GSE131928_integrated, cells = rownames(GSE131928_integrated@meta.data)[which(!GSE131928_integrated[['integrated_snn_res.0.5']][,1] %in% clustsToRemove)])
GSE131928_integrated_tumor = ScaleData(GSE131928_integrated_tumor)
GSE131928_integrated_tumor = RunPCA(GSE131928_integrated_tumor, npcs=30)
GSE131928_integrated_tumor = RunUMAP(GSE131928_integrated_tumor, reduction='pca', dims=1:30)
as.loom(GSE131928_integrated_tumor, assay='RNA', filename='GSE131928.loom')

# Subset cells
cells1 = colnames(GSE131928_integrated_tumor)
GSE131928_10X = GSE131928_10X[,intersect(cells1, colnames(GSE131928_10X))]
GSE131928_10X[['subtype']] = calls1[intersect(cells1, colnames(GSE131928_10X))]
as.loom(GSE131928_10X, assay='RNA', filename='GSE131928_10X.loom')
GSE131928_Smartseq2 = GSE131928_Smartseq2[,intersect(cells1, colnames(GSE131928_Smartseq2))]
GSE131928_Smartseq2[['cell_cycle']] = cc_clusts2
GSE131928_Smartseq2[['subtype']] = calls1[intersect(cells1, colnames(GSE131928_Smartseq2))]
as.loom(GSE131928_Smartseq2, assay='RNA', filename='GSE131928_Smartseq2.loom')


# Final plots
colors1 = c(rgb(123/255,175/255,65/255),rgb(201/255,149/255,43/255),rgb(243/255,118/255,110/255),rgb(31/255,189/255,194/255),rgb(166/255,129/255,186/255),rgb(225/255,109/255,170/255),rgb(43/255,181/255,103/255),rgb(74/255,161/255,217/255))
colors2 = c(rgb(243/255,118/255,110/255),rgb(166/255,129/255,186/255),rgb(43/255,181/255,103/255))
GSE131928_integrated_tumor@meta.data$ccAF = factor(GSE131928_integrated_tumor@meta.data$ccAF, levels=c('G1/other','Neural G0','G1','Late G1','S','S/G2','G2/M','M/Early G1'))
pdf('GSE131928_cell_cycle_UMAP_7_24_2020.pdf',width = 8.5, height=3)
p1 = DimPlot(GSE131928_integrated_tumor, reduction = 'umap', group.by='ccAF', cols=colors1)+theme(legend.position='None')+xlab('UMAP1')+ylab('UMAP2')
p2 = DimPlot(GSE131928_integrated_tumor, reduction = 'umap', group.by='Phase', cols=colors2)+theme(legend.position='None')+xlab('UMAP1')+ylab('UMAP2')
p4 = DimPlot(GSE131928_integrated_tumor, reduction = 'umap', group.by='subtype', label=T)+theme(legend.position='None')+xlab('UMAP1')+ylab('UMAP2')
CombinePlots(plots = list(p1, p2, p4), legend=NULL, ncol=3)
dev.off()

pdf('barplots.pdf', width=8.5, height=3)
p1 = ggplot(data=data.frame(table(GSE131928_integrated_tumor[['ccAF']][,1])), aes(x=Var1, y=Freq,fill=Var1))+geom_bar(stat="identity", position=position_dodge())+scale_fill_manual(values=colors1)+theme(legend.position='None',axis.title.x = element_blank(),axis.text.x = element_text(angle = 45))+ylab('Cells')
p2 = ggplot(data=data.frame(table(GSE131928_integrated_tumor[['Phase']][,1])), aes(x=Var1, y=Freq, fill=Var1))+geom_bar(stat="identity", position=position_dodge())+scale_fill_manual(values=colors2)+theme(legend.position='None',axis.title.x = element_blank())+ylab('Cells')
p3 = ggplot(data=data.frame(table(GSE131928_integrated_tumor[['ccAF']][,1],GSE131928_integrated_tumor[['subtype']][,1])), aes(x=Var2, y=Freq, fill=Var1))+geom_bar(stat="identity", position=position_dodge())+scale_fill_manual(values=colors1)+theme(legend.position='None',axis.title.x = element_blank())+ylab('Cells')
CombinePlots(plots = list(p1, p2, p3), legend=NULL, ncol=3)
dev.off()

## GSE131928_10X
# Bhaduri putative stem cells
res1 = matrix(nrow=7,ncol=5, dimnames=list(c('Neural G0','G1','Late G1','S','S/G2','G2/M','M/Early G1'),c('q','m','n','k','pv')))
for(clust1 in c('Neural G0','G1','Late G1','S','S/G2','G2/M','M/Early G1')) {
        q = sum(GSE131928_10X_clean@meta.data$ccAF==clust1 & (((GSE131928_10X_clean@assays$RNA@data['PROM1',]>0 | GSE131928_10X_clean@assays$RNA@data['FUT4',]>0 | GSE131928_10X_clean@assays$RNA@data['L1CAM',]>0) & GSE131928_10X_clean@assays$RNA@data['SOX2',]>0) & GSE131928_10X_clean@assays$RNA@data['TLR4',]==0),na.rm=T)
        m = sum((((GSE131928_10X_clean@assays$RNA@data['PROM1',]>0 | GSE131928_10X_clean@assays$RNA@data['FUT4',]>0 | GSE131928_10X_clean@assays$RNA@data['L1CAM',]>0) & GSE131928_10X_clean@assays$RNA@data['SOX2',]>0) & GSE131928_10X_clean@assays$RNA@data['TLR4',]==0),na.rm=T)
        n = dim(GSE131928_10X_clean)[2]-m
        k = sum(GSE131928_10X_clean@meta.data$ccAF==clust1,na.rm=T)
        pv = phyper(q, m, n, k, lower.tail=F)
        res1[clust1,] = c(q,m,n,k,pv)
}
write.csv(res1,'results/Bhaduri_putative_stem_cells_GSE131928_10X.csv')


## GSE131928_Smartseq2
# Bhaduri putative stem cells
res1 = matrix(nrow=7,ncol=5, dimnames=list(c('Neural G0','G1','Late G1','S','S/G2','G2/M','M/Early G1'),c('q','m','n','k','pv')))
for(clust1 in c('Neural G0','G1','Late G1','S','S/G2','G2/M','M/Early G1')) {
        q = sum(GSE131928_Smartseq2@meta.data$ccAF==clust1 & (((GSE131928_Smartseq2@assays$RNA@data['PROM1',]>1 | GSE131928_Smartseq2@assays$RNA@data['FUT4',]>1 | GSE131928_Smartseq2@assays$RNA@data['L1CAM',]>1) & GSE131928_Smartseq2@assays$RNA@data['SOX2',]>1) & GSE131928_Smartseq2@assays$RNA@data['TLR4',]==0),na.rm=T)
        m = sum((((GSE131928_Smartseq2@assays$RNA@data['PROM1',]>1 | GSE131928_Smartseq2@assays$RNA@data['FUT4',]>1 | GSE131928_Smartseq2@assays$RNA@data['L1CAM',]>1) & GSE131928_Smartseq2@assays$RNA@data['SOX2',]>1) & GSE131928_Smartseq2@assays$RNA@data['TLR4',]==0),na.rm=T)
        n = dim(GSE131928_Smartseq2)[2]-m
        k = sum(GSE131928_Smartseq2@meta.data$ccAF==clust1,na.rm=T)
        pv = phyper(q, m, n, k, lower.tail=F)
        res1[clust1,] = c(q,m,n,k,pv)
}
write.csv(res1,'results/Bhaduri_putative_stem_cells_GSE131928_Smartseq2.csv')


############################
### GSE102130 DMG H3K27M ###
############################
d1 = read.table(gzfile('GSE102130_K27Mproject.RSEM.vh20170621.txt.gz'),sep='\t',row.names=1,header=T)
dim(d1)
d1 = d1[which(apply(d1,1,sum)!=0),]
dim(d1)
GSE102130 = CreateSeuratObject(counts=d1)
GSE102130[['patient']] = as.factor(sapply(sapply(colnames(d1), function(x) { strsplit(x,"_")[[1]][1] }), function(x) { gsub('X','MGH',toupper(strsplit(x,"\\.")[[1]][1])) }))
GSE102130 = subset(GSE102130, cells = rownames(GSE102130@meta.data)[which(!GSE102130[['patient']][,1] %in% c('MGH101','MGH104','MGH66','OLIGO'))])

GSE102130 = NormalizeData(object = GSE102130)
GSE102130 = FindVariableFeatures(object = GSE102130)
GSE102130 = ScaleData(object = GSE102130)
GSE102130 = RunPCA(object = GSE102130, features = VariableFeatures(GSE102130))
pdf('GSE102130_PCA.pdf')
DimHeatmap(GSE102130, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(GSE102130, dims = 16:30, cells = 500, balanced = TRUE)
dev.off()

GSE102130 = RunUMAP(GSE102130, dims=1:15)
GSE102130 = FindNeighbors(GSE102130, dims=1:15)
GSE102130 = FindClusters(GSE102130, resolution = 0.5)

# Predict our cell_cycle
clusts.WT.GSE102130 = ClassifyCells3(GSE102130,classifier=ccAF,new.data=as.matrix(GSE102130@assays$RNA@data))
clusts.WT.GSE102130 = revalue(clusts.WT.GSE102130, c('G1/Differentiation'='Neural G0'))
names(clusts.WT.GSE102130) = colnames(GSE102130@assays$RNA@data)
GSE102130[['cell_cycle']] = clusts.WT.GSE102130
table(GSE102130[['cell_cycle']])

# Classify using Wang et al., 2017 (Verhaak) classifier
#res1 = runSsGSEAwithPermutation('GSE102130.gct',1000)
res1 = read.table('data/ssGSEA.GBM.classifications/p_result_GSE102130.gct.txt', sep='\t', header=T, row.names=1)
min1 = apply(res1[,4:6],1,min)
calls1 = c('PN','CL','MS')[sapply(rownames(res1), function(x){ which(res1[x,4:6]==min1[x])[[1]] })]
names(calls1) = rownames(res1)
GSE102130[['subtype']] = calls1[colnames(GSE102130)]
table(GSE102130[[c('subtype','cell_cycle')]])
table(GSE102130[[c('subtype','patient')]])

# Run the Seurat cell_cycle classifeir against the data
GSE102130 = CellCycleScoring(GSE102130, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

GSE102130 = FindClusters(GSE102130, resolution = 0.5)
pdf('GSE102130_allCells.pdf', width=10,height=10)
DimPlot(GSE102130, reduction='umap', label=T)
DimPlot(GSE102130, group.by='patient', reduction='umap', label=T)
DimPlot(GSE102130, group.by='cell_cycle', reduction='umap', label=T)
FeaturePlot(GSE102130, features = c('CD3','PTPRC','CD14','CX3CR1','AIF1','MBP','PLP1','GFAP','PRC1'))
dev.off()

# Removed: 6, 12, 13, 14, 15
GSE102130_clean = subset(GSE102130, cells = rownames(GSE102130@meta.data)[which(!GSE102130[['RNA_snn_res.0.5']][,1] %in% c(6,12,13,14,15))])
GSE102130_clean = RunUMAP(GSE102130_clean, dims=1:15)
as.loom(GSE102130_clean, assay='RNA', filename='GSE102130.loom')

pdf('GSE102130_clean.pdf')
DimPlot(GSE102130_clean, reduction='umap', label=T)
DimPlot(GSE102130_clean, group.by='patient', reduction='umap', label=T)
DimPlot(GSE102130_clean, group.by='cell_cycle', reduction='umap', label=T)
DimPlot(GSE102130_clean, group.by='subtype', reduction='umap', label=T)
FeaturePlot(GSE102130_clean, features = c('CD3','PTPRC','CD14','CX3CR1','AIF1','MBP','PLP1','GFAP','PRC1'))
dev.off()


##############################
### GSE84465 GBM => Counts ###
##############################

# Load up count data and start Seurat
d1 = read.table(gzfile('GSE84465_GBM_All_data.csv.gz'),sep=' ',row.names=1,header=T)
dim(d1)
d1 = d1[which(apply(d1,1,sum)!=0),]
dim(d1)
n1 = c(rep('BT_S2',1169),rep('BT_S1',489),rep('BT_S4',1542),rep('BT_S6',389))

# Load phenotypic information
p1 = read.csv('phenotypes_GSE84465.csv',row.names=1,header=T)
rownames(p1) = paste("X",rownames(p1),sep='')

# Subset to tumor cells only
tumor_cells = rownames(p1)[which(p1[,'neoplastic']=='Neoplastic')]
d1 = d1[,tumor_cells]
p1 = p1[tumor_cells,]
GSE84465 = CreateSeuratObject(counts=d1)
GSE84465[["percent.mt"]] = PercentageFeatureSet(GSE84465, pattern = "^MTRNR")

GSE84465[['patient']] = factor(p1[colnames(d1),'patient_ID'])
GSE84465[['tissue']] = factor(p1[colnames(d1),'tissue'])
GSE84465[['cell_type']] = factor(p1[colnames(d1),'cell_type'])
GSE84465[['neoplastic']] = factor(p1[colnames(d1),'neoplastic'])
GSE84465[['cellsLoc']] = factor(paste(p1[colnames(d1),'tissue'],p1[colnames(d1),'cell_type'],sep='_'))

d1 = read.csv('/files/ccAF/classifyPrimaryCells/results/ccAF_results_GSE84465.csv',header=T,row.names=1)

GSE84465[['ccAF']] = d1[,'Predictions',drop=F]
GSE84465[['subtype']] = d1[,'subtype',drop=F]


# Normalize and clean up data
GSE84465 = NormalizeData(object = GSE84465)
GSE84465 = FindVariableFeatures(object = GSE84465)
GSE84465 = ScaleData(object = GSE84465, vars.to.regress=c('nCount_RNA','patient'))
GSE84465 = RunPCA(object = GSE84465, features = VariableFeatures(GSE84465), npcs = 50, weigth.by.var=FALSE)

pdf('GSE84465_plots.pdf')
plotMe = rownames(GSE84465@meta.data)[which(GSE84465@meta.data$cellsLoc=='Tumor_Neoplastic')]
GSE84465 = RunTSNE(GSE84465, dims.use=1:10, do.fast=T, perplexity = 15, cells.use=plotMe)
TSNEPlot(GSE84465, cells.use=plotMe)
dev.off()

# Classifying GSE84465
clusts.WT.GSE84465 = ClassifyCells3(GSE84465,classifier=ccAF,new.data=as.matrix(GSE84465@assays$RNA@data))
names(clusts.WT.GSE84465) = colnames(GSE84465)
GSE84465[['cell_cycle']] = clusts.WT.GSE84465
table(GSE84465[['cell_cycle']])

# Classify using Wang et al., 2017 (Verhaak) classifier
#res1 = runSsGSEAwithPermutation('GSE84465.gct',1000)
res1 = read.table('data/ssGSEA.GBM.classification/p_result_GSE84465.gct.txt', sep='\t', header=T, row.names=1)
min1 = apply(res1[,4:6],1,min)
calls1 = c('PN','CL','MS')[sapply(rownames(res1), function(x){ which(res1[x,4:6]==min1[x])[[1]] })]
names(calls1) = rownames(res1)
GSE84465[['subtype']] = calls1[colnames(GSE84465)]

as.loom(GSE84465, assay='RNA', filename='GSE84465.loom')

# Bhaduri putative stem cells
res1 = matrix(nrow=7,ncol=5, dimnames=list(c('Neural G0','G1','Late G1','S','S/G2','G2/M','M/Early G1'),c('q','m','n','k','pv')))
for(clust1 in c('Neural G0','G1','Late G1','S','S/G2','G2/M','M/Early G1')) {
        q = sum(GSE84465@meta.data$ccAF==clust1 & (((GSE84465@assays$RNA@data['PROM1',]>0 | GSE84465@assays$RNA@data['FUT4',]>0 | GSE84465@assays$RNA@data['L1CAM',]>0) & GSE84465@assays$RNA@data['SOX2',]>0) & GSE84465@assays$RNA@data['TLR4',]==0),na.rm=T)
        m = sum((((GSE84465@assays$RNA@data['PROM1',]>0 | GSE84465@assays$RNA@data['FUT4',]>0 | GSE84465@assays$RNA@data['L1CAM',]>0) & GSE84465@assays$RNA@data['SOX2',]>0) & GSE84465@assays$RNA@data['TLR4',]==0),na.rm=T)
        n = dim(GSE84465)[2]-m
        k = sum(GSE84465@meta.data$ccAF==clust1,na.rm=T)
        pv = phyper(q, m, n, k, lower.tail=F)
        res1[clust1,] = c(q,m,n,k,pv)
}
write.csv(res1, 'results/Bhaduri_putative_stem_cells_GSE84465.csv')


####################################
### GSE84465 All cells => Counts ###
####################################

# Load up count data and start Seurat
setwd('/files/scGlioma/GSE84465')
d1 = read.table(gzfile('GSE84465_GBM_All_data.csv.gz'),sep=' ',row.names=1,header=T)
dim(d1)
d1 = d1[which(apply(d1,1,sum)!=0),]
dim(d1)
n1 = c(rep('BT_S2',1169),rep('BT_S1',489),rep('BT_S4',1542),rep('BT_S6',389))

# Load phenotypic information
p1 = read.csv('phenotypes_GSE84465.csv',row.names=1,header=T)
rownames(p1) = paste("X",rownames(p1),sep='')

# Subset to tumor cells only
#tumor_cells = rownames(p1)[which(p1[,'neoplastic']=='Neoplastic')]
#d1 = d1[,tumor_cells]
#p1 = p1[tumor_cells,]
GSE84465 = CreateSeuratObject(counts=d1)
GSE84465[["percent.mt"]] = PercentageFeatureSet(GSE84465, pattern = "^MTRNR")

GSE84465[['patient']] = factor(p1[colnames(d1),'patient_ID'])
GSE84465[['tissue']] = factor(p1[colnames(d1),'tissue'])
GSE84465[['cell_type']] = factor(p1[colnames(d1),'cell_type'])
GSE84465[['neoplastic']] = factor(p1[colnames(d1),'neoplastic'])
GSE84465[['cellsLoc']] = factor(paste(p1[colnames(d1),'tissue'],p1[colnames(d1),'cell_type'],sep='_'))

# Normalize and clean up data
GSE84465 = NormalizeData(object = GSE84465)
GSE84465 = FindVariableFeatures(object = GSE84465)
GSE84465 = ScaleData(object = GSE84465, vars.to.regress=c('nCount_RNA','patient'))
GSE84465 = RunPCA(object = GSE84465, features = VariableFeatures(GSE84465), npcs = 50, weigth.by.var=FALSE)

as.loom(GSE84465, assay='RNA', filename='GSE84465_all.loom')


###############################
### GSE139448 GBM => Counts ###
###############################

# Load up count data and start Seurat
setwd('/files/scGlioma/GSE139448')
d1 = read.table('GSM4141788_92017_dense.csv',sep=',',row.names=1,header=T)
rownames(d1) = paste('GBM27_',rownames(d1), sep='')
d2 = read.table('GSM4141789_92217_dense.csv',sep=',',row.names=1,header=T)
rownames(d2) = paste('GBM28_',rownames(d2), sep='')
d3 = read.table('GSM4141790_92717_dense.csv',sep=',',row.names=1,header=T)
rownames(d3) = paste('GBM29_',rownames(d3), sep='')
commonGenes = intersect(colnames(d1),intersect(colnames(d2),colnames(d3)))
d0 = t(rbind(d1[,commonGenes],d2[,commonGenes],d3[,commonGenes]))
dim(d0)
d0 = d0[which(apply(d0,1,sum)!=0),]
dim(d0)

GSE139448 = CreateSeuratObject(counts=d0)
GSE139448[['patient']] = c(rep('GBM27',nrow(d1)), rep('GBM28',nrow(d2)), rep('GBM29',nrow(d3)))

GSE139448 = NormalizeData(object = GSE139448, )
GSE139448 = FindVariableFeatures(object = GSE139448)
GSE139448 = ScaleData(object = GSE139448)
GSE139448 = RunPCA(object = GSE139448, features = VariableFeatures(GSE139448), npcs = 50, weigth.by.var=FALSE)

GSE139448 = RunUMAP(GSE139448, dims=1:16)
GSE139448 = FindNeighbors(GSE139448, dims=1:16)
GSE139448 = FindClusters(GSE139448, resolution = 0.5)

GSE139448_clean[['NANOG_G0']] = GSE139448_clean@assays$RNA@data['NANOG',]>0 & GSE139448_clean@meta.data$ccAF=='Neural G0'
nanog_g0 = FindMarkers(GSE139448_clean, ident.1=T, group.by='NANOG_G0')
write.csv(nanog_g0,'GSE139448.markers_NANOG_G0_SeuratV3.0.csv')

# Predict our cell_cycle
clusts.WT.GSE139448 = ClassifyCells3(GSE139448,classifier=ccAF,new.data=as.matrix(GSE139448@assays$RNA@data))
names(clusts.WT.GSE139448) = colnames(GSE139448@assays$RNA@data)
GSE139448[['cell_cycle']] = clusts.WT.GSE139448
table(GSE139448[['cell_cycle']])

# Classify using Wang et al., 2017 (Verhaak) classifier
#write.table(GSE139448@assays$RNA@data, 'GSE139448_forGCT.txt', sep='\t')
res1 = read.table('data/ssGSEA.GBM.classification/p_result_GSE139448.gct.txt', sep='\t', header=T, row.names=1)
min1 = apply(res1[,4:6],1,min)
calls1 = c('PN','CL','MS')[sapply(rownames(res1), function(x){ which(res1[x,4:6]==min1[x])[[1]] })]
names(calls1) = rownames(res1)
GSE139448[['subtype']] = calls1[colnames(GSE139448)]
table(GSE139448[[c('subtype','cell_cycle')]])
table(GSE139448[[c('subtype','patient')]])

# Run the Seurat cell_cycle classifeir against the data
GSE139448 = CellCycleScoring(GSE139448, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

pdf('GSE139448_allCells.pdf')
DimPlot(GSE139448, group.by='seurat_clusters', reduction='umap', label=T)
DimPlot(GSE139448, group.by='patient', reduction='umap', label=T)
DimPlot(GSE139448, group.by='cell_cycle', reduction='umap', label=T)
FeaturePlot(GSE139448, features = c('CD3','PTPRC','CD14','CX3CR1','AIF1','MBP','PLP1','GFAP','PRC1'))
dev.off()

# Did remove clusters 9 and 16 which correpsonded to:
# Immune cells = 9, 16
GSE139448_clean = subset(GSE139448, cells = rownames(GSE139448@meta.data)[which(!GSE139448[['RNA_snn_res.0.5']][,1] %in% c(9,16))])
GSE139448_clean = RunUMAP(GSE139448_clean, dims=1:15)
as.loom(GSE139448_clean, assay='RNA', filename='GSE139448.loom')

pdf('GSE139448_clean.pdf')
DimPlot(GSE139448_clean, group.by='seurat_clusters', reduction='umap', label=T)
DimPlot(GSE139448_clean, group.by='patient', reduction='umap', label=T)
DimPlot(GSE139448_clean, group.by='cell_cycle', reduction='umap', label=T)
DimPlot(GSE139448_clean, group.by='subtype', reduction='umap', label=T)
FeaturePlot(GSE139448_clean, features = c('CD3','PTPRC','CD14','CX3CR1','AIF1','MBP','PLP1','GFAP','PRC1'))
dev.off()

##################
## All Patients ##
##################
#patients = c("GBM27", "GBM28", "GBM29")
patients = c("GBM27")
pcs = list()
pcs$GBM27 = 9
pcs$GBM28 = 7
pcs$GBM29 = 9

res = list()
res$GBM27 = 0.2
res$GBM28 = 0.5
res$GBM29 = 0.5

for(pat1 in patients) {
    #for(sub1 in c('CL','PN','MS')) {
    for(sub1 in c('PN')) {
        patSeurat = subset(GSE139448, cells = rownames(GSE139448@meta.data)[which(GSE139448[['patient']][,1] == pat1)])
        patSeurat = subset(patSeurat, cells = rownames(patSeurat@meta.data)[which(!is.na(patSeurat[['cell_cycle']][,1]))])
        patSeurat = subset(patSeurat, cells = rownames(patSeurat@meta.data)[which(patSeurat[['subtype']][,1]==sub1)])
        
        patSeurat = ScaleData(patSeurat)
        patSeurat = RunPCA(patSeurat, features = VariableFeatures(object = patSeurat))
        
        pdf(paste(pat1,'_',sub1,'_PCA.pdf',sep=''),width = 8.5, height=11)
        DimHeatmap(patSeurat, dims = 1:15, cells = 500, balanced = TRUE)
        DimHeatmap(patSeurat, dims = 16:30, cells = 500, balanced = TRUE)
        dev.off()
        
        patSeurat = FindNeighbors(patSeurat, dims = 1:pcs[[pat1]])
        patSeurat = FindClusters(patSeurat, resolution = 0.1)
        patSeurat = FindClusters(patSeurat, resolution = 0.2)
        patSeurat = FindClusters(patSeurat, resolution = 0.3)
        patSeurat = FindClusters(patSeurat, resolution = 0.35)
        patSeurat = FindClusters(patSeurat, resolution = 0.4)
        patSeurat = FindClusters(patSeurat, resolution = 0.5)
        patSeurat = FindClusters(patSeurat, resolution = 0.6)
        patSeurat = FindClusters(patSeurat, resolution = 0.7)
        patSeurat = FindClusters(patSeurat, resolution = 0.8)
        patSeurat = RunUMAP(patSeurat, dims = 1:pcs[[pat1]])
        
        pdf(paste(pat1,'_',sub1,'_UMAP_res_test.pdf',sep=''),width = 11, height=8)
        a = DimPlot(patSeurat, reduction = "umap", group.by=c('RNA_snn_res.0.1','RNA_snn_res.0.2','RNA_snn_res.0.3','RNA_snn_res.0.4','RNA_snn_res.0.5','RNA_snn_res.0.6','RNA_snn_res.0.7','RNA_snn_res.0.8','cell_cycle','subtype'))
        print(a)
        dev.off()
        
        # Get marker genes
        patSeurat = FindClusters(patSeurat, resolution = res[[pat1]])
        png(paste(pat1,'_',sub1,'_UMAP_res_Final.png',sep=''),width =6, height=2.5, units='in',res=300)
        a = DimPlot(patSeurat, reduction = "umap", group.by=c('cell_cycle',paste('RNA_snn_res.',res[[pat1]],sep='')))+theme(legend.position='none')
        print(a)
        dev.off()

        png(paste(pat1,'_',sub1,'_featurePlot_res_Final.png',sep=''),width =6, height=2.5, units='in',res=300)
        a = FeaturePlot(patSeurat, features = c('MBP','ETNPPL'), reduction = "umap")
        print(a)
        dev.off()
    }
}
write.csv(patSeurat@meta.data,'GBM27_meta_data.csv')

GSE139448@meta.data[,'cell_cycle'] = factor(GSE139448@meta.data[,'cell_cycle'], levels=c('G1/other','Neural G0','G1','Late G1','S','S/G2','G2/M','M/Early G1'))
GSE139448_clean = subset(GSE139448, cells = rownames(GSE139448@meta.data)[which(GSE139448[['patient']][,1] %in% patients)])
GSE139448_clean = subset(GSE139448_clean, cells = rownames(GSE139448_clean@meta.data)[which(!is.na(GSE139448_clean[['cell_cycle']][,1]))])

# Bhaduri putative stem cells
res1 = matrix(nrow=7,ncol=5, dimnames=list(c('Neural G0','G1','Late G1','S','S/G2','G2/M','M/Early G1'),c('q','m','n','k','pv')))
for(clust1 in c('Neural G0','G1','Late G1','S','S/G2','G2/M','M/Early G1')) {
        q = sum(GSE139448_clean@meta.data$ccAF==clust1 & (((GSE139448_clean@assays$RNA@data['PROM1',]>0 | GSE139448_clean@assays$RNA@data['FUT4',]>0 | GSE139448_clean@assays$RNA@data['L1CAM',]>0) & GSE139448_clean@assays$RNA@data['SOX2',]>0) & GSE139448_clean@assays$RNA@data['TLR4',]==0),na.rm=T)
        m = sum((((GSE139448_clean@assays$RNA@data['PROM1',]>0 | GSE139448_clean@assays$RNA@data['FUT4',]>0 | GSE139448_clean@assays$RNA@data['L1CAM',]>0) & GSE139448_clean@assays$RNA@data['SOX2',]>0) & GSE139448_clean@assays$RNA@data['TLR4',]==0),na.rm=T)
        n = dim(GSE139448_clean)[2]-m
        k = sum(GSE139448_clean@meta.data$ccAF==clust1,na.rm=T)
        pv = phyper(q, m, n, k, lower.tail=F)
        res1[clust1,] = c(q,m,n,k,pv)
}
write.csv(res1, 'results/Bhaduri_putative_stem_cells_GSE139448.csv')


##################################
### Bhaduri_2019 GBM => Counts ###
##################################

# Load up count data and start Seurat
d1 = read.table(gzfile('exprMatrix.tsv.gz'),sep='\t',row.names=1,header=T)
dim(d1)
d1 = d1[which(apply(d1,1,sum)!=0),]
dim(d1)

# Make Seurat 3 object
Bhaduri_2019 = CreateSeuratObject(counts=d1)

# Load meta-data
meta = read.table('meta.tsv', header=T, row.names=1, sep='\t')
rownames(meta) = gsub('-','.', rownames(meta))

# Add meta data variables of interest
for(pheno1 in c('Study_ID','Celltype','Tumor.Normal.Classification','CNV.Status','Sample_ID','Highest.Correlated.Cell.Type','Subtype','Cell.Type.Assignment')) {
    Bhaduri_2019[[pheno1]] = meta[,pheno1]
}

# Subset based on Tumor.Normal.Classification
clustsToRemove = c('B Cells','Dividing B Cells','Endothelial','Microglia','Pericyte','Red blood cells','Tumor Associated Macrophage')
Bhaduri_2019 = subset(Bhaduri_2019, cells = rownames(Bhaduri_2019@meta.data)[which(!Bhaduri_2019[['Cell.Type.Assignment']][,1] %in% clustsToRemove)])

# Find variable genes and clusters
Bhaduri_2019 = FindVariableFeatures(object = Bhaduri_2019)
Bhaduri_2019 = ScaleData(object = Bhaduri_2019)
Bhaduri_2019 = RunPCA(object = Bhaduri_2019, features = VariableFeatures(Bhaduri_2019), npcs = 50, weigth.by.var=FALSE)
Bhaduri_2019 = RunUMAP(Bhaduri_2019, dims=1:16)
Bhaduri_2019 = FindNeighbors(Bhaduri_2019, dims=1:16)
Bhaduri_2019 = FindClusters(Bhaduri_2019, resolution = 0.5)

# Predict our cell_cycle
clusts.WT.Bhaduri_2019 = ClassifyCells3(Bhaduri_2019,classifier=ccAF,new.data=as.matrix(Bhaduri_2019@assays$RNA@data))
names(clusts.WT.Bhaduri_2019) = colnames(Bhaduri_2019@assays$RNA@data)
Bhaduri_2019[['cell_cycle']] = clusts.WT.Bhaduri_2019
table(Bhaduri_2019[['cell_cycle']])

# Classify using Wang et al., 2017 (Verhaak) classifier
#res1 = runSsGSEAwithPermutation('Bhaduri_2019.gct',1000)
res1 = read.table('data/ssGSEA.GBM.classification/p_result_Bhaduri_2019.gct.txt', sep='\t', header=T, row.names=1)
res1 = res1[colnames(Bhaduri_2019),]
min1 = apply(res1[,4:6],1,min)
calls1 = c('PN','CL','MS')[sapply(rownames(res1), function(x){ which(res1[x,4:6]==min1[x])[[1]] })]
names(calls1) = rownames(res1)
Bhaduri_2019[['subtype']] = calls1[colnames(Bhaduri_2019)]
table(Bhaduri_2019[[c('subtype','cell_cycle')]])
table(Bhaduri_2019[[c('subtype','Study_ID')]])

# Run the Seurat cell_cycle classifeir against the data
Bhaduri_2019 = CellCycleScoring(Bhaduri_2019, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

pdf('Bhaduri_2019_allCells.pdf')
DimPlot(Bhaduri_2019, group.by='seurat_clusters', reduction='umap', label=T)
DimPlot(Bhaduri_2019, group.by='Study_ID', reduction='umap', label=T)
DimPlot(Bhaduri_2019, group.by='cell_cycle', reduction='umap', label=T)
FeaturePlot(Bhaduri_2019, features = c('CD3','PTPRC','CD14','CX3CR1','AIF1','MBP','PLP1','GFAP','PRC1'))
dev.off()

# Removed:
# Immune cells = 13
# Oligodendrocyte = 12
Bhaduri_2019_clean = subset(Bhaduri_2019, cells = rownames(Bhaduri_2019@meta.data)[which(!Bhaduri_2019[['RNA_snn_res.0.5']][,1] %in% c(12,13))])
Bhaduri_2019_clean = RunUMAP(Bhaduri_2019_clean, dims=1:15)
as.loom(Bhaduri_2019_clean, assay='RNA', filename='Bhaduri.loom')

pdf('Bhaduri_2019_clean.pdf')
DimPlot(Bhaduri_2019_clean, group.by='seurat_clusters', reduction='umap', label=T)
DimPlot(Bhaduri_2019_clean, group.by='patient', reduction='umap', label=T)
DimPlot(Bhaduri_2019_clean, group.by='cell_cycle', reduction='umap', label=T)
DimPlot(Bhaduri_2019_clean, group.by='subtype', reduction='umap', label=T)
FeaturePlot(Bhaduri_2019_clean, features = c('CD3','PTPRC','CD14','CX3CR1','AIF1','MBP','PLP1','GFAP','PRC1'))
dev.off()

d1 = read.csv('/files/ccAF/classifyPrimaryCells/results/ccAF_results_Bhaduri.csv',header=T,row.names=1)
Bhaduri_2019_clean[['ccAF']] = d1[,'Predictions',drop=F]
Bhaduri_2019_clean[['subtype']] = d1[,'subtype',drop=F]

# Bhaduri putative stem cells
res1 = matrix(nrow=7,ncol=5, dimnames=list(c('Neural G0','G1','Late G1','S','S/G2','G2/M','M/Early G1'),c('q','m','n','k','pv')))
for(clust1 in c('Neural G0','G1','Late G1','S','S/G2','G2/M','M/Early G1')) {
        q = sum(Bhaduri_2019_clean@meta.data$ccAF==clust1 & (((Bhaduri_2019_clean@assays$RNA@data['PROM1',]>0 | Bhaduri_2019_clean@assays$RNA@data['FUT4',]>0 | Bhaduri_2019_clean@assays$RNA@data['L1CAM',]>0) & Bhaduri_2019_clean@assays$RNA@data['SOX2',]>0) & Bhaduri_2019_clean@assays$RNA@data['TLR4',]==0),na.rm=T)
        m = sum((((Bhaduri_2019_clean@assays$RNA@data['PROM1',]>0 | Bhaduri_2019_clean@assays$RNA@data['FUT4',]>0 | Bhaduri_2019_clean@assays$RNA@data['L1CAM',]>0) & Bhaduri_2019_clean@assays$RNA@data['SOX2',]>0) & Bhaduri_2019_clean@assays$RNA@data['TLR4',]==0),na.rm=T)
        n = dim(Bhaduri_2019_clean)[2]-m
        k = sum(Bhaduri_2019_clean@meta.data$ccAF==clust1,na.rm=T)
        pv = phyper(q, m, n, k, lower.tail=F)
        res1[clust1,] = c(q,m,n,k,pv)
}
write.csv(res1, 'results/Bhaduri_putative_stem_cells_Bhaduri.csv')


#################################
### GSE103322 HNSCC => Counts ###
#################################

# Load up count data and start Seurat
setwd('/files/scGlioma/GSE103322')
d1 = read.table(gzfile('GSE103322_HNSCC_all_data_noPhenos.txt.gz'),sep='\t',row.names=1,header=T)
dim(d1)
d1 = d1[which(apply(d1,1,sum)!=0),]
dim(d1)
# Load phenotypic information
p1 = read.csv('pheno1.csv',row.names=1,header=T)
# Keep only cancer cells
d1 = d1[,names(p1)[which(p1==1)]]
GSE103322 = CreateSeuratObject(counts=d1)
GSE103322 = FindVariableFeatures(object = GSE103322)

# Classifying GSE103322
clusts.WT.GSE103322 = ClassifyCells3(GSE103322,classifier=ccAF,new.data=as.matrix(GSE103322@assays$RNA@data))
names(clusts.WT.GSE103322) = colnames(GSE103322)
GSE103322[['cell_cycle']] = clusts.WT.GSE103322

# Write out loom file
as.loom(GSE103322, assay='RNA', filename='GSE103322.loom')


#############################################
### HEK293T cells from barnyard assay 10X ###
#############################################

setwd('/files/scGlioma/hgmm')
exp_mat = Read10X('filtered_feature_bc_matrix')

# Subset to only human genes
exp_mat = exp_mat[grep('hg19_',rownames(exp_mat)),] # Only human genes
rownames(exp_mat) = gsub('hg19_', '', rownames(exp_mat)) # Rename to remove hg19_

# Make Seruat object
HEK293T = CreateSeuratObject(counts=exp_mat)
HEK293T[["percent.mt"]] = PercentageFeatureSet(HEK293T, pattern = "^MT-")
pdf('HEK293T_vlnPlot.pdf', width=11, height=8.5)
VlnPlot(HEK293T, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 = FeatureScatter(HEK293T, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(HEK293T, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

# Subset out good HEK293T cells  (ncol = 3468, nrow = 57,905)
HEK293T = subset(HEK293T, subset = nFeature_RNA > 4500 & nFeature_RNA < 100000 & percent.mt < 30)
dim(HEK293T)

# Find variable features, scale data, and do dimensionality reduction
HEK293T = NormalizeData(object = HEK293T)
HEK293T = FindVariableFeatures(object = HEK293T)
HEK293T = ScaleData(object = HEK293T, vars.to.regress=c('percent.mt','nUMI'))
HEK293T = RunPCA(object = HEK293T, features = VariableFeatures(HEK293T))

# Determined that 14 PCs looked reasonable
pdf('HEK293T_PCA.pdf')
DimHeatmap(HEK293T, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(HEK293T, dims = 16:30, cells = 500, balanced = TRUE)
dev.off()

# Predict our cell_cycle
clusts.WT.HEK293T = ClassifyCells3(HEK293T,classifier=ccAF,new.data=as.matrix(HEK293T@assays$RNA@data))
clusts.WT.HEK293T = revalue(clusts.WT.HEK293T, c('G1/Differentiation'='Neural G0'))
names(clusts.WT.HEK293T) = colnames(HEK293T@assays$RNA@data)
HEK293T[['cell_cycle']] = clusts.WT.HEK293T
table(HEK293T[['cell_cycle']])

# Run the Seurat cell_cycle classifeir against the data
HEK293T = CellCycleScoring(HEK293T, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
table(HEK293T[[c('Phase','cell_cycle')]])

as.loom(HEK293T, assay='RNA', filename='HEK293.loom')


###########################
### Nowakowski_2017 GBM ###
###########################

# Load up count data and start Seurat
setwd('/files/scGlioma/Nowakowski')
d1 = read.table(gzfile('exprMatrix.tsv.gz'),sep='\t',row.names=1,header=T)
dim(d1)
d1 = d1[which(apply(d1,1,sum)!=0),]
dim(d1)

# Make Seurat 3 object
Nowakowski_2017 = CreateSeuratObject(counts=d1)

# Load meta-data
meta = read.table('meta.tsv', header=T, row.names=1, sep='\t', as.is=T)
meta[is.na(meta)] = 'none'

# Make a cell type conversion list
conversion = list()
conversion[['Astrocyte']] = 'Astrocyte'
conversion[['Choroid']] = 'Choroid'
conversion[['Endothelial']] = 'Endothelial'
conversion[['EN-PFC1']] = 'EN-PFC'
conversion[['EN-PFC2']] = 'EN-PFC'
conversion[['EN-PFC3']] = 'EN-PFC'
conversion[['EN-V1-1']] = 'EN-V1'
conversion[['EN-V1-2']] = 'EN-V1'
conversion[['EN-V1-3']] = 'EN-V1'
conversion[['Glyc']] = 'Glyc'
conversion[['IN-CTX-CGE1']] = 'IN-CTX-CGE'
conversion[['IN-CTX-CGE2']] = 'IN-CTX-CGE'
conversion[['IN-CTX-MGE1']] = 'IN-CTX-MGE'
conversion[['IN-CTX-MGE2']] = 'IN-CTX-MGE'
conversion[['IN-STR']] = 'IN-STR'
conversion[['IPC-div1']] = 'IPC-div'
conversion[['IPC-div2']] = 'IPC-div'
conversion[['IPC-nEN1']] = 'IPC-nEN'
conversion[['IPC-nEN2']] = 'IPC-nEN'
conversion[['IPC-nEN3']] = 'IPC-nEN'
conversion[['MGE-div']] = 'MGE-div'
conversion[['MGE-IPC1']] = 'MGE-IPC'
conversion[['MGE-IPC2']] = 'MGE-IPC'
conversion[['MGE-IPC3']] = 'MGE-IPC'
conversion[['MGE-RG1']] = 'MGE-RG'
conversion[['MGE-RG2']] = 'MGE-RG'
conversion[['Microglia']] = 'Microglia'
conversion[['Mural']] = 'Mural'
conversion[['none']] = 'none'
conversion[['nEN-early1']] = 'nEN-early'
conversion[['nEN-early2']] = 'nEN-early'
conversion[['nEN-late']] = 'nEN-late'
conversion[['nIN1']] = 'nIN'
conversion[['nIN2']] = 'nIN'
conversion[['nIN3']] = 'nIN'
conversion[['nIN4']] = 'nIN'
conversion[['nIN5']] = 'nIN'
conversion[['OPC']] = 'OPC'
conversion[['oRG']] = 'oRG'
conversion[['RG-div1']] = 'RG-div'
conversion[['RG-div2']] = 'RG-div'
conversion[['RG-early']] = 'RG-early'
conversion[['tRG']] = 'tRG'
conversion[['U1']] = 'U'
conversion[['U2']] = 'U'
conversion[['U3']] = 'U'
conversion[['U4']] = 'U'
conversion[['vRG']] = 'vRG'
meta[,'WGCNAcluster_restricted'] = unlist(sapply(meta[,'WGCNAcluster'], function(x) { ifelse(is.na(x), 'NA', conversion[[x]]) }))

# Add meta data variables of interest
for(pheno1 in colnames(meta)) {
    Nowakowski_2017[[pheno1]] = meta[,pheno1]
}

# Find variable genes and clusters
Nowakowski_2017 = NormalizeData(object = Nowakowski_2017)
Nowakowski_2017 = FindVariableFeatures(object = Nowakowski_2017)

# Predict our cell_cycle
clusts.WT.Nowakowski_2017 = ClassifyCells3(Nowakowski_2017,classifier=ccAF,new.data=as.matrix(Nowakowski_2017@assays$RNA@data))
#clusts.WT.Nowakowski_2017 = revalue(clusts.WT.Nowakowski_2017, c('G1/Differentiation'='Neural G0'))
names(clusts.WT.Nowakowski_2017) = colnames(Nowakowski_2017@assays$RNA@data)
Nowakowski_2017[['cell_cycle']] = clusts.WT.Nowakowski_2017
table(Nowakowski_2017[['cell_cycle']])

# Run the Seurat cell_cycle classifeir against the data
Nowakowski_2017 = CellCycleScoring(Nowakowski_2017, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

# Write out loom file

as.loom(Nowakowski_2017, assay='RNA', filename='Nowakowski_norm.loom')


#################################
## GSE67833: Llorens-Bobadilla ##
#################################

setwd('/files/scGlioma/GSE67833')
exp_mat = read.csv('GSE67833_Gene_expression_matrix.csv',header=T,row.names=1)

GSE67833 = CreateSeuratObject(counts = exp_mat)
dim(GSE67833)
GSE67833 = NormalizeData(object = GSE67833)
GSE67833 = FindVariableFeatures(object = GSE67833, selection.method='vst', nfeatures=2000)
GSE67833 = ScaleData(GSE67833)
GSE67833 = RunPCA(GSE67833, npcs = 30)

meta_data = read.csv('meta_data.csv', header=T, row.names=2)
GSE67833[['cell.type']] = meta_data[colnames(exp_mat), 2]

GSE67833 = RunUMAP(GSE67833, reduction = 'pca', dims = 1:30)
pdf('GSE67833_UMAP_cell.type.pdf')
DimPlot(GSE67833, reduction = 'umap', group.by = 'cell.type')
dev.off()

# For classification
GSE67833.NI = subset(GSE67833, cells = colnames(GSE67833)[-grep('Injured',GSE67833[['cell.type']][,1])])
GSE67833.NI = subset(GSE67833, cells = colnames(GSE67833.NI)[!is.na(GSE67833.NI[['cell.type']])])

as.loom(GSE67833.NI, assay='RNA', filename='GSE67833.loom')


###################
### PRJNA324289 ###
###################
conv1 = read.csv('data/geneConversions/human_hgnc_mouse_mgi.csv',header=1)
# Make 1:1 mouse to human conversion
conv2 = conv1[which(conv1[,'MGI.symbol'] %in% names(which(table(as.character(conv1[,'MGI.symbol']))==1))),]
conv2 = conv2[which(conv2[,'Gene.stable.ID'] %in% names(which(table(as.character(conv2[,'Gene.stable.ID']))==1))),]
d1 = read.csv('Counts_AllLiveSingleCells_IN_VIVO_ONLY.csv',header=T,row.names=1)
d2 = d1[which(rownames(d1) %in% conv2[,'MGI.symbol']),]
annotLookup = as.character(conv2[which(conv2[,'MGI.symbol'] %in% rownames(d2)),'Gene.stable.ID.Mus'])
names(annotLookup) = as.character(conv2[which(conv2[,'MGI.symbol'] %in% rownames(d2)),'MGI.symbol'])
rownames(d2) = annotLookup[rownames(d2)]
ct1 = sapply(colnames(d1), function(x) { strsplit(x,'_')[[1]][1] })
PRJNA324289 = CreateSeuratObject(counts = d2)
PRJNA324289[['cell_type']] = ct1
pdf('qcPlots_PRJNA324289.pdf')
VlnPlot(object = PRJNA324289, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)
FeatureScatter(PRJNA324289, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
dev.off()
dim(PRJNA324289@data)
PRJNA324289 = NormalizeData(object = PRJNA324289)
PRJNA324289 = FindVariableFeatures(object = PRJNA324289, selection.method='vst', nfeatures=2000)
PRJNA324289 = ScaleData(PRJNA324289)
PRJNA324289 = RunPCA(PRJNA324289, npcs = 30)
# Choose PCs to use moving forward

PRJNA324289 = RunUMAP(PRJNA324289, reduction = 'pca', dims = 1:30)
pdf('PRJNA324289_UMAP_cell.type.pdf')
DimPlot(PRJNA324289, reduction = 'umap', group.by = 'cell_type')
dev.off()

as.loom(PRJNA324289, assay='RNA', filename='PRJNA324289.loom')

