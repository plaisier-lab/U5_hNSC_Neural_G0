###########################################################
## ccAF:  U5-hNSC analysis to discover cell cycle phases ##
##  ______     ______     __  __                         ##
## /\  __ \   /\  ___\   /\ \/\ \                        ##
## \ \  __ \  \ \___  \  \ \ \_\ \                       ##
##  \ \_\ \_\  \/\_____\  \ \_____\                      ##
##   \/_/\/_/   \/_____/   \/_____/                      ##
## @Developed by: Plaisier Lab                           ##
##   (https://plaisierlab.engineering.asu.edu/)          ##
##   Arizona State University                            ##
##   242 ISTB1, 550 E Orange St                          ##
##   Tempe, AZ  85281                                    ##
## @Author:  Chris Plaisier, Samantha O'Connor           ##
## @License:  GNU GPLv3                                  ##
##                                                       ##
## If this program is used in your analysis please       ##
## mention who built it. Thanks. :-)                     ##
###########################################################

# To install Seurat V2.3.4 which is required for this script
# source("https://z.umn.edu/archived-seurat")

# Preparing for analysis
library(Seurat)
library(dplyr)
genome = 'hg19'

setwd('data/U5_hNSC')

##########
### WT ###
##########
# Setup Seurat object
exp_mat = Read10X(data.dir='WT/filtered_gene_bc_matrices/hg19')
WT = CreateSeuratObject(raw.data=exp_mat)
WT@meta.data$pert = "WT"
mito.genes = grep(pattern = "^MT-", x = rownames(x = WT@data), value = TRUE)
percent.mito = colSums(WT@data[mito.genes, ]) / colSums(WT@data)
WT = AddMetaData(object = WT, metadata = percent.mito, col.name = "percent.mito")
WT = NormalizeData(object = WT)
WT = FindVariableGenes(object = WT)
WT = ScaleData(object = WT, vars.to.regress = c('percent.mito','nUMI'), genes.use = WT@var.genes, model.use = 'negbinom')

###############
### WT CDT+ ###
###############
# Setup Seurat object
exp_mat = Read10X(data.dir='WT_CDTplus/filtered_gene_bc_matrices/hg19')
WT_CDTplus = CreateSeuratObject(raw.data=exp_mat)
WT_CDTplus@meta.data$pert = "WT_CDTplus"
mito.genes = grep(pattern = "^MT-", x = rownames(x = WT_CDTplus@data), value = TRUE)
percent.mito = colSums(WT_CDTplus@data[mito.genes, ]) / colSums(WT_CDTplus@data)
WT_CDTplus = AddMetaData(object = WT_CDTplus, metadata = percent.mito, col.name = "percent.mito")
WT_CDTplus = NormalizeData(object = WT_CDTplus)
WT_CDTplus = FindVariableGenes(object = WT_CDTplus)
WT_CDTplus = ScaleData(object = WT_CDTplus, vars.to.regress = c('percent.mito','nUMI'), genes.use = WT_CDTplus@var.genes, model.use = 'negbinom')

###############
### sgTAOK1 ###
###############
exp_mat = Read10X(data.dir='sGTAOK1/filtered_gene_bc_matrices/hg19')
taok1 = CreateSeuratObject(raw.data=exp_mat)
taok1@meta.data$pert = "sgTAOK1"
mito.genes = grep(pattern = "^MT-", x = rownames(x = taok1@data), value = TRUE
percent.mito = colSums(taok1@data[mito.genes, ]) / colSums(taok1@data)
taok1 = AddMetaData(object = taok1, metadata = percent.mito, col.name = "percent.mito")
taok1 = NormalizeData(object = taok1)
taok1 = FindVariableGenes(object = taok1)
taok1 = ScaleData(object = taok1, vars.to.regress = c('percent.mito','nUMI'), genes.use = taok1@var.genes, model.use = 'negbinom')

# Overlap between gene lists
g.1 = head(rownames(WT@hvg.info), 1000)
g.2 = head(rownames(taok1@hvg.info), 1000)
genes.use = unique(c(g.1, g.2))
genes.use = intersect(genes.use, rownames(WT@scale.data))
genes.use = intersect(genes.use, rownames(taok1@scale.data))

# Align datasets using CCA
cellcycle.combined = RunCCA(WT, taok1, genes.use = genes.use, num.cc = 30, add.cell.id1='WT', add.cell.id2='sgTAOK1')
cellcycle.combined = AlignSubspace(cellcycle.combined, reduction.type = "cca", grouping.var = "pert", dims.align = 1:20)

# Clustering
cellcycle.combined = FindClusters(cellcycle.combined, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:20)
cellcycle.combined = SetAllIdent(object = cellcycle.combined, id='res.0.6')
current.cluster.ids = c(0,1,2,3,4,5,6,7)
new.cluster.ids = c('G1','G1/Differentiation','G2/M','S/G2','M/Early G1','S','Late G1','G1/other')
tmp = plyr::mapvalues(cellcycle.combined@ident, from=current.cluster.ids, to=new.cluster.ids)
cellcycle.combined = AddMetaData(object = cellcycle.combined, metadata = tmp, col.name = "clusts_named")
cellcycle.combined = SetAllIdent(object = cellcycle.combined, id='clusts_named')

# Make tSNE embeddings
cellcycle.combined = RunTSNE(cellcycle.combined, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T, perplexity=26)
pdf(paste('results/TSNE_perpelexity_26.pdf',sep=''))
TSNEPlot(cellcycle.combined, cells.use=rownames(cellcycle.combined@meta.data)[cellcycle.combined@meta.data[,'pert']=='WT'], do.label=T, no.legend=T)
dev.off()
write.csv(cellcycle.combined@dr$tsne@cell.embeddings, 'results/tsne_cell_embeddings_Perplexity_26.csv')

# Identify Markers
cellcycle.markers = FindAllMarkers(cellcycle.combined)
write.csv(cellcycle.markers,'results/eightClusters_WT_sgTAOK1.csv')


######################
## Network analysis ##
######################
# Make medioids for clusters for genes.use
library(iGraph)
tmp = sapply(c('G1','G1/Differentiation','G2/M','S/G2','M/Early G1','S','Late G1','G1/other'), function(x) { apply(cellcycle.combined@data[unique(cellcycle.markers[,'gene']),intersect(which(cellcycle.combined@meta.data[,'pert']=='WT'),which(cellcycle.combined@meta.data[,'clusts_named']==x))], 1, mean) })
tmp2 = tmp[,-8]
#tmp2 = tmp
dist1 = as.matrix(dist(t(tmp2), method='canberra'))
diag(dist1) = NA
g = graph.adjacency(
  dist1,
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)
rowSums(dist1<=240,na.rm=T)

g = delete_edges(g, E(g)[which(E(g)$weight>240)])
edgeweights = (1/(scale(E(g)$weight,center=F)^6))*3
pdf('results/cellCycleNetwork.pdf')
plot(g, layout=layout.fruchterman.reingold, edge.width=edgeweights)
dev.off()

#######################
## Enrichment of YAP ##
#######################
induced_by_YAP_Zhang = c('SV2B', 'ALPP', 'SPARCL1', 'DMC1', 'KRT34', 'TM4SF20', 'ZNF385B', 'ITGB2', 'OLR1', 'ATRNL1', 'KBTBD10', 'ARHGAP29', 'FBLN5', 'THBS1', 'RNF165', 'CLIC5', 'GRIN2A', 'SIDT1', 'PRSS7', 'STXBP6', 'ARGLU1', 'ANKRD1', 'PCTK3', 'CBLB', 'CTGF', 'KCNT2', 'TOX2', 'FANCB', 'CENPM', 'ADAMTS1', 'MIER3', 'TUBE1', 'FANCD2', 'PLCB4', 'INHBA', 'SFTA1P', 'SERTAD4', 'LOC100129406', 'MYOCD', 'DENND5B', 'MATN3', 'KNG1', 'TRIO', 'CCDC18', 'NRAP', 'ORC1L', 'ZBED2', 'DAB2', 'BCAR4', 'GJA5', 'HSD3B1', 'FRY', 'TGM2', 'TMEM171', 'LOC199899', 'NT5DC4', 'SRGN', 'WTIP', 'SPATA22', 'IGFBP3', 'PRR16', 'NBEAL1', 'STARD8', 'CLIC5', 'LOC100131610', 'CNOT2', 'CLIC3', 'DRAM', 'SPIN3', 'PHLDB2', 'LRRN1', 'ERCC6L', 'KIAA2018', 'WDR69', 'BCHE', 'C16orf88', 'TBXA2R', 'SERPINE2', 'DNAJC6', 'LOC100128329', 'E2F1', 'HJURP', 'ENTPD4', 'LOXL2', 'IQCF1', 'SLIT2', 'SCARA3', 'ZCCHC24', 'ARHGAP28', 'ZNF546', 'LOC728987', 'HEG1', 'SFRS11', 'SNHG4', 'VLDLR', 'ADAMTSL4', 'OSGIN2', 'MYBL2', 'TNS1', 'COL8A1', 'MYCL1', 'TMEM56', 'KLRC4', 'TNNT2', 'TMEM119', 'GLDC', 'FMN2', 'DEPDC1B', 'BDNF', 'ZWINT', 'ZNF169', 'TSPAN15', 'CDC20', 'RIMS2', 'RAB5A', 'PTGS2', 'ZDHHC8P', 'ZNF335', 'RNF34', 'SPARC', 'HIVEP3', 'PKMYT1', 'FBN2', 'NUF2', 'COL4A4', 'PIK3C2B', 'MTIF3', 'FOXM1', 'PMP22', 'PLEC1', 'FGF9', 'EZR', 'CCNA2', 'ASF1B', 'AREG', 'AREGB', 'ANLN', 'RAD51', 'SHCBP1', 'C18orf24', 'RTKN2', 'SERPINE1', 'BDP1', 'TOP2A', 'UBASH3B', 'PLK4', 'COL4A1', 'SLC44A2', 'POLQ', 'RAD51AP1', 'OR1E1', 'TACC3', 'ST6GAL2', 'PRSS23', 'C9orf25', 'GINS2', 'ATP13A3', 'MCM10', 'GNRHR', 'HIC2', 'CDC6', 'DKFZp686L14188', 'EXPH5', 'LOC284513', 'SERPINH1', 'SYTL4', 'ERG', 'APOL2', 'NUDT10', 'DNAH11', 'RRM2', 'TBC1D24', 'IGHG1', 'STMN1', 'DDAH1', 'RHOF', 'ATP6V0A4', 'LOC399959', 'TCEAL2', 'SKP2', 'NCAPH', 'CCNB1', 'TYMS', 'KIF14', 'CENPA', 'MCM5', 'FGFR1', 'PRC1', 'COL4A3', 'PDLIM2', 'TMEM158', 'FSCN1', 'NBR1', 'GINS1', 'GPR98', 'TRPC3', 'KRT80', 'TMEM97', 'DIAPH3', 'FYN', 'NAV3', 'TSPAN2', 'HDAC7', 'FSTL1', 'KIF4A', 'PLK1', 'CCDC150', 'KIF18B', 'CASC5', 'CHST7', 'C2orf59', 'LOC541471', 'E2F7', 'LOC643650', 'POLA2', 'MKI67', 'DLGAP5', 'CCNB2', 'CDC25C', 'KIF23', 'TAGLN', 'ORC6L', 'DUSP6', 'SLC7A1', 'FAM111B', 'CENPF', 'DLC1', 'REEP1', 'NDC80', 'CDH4', 'CLDN19', 'VWA3B', 'RASSF7', 'PPM1H', 'WFS1', 'FAM158A', 'LAMA4', 'ANKRD28', 'SPTA1', 'UHRF1', 'PEG10', 'MYBL1', 'MCAM', 'LHFP', 'TSPAN32', 'PCDH7', 'DNM3', 'CNN3', 'DHRS13', 'L1CAM', 'RHOT2', 'CHORDC1', 'CAPN6', 'FLJ37035', 'LPAR1', 'NCAPG', 'KRT4', 'ZNF542', 'DTYMK', 'LOC727761', 'KIAA0101', 'AKAP12', 'RECQL4', 'SMG7', 'NPHP1', 'NBPF5', 'NBPF6', 'LRCH3', 'COL12A1', 'WNT8A', 'DEPDC1', 'FAM80B', 'ZNF367', 'TUFT1', 'CDCA7', 'PSG5', 'AXL', 'TMEM139', 'LRP8', 'PRAGMIN', 'COL4A2', 'CEP55', 'BIRC5', 'NCAPG2', 'TMEFF2', 'ETV5', 'LETM2', 'SFRS5', 'BRCA1', 'DCN', 'BANP', 'EXO1', 'LTBP2', 'ACSL6', 'CYR61', 'BRIP1', 'TK1', 'GLS', 'WHSC1', 'LOC100127998', 'SNTG1', 'FAM54A', 'QSOX1', 'BRCA2', 'EVI2B', 'KCNK1', 'CDC42EP3', 'DACH2', 'ALPK2', 'TM7SF2', 'NEK2', 'HMMR', 'MELK', 'VCAN', 'SLC25A29', 'MFAP5', 'SCHIP1', 'PXDN', 'KCTD15', 'CDCA2', 'DTL', 'CEP152', 'CDC2', 'C8orf55', 'NR3C2', 'TINAGL1', 'ZNF608', 'PARVA', 'AFAP1L1', 'CFLAR', 'CGN', 'BLMH', 'RFTN1', 'ALPP', 'ALPPL2', 'TMEM54', 'SORCS1', 'TMEM120B', 'WSB1', 'CDKN3', 'WDR51A', 'FRAS1', 'NT5DC3', 'CDC14B', 'FAM64A', 'GAS6', 'LOC100133684', 'KDELR3', 'STX11', 'C9orf140', 'FAM72A', 'FAM72B', 'GCUD2', 'KIF20A', 'TELO2', 'FAM107B', 'C13orf3', 'C6orf173', 'MCM6', 'RNASEH2A', 'IFNA7', 'BCL2L1', 'LOC152586', 'HMGB2', 'NALCN', 'CLCN5', 'TMEFF1', 'THSD4', 'ACAT2', 'TRIM15', 'PHF19', 'CCDC99', 'PHGDH', 'TNNT3', 'TBC1D12', 'RASAL2', 'C21orf23', 'ZNF273', 'TCTE3', 'MCM2', 'SPAG5', 'CPNE2', 'NUP210', 'C4orf27', 'C14orf65', 'C1orf118', 'LOC440338', 'SGOL2', 'PDK2', 'TMEM201', 'LOC100133229', 'TCF7L1', 'AREGB', 'KLHDC1', 'RABEPK', 'CLEC4M', 'SPRED1', 'CR1', 'PTGER4', 'HES2', 'SLC25A15', 'TGFB2', 'TPX2', 'DOT1L', 'SH2D5', 'FGF2', 'QSER1', 'TXNDC11', 'NKX3-1', 'SIM2', 'FZD1', 'RFC3', 'OSBP2', 'SYDE2', 'KIF2C', 'BUB1B', 'PBK', 'SAMD4A', 'LOC647190', 'LOC202051', 'HMGA1', 'LOC283177', 'JAK2', 'MLF1IP', 'C15orf50', 'LOC100132077', 'ESRRG', 'SGK1', 'STRA6', 'ESPL1', 'NEGR1', 'PLA2G4A', 'BUB1', 'AMIGO3', 'GMPPB', 'FZD7', 'NDRG1', 'TMEM170B', 'FAM110B', 'SFRS10', 'SDK1', 'HLA-DRB1', 'EGR1', 'GLIS2', 'SPC24', 'EMP2', 'LRRC3B', 'LOC91431', 'AURKA', 'SPC25', 'ITPR1', 'F3', 'OCC-1', 'PCOLCE2', 'TRMT5', 'BTBD14A', 'CKAP2L', 'ASPM', 'PSMD9', 'CDK2', 'RP1-27O5.1', 'PSG1', 'CDCA3', 'GPR107', 'C17orf69', 'PIN1', 'KRT18', 'RHOBTB3', 'DIS3L2', 'WTAP', 'HIST1H4A', 'HIST1H4B', 'HIST1H4C', 'HIST1H4D', 'HIST1H4E', 'HIST1H4F', 'HIST1H4H', 'HIST1H4I', 'HIST1H4J', 'HIST1H4K', 'HIST1H4L', 'HIST2H4A', 'HIST2H4B', 'HIST4H4', 'FADS1', 'FADS3', 'FBXO5', 'MSRB3', 'MCTP2', 'LTBP1', 'MYADM', 'SH2D4A', 'ATP8B1', 'FANCI', 'RHOBTB1', 'TTK', 'UBE2C', 'CPVL', 'SRrp35', 'CEP78', 'MCPH1', 'AURKB', 'OIP5', 'TMPO', 'ABCA8', 'ANKRD2', 'DGKZ', 'CITED4', 'LIMCH1', 'WWC2', 'C9orf64', 'MGAT4A', 'MOBKL2B', 'PSG1', 'PSG3', 'SVEP1', 'LOC646870', 'PTPN21', 'LOC728705', 'KLK2', 'FBN1', 'LOC731049', 'UBE2S', 'CHST3', 'H2AFX', 'CENPN', 'GPR126', 'PNMA2', 'NEIL3', 'MARCKS', 'CDC7', 'SFRP1', 'PMEPA1', 'UGCGL2', 'KIAA0802', 'WDHD1', 'PLAC8', 'SPRY4', 'TIMM10', 'ST3GAL3', 'C6orf105', 'LOC643287', 'HABP4', 'SEC24A', 'SLC2A3', 'CCL20', 'BARD1', 'DHCR7', 'EGFLAM', 'S1PR1', 'DUSP4', 'OR7E24', 'SHB', 'PDXK', 'KCTD14', 'SMARCA1', 'ACTR2', 'SLC25A25', 'PRR6', 'PLEKHG2', 'GCNT1', 'NF2', 'WWC1', 'PFTK1', 'ARSJ', 'LMNB1', 'TMEM194A', 'SH3TC1', 'NCAPD2', 'RTN4R', 'DYRK2', 'KPNA2', 'LOC728860', 'LOC284952', 'GLI2', 'FLJ42709', 'LY6E', 'ODC1', 'MCM7', 'CDK6', 'C1orf2', 'CDH2', 'SLC35B4', 'WNT5A', 'BLM', 'FLJ14213', 'NUAK1', 'EPB41L2', 'CCR6', 'ELAVL1', 'C10orf35', 'C10orf93', 'DHX37', 'ZNF783', 'SOX13', 'KAP2.1B', 'KRTAP2-1', 'KRTAP2-4', 'LOC644350', 'LOC728285', 'LOC728934', 'LOC730755', 'ASTN2', 'SEC22C', 'RAD51L3', 'GJC1', 'KIF15', 'SLC10A2', 'LRRC14', 'CDCA8', 'DSCC1', 'RNF207', 'MRPL12', 'SLC1A4', 'HELLS', 'BTC', 'TCF19', 'SLC7A11', 'MAD2L1', 'AKAP10', 'CNPY4', 'MVD', 'WDR4', 'SLC37A4', 'PARVB', 'FAM83D', 'POLD1', 'MOSPD2', 'FAM3C', 'TTC9', 'SMC4', 'PSG9', 'PASK', 'COTL1', 'ZSCAN1', 'ACTN1', 'TUBB3', 'CENPK', 'NUSAP1', 'LAT2', 'PCDH24', 'MDFIC', 'GABRR1', 'CLU', 'EXTL3', 'PRIM2', 'AKAP6', 'SLC2A14', 'SLC2A3', 'TDRKH', 'C1orf106', 'TPM2', 'MGLL', 'PTPRJ', 'RRBP1', 'RAD54B', 'HMGB3', 'KHDRBS3', 'IKZF4', 'FZD5', 'RFC4', 'EHD1', 'LRRC8C', 'AMOTL2', 'ARSD', 'CRY1', 'NASP', 'C18orf25', 'B4GALT4', 'FLNA', 'CYP2U1', 'ATG7', 'NR2F1', 'NDP', 'METTL11A', 'CHST14', 'HSPG2', 'CSRP1', 'C19orf42', 'AP4E1', 'SLC4A4', 'RNF150', 'ZEB1', 'SUV39H1', 'CBX5', 'FLJ16734', 'GNG11', 'SAMD1', 'ZNF398', 'ATAD2', 'ODZ3', 'ARL17P1', 'TEAD1', 'CCNF', 'IDH2', 'C10orf110', 'INHA', 'PCNA', 'DNMT1', 'LOC153577', 'TRIB1', 'TICAM2', 'TMED7', 'PVRL2', 'C1orf201', 'CST6', 'BET1L', 'TAL1', 'TRIM42', 'SYNJ2', 'LOC284100', 'PHTF1', 'GKAP1', 'KIF11', 'KRT7', 'GRAMD3', 'EXOSC9', 'ARMCX4', 'ASB14', 'TROAP', 'GPR161', 'SLFN13', 'NSD1', 'MITF', 'SIRT6', 'KIF20B', 'HMGCLL1', 'CREB3L2', 'THEX1', 'HEXB', 'GRK5', 'SPAG1', 'SSX2IP', 'IGF2BP2', 'RNF182', 'TMEM2', 'hCG_2024094', 'WDR18', 'FBXO4', 'CKS2', 'TRAM2', 'FZD2', 'PIAS2', 'PSIP1', 'OBSL1', 'ELOVL6', 'GGH', 'C19orf55', 'RAB23', 'RBM28', 'PDE5A', 'LOC100131993', 'FAM176A', 'PGP', 'NTN4', 'DNAJC21', 'SFXN2', 'ALCAM', 'WNT10A', 'IVD', 'SLC45A3', 'LOC100129762', 'PRUNE2', 'IDI1', 'ANXA9', 'NMNAT2', 'KIAA0355', 'CENPJ', 'CDCA7L', 'ARHGEF17', 'COL1A1', 'SH3TC2', 'LOC54492', 'CLN8', 'WDR62', 'NRCAM', 'PSG1', 'PSG4', 'RNF167', 'GAS2L3', 'ACOT7', 'ELOVL5', 'TRIP13', 'AGPAT3', 'CCNE2', 'PGAP1', 'AMIGO2', 'SETD3', 'ZNF738', 'MGC39900', 'TMSL8', 'NEURL', 'STC1', 'FAM92A1', 'DMD', 'ZNF626', 'GRK6', 'FAM169A', 'PROS1', 'DDX39', 'GALNT14', 'CCND1', 'MCM3', 'LOC131185', 'RAD23B', 'CEP164', 'MAGED1', 'LYN', 'EXT1', 'TRAF7', 'WEE1', 'FSTL3', 'KIAA1524', 'DIP2C', 'ZNF672', 'SQLE', 'PLOD2', 'MPP5', 'PENK', 'NACAP1', 'SDPR', 'ANKFN1', 'PTRF', 'SLC16A6', 'FTO', 'IFRD2', 'PANX2', 'RAPH1', 'PEAR1', 'FN1', 'PLUNC', 'MACF1', 'PWWP2A', 'PSRC1', 'SLC20A1', 'TBC1D1', 'CCRN4L', 'C1orf159', 'SYDE1', 'MID1', 'HCK', 'CENPL', 'PRDX2', 'HMBS', 'WDR61', 'C13orf34', 'TMEM160', 'SSBP1', 'ZNF496', 'LIMA1', 'NT5E', 'CDT1', 'GPR64', 'UBE2T', 'C10orf47', 'BCL10', 'CA8', 'FLRT2', 'GLMN', 'SEMA5A', 'CCDC21', 'PVR', 'TAP2', 'C21orf45', 'GATA3', 'CD302', 'PHLDA1', 'NFYB', 'TPTE', 'MSTO1', 'MSTO2P', 'SH3RF1', 'AKT3', 'C15orf23', 'PSG6', 'GALNT3', 'SHC4', 'MALT1', 'KRT8', 'CCDC80', 'PHYHIPL', 'CAMK2N1', 'EXOC7', 'MAP1B', 'TOP3B', 'MNS1', 'H2AFZ', 'RFC5', 'ING3', 'GLS2', 'PKP4', 'CLSPN', 'GLT25D1', 'KRTAP4-8', 'TSPAN3', 'LOX', 'CENPE', 'LMNB2', 'BIN1', 'PRKACB', 'RAD1', 'IL6ST', 'ITGB5', 'MCF2', 'KLHDC10', 'C6orf125', 'DUSP1', 'DDEF1', 'OBFC2A', 'ZFYVE9', 'RP4-662A9.2', 'STC2', 'TIMELESS', 'C1orf61', 'ZNF655', 'UBE2J1', 'TMEM48', 'C10orf125', 'ST3GAL6', 'CXADR', 'SAMHD1', 'SORBS3', 'CTPS', 'ID4', 'GMNN', 'PARD3', 'FOXRED1', 'ULK4', 'ASPHD2', 'ATP1B1', 'PIF1', 'CD1D', 'MYLK4', 'CDC2L1', 'CDC2L2', 'HMGA2', 'LIG1', 'SCLY', 'SNX24', 'TEAD4', 'AFP', 'FAM122B', 'LOC283713', 'C1orf79', 'PTTG3', 'ACSF3', 'SGEF', 'STIL', 'FXYD5', 'SNX5', 'HBB', 'USP13', 'SPTLC2', 'MPRIP', 'CLDN18', 'TYRO3', 'CLN6', 'HYOU1', 'LMNA', 'GPR108', 'C11orf82', 'BRAP', 'FANCG', 'VCP', 'MED16', 'SCAMP4', 'CHD6', 'C11orf75', 'DENND4B', 'JPH1', 'CHEK1', 'PRR11', 'HSP90B1', 'RAD18', 'GSTM5', 'NFE2L3', 'SPATA6', 'FLJ11151', 'DDEFL1', 'TUBG1', 'B4GALT6', 'FAM130A1', 'POLE2', 'GANC', 'ILVBL', 'CRIM1', 'NUDT15', 'SLC1A3', 'SFRS2IP', 'MRPL43', 'ATP8B3', 'CREB1', 'ECT2', 'ARHGAP30', 'NADSYN1', 'CRELD1', 'NAT11', 'ZNF680', 'BAIAP2L1', 'LOC100128461', 'C8orf30A', 'ACYP1', 'KAT2B', 'PAPPA', 'SEC61B', 'LOC170082', 'CALU', 'KIAA1333', 'CSE1L', 'GCAT', 'ACTBL3', 'MND1', 'CDC25A', 'LONRF3', 'CIT', 'LARP7', 'USP1', 'SAE1', 'KIAA0040', 'SNRPB', 'RACGAP1', 'SLC7A6', 'C12orf31', 'BRPF3', 'ACSL3', 'SMA4', 'CPNE1', 'RBM12', 'HS3ST3B1', 'HYI', 'DBF4', 'RRP7A', 'PDGFRL', 'MET', 'NOS1AP', 'LIFR', 'TCOF1', 'TMEM47', 'ZBED4', 'PTGDR', 'DUT', 'PTPN1', 'RGS16', 'FIGNL1', 'ZDHHC18', 'PPIL5', 'PATZ1', 'IRAK3', 'LOC151162', 'B3GALNT1', 'ASPH', 'DDX27', 'ELL2', 'RNH1', 'UNQ338', 'UGT8', 'PCTK1', 'UGCG', 'PARP16', 'IQGAP3', 'C4orf46', 'HRBL', 'HMGN2', 'KLK3', 'SMC2', 'BAIAP2L1', 'SVIP', 'LPAR5', 'ZYX', 'PDLIM7', 'FRMD6', 'C9orf86', 'SLC4A10', 'PPP2R2C', 'MYPN', 'SEC61A2', 'ADM', 'SSR3', 'TNFAIP8L1', 'MAN1A1', 'NAV1', 'MYL9', 'RTN1', 'SLC46A1', 'RIMS3', 'DAP', 'B4GALT5', 'CENPH', 'CCDC124', 'CITED2', 'GNAQ', 'USP5', 'HSD17B6', 'TAZ', 'BAK1', 'MRS2', 'CDH11', 'SGCD', 'RNF41', 'ZNF843', 'PPAPDC1B', 'NLK', 'DAAM1', 'CAP2', 'MYO1C', 'SPECC1L', 'MED8', 'GTF2H2', 'SFXN1', 'LRCH4', 'MICB', 'STAMBPL1', 'AASS', 'CYB5B', 'TACSTD1', 'FAM111A', 'POLR2A', 'LOC100128918', 'LOC92497', 'RCCD1', 'RGNEF', 'WDR5B', 'PARD6A', 'GYS1', 'FASTKD2', 'C1orf69', 'TCEAL3', 'ABCC9', 'CENTG3', 'ZNRF3', 'KIAA1161', 'MAP2K3', 'DHFR', 'EHD4', 'RBM10', 'FUT8', 'BCAP29', 'CALM1', 'CALM2', 'CALM3', 'CCDC88A', 'CCDC49', 'FAM178A', 'GK2', 'BPGM', 'PXMP2', 'ZFP36L2', 'MRPS34', 'USP37', 'SCAPER', 'C14orf145', 'PCSK6', 'LOC100133781', 'TUBA4A', 'HOXA1', 'MAN2A2', 'GADD45B', 'EFEMP1', 'ADAMTS17', 'TMEM206', 'FUT10', 'INSIG1', 'DCK', 'TREX2', 'UCHL5IP', 'MED22', 'PRMT5', 'ARHGEF7', 'C1orf144', 'PTH2', 'USP39', 'LDLRAP1', 'ATXN1', 'EPOR', 'PTPLAD1', 'LYPLA3', 'GTSE1', 'ERBB2', 'CETN3', 'LOC338653', 'MPEG1', 'MT1F', 'KIAA0564', 'ABBA-1', 'ATP5SL', 'BACE1', 'EZH1', 'OPN3', 'MAN2A1', 'C13orf15', 'DYRK3', 'RRM1', 'C9orf40', 'TWSG1', 'NARS2', 'PELO', 'DOCK5', 'LOC400655', 'RCAN1', 'TMEM173', 'NPAL3', 'C14orf82', 'FEN1', 'DPH2', 'PM20D2', 'ADCY6', 'LOC100127955', 'LOC100128374', 'CHML', 'TMEM127', 'GTPBP6', 'ARMC8', 'PRICKLE1', 'SHMT2', 'ZNF107', 'F12', 'CTNNAL1', 'CCBL1', 'HUS1', 'OGFOD2', 'IVNS1ABP', 'TUBB6', 'COQ7', 'NHSL1', 'WDR66', 'DPF3', 'TFDP1', 'LPIN2', 'MUTED', 'TXNDC5', 'CCDC97', 'WDR34', 'SEMA4F', 'RPL22L1', 'UBE2D4', 'FDPS', 'NUP62', 'AFG3L2', 'KIAA1432', 'PA2G4', 'HNRNPA1', 'HNRPA1L-2', 'LOC100128836', 'LOC120364', 'LOC391670', 'LOC402112', 'LOC440125', 'LOC642817', 'LOC643033', 'LOC644037', 'LOC645001', 'LOC728170', 'LOC728643', 'LOC728732', 'LOC729102', 'LOC729366', 'LOC730246', 'RP11-78J21.1', 'C20orf27', 'KLHL7', 'C9orf91', 'RBBP7', 'STK10', 'C14orf139', 'ABHD2', 'DDX11', 'PRKCDBP', 'GNPTAB', 'KDELR2', 'TPM1', 'TRAF4', 'METT5D1', 'ILK', 'CMTM7', 'RPL39L', 'GTF3C5', 'KDELR1', 'C10orf85', 'MEGF6', 'PRKD1', 'CSGALNACT1', 'TUBB2A', 'GPSM1', 'EGFR', 'ZAK', 'NDE1', 'AAAS', 'NAP1L5', 'TMEM55B', 'C1orf52', 'ASXL1', 'SPCS3', 'DFFA', 'RFT1', 'TSEN15', 'CDKN2C', 'TPD52', 'PPARA', 'RFPL2', 'CHST13', 'CTH', 'MAST2', 'MMEL1', 'RAD51C', 'USP28', 'AK2', 'C11orf44', 'NIN')

induced_by_YAP_and_TAZ_Zhang = c('KRT34', 'STXBP6', 'NAV3', 'DIAPH3', 'INHBA', 'BDNF', 'SERTAD4', 'OLR1', 'PTGS2', 'WDR69', 'HSD3B1', 'CTGF', 'LRRC8C', 'THBS1', 'ETV5', 'IGFBP3', 'ORC1L', 'TGM2', 'ANKRD1', 'RECQL4', 'SERPINE1', 'DEPDC1B', 'CCNA2', 'C14orf145', 'CENPM', 'SERPINE2', 'TMEFF2', 'DDAH1', 'TMEM171', 'CDC25C', 'MCM10', 'CDC2', 'CST6', 'ADAMTS1', 'RRM2', 'AKAP12', 'NT5E', 'CDC25A', 'RAD51AP1', 'ITGB2', 'RAD51', 'CDC6', 'CCDC18', 'C18orf24', 'SHCBP1', 'KIF14', 'PLCB4', 'CNN3', 'UHRF1', 'TBXA2R', 'ZWINT', 'ZBED2', 'EXO1', 'SLIT2', 'KIF23', 'ADM', 'KIAA0101', 'BIRC5', 'CCDC99', 'PLK4', 'MATN3', 'DTL', 'CCNB1', 'CEP55', 'CDCA7', 'AMIGO2', 'RIMS2', 'UGCGL2', 'FBXO5', 'AURKA', 'KDELR3', 'GINS1', 'GINS2', 'ZNF367', 'BCAR4', 'CENPA', 'ASF1B', 'DDX11', 'NAV1', 'CCNE2', 'BRCA2', 'FAM111B', 'EGR1', 'TBC1D24', 'TYMS', 'LMNB1', 'CLSPN', 'HMGA2', 'BUB1', 'PRR16', 'SCLY', 'PSRC1', 'STAMBPL1', 'ORC6L', 'AURKB', 'BARD1', 'FMN2', 'ATAD2', 'KIF2C', 'FGF2', 'STMN1', 'DAB2', 'POLE2', 'CDCA3', 'DUSP6', 'MAD2L1', 'MCM5', 'TOP2A', 'FAM54A', 'C13orf3', 'DLC1', 'HELLS', 'TMEM56', 'F3', 'CASC5', 'ANKRD28', 'ANLN', 'HIC2', 'MCM6', 'SLC20A1', 'MELK', 'FOXM1', 'WDHD1', 'POLQ', 'CCNF', 'TMPO', 'COL8A1', 'HMGB2', 'POLA2', 'RFC5', 'E2F7', 'MKI67', 'SSR3', 'PRC1', 'ARSJ', 'KIF4A', 'FBN2', 'PHLDA1', 'TACC3', 'CDCA8', 'EXOSC9', 'CDCA2', 'CLIC5', 'MCM7', 'NEK2', 'KIAA1524', 'MLF1IP', 'WHSC1', 'CENPN', 'SLC25A15', 'CENPK', 'DEPDC1', 'SLC16A6', 'SGOL2', 'LRP8', 'CDKN3', 'CEP152', 'HMMR', 'C9orf25', 'RFC4', 'CDT1', 'RFC3', 'GJA5', 'NEIL3', 'FYN', 'BRCA1', 'LOC440338', 'GMNN', 'ING3', 'C6orf173', 'TTK', 'C13orf34', 'CENPE', 'SPAG1', 'CDC7', 'HSP90B1', 'DNAJC6', 'SPAG5', 'CDCA7L', 'GAS2L3', 'KIF11', 'LHFP', 'PCDH7', 'BCAP29', 'MCAM', 'C9orf140', 'CENPF', 'CCL20', 'BUB1B', 'WDR51A', 'BANP', 'UBE2T', 'SERPINH1', 'KRT7', 'DBF4', 'TK1', 'MYBL2', 'OIP5', 'RAD18', 'WDR4', 'MYBL1', 'TPX2', 'ALPK2', 'SPRED1', 'UGCG', 'ODC1', 'CCNB2', 'KCTD14', 'AXL', 'CEP78', 'DGKZ', 'CKS2', 'FAM64A', 'CTPS', 'EPB41L2', 'CYR61', 'SPARC', 'ESPL1', 'CENPL', 'FANCD2', 'FIGNL1', 'SH3RF1', 'CCND1', 'KCNT2', 'FEN1', 'MND1', 'COL4A1', 'C6orf105', 'TRIP13', 'C21orf45', 'CHEK1', 'TEAD4', 'RNASEH2A', 'PBK', 'PCNA', 'ZAK', 'PHGDH', 'C15orf23', 'C19orf42', 'MFAP5', 'PPAPDC1B', 'TAGLN', 'GCNT1', 'GPR107', 'RAD54B', 'CDH2', 'CDK6', 'SYTL4', 'MT1F', 'PKMYT1', 'SAMD1', 'WWC1', 'PPIL5', 'GTSE1', 'DNMT1', 'CKAP2L', 'DDX39', 'MCM3', 'USP1', 'FGFR1', 'LETM2', 'BLM', 'SFXN2', 'FSCN1', 'SH2D4A', 'SMC4', 'NASP', 'RAPH1', 'KLHL7', 'SFRS10', 'TMEM97', 'ASPM', 'ATP13A3', 'TACSTD1', 'C11orf75', 'SLC44A2', 'GNG11', 'CALU', 'MGLL', 'CENPJ', 'CREB3L2', 'PTPLAD1', 'MCM2', 'NDE1', 'C1orf69', 'FAM83D', 'MAN2A1', 'UGT8', 'C10orf110', 'GALNT14', 'HMBS', 'TMEFF1', 'CPNE2', 'PHF19', 'ACSL3', 'PEG10', 'SPCS3', 'DUSP4', 'KIAA1333', 'KCTD15', 'APOL2', 'TMEM48', 'BRIP1', 'HMGB3', 'QSER1', 'FZD7', 'PKP4', 'TIMELESS', 'CTNNAL1', 'PRR11', 'SGEF', 'BLMH', 'MDFIC', 'ASPH', 'RAD51L3', 'RRM1', 'C14orf65', 'HYI', 'SEC22C', 'RFPL2', 'SCARA3', 'SUV39H1', 'SKP2', 'COL4A2', 'MCPH1', 'COL4A3', 'TNNT2', 'KCNK1', 'MED8', 'RHOT2', 'TMEM2', 'TUBB3', 'FZD5', 'PDGFRL', 'DCK', 'SLC25A25', 'STIL', 'TGFB2', 'KIF20A', 'RBM28', 'ANKRD2', 'PMP22', 'SLC45A3', 'TUBB2A', 'THSD4', 'PRMT5', 'UBE2J1', 'FAM3C', 'JPH1', 'KRT18', 'PSG9', 'COTL1', 'SLC7A11', 'LTBP1', 'PSG1', 'GLS', 'C9orf40', 'SETD3', 'RAD51C', 'C10orf47', 'SMARCA1', 'PRR6', 'E2F1', 'NAP1L5', 'NUSAP1', 'PTGER4', 'FAM111A', 'FLJ42709', 'CDH4', 'PLK1', 'TIMM10', 'TUFT1', 'H2AFZ', 'SHMT2', 'EXT1', 'IFRD2', 'PDXK', 'USP13', 'BCL2L1', 'DDEF1', 'FAM107B', 'HEXB', 'RFT1', 'LIFR', 'RACGAP1', 'TFDP1', 'USP39', 'H2AFX', 'EGFLAM', 'GALNT3', 'HMGA1', 'PLOD2', 'DHRS13', 'RASAL2', 'RRBP1', 'SYDE1', 'EXPH5', 'FOXRED1', 'NUDT15', 'RHOF', 'PASK', 'TPD52', 'IQGAP3', 'SFXN1', 'C4orf27', 'CSRP1', 'TMEM54', 'PLEC1', 'SNX24', 'DOCK5', 'KDELR2', 'NF2', 'OSBP2', 'PSIP1', 'GRAMD3', 'MNS1', 'NTN4', 'RBM10', 'WEE1', 'BAIAP2L1', 'AK2', 'TEAD1', 'CAMK2N1', 'PIN1', 'CENPH', 'ATG7', 'FSTL1', 'CFLAR', 'SYDE2', 'CHORDC1', 'FRAS1', 'SHB', 'ACAT2', 'CEP164', 'FAM92A1', 'MPP5', 'PSG6', 'ELL2', 'CDKN2C', 'SPTLC2', 'POLD1', 'TROAP', 'ELAVL1', 'TRAM2', 'SCHIP1', 'KHDRBS3', 'ZNF680', 'PXDN', 'ACYP1', 'CRIM1', 'FN1', 'HMGN2', 'MRPL12', 'ECT2', 'FRMD6', 'PRKACB', 'ACTN1', 'C1orf144', 'CSE1L', 'NMNAT2', 'USP37', 'HYOU1', 'LOC643287', 'DUT', 'HRBL', 'LOC151162', 'MALT1', 'PLA2G4A', 'SYNJ2', 'NARS2', 'SH3TC1', 'TRIM42', 'ATP6V0A4', 'COL12A1', 'GPR126', 'SLC1A4', 'TAP2', 'AP4E1', 'SAMD4A', 'TCOF1', 'LMNB2', 'NUP62', 'SMG7', 'PSG5', 'THEX1', 'C9orf91', 'CCDC49', 'CYB5B', 'EHD1', 'MAP1B', 'NLK', 'ATP1B1', 'WWC2', 'ABHD2', 'COQ7', 'PDLIM2', 'HOXA1', 'LONRF3', 'NUP210', 'SNRPB', 'ACOT7', 'GLDC', 'C10orf35', 'PARD6A', 'ODZ3', 'SIM2', 'RTN1', 'SNX5', 'PVR', 'ZNF738', 'STRA6', 'FRY', 'SSX2IP', 'CCRN4L', 'DFFA', 'WDR62', 'C1orf52', 'SVEP1', 'CREB1', 'DPH2', 'MSRB3', 'CDC42EP3', 'ELOVL6', 'FAM80B', 'SGCD', 'CIT', 'KRT8', 'PPARA', 'SNHG4', 'TCF7L1', 'HEG1', 'HS3ST3B1', 'PRSS23', 'PARVA', 'PIAS2', 'LOC202051', 'LYN', 'SH2D5', 'NKX3-1', 'CETN3', 'CMTM7', 'RBBP7', 'TUBB6', 'LDLRAP1', 'ITGB5', 'SMC2', 'OBFC2A', 'TRIB1', 'FSTL3', 'NIN', 'NRCAM', 'CDC20', 'EHD4', 'NPAL3', 'FBN1', 'GANC', 'GK2', 'NR3C2', 'ASPHD2', 'C1orf79', 'GGH', 'IVD', 'PSMD9', 'SFRP1', 'SLC7A1', 'TPM1', 'ZFP36L2', 'SAMHD1', 'TUBG1', 'GATA3', 'PHTF1', 'SLC35B4', 'CBX5', 'GNAQ', 'AFG3L2', 'AKT3', 'CLN8', 'CRY1', 'HABP4', 'MCTP2', 'RAB23', 'SDPR', 'ZNF496', 'BAK1', 'OGFOD2', 'ARHGAP29', 'MYADM', 'ZNF655', 'HBB', 'MRPS34', 'PARVB', 'SAE1', 'TOP3B', 'FLJ14213', 'GLT25D1', 'RPL39L', 'USP28', 'EMP2', 'SEC24A', 'WDR18', 'ZDHHC18', 'C6orf125', 'IDH2', 'MMEL1', 'FZD1', 'PTTG3', 'CHST3', 'IGF2BP2', 'AMOTL2', 'GCAT', 'GRK6', 'LY6E', 'CLN6', 'PENK', 'CHML', 'NSD1', 'SFRS2IP', 'RASSF7', 'TRAF7')

tead_dependent_YAP_and_TAZ_Zhang = c('ABHD2', 'BCAR4', 'BDNF', 'CHST3', 'CTGF', 'CYR61', 'DAB2', 'FMN2', 'FRY', 'GGH', 'GJA5', 'ITGB2', 'LHFP', 'LIFR', 'MFAP5', 'OLR1', 'PARVA', 'PDGFRL', 'PMP22', 'PRR16', 'PRSS23', 'PTGS2', 'PXDN', 'RASAL2', 'SCARA3', 'SMARCA1', 'SPARC', 'SPRED1', 'STXBP6', 'SYDE2', 'TGFB2', 'TMEM56', 'ZNF738')

forOverlap = list()
for(cluster in levels(cellcycle.markers[,'cluster'])) {
    forOverlap[[cluster]] = cellcycle.markers[which(cellcycle.markers[,'cluster']==cluster & cellcycle.markers[,'avg_logFC']>0 & cellcycle.markers[,'p_val_adj']<=0.05),'gene']
}

m2 = matrix(nrow=nrow(cellcycle.combined@raw.data),ncol=ncol(cellcycle.combined@raw.data))
rownames(m2) = rownames(cellcycle.combined@raw.data)
m2[,] = 0
m2[as.matrix(cellcycle.combined@raw.data)>3] = 1
rs1 = rowSums(m2)
rs1[c('FAM181B','MAP2','TAGLN3','LAMA4','SOX6','EDNRB','RND2','RGS6','LRRC4C','HIST1H2AC','COL11A1','PLEKHB1','NLRP1','PNRC1','LIMCH1','BCHE','EFHC1','SMIM14','SLC44A1')]
allGenes = rownames(cellcycle.combined@raw.data)[which(rowSums(m2)>10)]
lateG1 = forOverlap[['Late G1']]

c('FAM181B','MAP2','TAGLN3','LAMA4','SOX6','EDNRB','RND2','RGS6','LRRC4C','HIST1H2AC','COL11A1','PLEKHB1','NLRP1','PNRC1','LIMCH1','BCHE','EFHC1','SMIM14','SLC44A1')

#induced_by_YAP_and_TAZ_Zhang
intersect(intersect(lateG1,induced_by_YAP_and_TAZ_Zhang),allGenes)
q = length(intersect(intersect(lateG1,induced_by_YAP_and_TAZ_Zhang),allGenes))
m = length(intersect(lateG1,allGenes))
n = length(allGenes) - m
k = length(intersect(induced_by_YAP_and_TAZ_Zhang,allGenes))
phyper(q,m,n,k,lower.tail=F)

#induced_by_YAP_Zhang
intersect(intersect(lateG1,induced_by_YAP_Zhang),allGenes)
q = length(intersect(intersect(lateG1,induced_by_YAP_Zhang),allGenes))
m = length(intersect(lateG1,allGenes))
n = length(allGenes) - m
k = length(intersect(induced_by_YAP_Zhang,allGenes))
phyper(q,m,n,k,lower.tail=F)

#tead_dependent_YAP_and_TAZ_Zhang
intersect(intersect(lateG1,tead_dependent_YAP_and_TAZ_Zhang),allGenes)
q = length(intersect(intersect(lateG1,tead_dependent_YAP_and_TAZ_Zhang),allGenes))
m = length(intersect(lateG1,allGenes))
n = length(allGenes) - m
k = length(intersect(tead_dependent_YAP_and_TAZ_Zhang,allGenes))
phyper(q,m,n,k,lower.tail=F)


