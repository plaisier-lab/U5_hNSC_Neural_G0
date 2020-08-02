## [Neural G0:  a quiescent-like state found in neuroepithelial-derived cells and glioma](https://www.biorxiv.org/content/10.1101/446344v3)

![UMAPs](umaps.gif)

### Table of Contents
- **[Abstract](#abstract)**
- **[Data and code availability](#data-and-code-availability)**
- **[Instructions to setup data and code for recreating analyses](#instructions-to-setup-data-and-code-for-recreating-analyses)**
- **[Directory structure](#directory-structure)**
- **[Order of operations](#order-of-operations)**
- **[Contact](#contact)**
- **[Citation](#citation)**

### Abstract
In depth knowledge of the cellular states associated with normal and disease tissue homeostasis is critical for understanding disease etiology and uncovering therapeutic opportunities. Here, we used single cell RNA-seq to survey the cellular states of neuroepithelial-derived cells in cortical and neurogenic regions of developing and adult mammalian brain to compare with 38,474 cells obtained from 59 human gliomas, as well as pluripotent ESCs, endothelial cells, CD45+ immune cells, and non-CNS cancers. This analysis suggests that a significant portion of neuroepithelial-derived stem and progenitor cells and glioma cells that are not in G2/M or S phase exist in two states: G1 or Neural G0, defined by expression of certain neuro-developmental genes. In gliomas, higher overall Neural G0 gene expression is significantly associated with less aggressive gliomas, IDH1 mutation, and extended patient survival, while also anti-correlated with cell cycle gene expression. Knockout of genes associated with the Hippo/Yap and p53 pathways diminished Neural G0 in vitro, resulting in faster G1 transit, down regulation of quiescence-associated markers, and loss of Neural G0 gene expression. Thus, Neural G0 is a dynamic cellular state required for indolent cell cycles in neural-specified stem and progenitors poised for cell division. As a result, Neural G0 occupancy may be an important determinant of glioma tumor progression.

### Data and code availability
- All datafiles need to run this code can be found on [figshare](https://figshare.com/projects/Neural_G0_a_quiescent-like_state_found_in_neuroepithelial-derived_cells_and_glioma/86939).
- The analysis software and scripts are available on [github](https://github.com/plaisier-lab/U5_hNSC_Neural_G0/).

### Instructions to setup data and code for recreating analyses
In order to run the software and scripts you will need to setup a specific directory structure and download all the data and scripts. Here are the instructions to setup things up:
1. Clone the [U5_hNSC_Neural_G0 repository](https://github.com/plaisier-lab/U5_hNSC_Neural_G0/)
2. Make a "data" folder inside the U5_hNSC_Neural_G0 folder
3. Download (and unzip for zip files) all files from [figshare](https://figshare.com/projects/Neural_G0_a_quiescent-like_state_found_in_neuroepithelial-derived_cells_and_glioma/86939):
    - [U5_hNSC.zip](https://figshare.com/articles/dataset/U5_hNSC_zip/12751082) (needs to be unzipped) - contains all the U5 hNSC scRNA-seq datasets as output from cellranger.
    - [ccAF_1536_smaller.pkl](https://figshare.com/articles/software/ccAF_1536_smaller_pkl/12751058) (does not need to be unzipped) - the ccAF ACTINN loadings for classification of cell cycle phases for cells or transcriptome profiles.
    - [geneConversions.zip](https://figshare.com/articles/dataset/geneConversions_zip/12751073) (needs to be unzipped) - helpful gene ID conversion files.
    - [forClassification.zip](https://figshare.com/articles/dataset/forClassification_zip/12751079) (needs to be unzipped) - loom data files that were classified using ccAF.
    - [ssGSEA.GBM.classification.zip](https://figshare.com/articles/dataset/ssGSEA_GBM_classification_zip/12751076) (needs to be unzipped) - subtype classification results for all glioma datasets.
    - [cellcycle_int_integrated.loom](https://figshare.com/articles/dataset/cellcycle_int_integrated_loom/12751055) (does not need to be unzipped) - U5 hNSC data as a loom file that was used to build the ccAF classifier.


#### Directory structure
After downloading and unzipping the files the directory structure should look like this:

```
.
+-- U5_hNSC_Neural_G0
|   +-- actinn.py
|   +-- calculatingErrors.py
|   +-- calculatingErrors_CV.py
|   +-- calculatingErrors_Whitfield.py
|   +-- classifiersV3.py
|   +-- classifyPrimaryCells_gliomas.py
|   +-- classifyPrimaryCells_homoSapeins.py
|   +-- classifyPrimaryCells_musMusculus.py
|   +-- converting_to_loom.R
|   +-- cvClassification_FullAnalysis.py
|   +-- data
|   |   +-- ccAF_1536_smaller.pkl
|   |   +-- cellcycle_int_integrated.loom
|   |   +-- forClassification
|   |   |   +-- gliomas
|   |   |   |   +-- Bhaduri.loom
|   |   |   |   +-- GSE102130.loom
|   |   |   |   +-- GSE131928_10X.loom
|   |   |   |   +-- GSE131928_Smartseq2.loom
|   |   |   |   +-- GSE139448.loom
|   |   |   |   +-- GSE70630.loom
|   |   |   |   +-- GSE84465.loom
|   |   |   |   +-- GSE84465_all.loom
|   |   |   |   +-- GSE89567.loom
|   |   |   +-- GSE103322.loom
|   |   |   +-- GSE67833.loom
|   |   |   +-- HEK293T.loom
|   |   |   +-- Nowakowski_norm.loom
|   |   |   +-- PRJNA324289.loom
|   |   +-- geneConversions
|   |   |   +-- ensembl_entrez.csv
|   |   |   +-- hgnc_geneSymbols.txt
|   |   |   +-- hgnc_geneSymbols_ensmbl.txt
|   |   |   +-- human_hgnc_mouse_mgi.csv
|   |   |   +-- mart_export.txt
|   |   +-- highlyVarGenes_WT_sgTAOK1_1584.csv
|   |   +-- ssGSEA.GBM.classification
|   |   |   +-- p_result_Bhaduri_2019.gct.txt
|   |   |   +-- p_result_GSE102130.gct.txt
|   |   |   +-- p_result_GSE131928_GSM3828672.gct.txt
|   |   |   +-- p_result_GSE131928_GSM3828673_1.gct.txt
|   |   |   +-- p_result_GSE131928_GSM3828673_2.gct.txt
|   |   |   +-- p_result_GSE139448.gct.txt
|   |   |   +-- p_result_GSE70630.gct.txt
|   |   |   +-- p_result_GSE84465.gct.txt
|   |   |   +-- p_result_GSE89567.gct.txt
|   |   |   +-- res_All_GSE131928.csv
|   |   +-- U5_hNSC
|   |   |   +-- sgTAOK1
|   |   |   |   +-- filtered_gene_bc_matrices
|   |   |   |   |   +-- hg19
|   |   |   |   |   |   +-- barcodes.tsv
|   |   |   |   |   |   +-- genes.tsv
|   |   |   |   |   |   +-- matrix.mtx
|   |   |   +-- WT
|   |   |   |   +-- filtered_gene_bc_matrices
|   |   |   |   |   +-- hg19
|   |   |   |   |   |   +-- barcodes.tsv
|   |   |   |   |   |   +-- genes.tsv
|   |   |   |   |   |   +-- matrix.mtx
|   |   |   +-- WT_CDTplus
|   |   |   |   +-- filtered_gene_bc_matrices
|   |   |   |   |   +-- hg19
|   |   |   |   |   |   +-- barcodes.tsv
|   |   |   |   |   |   +-- genes.tsv
|   |   |   |   |   |   +-- matrix.mtx
|   |   +-- Whitfield
|   |   |   +-- data
|   |   |   |   +-- whitfield_dataPlusScores_6_30_2020_SHAKE.T_1334.csv
|   |   |   |   +-- whitfield_dataPlusScores_6_30_2020_SHAKE_1334.csv
|   |   |   |   +-- whitfield_dataPlusScores_6_30_2020_TN.T_1334.csv
|   |   |   |   +-- whitfield_dataPlusScores_6_30_2020_TN_1334.csv
|   |   |   |   +-- whitfield_dataPlusScores_6_30_2020_TT1.T_1334.csv
|   |   |   |   +-- whitfield_dataPlusScores_6_30_2020_TT1_1334.csv
|   |   |   |   +-- whitfield_dataPlusScores_6_30_2020_TT2.T_1334.csv
|   |   |   |   +-- whitfield_dataPlusScores_6_30_2020_TT2_1334.csv
|   |   |   |   +-- whitfield_dataPlusScores_6_30_2020_TT3.T_1334.csv
|   |   |   |   +-- whitfield_dataPlusScores_6_30_2020_TT3_1334.csv
|   |   |   +-- markergenes_ForPlotting.csv
|   |   |   +-- metaInformation.csv
|   +-- LICENSE
|   +-- makeDatasetsForClassification.R
|   +-- plotNowakowski.py
|   +-- plottingClassifiers.py
|   +-- README.md
|   +-- results
|   +-- sensitivityAnalysis_plot.py
|   +-- sensitivityAnalysis_run.py
|   +-- ssgse.GBM.classification
|   +-- U5_hNSC_scRNA_seq_Analysis.R
|   +-- Whitfield_classification_ACTINN_analysis.py
```

The results directory will hold the output from analysis scripts.

### Docker container
We facilitate the use of our code and data by providing a Docker Hub container [cplaisier/scrna_seq_velocity](https://hub.docker.com/r/cplaisier/scrna_seq_velocity) which has all the dependencies and libraries to run the scripts. Please install Docker and then from the command line run:

```shell
docker pull cplaisier/scrna_seq_velocity
```

### Order of operations
The order of analyses in this study and the details of each analysis are described below:

#### 1. Identification of cell cycle phases and candidate G0/G1 subpopulations in human NSCs


### Contact
For issues or comments please contact:  [Chris Plaisier](mailto:plaisier@asu.edu)

### Citation
(Neural G0: a quiescent-like state found in neuroepithelial-derived cells and glioma.)[https://doi.org/10.1101/446344] Samantha A. O'Connor, Heather M. Feldman, Chad M. Toledo, Sonali Arora, Pia Hoellerbauer, Philip Corrin, Lucas Carter, Megan Kufeld, Hamid Bolouri, Ryan Basom, Jeffrey Delrow, José L. McFaline-Figueroa, Cole Trapnell, Steven M. Pollard, Anoop Patel, Patrick J. Paddison, Christopher L. Plaisier. bioRxiv 446344; doi: (https://doi.org/10.1101/446344)[https://doi.org/10.1101/446344]
