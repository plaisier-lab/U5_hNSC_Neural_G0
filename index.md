## [Neural G0:  a quiescent-like state found in neuroepithelial-derived cells and glioma](https://www.biorxiv.org/content/10.1101/446344v3)

![UMAPs](umaps.gif)

### Abstract
In depth knowledge of the cellular states associated with normal and disease tissue homeostasis is critical for understanding disease etiology and uncovering therapeutic opportunities. Here, we used single cell RNA-seq to survey the cellular states of neuroepithelial-derived cells in cortical and neurogenic regions of developing and adult mammalian brain to compare with 38,474 cells obtained from 59 human gliomas, as well as pluripotent ESCs, endothelial cells, CD45+ immune cells, and non-CNS cancers. This analysis suggests that a significant portion of neuroepithelial-derived stem and progenitor cells and glioma cells that are not in G2/M or S phase exist in two states: G1 or Neural G0, defined by expression of certain neuro-developmental genes. In gliomas, higher overall Neural G0 gene expression is significantly associated with less aggressive gliomas, IDH1 mutation, and extended patient survival, while also anti-correlated with cell cycle gene expression. Knockout of genes associated with the Hippo/Yap and p53 pathways diminished Neural G0 in vitro, resulting in faster G1 transit, down regulation of quiescence-associated markers, and loss of Neural G0 gene expression. Thus, Neural G0 is a dynamic cellular state required for indolent cell cycles in neural-specified stem and progenitors poised for cell division. As a result, Neural G0 occupancy may be an important determinant of glioma tumor progression.

### Data and code availability
All datafiles need to run this code can be found on [figshare](https://figshare.com/projects/Neural_G0_a_quiescent-like_state_found_in_neuroepithelial-derived_cells_and_glioma/86939). The analysis software and scripts are available on [github](https://github.com/plaisier-lab/U5_hNSC_Neural_G0/).

### Instructions to setup data and code for recreating analyses
In order to run the software and scripts you will need to setup a specific directory structure and download all the data and scripts. Here are the instructions to setup things up:
1. Clone the [U5_hNSC_Neural_G0 repository](https://github.com/plaisier-lab/U5_hNSC_Neural_G0/)
2. Make a data folder inside the U5_hNSC_Neural_G0 folder
3. Download (and unzip for zip files) all files from [figshare](https://figshare.com/projects/Neural_G0_a_quiescent-like_state_found_in_neuroepithelial-derived_cells_and_glioma/86939)

Here is what the fila directory structure should look like:
'''
.
+-- U5_hNSC_Neural_G0
    +-- actinn.py
    +-- calculatingErrors.py
    +-- calculatingErrors_CV.py
    +-- calculatingErrors_Whitfield.py
    +-- classifiersV3.py
    +-- classifyPrimaryCells_gliomas.py
    +-- classifyPrimaryCells_homoSapeins.py
    +-- classifyPrimaryCells_musMusculus.py
    +-- converting_to_loom.R
    +-- cvClassification_FullAnalysis.py
    +-- data
        +-- ccAF_1536_smaller.pkl
        +-- cellcycle_int_integrated.loom
        +-- forClassification
            +-- gliomas
                +-- Bhaduri.loom
                +-- GSE102130.loom
                +-- GSE131928_10X.loom
                +-- GSE131928_Smartseq2.loom
                +-- GSE139448.loom
                +-- GSE70630.loom
                +-- GSE84465.loom
                +-- GSE84465_all.loom
                +-- GSE89567.loom
            +-- GSE103322.loom
            +-- GSE67833.loom
            +-- HEK293T.loom
            +-- Nowakowski_norm.loom
            +-- PRJNA324289.loom
        +-- geneConversions
            +-- ensembl_entrez.csv
            +-- hgnc_geneSymbols.txt
            +-- hgnc_geneSymbols_ensmbl.txt
            +-- human_hgnc_mouse_mgi.csv
            +-- mart_export.txt
        +-- highlyVarGenes_WT_sgTAOK1_1584.csv
        +-- ssGSEA.GBM.classification
            +-- p_result_Bhaduri_2019.gct.txt
            +-- p_result_GSE102130.gct.txt
            +-- p_result_GSE131928_GSM3828672.gct.txt
            +-- p_result_GSE131928_GSM3828673_1.gct.txt
            +-- p_result_GSE131928_GSM3828673_2.gct.txt
            +-- p_result_GSE139448.gct.txt
            +-- p_result_GSE70630.gct.txt
            +-- p_result_GSE84465.gct.txt
            +-- p_result_GSE89567.gct.txt
            +-- res_All_GSE131928.csv
        +-- U5_hNSC
            +-- 
        +-- Whitfield
    +-- LICENSE
    +-- makeDatasetsForClassification.R
    +-- makeDatasetsForClassification.R~
    +-- plotNowakowski.py
    +-- plottingClassifiers.py
    +-- README.md
    +-- results
    +-- sensitivityAnalysis_plot.py
    +-- sensitivityAnalysis_run.py
    +-- ssgse.GBM.classification
    +-- U5_hNSC_scRNA_seq_Analysis.R
    +-- Whitfield_classification_ACTINN_analysis.py
'''

#### Description of datafiles



You can use the [editor on GitHub](https://github.com/plaisier-lab/U5_hNSC_Neural_G0/edit/gh-pages/index.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/plaisier-lab/U5_hNSC_Neural_G0/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
