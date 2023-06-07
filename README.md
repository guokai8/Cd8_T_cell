# The chemokine receptor CXCR3 promotes CD8 T cell dependent lung pathology during influenza pathogenesis
### Description
Data were collected after mapping all fastq files to the mm10 genome with the cellranger v5.0 software with default paramaters by using "cellranger count". Details for alignment will be found in the "web_summary.html" in the each data file
### Required packages:
Seurat v3, CellChat, tidyverse, richR, scGSVA, VennDetail, SCENIC(v1.2.2), ggrepel
```{}
install.packages(c('Seurat','tidyverse'))
BicManager::install(c("CellChat","VennDetail","SCENIC))
devtools::install_github('guokai8/richR')
devtools::install_github('guokai8/scGSVA')
```
### Download data 
```{}
git clone https://github.com/guokai8/Cd8_T_cell/
cd Cd8_T_cell
```
Then load data into R 4.0.2 and do the analysis by following below steps
#### Step1 (Fig.2, Fig.S2, Fig.S3, Data S1-S9)
run the Tcell_clustering.R for clustering, cell annotation, Differential expressed analysis and GSEA.

#### Step2 (Fig.4, Fig.S3)
run the tcell_cell.r for the cell cell communication analysis

#### Step3 (Fig.4)
run the tcell_scenic.r for the TF network analysis
