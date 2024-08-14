# HTOreader
An R toolkit for cell hashing (HTO) demultiplexing

## How to install
make sure you installed these two dependencies:<br>
`install.packages("flexmix")` <br> `install.packages("hash")`<br>
then:<br>
`if (!require("devtools")) {install.packages("devtools")}`<br>
`devtools::install_github("WilsonImmunologyLab/HTOreader")`

## Hnads on Example of Hybrid demultiplexing using real-world data
Please click this link for a detailed example<br>
<img src="https://github.com/WilsonImmunologyLab/HTOreader/blob/main/docs/workflow.png">

<img src="https://github.com/WilsonImmunologyLab/HTOreader/blob/main/docs/arrow_down.png" width="100">

>[Hybrid demultiplexing using real-world data](https://wilsonimmunologylab.github.io/Librator/Hybrid)

<img src="https://github.com/WilsonImmunologyLab/HTOreader/blob/main/docs/arrow_up.png" width="100">

<br>
Please refer to /docs folder of this repository. We have provided the HTML version and the R markdown file. 
Data for this sample under /data folder. 

## Quick start
The `pbmc.hashtag` in the following example is a Seurat object. Please load your "HTO" assay using the Seurat package (More details in our example). 

There are two main functions, `HTOClassification` for Hashtag demultiplexing, `PlotHTO` for plotting distributions of hashtags <br>
`pbmc.hashtag <- HTOClassification(pbmc.hashtag, assay = 'HTO', method = 'CLR')` <br>
`VlnPlot(pbmc.hashtag,features = rownames(pbmc.hashtag@assays[["HTO"]]@data), group.by = 'HTOid')` <br>
![000009](https://user-images.githubusercontent.com/4589583/161608424-2a748fdf-5872-49fa-b519-ef0519b30b48.png)

Users can use PlotHTO function to check the quality of the cutoff of each HTO <br>
`plots1 <- PlotHTO(pbmc.hashtag, assay = 'HTO', method = 'CLR')` <br>
this function returns a list of ggplot objects, use plot_grid function to plot them: <br>
`plot_grid(plotlist = plots1)` <br>
![000009 (1)](https://user-images.githubusercontent.com/4589583/161609949-0599145c-a03b-466c-bf98-1c652cd4ce83.png)


Users can also specify the cutoff manually:  <br>
`pbmc.hashtag <- HTOClassification(pbmc.hashtag, assay = 'HTO', method = 'CLR', specify_cutoff = c(2.5,2,2,2,2,2,2,2))`  <br>
or use a different normalization method (CLR or log):  <br>
`pbmc.hashtag <- HTOClassification(pbmc.hashtag, assay = 'HTO', method = 'log')`  <br>
users also can specify the cutoff manually for `PlotHTO` function to find out the best cutoff: <br>
`plots1 <- PlotHTO(pbmc.hashtag, assay = 'HTO', method = 'CLR', specify_cutoff = c(2.5,2,2,2,2,2,2,2)) ` <br>
