# Final Project
TRGN 510 Final Project: Colon Cancer

## Section 1: Objective
For this project, I will be showing the log2(fold-change) in differential gene expression for RNA-seq data in two subsets of colon cancer: adenocarcinoma and cystic, mucinous and serous neoplasms. The final deliverable will be the the Glimma Vignette. The input files for the Glimma vignette will be HTSeq-count data obtained from The Cancer Genome Atlas Program (TCGA). I will be turning in an HTML containing my R Notebook.

## Section 2: Datasets
### Loading Data

Data is obtained from TCGA. I filtered for RNA-Seq experimental strategy, TXT data format and HTSeq-counts workflow type. HTSeq-counts is a tool that quantifies the aligned reads overlapping a gene's exons. HTSeq data does not have a header, is tab-delimited, the first column is the Ensembl gene ID and the second column is the number of mapped reads of the gene. The counts will be used in differential gene expression analysis using edgeR as the method. To look at the differential gene expression, the counts will be normalized using the calcNormFactors in edgeR and only reads that unambigously map to one gene are used. 
  
  
[Cystic, mucinous, and serous neoplasms](https://portal.gdc.cancer.gov/repository?filters=%7B"op"%3A"and"%2C"content"%3A%5B%7B"content"%3A%7B"field"%3A"cases.case_id"%2C"value"%3A%5B"set_id%3AAW45M6AKTt_rMbGdDakT"%5D%7D%2C"op"%3A"IN"%7D%2C%7B"op"%3A"in"%2C"content"%3A%7B"field"%3A"cases.disease_type"%2C"value"%3A%5B"Cystic%2C%20Mucinous%20and%20Serous%20Neoplasms"%5D%7D%7D%2C%7B"op"%3A"in"%2C"content"%3A%7B"field"%3A"files.analysis.workflow_type"%2C"value"%3A%5B"HTSeq%20-%20Counts"%5D%7D%7D%2C%7B"op"%3A"in"%2C"content"%3A%7B"field"%3A"files.data_format"%2C"value"%3A%5B"TXT"%5D%7D%7D%2C%7B"op"%3A"in"%2C"content"%3A%7B"field"%3A"files.experimental_strategy"%2C"value"%3A%5B"RNA-Seq"%5D%7D%7D%5D%7D&searchTableTab=cases)
, I will choose 30.

[Adenocarcinoma](https://portal.gdc.cancer.gov/repository?filters=%7B"op"%3A"and"%2C"content"%3A%5B%7B"content"%3A%7B"field"%3A"cases.case_id"%2C"value"%3A%5B"set_id%3AAW45M6AKTt_rMbGdDakT"%5D%7D%2C"op"%3A"IN"%7D%2C%7B"op"%3A"in"%2C"content"%3A%7B"field"%3A"cases.disease_type"%2C"value"%3A%5B"Adenomas%20and%20Adenocarcinomas"%5D%7D%7D%2C%7B"op"%3A"in"%2C"content"%3A%7B"field"%3A"files.analysis.workflow_type"%2C"value"%3A%5B"HTSeq%20-%20Counts"%5D%7D%7D%2C%7B"op"%3A"in"%2C"content"%3A%7B"field"%3A"files.data_format"%2C"value"%3A%5B"TXT"%5D%7D%7D%2C%7B"op"%3A"in"%2C"content"%3A%7B"field"%3A"files.experimental_strategy"%2C"value"%3A%5B"RNA-Seq"%5D%7D%7D%5D%7D&searchTableTab=cases)
, I will choose 30. 
  
### Unit test
When data is read in: 
````
ENSG00000000003.13
<fctr>
X5290
<int>
ENSG00000000005.5	47			
ENSG00000000419.11	1212			
ENSG00000000457.12	1176			
ENSG00000000460.15	121			
ENSG00000000938.11	166			
ENSG00000000971.14	1012			
ENSG00000001036.12	4401			
ENSG00000001084.9	1977			
ENSG00000001167.13	976			
ENSG00000001460.16	1638			
1-10 of 60 rows
````
When dataset is made using the 60 text files:
````
[1] "DGEList"
attr(,"package")
[1] "edgeR"
[1] 60487    60
````

### Unit test for normalized data
Make a box plots of unnormalized and normalized data in R
```{r}
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors
```
## Proposed Analysis

This project will compute and analyze the logarithmic ratio of differential gene expression of two subtypes of colon cancer. EdgeR will be used to import, organize, and normalize the data,  Mus.musculus will be used for gene annotions, limma will be used to examine the gene expression anaylsis and make exploratory plots, Glimma will be used to make these plots interactive. RColorBrewer and gplots will be used to make heatmaps.

### Flowchart
```{r}
library(DiagrammeR)
grViz("digraph flowchart {
      # node definitions with substituted label text
      node [fontname = Helvetica, shape = rectangle]        
      tab1 [label = '@@1']
      tab2 [label = '@@2']
      tab3 [label = '@@3']
      tab4 [label = '@@4']
      tab5 [label = '@@5']
      tab6 [label = '@@6']
      tab7 [label = '@@7']
      tab8 [label = '@@8']
      tab9 [label = '@@9']
      tab10 [label = '@@10']

      # edge definitions with the node IDs
      tab1 -> tab2;
      tab2 -> tab3;
      tab3 -> tab4;
      tab4 -> tab8 -> tab5;
      tab4 -> tab5 -> tab6 -> tab7 -> tab9 -> tab10
      }

      [1]: 'Download necessary libraries'
      [2]: 'Load and read datasets'
      [3]: 'Join datasets'
      [4]: 'Unit test: Data was properly loaded?'
      [5]: 'Normalize data'
      [6]: 'Find Mean Varience Trend'
      [7]: 'Analyze DE genes'
      [8]: 'Troubleshoot and fix errors'
      [9]: 'Make Interactive MDS plot'
      [10]: 'Make HeatMap of log-CPM data'
      ")

```
![FLowchart](https://github.com/Aenorieg/FinalProject/blob/master/TRGN%20510%20Final%20Project%20Flowchart.png)


### Loading libraries
```{r}
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
library(gplots)
library(RColorBrewer)
```

## Proposed Timeline and Milestones
Week 1: Run the Glimma vignette. I will install the necessary packages in R and understand each step in the vignette.


Week 2: Load in the data (joins, creating datasets) and do a simple, 1 line unit test to look at the data. I will download 60 datasets (30 from each subtype) and join multiple datasets. Emailed Dr. Craig on 11/19 and agreed on turning this milestone in on Sat Nov 23, 2019.


Week 3: Confirm that the data was loaded in correctly and analyze data using the Glimma vignette.Emailed Dr. Craig on 11/26 and agreed on turning this milestone in on 12/1. 


Week 4: Troubleshoot for more errors and enhance the user interface.

## User Interface
I anticipate having boxplots, heatmaps, and interactive multi-dimensional scaling (MDS) plots done in an R Notebook. I will submit an HTML page of my completed R Notebook.

### MDS plot
```{r}
library(limma)
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")
```

### Make MDS plot interactive
```{r}
library(glimma)
glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), 
          groups=x$samples[,c(2,5)], launch=FALSE)
```

### HeatMap of log-CPM data
```{r}
library(gplots)
basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",
   labRow=v$genes$SYMBOL[i], labCol=group, 
   col=mycol, trace="none", density.info="none", 
   margin=c(8,6), lhei=c(2,10), dendrogram="column")
```



