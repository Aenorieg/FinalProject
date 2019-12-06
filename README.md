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

# Ashley E Noriega,
# Nov 13, 2019
# TRGN 510 Final Project: Milestone 1
# A script for setting up the Glimma Vignette
{
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")
library(limma)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Glimma") 
library(Glimma)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("edgeR")
library(edgeR)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Mus.musculus")
library(Mus.musculus)

library(R.utils)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("CAMERA")
library(CAMERA)
}

# Data Packaging

## Read in count data
```{r}
library(R.utils)
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb") 
utils::untar("GSE63310_RAW.tar", exdir = ".")
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
  "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
  "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
for(i in paste(files, ".gz", sep=""))
  R.utils::gunzip(i, overwrite=TRUE)
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt",
"GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt", 
   "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt", 
   "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", 
   "GSM1545545_JMS9-P8c.txt")
read.delim(files[1], nrow=5)
```
EntrezID
<int>
GeneLength
<int>
Count
<int>
497097	3634	1		
100503874	3259	0		
100038431	1634	0		
19888	9747	0		
20671	3130	1		
5 rows
data.frame
5 x 3
EntrezID
<int>
GeneLength
<int>
Count
<int>
497097	3634	1		
100503874	3259	0		
100038431	1634	0		
19888	9747	0		
20671	3130	1		
5 rows
  
## Read the 9 text files into R and combining into a matrix of counts
```{r}
library(edgeR)
x <- readDGE(files, columns=c(1,3))
class(x)
dim(x)
```
[1] "DGEList"
attr(,"package")
[1] "edgeR"
[1] 27179     9  
```{r}
names(x)
str(x)
```
[1] "samples" "counts" 
Formal class 'DGEList' [package "edgeR"] with 1 slot
  ..@ .Data:List of 2
  .. ..$ :'data.frame':	9 obs. of  4 variables:
  .. .. ..$ files       : chr [1:9] "GSM1545535_10_6_5_11.txt" "GSM1545536_9_6_5_11.txt" "GSM1545538_purep53.txt" "GSM1545539_JMS8-2.txt" ...
  .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1
  .. .. ..$ lib.size    : num [1:9] 32863052 35335491 57160817 51368625 75795034 ...
  .. .. ..$ norm.factors: num [1:9] 1 1 1 1 1 1 1 1 1
  .. ..$ : num [1:27179, 1:9] 1 0 0 0 1 431 768 4 810 452 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ Tags   : chr [1:27179] "497097" "100503874" "100038431" "19888" ...
  .. .. .. ..$ Samples: chr [1:9] "GSM1545535_10_6_5_11" "GSM1545536_9_6_5_11" "GSM1545538_purep53" "GSM1545539_JMS8-2" ...

## Annotate the samples
```{r}
x$samples
```
files
<chr>
group
<fctr>
lib.size
<dbl>
norm.factors
<dbl>
GSM1545535_10_6_5_11	GSM1545535_10_6_5_11.txt	1	32863052	1
GSM1545536_9_6_5_11	GSM1545536_9_6_5_11.txt	1	35335491	1
GSM1545538_purep53	GSM1545538_purep53.txt	1	57160817	1
GSM1545539_JMS8-2	GSM1545539_JMS8-2.txt	1	51368625	1
GSM1545540_JMS8-3	GSM1545540_JMS8-3.txt	1	75795034	1
GSM1545541_JMS8-4	GSM1545541_JMS8-4.txt	1	60517657	1
GSM1545542_JMS8-5	GSM1545542_JMS8-5.txt	1	55086324	1
GSM1545544_JMS9-P7c	GSM1545544_JMS9-P7c.txt	1	21311068	1
GSM1545545_JMS9-P8c	GSM1545545_JMS9-P8c.txt	1	19958838	1
9 rows
  
## Organize sample information
```{r}
library(edgeR)
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames
colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", 
                     "Basal", "ML", "LP"))
x$samples$group <- group
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane
x$samples
```
files
<chr>
group
<fctr>
lib.size
<dbl>
norm.factors
<dbl>
lane
<fctr>
10_6_5_11	GSM1545535_10_6_5_11.txt	LP	32863052	1	L004
9_6_5_11	GSM1545536_9_6_5_11.txt	ML	35335491	1	L004
purep53	GSM1545538_purep53.txt	Basal	57160817	1	L004
JMS8-2	GSM1545539_JMS8-2.txt	Basal	51368625	1	L006
JMS8-3	GSM1545540_JMS8-3.txt	ML	75795034	1	L006
JMS8-4	GSM1545541_JMS8-4.txt	LP	60517657	1	L006
JMS8-5	GSM1545542_JMS8-5.txt	Basal	55086324	1	L006
JMS9-P7c	GSM1545544_JMS9-P7c.txt	ML	21311068	1	L008
JMS9-P8c	GSM1545545_JMS9-P8c.txt	LP	19958838	1	L008
9 rows
  
## Organize gene annotations
```{r}
library(Mus.musculus)
geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENTREZID")
head(genes)
```
ENTREZID
<chr>
SYMBOL
<chr>
TXCHROM
<chr>
1	497097	Xkr4	chr1	
2	100503874	Gm19938	NA	
3	100038431	Gm10568	NA	
4	19888	Rp1	chr1	
5	20671	Sox17	chr1	
6	27395	Mrpl15	chr1	
6 rows
data.frame
6 x 3
 
 
ENTREZID
<chr>
SYMBOL
<chr>
TXCHROM
<chr>
1	497097	Xkr4	chr1	
2	100503874	Gm19938	NA	
3	100038431	Gm10568	NA	
4	19888	Rp1	chr1	
5	20671	Sox17	chr1	
6	27395	Mrpl15	chr1	
6 rows
## Resolve duplicate gene IDs
```{r}
library(Mus.musculus)
genes <- genes[!duplicated(genes$ENTREZID),]
```

## Package in a DGEList-object containing raw count data with associated sample information and gene annotations
```{r}
library(Mus.musculus)
x$genes <- genes
x
```

# Data Pre-processing

## Transformations from the raw-scale: convert raw counts to counts per million (CPM) and log2-counts per million (log-CPM)
```{r}
library(edgeR)
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm)
```

## Remove lowly expressed genes
```{r}
library(edgeR)
table(rowSums(x$counts==0)==9)
```

## Filter genes while keeping as many genes as possible with worthwile counts
```{r}
library(edgeR)
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
```

## Plot the density of log-CPM values for raw and filtered data
```{r}
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
den <- density(lcpm[,i])
lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
den <- density(lcpm[,i])
lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
```

## Normalize gene expression distributions
```{r}
library(edgeR)
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
```

## Improve visualization by duplicating data then adjusting the counts
```{r}
x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5
```

## Boxplot expression distribution of samples for normalised and unnormalised data
```{r}
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")
```

## Unsupervised clustering of cells: make multi-dimensional scaling plot (MDS) to show simmilarities and dissimilarities between samples in an unsupervised manner
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

## Make interactive using Glimma 
HTML page will be generarted and opened in a browser if launch=TRUE
```{r}
library(Glimma)
glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), 
          groups=x$samples[,c(2,5)], launch=FALSE)

```

# Differential Expression Analysis

## Create a design matrix
```{r}
design <- model.matrix(~0+group+lane)
colnames(design) <- gsub("group", "", colnames(design))
design
```

## Contrasts for pairwise comparisons between cell populations
```{r}
library(Glimma)
contr.matrix <- makeContrasts(
   BasalvsLP = Basal-LP, 
   BasalvsML = Basal - ML, 
   LPvsML = LP - ML, 
   levels = colnames(design))
contr.matrix
```

## Remove heteroscedascity from count data
```{r}
library(Glimma)
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v
```

## Apply voom precision weights to data
```{r}
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
```

## Examine the number of DE genes
```{r}
summary(decideTests(efit))
```

## Set a minimum log-fold change(log-FC) of 1
```{r}
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
```

## Extract genes that are DE in multiple comparisons
```{r}
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
head(tfit$genes$SYMBOL[de.common], n=20)
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))
```

## Extract and write results for all 3 comparisons (basalvsLP, basalvsML, and LPvsML) to a single output file
```{r}
write.fit(tfit, dt, file="results.txt")
```

## Examine individual DE genes from top to bottom
```{r}
basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)
head(basal.vs.lp)
```

```{r}
head(basal.vs.ml)
```

## Summarize results for genes using mean-difference plots that highlight differentially expressed genes
```{r}
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))
```

## Make interactive mean-difference plot 
To open HTML page in a browser, make launch=TRUE
```{r}
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],         side.main="ENTREZID", counts=lcpm, groups=group, launch=TRUE)
```

## Make heatmap
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

# Gene set testing by applying the camera method on c2 gene signatures from the Broad Institute’s MSigDB c2 collection
```{r}
load("/Users/ashleynoriega/Downloads/mouse_c2_v5p2.rdata")
idx <- ids2indices(Mm.c2,id=rownames(v))
cam.BasalvsLP <- camera(v,idx,design,contrast=contr.matrix[,1])
head(cam.BasalvsLP,5)
cam.BasalvsML <- camera(v,idx,design,contrast=contr.matrix[,2])
head(cam.BasalvsML,5)
cam.LPvsML <- camera(v,idx,design,contrast=contr.matrix[,3])
head(cam.LPvsML,5)
```


```{r}
barcodeplot(efit$t[,3], index=idx$LIM_MAMMARY_LUMINAL_MATURE_UP, 
            index2=idx$LIM_MAMMARY_LUMINAL_MATURE_DN, main="LPvsML")
```

# Ashley E Noriega,
# Nov 20, 2019
# TRGN 510 Final Project: Milestone 2
# Loading in RNA seq colon cancer data from TCGA 
## RNASeq files for cystic, mucinous, and serous neoplasms (CMS) and adenocarcinoma
I created a new folder called COAD and stored all 60 files in it. I then changed them to TXT files and opened them.

# Data Packaging

## The urls for the 30 files I downloaded locally for cystic, mucinous and serous colon cancer (CMS)

1. https://portal.gdc.cancer.gov/files/536f5a77-0087-457d-ac95-6d1a9abad8cb, UUID 536f5a77-0087-457d-ac95-6d1a9abad8cb, case: TCGA-AA-3516 

2. https://portal.gdc.cancer.gov/files/ed52de66-66fa-44ce-b679-cf641b0d92cd, UUID ed52de66-66fa-44ce-b679-cf641b0d92cd, case: TCGA-AA-3516 

3. https://portal.gdc.cancer.gov/files/b28090c5-c42d-4836-9bb1-ce906d3ead95, UUID: b28090c5-c42d-4836-9bb1-ce906d3ead95, case TCGA-AA-3854 

4. https://portal.gdc.cancer.gov/cases/57cdaa1c-4e94-4a28-ab3b-300c0457555f, UUID: 49e29c69-d9d7-4496-9f24-26f42c8b6d8e, case: TCGA-A6-2674 

5. https://portal.gdc.cancer.gov/files/08ed32e4-fb94-4bc0-8715-83ee2143a13d, UUID: 08ed32e4-fb94-4bc0-8715-83ee2143a13d, case: TCGA-AA-A00J 

6. https://portal.gdc.cancer.gov/files/6e571f71-d5fb-42f3-a35b-554c5ab76587, UUID: 6e571f71-d5fb-42f3-a35b-554c5ab76587, case: TCGA-AA-A01G 

7. https://portal.gdc.cancer.gov/files/8b12a000-f588-4a78-a9eb-f06041a65789, UUID: 8b12a000-f588-4a78-a9eb-f06041a65789, case: TCGA-A6-6780 

8. https://portal.gdc.cancer.gov/files/02734d4d-fc8f-4ef7-ac82-1b4d7184cc5e, UUID: 02734d4d-fc8f-4ef7-ac82-1b4d7184cc5e, case: TCGA-CK-4950 

9. https://portal.gdc.cancer.gov/files/6466a8b1-d1e2-4195-a353-0800576c13c8, UUID: 6466a8b1-d1e2-4195-a353-0800576c13c8, case: TCGA-G4-6322 

10. https://portal.gdc.cancer.gov/files/bc47f01c-1994-4ff8-a356-94d9679b66ee, UUID: bc47f01c-1994-4ff8-a356-94d9679b66ee, case: TCGA-AA-3947 

11. https://portal.gdc.cancer.gov/files/b045ee79-82a6-4636-a875-1a58603d89ff, UUID: b045ee79-82a6-4636-a875-1a58603d89ff, case: TCGA-A6-A566 

12. https://portal.gdc.cancer.gov/files/c383ba2c-b00a-4bd2-82cb-b3f04c2a8172, UUID: c383ba2c-b00a-4bd2-82cb-b3f04c2a8172, case: TCGA-AA-3877 

13. https://portal.gdc.cancer.gov/files/b52775aa-273e-484e-82c7-c625f09415fa, UUID: b52775aa-273e-484e-82c7-c625f09415fa, case: TCGA-A6-3809 

14. https://portal.gdc.cancer.gov/files/7b15a87a-805c-4b8a-84de-549cec9c44e3, UUID: 7b15a87a-805c-4b8a-84de-549cec9c44e3, case: TCGA-AA-3684 

15. https://portal.gdc.cancer.gov/files/b4f3dbbb-2686-4896-9e60-5bef6c9150b4, UUID: b4f3dbbb-2686-4896-9e60-5bef6c9150b4, case: TCGA-AA-3692 

16. https://portal.gdc.cancer.gov/files/0b16e2bd-3ec7-4901-9ff0-a389670e5019, UUID: 0b16e2bd-3ec7-4901-9ff0-a389670e5019, case: TCGA-D5-6534 

17. https://portal.gdc.cancer.gov/files/a6690007-f347-49c3-a0ba-28e01d131971, UUID: a6690007-f347-49c3-a0ba-28e01d131971, case: TCGA-A6-3809 

18. https://portal.gdc.cancer.gov/files/a1742cf6-c3c5-43e7-879c-489494460e78, UUID: a1742cf6-c3c5-43e7-879c-489494460e78, case: TCGA-AA-A00N 

19. https://portal.gdc.cancer.gov/files/d5be795d-beb6-4def-bda8-f485ee45bfc1, UUID: d5be795d-beb6-4def-bda8-f485ee45bfc1, case: TCGA-A6-2674 

20. https://portal.gdc.cancer.gov/files/46306072-c59c-4b4b-963c-9c4e778ff34b, UUID: 46306072-c59c-4b4b-963c-9c4e778ff34b, case: TCGA-A6-6780 

21. https://portal.gdc.cancer.gov/files/a938cb2c-c8e8-4395-915b-37e1e279a4da, UUID: a938cb2c-c8e8-4395-915b-37e1e279a4da, case: TCGA-G4-6302 

22. https://portal.gdc.cancer.gov/files/7fec7c90-fd2e-4ee2-ba1a-77f85920771f, UUID: 7fec7c90-fd2e-4ee2-ba1a-77f85920771f, case: TCGA-DM-A282 

23. https://portal.gdc.cancer.gov/files/2c3fd34c-70d1-4331-9628-260b77329b53, UUID: 2c3fd34c-70d1-4331-9628-260b77329b53, case: TCGA-F4-6704 

24. https://portal.gdc.cancer.gov/files/4168a720-521e-47ff-afb5-4abe3e815490, UUID: 4168a720-521e-47ff-afb5-4abe3e815490, case: TCGA-AA-3950 

25. https://portal.gdc.cancer.gov/files/ecc90bd1-f594-41ea-ba4b-d42f4c64880b, UUID: ecc90bd1-f594-41ea-ba4b-d42f4c64880b, case: TCGA-A6-6781 

26. https://portal.gdc.cancer.gov/files/8736ed27-2141-48d9-b677-b1a0e14d4b50, UUID: 8736ed27-2141-48d9-b677-b1a0e14d4b50, case: TCGA-CA-6717 

27. https://portal.gdc.cancer.gov/files/3b8d04cd-d658-46ba-adca-079fee531e17, UUID: 3b8d04cd-d658-46ba-adca-079fee531e17, case: TCGA-AA-3821 

28. https://portal.gdc.cancer.gov/files/b27da518-d023-4f9c-a9ab-5cd68ee37870, UUID: b27da518-d023-4f9c-a9ab-5cd68ee37870, case: TCGA-CK-4951 

29. https://portal.gdc.cancer.gov/files/e7005df6-f78b-4e47-abe7-61ae6a2ee026, UUID: e7005df6-f78b-4e47-abe7-61ae6a2ee026, case: TCGA-AA-A01R 

30. https://portal.gdc.cancer.gov/files/e3598d14-292c-41cc-9b59-4497fa078272, UUID: e3598d14-292c-41cc-9b59-4497fa078272, case: TCGA-D5-6930 


## The urls for the 30 files adenocarcinoma I downloaded locally for adenocarcinoma

1. https://portal.gdc.cancer.gov/files/f1185347-ad15-43ae-9ef3-d5343b31a0fc, UUID: f1185347-ad15-43ae-9ef3-d5343b31a0fc, case: TCGA-A6-6654

2. https://portal.gdc.cancer.gov/files/0d53cb1c-97c4-4088-9e43-029de88fd66d, UUID: 0d53cb1c-97c4-4088-9e43-029de88fd66d, case: TCGA-DM-A1D4

3. https://portal.gdc.cancer.gov/files/a74bbce0-7f3d-434e-b294-7fa45e5b3a60, UUID: a74bbce0-7f3d-434e-b294-7fa45e5b3a60, case: TCGA-A6-2684

4. https://portal.gdc.cancer.gov/files/47554e4e-cd13-4b92-80be-e1940f9a950f, UUID: 47554e4e-cd13-4b92-80be-e1940f9a950f, case: TCGA-A6-5657

5. https://portal.gdc.cancer.gov/files/de60dbd7-8a93-47a5-b1ea-a3f95beade8a, UUID: de60dbd7-8a93-47a5-b1ea-a3f95beade8a, case: TCGA-F4-6854

6. https://portal.gdc.cancer.gov/files/70883b31-d130-4efd-a7c6-169c8d4a253d, UUID: 70883b31-d130-4efd-a7c6-169c8d4a253d, case: TCGA-AD-A5EJ

7. https://portal.gdc.cancer.gov/files/042bda3d-77aa-4522-8a97-c121711a760e, UUID: 042bda3d-77aa-4522-8a97-c121711a760e, case: TCGA-AG-3582

8. https://portal.gdc.cancer.gov/files/b6388e09-7ed5-4041-97bb-4427ba5571ba, UUID: b6388e09-7ed5-4041-97bb-4427ba5571ba, case: TCGA-AY-6197

9. https://portal.gdc.cancer.gov/files/54394c0b-6ae3-4b48-8e89-350ad5349611, UUID: 54394c0b-6ae3-4b48-8e89-350ad5349611, case: TCGA-AA-3554

10. https://portal.gdc.cancer.gov/files/f7e21d61-19b6-4e99-887f-463d4419628c, UUID: f7e21d61-19b6-4e99-887f-463d4419628c, case: TCGA-AG-4015

11. https://portal.gdc.cancer.gov/files/b4114885-38cd-4e8a-874b-b78da8d95e2c, UUID: b4114885-38cd-4e8a-874b-b78da8d95e2c, case: TCGA-CM-6171

12. https://portal.gdc.cancer.gov/files/f9fda40d-67e4-4cb9-859c-ddc2ea84b7e4, UUID: f9fda40d-67e4-4cb9-859c-ddc2ea84b7e4, case: TCGA-CM-6170

13. https://portal.gdc.cancer.gov/files/b4aebb2a-d0b8-43d8-bd1f-78af2065d8f9, UUID: b4aebb2a-d0b8-43d8-bd1f-78af2065d8f9, case: TCGA-AA-3846

14. https://portal.gdc.cancer.gov/files/6a750710-5ed9-4d24-b2bf-3a4e3211878f, UUID: 6a750710-5ed9-4d24-b2bf-3a4e3211878f, case: TCGA-CM-6677

15. https://portal.gdc.cancer.gov/files/93d1a78f-423e-4560-b4d3-ee4a89ac922b, UUID: 93d1a78f-423e-4560-b4d3-ee4a89ac922b, case: TCGA-RU-A8FL

16. https://portal.gdc.cancer.gov/files/2e632fd9-fa17-4290-9601-a5d462cf152c, UUID: 2e632fd9-fa17-4290-9601-a5d462cf152c, case: TCGA-AZ-4323

17. https://portal.gdc.cancer.gov/files/7239b026-2587-489d-81fe-7bc657b7523c, UUID: 7239b026-2587-489d-81fe-7bc657b7523c, case: TCGA-CM-6164

18. https://portal.gdc.cancer.gov/files/90e86a26-fffa-4c38-b2e0-bf0704ee3615, UUID: 90e86a26-fffa-4c38-b2e0-bf0704ee3615, case: TCGA-AZ-4315

19. https://portal.gdc.cancer.gov/files/9ff11fe0-037c-405e-95c3-dc4a15413db8, UUID: 9ff11fe0-037c-405e-95c3-dc4a15413db8, case: TCGA-G4-6311

20. https://portal.gdc.cancer.gov/files/b8eed826-6051-4358-9b3d-44d1553dd9ad, UUID: b8eed826-6051-4358-9b3d-44d1553dd9ad, case: TCGA-AA-3522

21. https://portal.gdc.cancer.gov/files/c172bc07-d4f0-41be-a558-49abc81065c2, UUID: c172bc07-d4f0-41be-a558-49abc81065c2, case: TCGA-AA-3667

22. https://portal.gdc.cancer.gov/files/260edc5e-1ca6-4b07-b96d-59594d03ac54, UUID: 260edc5e-1ca6-4b07-b96d-59594d03ac54, case: TCGA-AA-A00U

23. https://portal.gdc.cancer.gov/files/0c5c1a38-7e9c-4b43-810d-0761c3af49b1, UUID: 0c5c1a38-7e9c-4b43-810d-0761c3af49b1, case: TCGA-AA-3506

24. https://portal.gdc.cancer.gov/files/7024ba0c-be56-4907-9254-cdb2579e536e, UUID: 7024ba0c-be56-4907-9254-cdb2579e536e, case: TCGA-NH-A8F7

25. https://portal.gdc.cancer.gov/files/031cf2a5-74e0-4b5f-98bd-da60628c0854, UUID: 031cf2a5-74e0-4b5f-98bd-da60628c0854, case: TCGA-AA-3680

26. https://portal.gdc.cancer.gov/files/91991ecf-cc54-4110-8a4e-9236bf8aa072, UUID: 91991ecf-cc54-4110-8a4e-9236bf8aa072, case: TCGA-A6-4105

27. https://portal.gdc.cancer.gov/files/47aceec1-a01d-419f-9689-c46284c79bcb, UUID: 47aceec1-a01d-419f-9689-c46284c79bcb, case: TCGA-D5-6922

28. https://portal.gdc.cancer.gov/files/ce84c955-63db-473a-a6d7-0e3daad6efd4, UUID: ce84c955-63db-473a-a6d7-0e3daad6efd4, case: TCGA-AA-3524

29. https://portal.gdc.cancer.gov/files/a071fc45-61ea-4815-93bf-be34980e59ee, UUID: a071fc45-61ea-4815-93bf-be34980e59ee, case: TCGA-AA-3855

30. https://portal.gdc.cancer.gov/files/8b275144-b885-4fb0-af39-fea1e48a970a, UUID: 8b275144-b885-4fb0-af39-fea1e48a970a, case: TCGA-AA-A00Q

## Load in count data
First rename the ".count" files to ".txt" and unzip each one by opening each file.
```{r}
setwd('~/Desktop/COAD_Data/')
COAD_files <- c("9e8b528b-1172-4c07-a09b-ebb23cf2310c.htseq.txt", "bda1a9a4-a14f-4463-81d2-a4fcca65d6f1.htseq.txt", 
   "5697212f-b3fd-479f-84b0-ec0aae54534a.htseq.txt", "7f9a629b-12ed-48cc-8d5c-1c2f5db9cf1f.htseq.txt",
   "15864159-be88-41c8-bdef-c2c5927cb1a1.htseq.txt", "649b19e1-96e2-4b55-951d-3b6ee9f4b91f.htseq.txt",
   "86679663-dfc5-46ad-8cf9-c7954c4b339b.htseq.txt", "28004569-048d-4f8c-99aa-7a8c69a98dcc.htseq.txt",                 "911f6378-8a25-4570-9d3b-80f5b5bfc085.htseq.txt", "d5dca54e-d7e9-4328-b2ca-1d191a2b8b4c.htseq.txt",                 "f3895ae4-1228-49b3-9342-3c3b86cb5243.htseq.txt", "f590941d-19dc-427a-95b6-942c97ea8333.htseq.txt",                 "55aa6d16-3598-42ca-8844-0fe84739ef66.htseq.txt", "0e7094cf-4c79-43f4-8b72-9de259e5e18f.htseq.txt",                 "9a62fe1f-36ec-4e8e-b3d9-bdfc62f71905.htseq.txt", "d2587070-cb7d-440d-ae49-52f5077248e6.htseq.txt",                 "7800bdb2-aa8b-43e0-8e45-1b968872b34e.htseq.txt", "2bcd2efd-4fd6-40ee-86a4-867ae82711b0.htseq.txt",                 "424d8e5f-9fc6-470b-ad2c-b4447b0eb07e.htseq.txt", "934f9dc6-1260-4268-b022-870f1e37dd6f.htseq.txt",                 "0fa55c0e-6f8f-44a6-82fe-9a42495d3484.htseq.txt", "c8544a8a-4352-438d-94d4-3495af2e9a78.htseq.txt",                 "dade0b16-ecc3-43b3-b328-3819a8fc18c6.htseq.txt", "e875ae4e-4645-4e84-b0ff-9c9a694717a9.htseq.txt",                 "debd6982-7c27-42e8-b778-20afcc78a5f3.htseq.txt", "17c88994-9e8e-4f16-8c41-34e98a0d8c52.htseq.txt",                 "7f5a924a-ddf3-45ff-be1f-5b5909305f46.htseq.txt", "fa73bdce-67fb-42aa-883f-635f0e7bcdc6.htseq.txt",                 "abe20df7-6b97-4397-8864-881bac27e92c.htseq.txt", "62f84581-4c7d-4c8e-835c-9304bcec3106.htseq.txt", "3abbd2b5-04db-4fe0-8dd1-ea2b48caa4c1.htseq.txt", "087666cd-47ae-4f56-b947-d6aa1c25e8a7.htseq.txt", 
   "c14f98e2-8e9b-49f4-a244-3d06c6cb7126.htseq.txt", "13abc91e-fbfc-4c55-bf54-fbd134979ccc.htseq.txt",
   "6ae2dd6c-2a39-411f-a1fc-11e0e6e82165.htseq.txt", "8f77f4f4-b184-40c7-8ab8-2f95b13620b5.htseq.txt",
   "168e5cb2-7390-45ad-ad04-c9aa4416e950.htseq.txt", "0ed65bdf-cb92-47c1-8aeb-42518ce639b8.htseq.txt",                 "4e7c6811-88e4-4bb7-a88f-7491dfa6d072.htseq.txt", "7fb73a84-867a-4c28-aa02-93068efffb7b.htseq.txt",                 "b53f9a9d-b24d-410d-b3e9-f2a8bf22ca27.htseq.txt", "f7ce175f-763e-4a55-97e3-0381d889b0eb.htseq.txt",                 "f346f2d2-285c-455c-ba34-ea8eec3fa881.htseq.txt", "e53e1a83-1979-4e12-bbb7-79b37d0cfe03.htseq.txt",                 "a26d49db-2309-46a0-a3ed-275378d484e7.htseq.txt", "a3f88a5d-7169-465b-bb80-e5999590681c.htseq.txt",                 "c264fe3b-482b-44ec-83a4-73df565663ff.htseq.txt", "bd2dfab3-88a8-4673-ba36-3daf252d0b4d.htseq.txt",                 "7261b656-c79c-4581-a503-15b653e2b5d2.htseq.txt", "ee4dcccc-514b-4cc6-ae63-6ed3e7519a40.htseq.txt",                 "f596eabc-e39a-4e35-9fc6-edade04eb785.htseq.txt", "bf9c448b-bdc9-4f74-b13a-374e6add7939.htseq.txt",                 "564daa81-cfef-45b6-94a0-3249b2724d9b.htseq.txt", "82e00e45-734c-471f-ba97-79ec3b7e0baa.htseq.txt",                 "9c52ed00-325f-4664-8873-327bcaa5ea74.htseq.txt", "fabefb10-5546-4017-8ea1-29982a10fb3c.htseq.txt",                 "32a115cf-570f-4ad9-a123-8e1970062f51.htseq.txt", "05eef9f8-a246-403a-b0be-07d274b6f93a.htseq.txt",                 "5c18c6a8-9ad2-43a8-a3a0-83d8fc0cc257.htseq.txt", "43b292be-5d63-4523-a43f-666d20039208.htseq.txt")
read.delim(COAD_files[1], nrows = 60)
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/Screen%20Shot%202019-12-05%20at%202.15.07%20PM.png)

## Create dataset, join the 60 loaded txt files
Use edgeR to create a matrix of 60 text files.

#### Known issue: Working Directory
Spoke to professor Craig on 12/4 and it is ok to ignore the error message and not change root, just setwd to desktop since files were downloaded locally.
```{r}
setwd('~/Desktop/COAD_Data/')
library(edgeR)
x <- readDGE(COAD_files, columns=c(1,2)) #joins my 60 files and creates a dataset
class(x)
dim(x)
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/class(x).png)

```{r}
names(x) #accessor function  
str(x) #displays the structure of x in compact way, alternative to summary and best for displaying contents of lists
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/names(x).png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/str(x).png)

## Annotate the samples
```{r}
x$samples
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/x%24samples.png)

# Organize sample information

## Associate sample-level information with the columns of the counts matrix
```{r}
samplenames <- substring(colnames(x), 1, nchar(colnames(x)))
samplenames
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/samplenames.png)

## Specify which files are Cystic, Mucinous, and Serous (CMS) and which files are Adenocarcinoma
```{r}
colnames(x) <- samplenames
group <- as.factor(c("CMS", "CMS", "CMS", "CMS", "CMS", "CMS",
                     "CMS", "CMS", "CMS", "CMS", "CMS", "CMS",
                     "CMS", "CMS", "CMS", "CMS", "CMS", "CMS",
                     "CMS", "CMS", "CMS", "CMS", "CMS", "CMS",
                     "CMS", "CMS", "CMS", "CMS", "CMS", "CMS",
                     "ADENOCARCINOMA", "ADENOCARCINOMA", "ADENOCARCINOMA", "ADENOCARCINOMA",                                 "ADENOCARCINOMA", "ADENOCARCINOMA", "ADENOCARCINOMA", "ADENOCARCINOMA",
                     "ADENOCARCINOMA", "ADENOCARCINOMA", "ADENOCARCINOMA", "ADENOCARCINOMA",
                     "ADENOCARCINOMA", "ADENOCARCINOMA", "ADENOCARCINOMA", "ADENOCARCINOMA",
                     "ADENOCARCINOMA", "ADENOCARCINOMA", "ADENOCARCINOMA", "ADENOCARCINOMA",
                     "ADENOCARCINOMA", "ADENOCARCINOMA", "ADENOCARCINOMA", "ADENOCARCINOMA",
                     "ADENOCARCINOMA", "ADENOCARCINOMA", "ADENOCARCINOMA", "ADENOCARCINOMA",
                     "ADENOCARCINOMA", "ADENOCARCINOMA"))

x$samples$group <- group
x$samples
DF<-x$samples #for my own visualization purposes
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/specified%20samples.png)

# Script to organize gene annotations
{
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Homo.sapiens")
library(Homo.sapiens)
install.packages(gsubfn)
library(gsubfn)
}

## Annotate Genes
First install Homo.sapiens, then use a script remove the decimals and numbers after the decimal points in all 60487 ENSEMBL geneid elements.
```{r}
library(Homo.sapiens)
#library(stringr)
library(gsubfn)
geneid <- rownames(x)
#geneid_test <- c("ENSG00000000005", 
#	"ENSG00000000419",
#	"ENSG00000000457",
#	"ENSG00000000938") 
#geneid <- str_remove(geneid, "[.]") removes decimals only
geneid <- gsub("\\.[0-9]*$", "", geneid) #remove decimals and numbers after decimals
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENSEMBL")
head(genes)
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/genes.png)

## Remove duplicate genes
```{r}
genes <- genes[!duplicated(genes$ENSEMBL),]
```

## Package in a DGEList-object containing raw count data with associated sample information and gene annotations
```{r}
x$genes <- genes
x
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/R%20console.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/data.frame%205x4.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/data.frame%205x3.png)

# Ashley E Noriega
# Nov 30, 2019
# TRGN 510 Final Project: Milestone 3
# Running the Glimma Vignette

# Data Pre-processing

## Transformations from the raw-scale: convert raw counts to counts per million (CPM) and log2-counts per million (log-CPM)
```{r}
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm)
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/lcpm1.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/lcpm2.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/lcpm3.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/lcpm4.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/lcpm5.png)

## Remove lowly expressed genes
True signifies how many genes have counts equal to zero, meaning genes are unexpressed throughout all samples. 
```{r}
table(rowSums(x$counts==0)==9)
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/Counts%3D0.png)

## Filter genes while keeping as many genes as possible with worthwile counts
```{r}
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/dim(x).png)

## Plot the density of log-CPM values for raw and filtered data
There is a sample that is a potential outlier (green colored line), could remove the sample for future analysis but spoke to porfessor Craig on 12/4 and agreed to leave the sample in since the vignette has a normalisation step.

#### Known issue: color palette
Spoke to professor Craig on 12/4 and agreed to stop working on this issue. I understand that the "Paired" palatte only offers 12 colors so every 13th sample repeats color scheme. I tried increasing the number of colors available with colorRampPalatte but was unsuccesful.
```{r}
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
#library(colorRamps)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired") #results in the error message: n too large, allowed maximum for palette Paired is 12. Returning the palette you asked for with that many colors
#nb.cols = 60
#col <- colorRampPalette(brewer.pal(nsamples, "Paired"))(nb.cols) #colorRampPalette is a constructor function that builds palettes with arbitrary number of colors by interpolating existing palette 
par(mfrow=c(1,2)) #1 row, 2 columns
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
den <- density(lcpm[,i])
lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
den <- density(lcpm[,i])
lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/Raw%20and%20Filtered%20data.png)

## Normalising gene expression distributions
```{r}
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/normalize%20gene%20expression%20distribution.png)

## Improve visualization by duplicating data then adjusting the counts
```{r}
x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5
```

## Boxplot expression distribution of samples for unnormalised data
```{r}
par(mfrow=c(1,1)) #makes boxplot look less cramped 
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/Unnormalised%20boxplot.png) #what boxplot looked like in previous milestone

![image](https://github.com/Aenorieg/FinalProject/blob/master/new%20unormalised%20boxplot.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/x2.png)

## Boxplot expression distribution of samples for normalised data
This step forces the samples to even out, may not be a good thing since there is a potential outlier.
```{r}
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/normalised%20boxplot.png)

## Unsupervised clustering of cells: make multi-dimensional scaling plot (MDS) to show simmilarities and dissimilarities between samples in an unsupervised manner

#### Known issue: color palette 
I spoke to professor Craig on 12/4, ok to ignore error since I am only comparing 2 different subsets of colon cancer. To get rid of this error I would need to add an additional factor: lane.
```{r}
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,1)) #1 row, 1 column 
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1") #n= number of different colors in a palette with the min being 3 
col.group <- as.character(col.group)
#col.lane <- lane did not have lanes for my data
#levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
#col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
#plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
#title(main="B. Sequencing lanes")
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/MDS.png)

## Make interactive using Glimma
HTML page will be generarted and opened in a browser if launch=TRUE
```{r}
library(Glimma)
glMDSPlot(lcpm, labels=paste(group, sep="_"), 
          groups=x$samples[,c(1,2)], launch=TRUE)
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/Interactive%20MDS.png)

# Differential expression analysis

## Creating a design matrix
```{r}
design <- model.matrix(~0+group) #removes intercept from the factor group
#design <- model.matrix(~group) leaves intercept from factor group, but model contrasts are more straight forward without intercept
colnames(design) <- gsub("group", "", colnames(design))
design
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/Design%20Matrix%201.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/Design%20Matrix%202.png)

## Contrasts for pairwise comparisons between cell populations
Since I am only comparing CMS and Adenocarcinoma, I will only have 1 pairwise comparison.
```{r}
library(limma)
contr.matrix <- makeContrasts(
   ADENOCARCINOMAvsCMS = ADENOCARCINOMA-CMS, 
   levels = colnames(design))
contr.matrix
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/pairwise%20comparisons.png)

## Remove heteroscedascity from count data

### Voom plot
Each black dot represents a gene. The red curve is the estimated mean-varience trend used to compute the voom weights.
```{r}
 par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE) #voom converts raw counts to log-CPM values by extracting library sizes and normalisation factors from x
v
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/voom.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/v%201.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/v%202.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/v%203.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/v%204.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/v%205.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/v%206.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/v%207.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/v%208.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/v%209.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/v%2010.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/v%2011.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/v%2012.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/v%2013.png)

### Apply voom precision weights to data
Each black dot is a gene. The blue line is the average log2 residual standard deviation computed with the Bayes algorithm.
```{r}
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend") #plots log2 residual standard deviations against mean log-CPM values
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/voom%20precision.png)

## Examine the number of DE genes
Quick view at how many genes are down-regulated, up-regulated, and not statistically significant. The adjusted p-value cutoff is 5% by default.
```{r}
summary(decideTests(efit))
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/efit.png)

## Set a minimum log-fold change(log-FC) of 1
This is a stricter definition of significance and could be overcorrecting since now I don't have any down-regulated or up-regulated genes.
```{r}
tfit <- treat(vfit, lfc=1) #p-values calculated from empirical Bayes moderated t-statistics with a minimum log-FC requirement.
dt <- decideTests(tfit)
#dt <- decideTests(efit) #for testing purposes
summary(dt)
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/vfit.png)

## Extract genes that are DE in multiple comparisons
I don't have any DE genes if tfit is used. If efit is used, I have 3295 DE genes.
```{r}
de.common <- which(dt[,1]!=0)
length(de.common)
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/de.common.png)

## The first 20 DE genes
If efit is used the genes are: "DPM1", "CFH", "LAS1L", "CFTR", "TMEM176A", "DBNDD1", "TFPI", "SLC7A2", "ARF5", "POLDIP2", "ARHGAP33", "UPF1", "MCUB", "POLR2J", "THSD7A", "LIG3", "SPPL2B", "IBTK", "PDK2", "REX1BD" 
```{r}
head(tfit$genes$SYMBOL[de.common], n=20)
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/head%20de.common.png)

## Make Venn Diagram
My diagram only has 1 circle because I only have 1 pairwise comparison.
```{r}
vennDiagram(dt[,1], circle.col=c("turquoise", "salmon"))
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/Venn%20Diagram.png)

## Extract and write results for comparisons of ADENOCARCINOMAvsCMS to a single output file
```{r}
write.fit(tfit, dt, file="results.txt")
```

## Examining individual DE genes from top to bottom
```{r}
ADENOCARCINOMA.vs.CMS <- topTreat(tfit, coef=1, n=Inf)
head(ADENOCARCINOMA.vs.CMS)
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/DE%20genes%20top%20to%20bottom.png)

## Summarize results for genes using mean-difference plots that highlight differentially expressed genes
If efit is used, will have read, black and blue genes. Since tfit is used, all genes are black.
```{r}
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/mean-difference%20plot%20w:%20vfit.png)

## Make interactive mean-difference plot 
To open html page in a browser make launch=TRUE
```{r}
library(Glimma)
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         side.main="ENSEMBL", counts=lcpm, groups=group, launch=TRUE)
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/interactive%20mean-difference%20plot.png)

![image](https://github.com/Aenorieg/FinalProject/blob/master/interactive%20mean-difference%20if%20efit%20used.png) #if tfit was not applied

## Make heatmap
Install heatmap.plus beacuse heatmap.2 did not work for my data.
```{r}
library(gplots)
library(heatmap.plus)
ADENOCARCINOMA.vs.CMS.topgenes <- ADENOCARCINOMA.vs.CMS$ENSEMBL[1:100]
i <- which(v$genes$ENSEMBL %in% ADENOCARCINOMA.vs.CMS.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
#par("mar") OUTPUT SHOULD BE [1] 5.1 4.1 4.1 2.1
par(cex.main=0.8,mar=c(1,1,1,1)) #mar=c(1,1,1,1) ensures margins are large enough
heatmap.plus(lcpm[i,], col=bluered(20),cexRow=1,cexCol=0.2, margins = c(10,10), main = "HeatMap") #changed the margins to have a more legible heatmap
```
![image](https://github.com/Aenorieg/FinalProject/blob/master/heatmap.png) #heatmap in prior milestone

![image](https://github.com/Aenorieg/FinalProject/blob/master/new%20heatmap.png)

## Gene set testing with Camera
I spoke to professor Craig on 12/4 and agreed that this step would not work for me becuase I did not have differentlially expressed genes.
