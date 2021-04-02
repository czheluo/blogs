---
layout: single
author_profile: false
comments: true
enable_mathjax: true
output: html_document
title : "scRNA-seq: Quality Control"
toc: true
tags: [Quality Control, QC, scRNA-seq, data sciences]
---

<iframe width="560" height="315" src="https://www.youtube.com/embed/RS5xoGE2yZ4" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/latest.js?config=TeX-MML-AM_CHTML">
</script>

##  Intro Single-cell RNA-seq
### Why Single-cell RNA-seq？
- 为了更好的了解组织和存在的细胞类型，需要更高分辨率的技术
- scRNA-seq提供了在单个细胞水平上表达哪些基因的信息
- 探索组织中存在哪些细胞类型
- 识别未知/稀有的细胞类型或状态
- 阐明分化过程中或跨时间或不同状态下的基因表达变化
- 识别在特定条件下（例如，治疗或疾病）在特定细胞类型中差异表达的基因

### Common applications of single-cell RNA sequencing

- 细胞异质性研究：能够鉴定细胞亚型和稀有细胞类型
- 细胞状态转变的轨迹分析：鉴定谱系特异性基因表达和驱动分支的关键基因
- 解剖转录动力学：转录爆发，基因在每个细胞中的打开和关闭
- 网络推断：推断模块，共同调节的基因-推断基因调节网络


<p align="center">
<img src="/images/Blog/QC/application.png" width="600">
</p>


ref: *Liu S, Trapnell C. Single-cell transcriptome sequencing: recent advances and remaining challenges. F1000Res. 2016 Feb 17;5:F1000 Faculty Rev-182. doi: 10.12688/f1000research.7223.1. PMID: 26949524; PMCID: PMC4758375.*

### Challenges(complexities) of scRNA-seq analysis
#### Large volume of data (high dimension)
- Expression data from scRNA-seq experiments represent ten or hundreds of thousands of reads for thousandsof cells. The data output is much larger, requiring higher amounts of memory to analyze, larger storage requirements, and more time to run the analyses

#### Low depth of sequencing per cell
- For the droplet-based methods of scRNA-seq, the depth of sequencing is shallow, often detecting only 10- 50% of the transcriptome per cell. This results in cells showing zero counts for many of the genes. However, in a particular cell, a zero count for a gene could either mean that the gene was not being expressed or the transcripts were just not detected. Across cells, genes with higher levels of expression tend to have fewer zeros. Due to this feature, many genes will not be detected in any cell and gene expression will be highly variable between cells.

#### Biological variability across cells/samples
- Transcriptional bursting: Gene transcription is not turned on all of the time for all genes. Time of harvest will determine whether gene is on or off in each cell.
- Varying rates of RNA processing: Different RNAs are processed at different rates.
- Continuous or discrete cell identities : Continuous phenotypes are by definitition variable in gene expression, and separating the continuous from the discrete can sometimes be difficult.
- Environmental stimuli: The local environment of the cell can influence the gene expression depending on spatial position, signaling molecules, etc.
- Temporal changes: Fundamental fluxuating cellular processes, such as cell cycle, can affect the gene expression profiles of individual cells.




<p align="center">
<img src="/images/Blog/QC/challenge.png" width="600">
</p>
ref: *Wagner A, Regev A, Yosef N. Revealing the vectors of cellular identity with single-cell genomics. Nat Biotechnol. 2016 Nov 8;34(11):1145-1160. doi: 10.1038/nbt.3711. PMID: 27824854; PMCID: PMC5465644.*


#### Technical variability across cells/samples
- Cell-specific capture efficiency: Different cells will have differing numbers of transcripts captured resulting in differences in sequencing depth (e.g. 10-50% of transcriptome).
- Library quality: Degraded RNA, low viability/dying cells, lots of free floating RNA, poorly dissociated cells, and inaccurate quantitation of cells can result in low quality metrics
- Amplification bias: During the amplification step of library preparation, not all transcripts are amplified to the same level.
- Batch effects: Batch effects are a significant issue for scRNA-Seq analyses, since you can see significant differences in expression due solely to the batch effect.

<p align="center">
<img src="/images/Blog/QC/batch.png" width="600">
</p>

ref:*Stephanie C Hicks, F William Townes, Mingxiang Teng, Rafael A Irizarry, Missing data and technical variability in single-cell RNA-sequencing experiments, Biostatistics,October 2018.*

#### While scRNA-seq is a powerful and insightful method for the analysis of gene expression with single-cell resolution, there are many challenges and sources of variation that can make the analysis of the data complex or limited.


## Workflows

<p align="center">
<img src="/images/Blog/QC/scRNA.png" width="600">
</p>

### 10x Genomics
- 10X Chromium 单细胞转录组测序可分为3'端polyA 附近区域捕获和5'端转录起始位置附近捕获建库测序。3'端转录本测序适用于各种类型的细胞，对于10X 单细胞3'的V2 试剂盒处理的样本，Read 1 由16 bp 的10X 细胞Barcode 和10 bp 的UMI 序列组成；而V3 试剂盒处理的样本，Read 1 由16 bp 的10X 细胞Barcode 和12 bp 的UMI 序列组成。其中，10X chromium 的Barcode 用于标记单个细胞，存在于逆转录引物上的随机核苷酸序列上。Read 2 是151 bp 的cDNA 序列, 一般只将前98 bp 用于下游分析。



<p align="center">
<img src="/images/Blog/QC/10x.png" width="600">
</p>
<p align="center">
<img src="/images/Blog/QC/10xAdapter.png" width="600">
</p>
[10x Genomics methods](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html)

- https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3fb.htm

<p align="center">
<img src="/images/Blog/QC/UMI.png" width="600">
</p>

#### 液滴方法:
Sample index(样本索引)：确定read来自哪个样本(在库准备过程中添加—需要记录)
Cellular barcode：确定read来自哪个细胞(每种库制备方法都有在库制备过程中使用的细胞条形码的库)
UMI(唯一分子标识符)：确定read来自哪个转录分子
Sequencing read1：Read1序列
Sequencing read2：Read2序列

### Assessing the quality metrics

- Cell counts
- UMI counts per cell
- Genes detected per cell
- UMIs vs. genes detected
- Mitochondrial counts ratio
- doublets: doublets are generated from two cells. They typically arise due to errors in cell sorting or capture, especially in droplet-based protocols involving thousands of cells. Doublets are obviously undesirable when the aim is to characterize populations at the single-cell level. 
*对于大多数scRNA­seq方法，doublets的产生，即当两个或两个以上的细胞被分配到相同的细胞条形码，会在下游分析中产生假的cluster，因为合并两种不同细胞类型的基因表达模式可能会产生一种独特的表达特征，这在任何真正的细胞类型中都找不到。然而，手动的从真实的cluster中区分出doublet cluster会有很大的挑战，特别是对有多种细胞类型的大数据集。识别doublets的一种常见策略是通过组合数据集中不同cluster中的细胞来生成模拟doublets，并评估哪些细胞具有与模拟doublets细胞相似的表达谱。然而，这种策略只有在数据集包含离散的细胞类型而不是连续的细胞轨迹时才可行。*
*在单细胞RNA测序实验中，双胞体是由两个细胞产生的。它们通常是由于细胞分选或捕获中的错误引起的，特别是在涉及数千个细胞的基于液滴的协议中。当目标是描述单细胞水平的群体特征时，双峰显然是不可取的。具体地说，他们可能错误地暗示存在实际并不存在的中间群体或短暂状态。因此，需要移除双峰文库，以便它们不会影响结果的解释。*


<p align="center">
<img src="/images/Blog/QC/qm.png" width="600">
</p>


### filtering in different ways and exploring variablility data
#### [3 different PBMC datasets from the 10x Genomics](https://support.10xgenomics.com/single-cell-gene-expression/datasets) 	
- 1k PBMCs using 10x v2 chemistry
- 1k PBMCs using 10x v3 chemistry
- 1k PBMCs using 10x v3 chemistry in combination with cell surface proteins, but disregarding the protein data and only looking at gene expression.

#### QC-features :nFeature_RNA, nCount_RNA  
```r
suppressMessages(require(Seurat))
suppressMessages(require(scater))
suppressMessages(require(Matrix))
#Download data
#curl -O http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_filtered_feature_bc_matrix.h5
#curl -O http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5
#curl -O http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5
v3.1k <- Read10X_h5("pbmc_1k_v3_filtered_feature_bc_matrix.h5", use.names = T)
v2.1k <- Read10X_h5("pbmc_1k_v2_filtered_feature_bc_matrix.h5", use.names = T)
p3.1k <- Read10X_h5("pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5", use.names = T)
#First, create Seurat objects for each of the datasets, and then merge into one large seurat object.

sdata.v2.1k <- CreateSeuratObject(v2.1k, project = "v2.1k")
sdata.v3.1k <- CreateSeuratObject(v3.1k, project = "v3.1k")
sdata.p3.1k <- CreateSeuratObject(p3.1k, project = "p3.1k")

# merge into one single seurat object. Add cell ids just in case you have overlapping barcodes between the datasets.
alldata <- merge(sdata.v2.1k, c(sdata.v3.1k,sdata.p3.1k), add.cell.ids=c("v2.1k","v3.1k","p3.1k"))


# also add in a metadata column that indicates v2 vs v3 chemistry
chemistry <- rep("v3",ncol(alldata))
chemistry[Idents(alldata) == "v2.1k"] <- "v2"
alldata <- AddMetaData(alldata, chemistry, col.name = "Chemistry")
alldata


## An object of class Seurat
## 33538 features across 2931 samples within 1 assay
## Active assay: RNA (33538 features)
# check number of cells from each sample, is stored in the orig.ident slot of metadata and is autmatically set as active ident.
table(Idents(alldata))

##
## p3.1k v2.1k v3.1k
##   713   996  1222


##Calculate mitochondrial proportion

#Seurat automatically calculates some QC-stats,
#like number of UMIs and features per cell. Stored in columns nCount_RNA & nFeature_RNA of the metadata.

head(alldata@meta.data)


#Manually calculate the proportion of mitochondrial reads and add to the metadata table.

mt.genes <- rownames(alldata)[grep("^MT-",rownames(alldata))]
C<-GetAssayData(object = alldata, slot = "counts")

percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
alldata <- AddMetaData(alldata, percent.mito, col.name = "percent.mito")


#plotQC
VlnPlot(alldata, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
VlnPlot(alldata, features = "nCount_RNA", pt.size = 0.1) + NoLegend()

VlnPlot(alldata, features = "percent.mito", pt.size = 0.1) + NoLegend()

VlnPlot(alldata, features = "percent.ribo", pt.size = 0.1) + NoLegend()

#And we can plot the different QC-measures as scatter plots


FeatureScatter(alldata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

FeatureScatter(alldata, feature1 = "nFeature_RNA", feature2 = "percent.mito")

FeatureScatter(alldata, feature1="percent.ribo", feature2="nFeature_RNA")


p<-FeatureScatter(alldata, feature1="percent.ribo", feature2="percent.mito",
                  cells = WhichCells(alldata, expression = orig.ident == "v3.1k"))
FeatureScatter(alldata, feature1="percent.ribo", feature2="percent.mito",
               cells = WhichCells(alldata, expression = orig.ident == "v2.1k"))
ggsave(p,filename ="ribo.mito.pdf" )

FeatureScatter(alldata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
               cells = WhichCells(alldata, expression = orig.ident == "v2.1k") )


```



<p align="center">
<img src="/images/Blog/QC/qcfc.png" width="600">
</p>

*The v2 chemistry gives lower gene detection, but higher detection of ribosomal proteins. As the ribosomal proteins are highly expressed they will make up a larger proportion of the transcriptional landscape when fewer of the lowly expressed genes are detected.*

- *V2给出了较低的基因检测，但较高的检测核糖体蛋白。由于核糖体蛋白是高表达的，当低表达的基因较少时，它们将在转录占更大的比例*

##### correlation

```r


```

<p align="center">
<img src="/images/Blog/QC/qcc.png" width="600">
</p>

##### Calculate mitochondrial, ribosomal proportion

```r

#select cells with percent.mito < 25
selected <- WhichCells(alldata, expression = percent.mito < 25)
length(selected)
# and subset the object to only keep those cells
data.filt <- subset(alldata, cells = selected)
data.filt <- subset(alldata, subset = percent.mito < 25)
# plot violins for new data
VlnPlot(data.filt, features = "percent.mito")

#start with cells with many genes detected.
high.det.v3 <- WhichCells(data.filt, expression = nFeature_RNA > 4100)
high.det.v2 <- WhichCells(data.filt, expression = nFeature_RNA > 2000 & orig.ident == "v2.1k")

# remove these cells
data.filt <- subset(data.filt, cells=setdiff(WhichCells(data.filt),c(high.det.v2,high.det.v3)))

# check number of cells
ncol(data.filt)


#Filter the cells with low gene detection (low quality libraries) with less than 1000 genes for v2 and < 500 for v2.

#start with cells with many genes detected.
low.det.v3 <- WhichCells(data.filt, expression = nFeature_RNA < 1000 & orig.ident != "v2.1k")
low.det.v2 <- WhichCells(data.filt, expression = nFeature_RNA < 500 & orig.ident == "v2.1k")

# remove these cells
data.filt <- subset(data.filt, cells=setdiff(WhichCells(data.filt),c(low.det.v2,low.det.v3)))

# check number of cells
ncol(data.filt)


#Plot QC-stats again
VlnPlot(data.filt, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()

VlnPlot(data.filt, features = "nCount_RNA", pt.size = 0.1) + NoLegend()

VlnPlot(data.filt, features = "percent.mito", pt.size = 0.1) + NoLegend()

VlnPlot(data.filt, features = "percent.ribo", pt.size = 0.1) + NoLegend()

# and check the number of cells per sample before and after filtering
table(Idents(alldata))

table(Idents(data.filt))


```

<p align="center">
<img src="/images/Blog/QC/qcmr.png" width="600">
</p>

<p align="center">
<img src="/images/Blog/QC/qcmrc.png" width="600">
</p>


##### Gene detection filtering
- Extremely high number of detected genes could indicate doublets. However, depending on the celltype composition in your sample, you may have cells with higher number of genes (and also higher counts) from one celltype.

- 被检测基因数量极高的可能表明doublets
- 基因检测方面，v2和v3也有明显的差异，数据的过滤不能用都采用相同的界限
- 有蛋白质分析数据中，有许多细胞几乎没有检测到的基因，但呈双峰分布。这种类型的分布在其他两个数据集中没有看到。考虑到它们都是PBMC数据集，把这个分布看作低质量的库是有意义的
- 过滤高基因检测的细胞(假定doublets)，v3的cutoff为4100，v2的cutoff为2000

##### Mitochondrial filtering
- We have quite a lot of cells with high proportion of mitochondrial reads. It could be wise to remove those cells, if we have enough cells left after filtering. Another option would be to either remove all mitochondrial reads from the dataset and hope that the remaining genes still have enough biological signal. A third option would be to just regress out the percent.mito variable during scaling.

- In this case we have as much as 99.7% mitochondrial reads in some of the cells, so it is quite unlikely that there is much celltype signature left in those.

- 有相当多的细胞具有高比例的线粒体reads。如果过滤后我们还有足够的细胞，最好去除这些细胞。另一种是从数据集中移除所有线粒体，剩余的基因仍然有足够的信号
- 以上数据分析中，有高达99.7%的线粒体存在细胞中，所以不太可能有很多细胞类型的标记留在这些细胞中
- 看图作出合理的决定，在哪里划出界限。看以上的数据中，大部分细胞的线粒体读数低于25%

##### Removal of cell cycle effect
```r

#Calculate cell-cycle scores


#Seurat has a function for calculating cell cycle scores based on a list of know S-phase and G2/M-phase genes.

data.filt <- CellCycleScoring(
  object = data.filt,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)

VlnPlot(data.filt, features = c("S.Score","G2M.Score"))


```

- Calculating cell cycle scores based on a list of know S-phase and G2/M-phase genes (In this case it looks like we only have a few cycling cells in the datasets)

<p align="center">
<img src="/images/Blog/QC/cce.png" width="600">
</p>



### [SCATER: pre-processing, quality control, normalisation and visualisation of single-cell RNA-seq data in R](https://bioconductor.org/packages/release/bioc/html/scater.html)

##### Most expressed features

```r

sce <- as.SingleCellExperiment(data.filt)


#Calculate QC-metrics

# calculate all qc-metrics
sce <- calculateQCMetrics(sce, feature_controls = list(mito = mt.genes))

# check what all entries are -
colnames(colData(sce))

colnames(rowData(sce))

#Most expressed features
plotHighestExprs(sce, exprs_values = "counts")



```
<p align="center">
<img src="/images/Blog/QC/mef.png" width="600">
</p>

##### Cumulative expression

```r
#Cumulative expression

# plot each sample separately
plotScater(sce, block1 = "ident", nfeatures = 1000)

#Plot gene stats
plotRowData(sce, x = "n_cells_by_counts", y = "mean_counts")

#Plot cell stats
#In the same manner plotColData can plot any of the qc-measures for cells.

p1 <- plotColData(sce, x = "total_counts",
    y = "total_features_by_counts", colour_by = "ident")
p2 <- plotColData(sce, x = "pct_counts_feature_control",
    y = "total_features_by_counts", colour_by = "ident")
p3 <- plotColData(sce, x = "pct_counts_feature_control",
    y = "pct_counts_in_top_50_features", colour_by = "ident")
multiplot(p1, p2, p3, cols = 2)


```
<p align="center">
<img src="/images/Blog/QC/ce.png" width="600">
</p>


##### Explanatory factors

```r

```
<p align="center">
<img src="/images/Blog/QC/ef.png" width="600">
</p>




## references

- Ding, J., et al. (2020). "Systematic comparison of single-cell and single-nucleus RNA-sequencing methods." Nature Biotechnology 38(6): 737-746.
- Hicks, S. C., et al. (2018). "Missing data and technical variability in single-cell RNA-sequencing experiments." Biostatistics 19(4): 562-578.
- Kang, H. M., et al. (2018). "Multiplexed droplet single-cell RNA-sequencing using natural genetic variation." Nat Biotechnol 36(1): 89-94.
- Klein, A. M., et al. (2015). "Droplet barcoding for single-cell transcriptomics applied to embryonic stem cells." Cell 161(5): 1187-1201.
- Liu, S. and C. Trapnell (2016). "Single-cell transcriptome sequencing: recent advances and remaining challenges." F1000Res 5.
- McCarthy, D. J., et al. (2017). "Scater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R." Bioinformatics 33(8): 1179-1186.
- Mereu, E., et al. (2020). "Benchmarking single-cell RNA-sequencing protocols for cell atlas projects." Nat Biotechnol 38(6): 747-755.
- Slyper, M., et al. (2020). "A single-cell and single-nucleus RNA-Seq toolbox for fresh and frozen human tumors." Nat Med 26(5): 792-802.
- Stuart, T. and R. Satija (2019). "Integrative single-cell analysis." Nat Rev Genet 20(5): 257-272.
Wagner, A., et al. (2016). "Revealing the vectors of cellular identity with single-cell genomics." Nature Biotechnology 34(11): 1145-1160.
- Wagner, D. E. and A. M. Klein (2020). "Lineage tracing meets single-cell omics: opportunities and challenges." Nat Rev Genet.
- Ziegenhain, C., et al. (2017). "Comparative Analysis of Single-Cell RNA Sequencing Methods." Mol Cell 65(4): 631-643.e634.
- Wu, Y., Zhang, K. Tools for the analysis of high-dimensional single-cell RNA sequencing data. Nat Rev Nephrol (2020).
- Lafzi, A., Moutinho, C., Picelli, S. et al. Tutorial: guidelines for the experimental design of single-cell RNA sequencing studies. Nat Protoc 13, 2742–2757 (2018). 
- Malte D Luecken;Fabian J Theis.Current best practices in single‐cell RNA‐seq analysis: a tutorial.Mol Syst Biol. (2019)




