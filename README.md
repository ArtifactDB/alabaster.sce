# Save `SingleCellExperiment`s to file

|Environment|Status|
|---|---|
|[BioC-release](https://bioconductor.org/packages/release/bioc/html/alabaster.sce.html)|[![Release OK](https://bioconductor.org/shields/build/release/bioc/alabaster.sce.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/alabaster.sce/)|
|[BioC-devel](https://bioconductor.org/packages/devel/bioc/html/alabaster.sce.html)|[![Devel OK](https://bioconductor.org/shields/build/devel/bioc/alabaster.sce.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/alabaster.sce/)|

The **alabaster.sce** package implements methods for saving and loading `SingleCellExperiment` objects under the **alabaster** framework.
It provides a language-agnostic method for serializing reduced dimensions and alternative experiments along with the usual bits and pieces from the `SummarizedExperiment` parent class.
To get started, install the package and its dependencies from Bioconductor:

```r
# install.packages("BiocManager")
BiocManager::install("alabaster.sce")
```

In the example below, we save a `SingleCellExperiment` object to file:

```r
library(SingleCellExperiment)
example(SingleCellExperiment, echo=FALSE) # can't be bothered to copy it here.
sce
## class: SingleCellExperiment
## dim: 200 100
## metadata(0):
## assays(2): counts logcounts
## rownames: NULL
## rowData names(0):
## colnames: NULL
## colData names(0):
## reducedDimNames(2): PCA tSNE
## mainExpName: NULL
## altExpNames(0):

library(alabaster.sce)
tmp <- tempfile()
saveObject(sce, tmp)

roundtrip <- readObject(tmp)
class(roundtrip)
## [1] "SingleCellExperiment"
## attr(,"package")
## [1] "SingleCellExperiment"
```
