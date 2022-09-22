# Save `SingleCellExperiment`s to file

The **alabaster.sce** package implements methods for saving and loading `SingleCellExperiment` objects under the **alabaster** framework.
It provides a language-agnostic method for serializing reduced dimensions and alternative experiments along with the usual bits and pieces from the `SummarizedExperiment` parent class.
To get started, install the package and its dependencies from GitHub:

```r
devtools::install_github("ArtifactDB/alabaster.schemas")
devtools::install_github("ArtifactDB/alabaster.base")
devtools::install_github("ArtifactDB/alabaster.ranges")
devtools::install_github("ArtifactDB/alabaster.matrix")
devtools::install_github("ArtifactDB/alabaster.se")
devtools::install_github("ArtifactDB/alabaster.sce")
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
dir.create(tmp)
meta <- stageObject(sce, tmp, "sce")
meta[["$schema"]]
## [1] "single_cell_experiment/v1.json"

roundtrip <- loadObject(meta, tmp)
class(roundtrip)
## [1] "SingleCellExperiment"
## attr(,"package")
## [1] "SingleCellExperiment"
```
