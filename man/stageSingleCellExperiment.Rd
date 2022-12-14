% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stageSingleCellExperiment.R
\name{stageObject,SingleCellExperiment-method}
\alias{stageObject,SingleCellExperiment-method}
\alias{stageObject,ReducedDimensions-method}
\title{Stage an experiment}
\usage{
\S4method{stageObject}{SingleCellExperiment}(x, dir, path, child = FALSE, rd.name = "dimreds", ...)
}
\arguments{
\item{x}{A \linkS4class{SingleCellExperiment} object or one of its subclasses.}

\item{dir}{String containing the path to the staging directory.}

\item{path}{String containing a prefix of the relative path inside \code{dir} where \code{x} is to be saved.
The actual path used to save \code{x} may include additional components, see Details.}

\item{child}{Logical scalar indicating whether \code{x} is a child of a larger object.}

\item{rd.name}{String containing the prefix of the file to save the reduced dimensions.}

\item{...}{Further arguments to pass to the RangedSummarizedExperiment method.}
}
\value{
A named list of metadata that follows the \code{single_cell_experiment} schema.
The contents of \code{x} are saved into a \code{path} subdirectory inside \code{dir}.
}
\description{
Save a \linkS4class{SingleCellExperiment} to file inside the staging directory.
}
\examples{
# Mocking up an SCE:
mat <- matrix(rpois(10000, 10), ncol=10)
colnames(mat) <- letters[1:10]
rownames(mat) <- sprintf("GENE_\%i", seq_len(nrow(mat)))

se <- SingleCellExperiment(list(counts=mat))
se$stuff <- LETTERS[1:10]
se$blah <- runif(10)
reducedDims(se) <- list(
    PCA=matrix(rnorm(ncol(se)*10), ncol=10),
    TSNE=matrix(rnorm(ncol(se)*2), ncol=2)
)
altExps(se) <- list(spikes=SummarizedExperiment(list(counts=mat[1:2,])))

# Staging it: 
tmp <- tempfile()
dir.create(tmp)
stageObject(se, dir=tmp, "rna-seq") 
list.files(file.path(tmp, "rna-seq"))

}
\author{
Aaron Lun
}
