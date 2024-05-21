# This tests the stageObject generic for SCE's.
# library(testthat); library(alabaster.sce); source("test-SingleCellExperiment.R")

# Making an SE and annotating it.
mat <- matrix(rpois(2000, 10), ncol=10)
colnames(mat) <- paste0("SAMPLE_", seq_len(ncol(mat)))

se <- SingleCellExperiment(list(counts=mat, cpm=mat/10))
se$stuff <- LETTERS[1:10]
se$blah <- runif(10)
rowData(se)$whee <- runif(nrow(se))
rownames(se) <- sprintf("GENE_%i", seq_len(nrow(se)))

test_that("stageObject works as expected for SCE objects", {
    tmp <- tempfile()
    dir.create(tmp)

    info <- stageObject(se, tmp, "rnaseq")
    alabaster.base::.writeMetadata(info, tmp)

    expect_identical(vapply(info$summarized_experiment$assays, function(x) x$name, ""), assayNames(se))
    cdf <- read.csv(file.path(tmp, info$summarized_experiment$column_data$resource$path), row.names=1)
    expect_equal(cdf, as.data.frame(colData(se)))
    expect_identical(info$single_cell_experiment, setNames(list(), character(0)))

    # Round trip works.
    out <- loadSingleCellExperiment(info, tmp)
    expect_s4_class(out, "SingleCellExperiment")
    expect_equal(rowData(out), rowData(se))
    expect_equal(colData(out), colData(se))

    # New world works.
    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_s4_class(out2, "SingleCellExperiment")
    expect_identical(rowData(out2), rowData(se))
    expect_identical(colData(out2), colData(se))
})

test_that("saveObject works with some non-trivial rowRanges", {
    all.ranges <- GRanges("chrX", IRanges(seq_len(nrow(mat) * 2), width=1))
    rowRanges(se) <- splitAsList(all.ranges, rep(seq_len(nrow(mat)), length.out=length(all.ranges)))

    tmp <- tempfile()
    saveObject(se, tmp, "rnaseq")
    out2 <- readObject(tmp)

    expect_identical(rowData(out2), rowData(se))
    expect_identical(rowRanges(out2), rowRanges(se))
    expect_identical(rownames(out2), rownames(se))
})

test_that("stageObject works as expected with reduced dims inside", {
    reducedDims(se) <- list(PCA=matrix(rnorm(ncol(mat)*50), ncol=50), TSNE=cbind(TSNE1=runif(ncol(mat)), TSNE2=runif(ncol(mat))))

    tmp <- tempfile()
    dir.create(tmp)

    info <- stageObject(se, tmp, "rnaseq")
    alabaster.base::.writeMetadata(info, tmp)

    rd.info <- info$single_cell_experiment$reduced_dimensions
    expect_identical(length(rd.info), 2L)
    expect_identical(vapply(rd.info, function(x) x$name, ""), c("PCA", "TSNE"))

    # Round trip works.
    out <- loadSingleCellExperiment(info, tmp)
    expect_identical(reducedDimNames(out), reducedDimNames(se))
    expect_identical(as.matrix(reducedDim(out, "PCA")), reducedDim(se, "PCA"))
    expect_identical(as.matrix(reducedDim(out, "TSNE")), reducedDim(se, "TSNE"))

    # New world works.
    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_identical(reducedDimNames(out2), reducedDimNames(se))
    expect_identical(as.matrix(reducedDim(out2, "PCA")), reducedDim(se, "PCA"))
    expect_identical(as.matrix(reducedDim(out2, "TSNE")), reducedDim(se, "TSNE"))
})

test_that("stageObject fails with non-unique, non-empty reduced dim names", {
    reducedDims(se) <- list(PCA=matrix(rnorm(ncol(mat)*50), ncol=50), TSNE=cbind(TSNE1=runif(ncol(mat)), TSNE2=runif(ncol(mat))))
    colnames(int_colData(se)$reducedDims) <- c("FOO", "FOO") # need to force the issue as the setters don't normally allow duplicates.

    tmp <- tempfile()
    dir.create(tmp)
    expect_error(info <- stageObject(se, tmp, "rnaseq"), "duplicate")

    colnames(int_colData(se)$reducedDims) <- c("", "FOO") 
    unlink(file.path(tmp, "rnaseq"), recursive=TRUE)
    expect_error(info <- stageObject(se, tmp, "rnaseq"), "empty")
})

test_that("stageObject works as expected with alternative experiments inside", {
    tmp <- tempfile()
    dir.create(tmp)

    copy <- se
    altExps(se) <- list(spikes=copy[1:2,], protein=copy[3:5,])

    info <- stageObject(se, tmp, "rnaseq")
    alabaster.base::.writeMetadata(info, tmp)

    expect_identical(vapply(info$single_cell_experiment$alternative_experiments, function(x) x$name, ""), altExpNames(se))

    round <- loadSingleCellExperiment(info, tmp)
    expect_identical(altExpNames(round), altExpNames(se))
    expect_identical(rownames(altExp(round, 1)), rownames(altExp(se, 1)))
    expect_identical(rownames(altExp(round, 2)), rownames(altExp(se, 2)))

    # New world works.
    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_identical(altExpNames(out2), altExpNames(se))
    expect_identical(rownames(altExp(out2, 1)), rownames(altExp(se, 1)))
    expect_identical(rownames(altExp(out2, 2)), rownames(altExp(se, 2)))
})

test_that("stageObject fails with non-unique, non-empty altexp names", {
    copy <- se
    altExps(se) <- list(spikes=copy[1:2,], protein=copy[3:5,])
    colnames(int_colData(se)$altExps) <- c("FOO", "FOO") # need to force the issue as the setters don't normally allow duplicates.

    tmp <- tempfile()
    dir.create(tmp)
    expect_error(info <- stageObject(se, tmp, "rnaseq"), "duplicate")

    colnames(int_colData(se)$altExps) <- c("", "FOO") 
    unlink(file.path(tmp, "rnaseq"), recursive=TRUE)
    expect_error(info <- stageObject(se, tmp, "rnaseq"), "empty")
})

test_that("stageObject works as expected when we slap in a main name", {
    tmp <- tempfile()
    dir.create(tmp)
    mainExpName(se) <- "FOO"

    info <- stageObject(se, tmp, "rnaseq")
    alabaster.base::.writeMetadata(info, tmp)

    expect_identical(info$single_cell_experiment$main_experiment_name, "FOO")

    round <- loadSingleCellExperiment(info, tmp)
    expect_identical(mainExpName(round), "FOO")

    # New world works.
    tmp2 <- tempfile()
    saveObject(se, tmp2)
    out2 <- readObject(tmp2)
    expect_identical(mainExpName(out2), "FOO")

    # Fails if there's an alternative experiment with the same name.
    altExp(se, "FOO") <- se[1:10,]
    unlink(file.path(tmp, "rnaseq"), recursive=TRUE)
    expect_error(stageObject(se, tmp, "rnaseq"), "conflicting name 'FOO'")

})
