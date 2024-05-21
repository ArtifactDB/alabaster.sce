#' Save a SingleCellExperiment to disk
#'
#' Save a \linkS4class{SingleCellExperiment} to its on-disk representation.
#' 
#' @param x A \linkS4class{SingleCellExperiment} object or one of its subclasses.
#' @inheritParams alabaster.base::saveObject
#' @param ... Further arguments to pass to the RangedSummarizedExperiment method.
#'
#' @author Aaron Lun
#' 
#' @return \code{x} is saved into \code{path} and \code{NULL} is invisibly returned.
#'
#' @seealso
#' \code{\link{readSingleCellExperiment}}, to read the SingleCellExperiment back into the R session.
#' 
#' @examples
#' # Mocking up an SCE:
#' mat <- matrix(rpois(10000, 10), ncol=10)
#' colnames(mat) <- letters[1:10]
#' rownames(mat) <- sprintf("GENE_%i", seq_len(nrow(mat)))
#'
#' se <- SingleCellExperiment(list(counts=mat))
#' se$stuff <- LETTERS[1:10]
#' se$blah <- runif(10)
#' reducedDims(se) <- list(
#'     PCA=matrix(rnorm(ncol(se)*10), ncol=10),
#'     TSNE=matrix(rnorm(ncol(se)*2), ncol=2)
#' )
#' altExps(se) <- list(spikes=SummarizedExperiment(list(counts=mat[1:2,])))
#'
#' # Staging it: 
#' tmp <- tempfile()
#' saveObject(se, tmp)
#' list.files(tmp, recursive=TRUE)
#' 
#' @name saveSingleCellExperiment
#' @aliases 
#' stageObject,SingleCellExperiment-method
NULL

#' @export
#' @rdname saveSingleCellExperiment
#' @import SingleCellExperiment alabaster.base methods
#' @importFrom jsonlite toJSON
setMethod("saveObject", "SingleCellExperiment", function(x, path, ...) {
    base <- as(x, "RangedSummarizedExperiment")
    altSaveObject(base, path, ...)

    red.nms <- reducedDimNames(x)
    if (length(red.nms)) {
        rddir <- file.path(path, "reduced_dimensions")
        dir.create(rddir, showWarnings=FALSE)
        for (i in seq_along(red.nms)) {
            red.path <- file.path(rddir, i - 1L)
            tryCatch({
                altSaveObject(reducedDim(x, red.nms[i], withDimnames=FALSE), red.path, ...)
            }, error=function(e) {
                stop("failed to stage 'reducedDim(<", class(x)[1], ">, \"", red.nms[i], "\")'\n  - ", e$message)
            })
        }
        write(toJSON(red.nms), file=file.path(rddir, "names.json"))
    }

    alt.nms <- altExpNames(x) 
    if (length(alt.nms)) {
        aedir <- file.path(path, "alternative_experiments")
        dir.create(aedir, showWarnings=FALSE)
        for (i in seq_along(alt.nms)) {
            alt.path <- file.path(aedir, i - 1L)
            tryCatch({ 
                altSaveObject(altExp(x, alt.nms[i], withDimnames=FALSE), alt.path, ...)
            }, error =function(e) {
                stop("failed to stage 'altExp(<", class(x)[1], ">, \"", alt.nms[i], "\")'\n  - ", e$message)
            })
        }
        write(toJSON(alt.nms), file=file.path(aedir, "names.json"))
    }

    info <- list(version="1.0")
    main.nm <- mainExpName(x)
    if (!is.null(main.nm)) {
        info$main_experiment_name <- main.nm
    }

    meta <- readObjectFile(path)
    meta$single_cell_experiment <- info
    saveObjectFile(path, "single_cell_experiment", meta)

    invisible(NULL)
})

##################################
######### OLD STUFF HERE #########
##################################

#' @export
setMethod("stageObject", "SingleCellExperiment", function(x, dir, path, child=FALSE, rd.name="dimreds", ...) {
    meta <- callNextMethod()
    sce.details <- list()
    names(sce.details) <- character(0)

    # Staging the reduced dimensions.
    red.nms <- reducedDimNames(x)
    if (length(red.nms)) {
        if (anyDuplicated(red.nms)) {
            stop("detected duplicate names for reduced dimensions in a ", class(x)[1], " object")
        }
        if (any(red.nms == "")) {
            stop("detected empty names for reduced dimensions in a ", class(x)[1], " object")
        }

        entries <- vector('list', length(red.nms))
        for (i in seq_along(red.nms)) {
            red.path <- file.path(path, paste0("reddim-", i))

            red.processed <- tryCatch({
                rd.meta <- .stageObject(reducedDim(x, red.nms[i]), dir, red.path, child=TRUE)
                .writeMetadata(dir=dir, rd.meta)
            }, error=function(e) {
                stop("failed to stage 'reducedDim(<", class(x)[1], ">, \"", red.nms[i], "\")'\n  - ", e$message)
            })

            entries[[i]] <- list(name = red.nms[i], resource = red.processed)
        }

        sce.details$reduced_dimensions <- entries
    }

    # Staging the alternative experiments.
    alt.nms <- altExpNames(x) 
    if (length(alt.nms)) {
        if (anyDuplicated(alt.nms)) {
            stop("detected duplicate names for alternative experiments in a ", class(x)[1], " object")
        }
        if (any(alt.nms == "")) {
            stop("detected empty names for alternative experiments in a ", class(x)[1], " object")
        }

        entries <- vector("list", length(alt.nms))
        for (i in seq_along(alt.nms)) {
            alt.path <- file.path(path, paste0("altexp-", i))

            alt.processed <- tryCatch({ 
                alt.meta <- .stageObject(altExp(x, alt.nms[i]), dir=dir, path = alt.path, child=TRUE)
                .writeMetadata(meta=alt.meta, dir=dir)
            }, error =function(e) {
                stop("failed to stage 'altExp(<", class(x)[1], ">, \"", alt.nms[i], "\")'\n  - ", e$message)
            })

            entries[[i]] <- list(name = alt.nms[i], resource = alt.processed)
        }

        sce.details$alternative_experiments <- entries
    }

    main.nm <- mainExpName(x)
    if (!is.null(main.nm)) {
        if (any(alt.nms == main.nm)) {
            stop("conflicting name '", main.nm, "' for main and alternative experiments in a ", class(x)[1], " object")
        }
        sce.details$main_experiment_name <- main.nm
    }

    meta$single_cell_experiment <- sce.details
    meta[["$schema"]] <- "single_cell_experiment/v1.json"

    meta
})
