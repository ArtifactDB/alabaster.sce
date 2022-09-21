#' Stage an experiment
#'
#' Save a \linkS4class{SingleCellExperiment} to file inside the staging directory.
#' 
#' @param x A \linkS4class{SingleCellExperiment} object or one of its subclasses.
#' @inheritParams alabaster.base::stageObject
#' @param rd.name String containing the prefix of the file to save the reduced dimensions.
#' @param ... Further arguments to pass to the RangedSummarizedExperiment method.
#'
#' @author Aaron Lun
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
#' dir.create(tmp)
#' stageObject(se, dir=tmp, "rna-seq") 
#' list.files(file.path(tmp, "rna-seq"))
#' 
#' @export
#' @rdname stageSingleCellExperiment
#' @aliases stageObject,ReducedDimensions-method
#' @import SingleCellExperiment alabaster.base methods
setMethod("stageObject", "SingleCellExperiment", function(x, dir, path, child=FALSE, rd.name="dimreds", ...) {
    meta <- callNextMethod()
    sce.details <- list()
    names(sce.details) <- character(0)

    # Staging the reduced dimensions.
    red.nms <- reducedDimNames(x)
    if (length(red.nms)) {
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
        sce.details$main_experiment_name <- main.nm
    }

    meta$single_cell_experiment <- sce.details
    meta[["$schema"]] <- "single_cell_experiment/v1.json"

    meta
})
