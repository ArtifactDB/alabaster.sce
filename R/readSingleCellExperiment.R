#' Read a SingleCellExperiment from disk
#'
#' Read a \linkS4class{SingleCellExperiment} object from its on-disk representation.
#' This is usually not directly called by users, but is instead called by dispatch in \code{\link{readObject}}.
#'
#' @param path String containing a path to a directory, itself created using the \code{\link{saveObject}} method for \linkS4class{SingleCellExperiment} objects.
#' @param metadata Named list of metadata for this object, see \code{\link{readObjectFile}} for details.
#' @param ... Further arguments passed to \code{\link{readRangedSummarizedExperiment}} and internal \code{\link{altReadObject}} calls.
#' 
#' @return A \linkS4class{SingleCellExperiment} object.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{"\link{saveObject,SingleCellExperiment-method}"}, to save the SingleCellExperiment to disk.
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
#' readObject(tmp)
#'
#' @export
#' @aliases loadSingleCellExperiment
#' @importFrom alabaster.se loadSummarizedExperiment
#' @import SingleCellExperiment alabaster.base 
#' @importFrom jsonlite fromJSON
readSingleCellExperiment <- function(path, metadata, ...) {
    metadata$type <- "ranged_summarized_experiment"
    se <- altReadObject(path, metadata, ...)
    se <- as(se, "SingleCellExperiment")

    rddir <- file.path(path, "reduced_dimensions")
    if (file.exists(rddir)) {
        reddim.names <- fromJSON(file.path(rddir, "names.json"))
        all.reddims <- vector("list", length(reddim.names))
        for (y in seq_along(reddim.names)) {
            all.reddims[[y]] <- altReadObject(file.path(rddir, y - 1L), ...)
        }
        names(all.reddims) <- reddim.names
        reducedDims(se, withDimnames=FALSE) <- all.reddims
    }

    aedir <- file.path(path, "alternative_experiments")
    if (file.exists(aedir)) {
        altexp.names <- fromJSON(file.path(aedir, "names.json"))
        all.altexps <- vector("list", length(altexp.names))
        for (y in seq_along(altexp.names)) {
            all.altexps[[y]] <- altReadObject(file.path(aedir, y - 1L), ...)
        }
        names(all.altexps) <- altexp.names 
        altExps(se, withDimnames=FALSE) <- all.altexps
    }

    info <- metadata$single_cell_experiment
    main.nm <- info$main_experiment_name
    if (!is.null(main.nm)) {
        mainExpName(se) <- main.nm
    }

    se
}

##################################
######### OLD STUFF HERE #########
##################################

#' @export
loadSingleCellExperiment <- function(exp.info, project, ...) {
    se <- loadSummarizedExperiment(exp.info, project, ...)
    se <- as(se, "SingleCellExperiment")

    all.reds <- exp.info$single_cell_experiment$reduced_dimensions
    if (length(all.reds)) {
        red.names <- character(length(all.reds))
        for (i in seq_along(all.reds)) {
            red.info <- all.reds[[i]]
            red.exp.info <- acquireMetadata(project, path=red.info$resource$path)
            all.reds[[i]] <- .loadObject(red.exp.info, project=project)
            red.names[i] <- red.info$name
        }
        names(all.reds) <- red.names
        reducedDims(se) <- all.reds
    }

    all.alts <- exp.info$single_cell_experiment$alternative_experiments
    if (length(all.alts)) {
        alt.names <- character(length(all.alts))
        for (i in seq_along(all.alts)) {
            alt.info <- all.alts[[i]]
            alt.exp.info <- acquireMetadata(project, path=alt.info$resource$path)
            all.alts[[i]] <- .loadObject(alt.exp.info, project=project)
            alt.names[i] <- alt.info$name
        }
        names(all.alts) <- alt.names
        altExps(se) <- all.alts
    }

    main.nm <- exp.info$single_cell_experiment$main_experiment_name
    if (!is.null(main.nm)) {
        mainExpName(se) <- main.nm
    }

    se
}
