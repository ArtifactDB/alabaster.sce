.onLoad <- function(libname, pkgname) {
    registerReadObjectFunction("single_cell_experiment", readSingleCellExperiment)
}

.onUnload <- function(libname, pkgname) {
    registerReadObjectFunction("single_cell_experiment", NULL)
}
