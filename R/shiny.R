#' Shiny application for PMD analysis
#' @export
runPMD <- function() {
    file <- system.file("shinyapp", "PMD.Rmd",
        package = "pmd")
    if (file == "") {
        stop("Could not find directory. Try re-installing `pmd`.",
            call. = FALSE)
    }
    rmarkdown::run(file)
}
#' Shiny application for PMD network analysis
#' @export
runPMDnet <- function() {
    file <- system.file("shinyapp", "pmdnet.Rmd",
                        package = "pmd")
    if (file == "") {
        stop("Could not find directory. Try re-installing `pmd`.",
             call. = FALSE)
    }
    rmarkdown::run(file)
}
