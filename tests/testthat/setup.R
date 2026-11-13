library(data.table)
library(dplyr)

resolve_pkg_root <- function() {
  candidates <- c(
    file.path("..", ".."),
    file.path("..", "..", "00_pkg_src", "medusa"),
    file.path("..", "..", ".."),
    file.path("..", "..", "..", "00_pkg_src", "medusa")
  )
  for (cand in candidates) {
    if (file.exists(file.path(cand, "DESCRIPTION")) && dir.exists(file.path(cand, "inst"))) {
      return(normalizePath(cand, mustWork = TRUE))
    }
  }
  installed_val <- system.file("validation", package = "medusa")
  if (nzchar(installed_val) && dir.exists(installed_val)) {
    staged <- file.path(tempdir(), "medusa_staged_root")
    dir.create(file.path(staged, "inst"), recursive = TRUE, showWarnings = FALSE)
    file.copy(installed_val, file.path(staged, "inst"), recursive = TRUE)
    return(normalizePath(staged, mustWork = TRUE))
  }
  normalizePath(file.path("..", ".."), mustWork = FALSE)
}
pkg_root <- resolve_pkg_root()

source_tree_file <- function(...) {
  parts <- c(...)
  path <- do.call(file.path, as.list(c(pkg_root, parts)))
  if (file.exists(path)) {
    return(normalizePath(path, mustWork = TRUE))
  }
  if (length(parts) > 1L && parts[[1L]] == "inst") {
    installed <- system.file(do.call(file.path, as.list(parts[-1L])), package = "medusa")
    if (nzchar(installed) && file.exists(installed)) {
      return(normalizePath(installed, mustWork = TRUE))
    }
  }
  testthat::skip(sprintf("Source-tree file not available in installed-package checks: %s", do.call(file.path, as.list(parts))))
}
