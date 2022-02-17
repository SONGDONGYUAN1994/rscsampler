scsampler_python <- NULL
.onLoad <- function(libname, pkgname) {
  scsampler_python <<- reticulate::import("scsampler", delay_load = TRUE)
}
