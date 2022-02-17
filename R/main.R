#' The scSampler function
#'
#' \code{scSampler} takes the input single cell data and return the selected diverse subsample.
#'
#' @param dat The input single cell dataset. Can be one of: a \code{SingleCellExperiment} object, a \code{SeuratObject} object, or a matrix.
#' @param use_pca A logical value. If TRUE, the algorithm will use the PC space for subsampling. We highly recommend setting it as TRUE.
#' @param run_pca A logical value. If TRUE, run PCA by \code{irlba::prcomp_irlba}. If FALSE, use the PCs from the input object.
#' @param use_dim An integer of how many PCs are used. Defalut is 50.
#' @param assay_use The input expression space.
#' @param fraction A numeric value from 0 to 1. The fraction of the subsample.
#' @param n_obs An integer of the size of the subsample.
#' @param random_state An integer of the random seed.
#' @param copy A logical value. If TRUE, return a subset of the origianl object, otherwise the indices of the subsample. Default is FALSE.
#' @param random_split An integer of how many splits for speeding up. Default is 1 (no split).
#'
#' @return A subset of the original object or the subsample indicies.
#'
#' @export scsampler
#'
#' @examples
#' data("pbmc68k")
#' res <- scsampler(dat = pbmc68k, use_pca = FALSE, run_pca = FALSE)

scsampler <- function(dat,
                      use_pca = TRUE,
                      run_pca = TRUE,
                      use_dim = 50,
                      assay_use = "logcounts",
                      fraction = 0.1,
                      n_obs = NULL,
                      random_state = 1,
                      copy = FALSE,
                      random_split = 1) {

  if(class(dat)[1] == "SingleCellExperiment") {
    mat <- SummarizedExperiment::assay(dat, assay_use)
    mat <- t(as.matrix(mat))
    if(use_pca) {
      if(run_pca) {
        mat <- irlba::prcomp_irlba(mat, n = use_dim)
      } else {
        mat <- SingleCellExperiment::reducedDim(dat, "PCA")
      }
    }
  } else if(class(dat)[1] == "SeuratObject") {
    mat <- Seurat::GetAssayData(object = dat, slot = assay_use)
    mat <- t(as.matrix(mat))
    if(use_pca) {
      if(run_pca) {
        mat <- irlba::prcomp_irlba(mat, n = use_dim)
      } else {
        mat <- Seurat::Embeddings(dat, "pca")
      }
    }
  } else {
    mat <- dat
    if(use_pca) {
      if(run_pca) {
        mat <- irlba::prcomp_irlba(mat, n = use_dim)
      }
    }
  }

  # use superassignment to update global reference to scipy
  res <- scsampler_python$scsampler(data = mat, fraction = fraction, n_obs = n_obs, random_state = as.integer(random_state), copy = FALSE, random_split = as.integer(random_split))
  ## Python index starts with 0. Add 1 for R.
  res <- as.vector(res + 1)

  if(copy) {
    res <- dat[, res]
  }

  return(res)
}
