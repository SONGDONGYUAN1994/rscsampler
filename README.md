# rscsampler

## Overview
An R interface of Python module [scsampler](https://github.com/SONGDONGYUAN1994/scsampler). It is designed for fast diversity-preserving subsampling of large-scale single-cell transcriptomic data.

## Installation
The package relies on the R package `reticulate` for calling `scsampler` in Python (which is designed for `scanpy`).
Please first install `scsampler` and `scanpy` from PyPI:
```{python}
pip install scanpy
pip install scsampler
```
Then install `rscsampler` in R:
```{r}
install.packages("devtools")
install.packages("reticulate")
library(devtools)
devtools::install_github("SONGDONGYUAN1994/rscsampler")
```
## Usage
The package can take the input dataset in the format of: `SingleCellExperiment`, `SeuratObject` or a matrix. By default, we use the top 50 PCs for cell selection. Here
we show an example of the top PCs examplary dataset:
```{r}
library(reticulate)
library(rscsampler)
data("pbmc68k")
subsample_index <- scsampler(dat = pbmc68k, use_pca = FALSE, run_pca = FALSE, fraction = 0.1) ## 10% subsamples
```
We still recommend users using the Python module `scsampler` since R is not very ideal for dealing with large datasets.

## Contact
Any questions or suggestions on `rscsampler` are welcomed! If you have any questions, please report it on [issues](https://github.com/SONGDONGYUAN1994/rscsampler/issues) or contact Dongyuan (<dongyuansong@ucla.edu>).
