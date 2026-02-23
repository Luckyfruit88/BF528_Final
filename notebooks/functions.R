load_required_packages <- function() {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(DiffBind)
    library(edgeR)
    library(ChIPseeker)
    library(motifmatchr)
    library(JASPAR2020)
  })
}

theme_publication <- function() {
  ggplot2::theme_bw(base_size = 12)
}
