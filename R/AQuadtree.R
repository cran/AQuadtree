#' AQuadtree: A package to anonymise spatial points data.
#'
#' @description
#' This package provides an S4 class for creating, manipulating
#' and exporting spatial quadtree varying size grids, and for methods
#' including print/show, plot, spplot, subset, [, [[, names,
#' dim, summary, write.
#'
#' @section Introduction:
#' The quadtree functions and class provide the tools to build a varying size
#' quadtree grid performing a bottom-up aggregation considering a minimum
#' threshold for each the cell.
#' The main goal of the package is the anonymization of a set of spatial
#' point data by an aggregation process as in a k-anonymity model. The grid
#' created follows the INSPIRE Specification on Geographical Grid Systems.
#'
#' @references
#'    D2.8.I.2 INSPIRE Specification on Geographical Grid Systems – Guidelines
#'    \url{https://inspire.ec.europa.eu/documents/Data_Specifications/INSPIRE_Specification_GGS_v3.0.1.pdf}
#'
#'    EEA reference grid dataset
#'    \url{https://data.europa.eu/euodp/data/dataset/data_eea-reference-grids-2}
#'
#' @docType package
#' @aliases AQuadtree-package
"_PACKAGE"
