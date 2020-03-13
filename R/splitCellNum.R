#' @title Split CellNum sequence into a vector
#' @description
#' \code{createQuadtree} returns a vector decomposing the given CellNum into the
#' sequence of the different cell numbers for each level.
#' @details
#' CellNum is an integer with the concatenated sequence of hierarchical cell positions
#' inside a main cell. \code{splitCellNum} splits that sequence into a vector.
#' For instance, the CellNum of the top right cell would be 416 (fourth
#' in first partition, sixteenth in second partition)\cr
#' \if{html}{\figure{CellNum.jpg}{options: width=200 alt="Hyerarchical CellNums"}}
#' \if{latex}{\figure{CellNum.jpg}{options: width=4cm}}
#' @seealso
#' \itemize{
#'   \item{
#'    D2.8.I.2 INSPIRE Specification on Geographical Grid Systems – Guidelines
#'    \url{https://inspire.ec.europa.eu/documents/Data_Specifications/INSPIRE_Specification_GGS_v3.0.1.pdf}
#'   }
#'   \item{
#'    EEA reference grid dataset
#'    \url{https://data.europa.eu/euodp/data/dataset/data_eea-reference-grids-2}
#'   }
#' }
#'
#' @param x a character or character vector containing a sequence of cell numbers or
#' an R object with a field named 'cellNum'
#' @return integer vector or list of integer vectors with the sequence
#' CellNums splitted
#' @export
#' @examples
#' data("CharlestonPop")
#' CharlestonPop.IDs<-spatialPointsCellCodes(CharlestonPop, layers=2)
#' splitCellNum(CharlestonPop.IDs)
#'
splitCellNum<-function(x){
  if ('cellNum' %in% names(x)) {
    x<-x$cellNum
  }
  if (length(x)>1){
    sapply(x, function(s){
      IDs<-c()
      t<-1
      while (nchar(s) > 0) {
        pos<-ceiling(log10(2^(2*t)))
        IDs<-c(IDs, as.integer(substr(s, 1, pos )))
        s<-substr(s, pos + 1 , nchar(s))
        t<-t+1
      }
      return(IDs)
    }, simplify = TRUE, USE.NAMES=FALSE)
  } else {
    IDs<-c()
    t<-1
    while (nchar(x) > 0) {
      pos<-ceiling(log10(2^(2*t)))
      IDs<-c(IDs, as.integer(substr(x, 1, pos )))
      x<-substr(x, pos + 1 , nchar(x))
      t<-t+1
    }
    return(IDs)
  }
}
