#' @title Add cell identifiers to SpatialPoints as in INSPIRE Specification
#' @description
#' \code{spatialPointsCellCodes} returns a SpatialPointsDataFrame with identifiers
#' (CellCode and CellNum) for each point as in INSPIRE Specification on Geographical
#' Grid Systems.
#' @details
#' INSPIRE coding system for point identifiers\cr
#' The objective of the coding system is to generate unique identifiers for each
#' point, for any of the recommended resolutions.\cr
#' The cellCode is a text string, composed of cell size and cell coordinates.
#' Cell codes start with a cell size prefix. The cell size is denoted in meter (m)
#' for cell sizes below 1000 m and kilometre (km) for cell sizes from 1000 m and
#' above.\cr
#' Examples: a 100 meter cell has an identifier starting with “100m”, the
#' identifier  of a 10000 meter cell starts with “10km”.\cr
#' The coordinate part of the cell code reflects the distance of the lower left
#' grid cell corner from the false origin of the CRS. In order to reduce the
#' length of the string, Easting (E) and Northing (N) values are divided by
#' 10n (n is the number of zeros in the cell size value). Example for a cell
#' size of 10000 meters: The number of zeros in the cell size value is 4.
#' The resulting divider for Easting and Northing values is 104 = 10000.\cr
#' The cellNum is a sequence of concatenated integers identifying all the
#' hierarchical partitions of the main cell in which the point resides.
#' For instance, the cellNum of the top right cell would be 416 (fourth
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
#' @param points object of class "SpatialPoints" or "SpatialPointsDataFrame".
#' @param dim a single integer specifying the initial cell sizes, defaults to
#' 1km.
#' @param layers a single integer specifying the number of divisions of the
#' initial cells, defaults to 1.
#' @return A "SpatialPointsDataFrame" containing all the points given. For
#' each point a cellCode and cellNum identify the cell to which the point belongs.\cr
#' CellCode is a text string, composed of cell size and cell coordinates.
#' Cell codes start with a cell size prefix. The cell size is denoted in meter (m)
#' for cell sizes below 1000 m and kilometer (km) for cell sizes from 1000 m and
#' above.\cr
#' The cellNum is a sequence identifying the different partitions of the main
#' cell in which the point resides.
#' @export
#' @examples
#' data("BarcelonaPop")
#' BarcelonaPop.IDs<-spatialPointsCellCodes(BarcelonaPop)
#' BarcelonaPop.IDs.10km<-spatialPointsCellCodes(BarcelonaPop, 10000, 3)
#'
spatialPointsCellCodes <- function(points, dim=1000, layers=1){
  ## aux function to determine the number of trailing zeros
  trailingZeros<-function(x){
    i<-0
    while (x %% 10 == 0 ) {
      x<-x%/%10
      i<-i+1
    }
    return(i)
  }
  if (missing(points)) stop("argument 'points' is missing missing, with no default", call.="FALSE")
  if (length(points)==0) stop("argument 'points' has length 0", call.="FALSE")
  #stopifnot(require("sp"))
  stopifnot(dim>0, (class(points) %in% c("SpatialPoints", "SpatialPointsDataFrame")))

  # calculate string CellCode of the form "1kmNyyyyExxxx"
  sizePrefix<-ifelse(dim>=1000, paste0(dim/1000, "km"), paste0(dim, "m"))

  # add x, y cell origin to each point
  points$CellOrigin.x<-points@coords[,1]%/%dim*dim
  points$CellOrigin.y<-points@coords[,2]%/%dim*dim
  # calculate string CellCode of the form "1kmNyyyyExxxx"
  zerosToRemove<-as.integer(10^trailingZeros(dim))
  cellCodeE<-points$CellOrigin.x/zerosToRemove
  cellCodeN<-points$CellOrigin.y/zerosToRemove
  lenCellCode<-nchar(max(cellCodeE, cellCodeN))
  points$cellCode<-paste0(sizePrefix, "N",formatC(cellCodeN, width=lenCellCode, flag=0, mode = 'integer'), "E",formatC(cellCodeE, width=lenCellCode, flag=0, mode = 'integer'))

  points$cellNum<-''
  if (layers>1) {
    # cellNumPos stores the number of digits for each subpart of the cell codes
    #   cellNumPosStart<-c( 1,14,15,17,19,22)
    #   cellNumPosStop <-c(13,14,16,18,21,25)
    cellNumPosStart<-c(1)
    cellNumPosStop<-c(nchar(points$CellCode[1]))
    # calculate CellNum of the form "nmmooopppp..." with n 1:4   mm 01:16  ooo 001:256   pppp 0001:1024 ...
    for (i in 2:layers){
      zeros<-ceiling(log10(2^(2*(i-1))))
      cellNumPosStop[i]<-cellNumPosStop[i-1]+zeros
      cellNumPosStart[i]<-cellNumPosStop[i-1]+1
      size<-dim/2^(i-1)
      points$cellNum<-paste0(points$cellNum, formatC((points$x-points$CellOrigin.x)%/%size + (2^(i-1))*(points$y-points$CellOrigin.y)%/%size+1, width=zeros, flag=0, mode = 'integer'))
    }
  }
  points$CellOrigin.x<-NULL
  points$CellOrigin.y<-NULL
  return(points)
}
