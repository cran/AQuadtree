#' @title Class "AQuadtree".
#' @description
#' An S4 class representing a Quadtree hierarchical geographic Grid
#' to anonymise spatial point data.
#'
#' Objects can be created by calls to the function \link{AQuadtree}
#'
#' @details
#' Given a set of points, the \code{AQuadtree class} represents a
#' varying size Quadtree grid created performing a
#' bottom-up aggregation considering a minimum threshold for each cell.
#' Cells with a value under the threshold for the \code{thresholdField} are
#' aggregated to the upper level in a quadtree manner.\cr
#' When no \code{thresholdField} is given, total number  of points in the cell
#' will be used, and so, given a threshold of k, none of the cells in the
#' resulting grid have a value less than k individuals as in a k-anonymity model.\cr
#' The Quadtree produced balances information loss and accuracy. For instance,
#' for the set of cells in the left image, where numbers in the cells represent
#' the values in the \code{thresholdField}, using a \code{threshold} value of 100,
#' the resulting Quadtree will be the one on the right. As we can see, some cells
#' will be discarded, and some aggregated to maintain as much information as
#' possible, keeping at the same time as much disaggregation as possible\cr
#' \if{html}{\figure{QTexampleA.png}{options: width=260 alt="62.5m2 cells"}}
#' \if{latex}{\figure{QTexampleA.png}{options: width=4.5cm}}
#' \if{html}{\figure{QTexampleB.png}{options: width=250 alt="resulting Quadtree"}}
#' \if{latex}{\figure{QTexampleB.png}{options: width=4.4cm}}\cr
#' The INSPIRE coding system for cell identifiers will be used to generate a
#' CellCode and CellNum for each cell in the Quadtree.
#' The objective of the coding system is to generate unique
#' identifiers for each cell, for any of the resolutions.\cr
#' The cellCode is a text string, composed of cell size and cell coordinates.
#' Cell codes start with a cell size prefix. The cell size is denoted in meter (m)
#' for cell sizes below 1000 m and kilometre (km) for cell sizes from 1000 m and
#' above.\cr
#' Example: a 100 meter cell has an identifier starting with “100m”, the
#' identifier  of a 10000 meter cell starts with “10km”.\cr
#' The coordinate part of the cell code reflects the distance of the lower left
#' grid cell corner from the false origin of the CRS. In order to reduce the
#' length of the string, Easting (E) and Northing (N) values are divided by
#' 10^n (n is the number of zeros in the cell size value). Example for a cell
#' size of 10000 meters: The number of zeros in the cell size value is 4.
#' The resulting divider for Easting and Northing values is 10^4 = 10000.\cr
#' The CellNum is a sequence of concatenated integers identifying all the
#' hierarchical partitions of the main cell in which the point resides.
#' For instance, the CellNum of the top right cell would be 416 (fourth
#' in first partition, sixteenth in second partition)\cr
#' The input object must be projected and units should be in 'meters'
#' because the system uses the INSPIRE coding system.
#'
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
#' @importFrom methods as as<- slot slot<- callGeneric
#' @name AQuadtree-class
#' @aliases AQuadtree-class
#' @slot dim = "numeric"
#' @slot layers = "numeric",
#' @slot colnames = "character",
#' @slot threshold = "numeric",
#' @slot thresholdField = "character",
#' @slot loss = "numeric"
#' @exportClass AQuadtree
#'
#' @examples
#' data("BarcelonaPop", "BarcelonaCensusTracts")
#' aquadtree.Barcelona<-AQuadtree(BarcelonaPop, layers = 3)
#' plot(aquadtree.Barcelona)
#'
#' aQuadtree.Charleston<-AQuadtree(CharlestonPop, colnames="sex", threshold=17,
#'   thresholdField=c("sex.male", "sex.female"))
#'
#' \dontrun{
#' ## spatial object not projected
#' sp.not.projected<-spTransform(CharlestonPop,CRS("+proj=longlat +datum=NAD27"))
#' is.projected(sp.not.projected)
#' aqt<-AQuadtree(sp.not.projected)
#'
#' ## not an SpatialPoints object
#' aqt<-AQuadtree(CharlestonCensusTracts)
#'
#' ## too many subdivisions
#' aqt<-AQuadtree(CharlestonPop, layers=15)
#'
#' }
setClass(
  Class="AQuadtree",
  contains = "SpatialPolygonsDataFrame",
  slots = c(
    dim = "numeric",
    layers = "numeric",
    colnames = "character",
    threshold = "numeric",
    thresholdField = "character",
    loss = "numeric"
  )
)
#' Wrapper function AQuadtree.
#'
#' @rdname AQuadtree-class
#' @title AQuadtree
#' @details function to create an object of class AQuadtree
#'
#' @param points object of class "SpatialPoints" or "SpatialPointsDataFrame".
#' @param dim a single integer specifying the initial cell sizes in meters, defaults to 1000.
#' @param layers a single integer specifying the number of divisions of the
#' initial cells, defaults to 5.
#' @param colnames character string or character string vector specifying the
#' columns to summarise in the resulting quadtree.
#' @param threshold number. The threshold minimum value each cell must have
#' in the column \code{thresholdField}.
#' @param thresholdField character string specifying the column at which the
#' \code{threshold} value will apply.
#' @param funs character string or character string vector specifying the summary
#' functions for each of the \code{colnames}. If vector, the size must be the
#' same as colnames.
#' @param ineq.threshold inequality threshold value to be considered on the
#' disaggregation process. Forces disaggregation under the given inequality
#' threshold.
#' @param loss.threshold loss threshold value to be considered on the
#' disaggregation process. Forces aggregation when there's much loss
#' (i.e loss rate > ineq.threshold ).
#' @return AQuadtree object representing a varying size Quadtree
#' aggregation for the given points.
#' @export
#'
AQuadtree<-function(points, dim=1000, layers=5, colnames=NULL, threshold=100, thresholdField=NULL, funs=NULL, ineq.threshold=0.25, loss.threshold=0.4) {
  aquadtree<-createAQuadtree(points, dim, layers, colnames, threshold, thresholdField, funs, as="AQuadtree", ineq.threshold, loss.threshold)
}
#'
#'
#' Method show.
#' @rdname show
#' @title show AQuadtree-method
#' @details Display the AQuadtree object
#' @aliases show
#'
#' @param object an object of class AQuadtree.
#' @return A data.frame showing the information attributes contained in the AQuadtree object.
#'
setMethod("show", "AQuadtree",
          function(object){
            cat("An object of class \"",class(object),"\" with ", sep="")
            cat(length(object),
                "grid cells with sizes between",
                ifelse(object@dim>=1000, paste0(object@dim/1000, "km"), paste0(object@dim, "m")), "and",
                ifelse(object@dim/(2^(object@layers-1))>=1000, paste0(object@dim/(2^(object@layers-1))/1000, "km"), paste0(object@dim/(2^(object@layers-1)), "m")),
                "\n")
            print(object@data)
          }
)
#'
#'
#' Method print.
#' @rdname print
#' @title print AQuadtree-method
#' @details Prints the AQuadtree object
#' @aliases print
#'
#' @param x an object of class AQuadtree.
#' @param ... passed through.
#' @return none
#'
setMethod("print", "AQuadtree",
          function(x, ...){
            cat("* dim: "); print(x@dim)
            cat("* layers: "); print(x@layers)
            cat("* threshold: "); print(x@threshold)
            cat("* colnames: "); print(x@colnames)
            cat("* loss: "); print(x@loss)
            cat("* residual cells: "); print(length(x$residual))
            cat("* grid cells: "); print(length(x))
            print(x@data)
          }
)
#'
#'
#' Method summary.
#' @rdname summary
#' @title summary AQuadtree-method
#' @details summarize information of an object of class AQuadtree
#' @aliases summary
#'
#' @param object an object of class AQuadtree.
#' @param ... passed through.
#' @return An object of class "table" with summarising information in the AQuadtree input object
#'
setMethod("summary", "AQuadtree",
  function(object, ...) {
    cat("Object of class \"", class(object), "\"\n", sep="")
    cat(length(object),
        "grid cells with sizes between",
        ifelse(object@dim>=1000, paste0(object@dim/1000, "km"), paste0(object@dim, "m")),
        "and",
        ifelse(object@dim/(2^(object@layers-1))>=1000, paste0(object@dim/(2^(object@layers-1))/1000, "km"), paste0(object@dim/(2^(object@layers-1)), "m")),
        "\n")
    cat("Coordinates:\n")
    print(bbox(object))
    if (!is.na(is.projected(object))) cat("Is projected:", is.projected(object), "\n")
    else cat("Is projected: NA", "\n")
    if (!is.na(proj4string(object))) cat("proj4string:\n", proj4string(object), "\n")
    cat("Initial Cell Size:", ifelse(object@dim>=1000, paste0(object@dim/1000, "km"), paste0(object@dim, "m")), "\n")
    if (!is.null(object@data$residual)){
      cat("Number of valid grid Cells:", length(object[!object@data$residual,]), "\n")
      cat("Number of residual grid Cells:", length(object[object@data$residual,]), "\n")
    }
    cat("Data attributes:\n")
    print(summary(object@data[object@colnames]))
  }
)
#'
#'
#' Method [
#' @rdname extract.aquadtree.data
#' @title [ AQuadtree-method
#' @details Extract a part of a AQuadtree object
#' @aliases [
#'
#' @param x an object of class AQuadtree.
#' @param i,j elements to extract.
#' @param ... passed through.
#' @param drop passed on to [ indexing operator.
#' @return An AQuadtree object with the selected subset of rows or columns from the input object.
#'
setMethod("[", "AQuadtree",
          function(x, i, j, ..., drop){
            if (missing(j)){
              if (length(i)==1 && is.character(i)) {
                if(i=="dim"){return(x@dim)}
                else if(i=="layers"){return(x@layers)}
                else if(i=="threshold"){return(x@threshold)}
                else if(i=="thresholdField"){return(x@thresholdField)}
                else if(i=="colnames"){return(x@colnames)}
                else if(i=="loss"){return(x@loss)}
              } else {
                x.SP<-as(x, "SpatialPolygonsDataFrame")
                x.SP<-x.SP[i, , ...]
              }
            } else if (missing(i)) {
              x.SP<-as(x, "SpatialPolygonsDataFrame")
              x.SP<-x.SP[, j, ...]
            } else{
              x.SP<-as(x, "SpatialPolygonsDataFrame")
              x.SP<-x.SP[i, j, ...]
            }
            return(
              new("AQuadtree",
                  x.SP,
                  dim=x@dim,
                  layers=x@layers,
                  threshold=x@threshold,
                  thresholdField = x@thresholdField,
                  colnames= x@colnames[x@colnames %in% names(x.SP)],
                  loss=x@loss
              )
            )
          }
)
#'
#'
#' Method [<-
#' @rdname replace.aquadtree.data
#' @title [<- AQuadtree-method
#' @details An AQuadtree object cannot be assigned directly
#' @aliases [<-
#'
#' @param x an object of class AQuadtree.
#' @param i,j elements to extract or replace.
#' @param ... passed through.
#' @param value value to set.
#' @return none
#'
setReplaceMethod("[", "AQuadtree",
                 function(x, ...){
                   stop("Error: quadtree slots cannot be changed", call.=FALSE)
                 }
)
#'
#'
#' Method plot
#' @rdname plot
#' @title plot AQuadtree-method
#' @details Plot an object of class AQuadtree.
#' @aliases plot
#' @importFrom sp plot
#' @export
#'
#' @param x an object of class AQuadtree.
#' @param residual logical; if TRUE cells marked as residual cells are included
#' @param add logical. TRUE to add plot to the current existing plot
#' @param col default plotting color
#' @param ... passed through.
#' @return none
#'
setMethod("plot", signature = c(x="AQuadtree", y="missing"),
          function(x, ..., residual=TRUE, add=FALSE, col){
            if (residual) {
              if (missing(col)) {
                callGeneric(obj.SP<-as(x, "SpatialPolygonsDataFrame")[(x$residual),], add=add, col="red", xlim=x@bbox[1,], ylim=x@bbox[2,], ...)
                callGeneric(obj.SP<-as(x, "SpatialPolygonsDataFrame")[!(x$residual),], add=TRUE, col="green", ...)
              } else {
                callGeneric(obj.SP<-as(x, "SpatialPolygonsDataFrame"), add=add, col=col, ...)
              }
            } else {
              if (missing(col)) col="green"
              callGeneric(obj.SP<-as(x, "SpatialPolygonsDataFrame")[!(x$residual),], add=add, col=col, ...)
            }
          }
)
#'
#'
#' Method spplot
#' @rdname spplot
#' @title spplot AQuadtree-method
#' @details Plots a AQuadtree object as a spatial object with its data
#' @aliases spplot
#' @importFrom sp spplot
#' @export
#'
#' @param obj an object of class AQuadtree.
#' @param by.density logical; if TRUE cell values specified in zcol are divided by cell areas
#' @param zcol character; attribute name(s) or column number(s) in attribute table
#' @param residual logical; if TRUE cells marked as residual cells are included
#' @param ... passed through.
#' @return Creates a lattice plot of class "trellis" created with the spplot method in the sp package
#'
setMethod("spplot", "AQuadtree",
          function(obj, zcol=NULL, by.density=TRUE, residual=TRUE, ...){
            if (is.null(zcol)) zcol=obj@colnames
            if (residual) obj.SP<-as(obj, "SpatialPolygonsDataFrame")
            else obj.SP<-as(obj[!obj$residual], "SpatialPolygonsDataFrame")
            if (by.density) {
              obj.SP@data[,zcol]<-obj.SP@data[,zcol]*(1000*2^(obj.SP@data$level-1)/obj@dim)^2
              callGeneric(obj.SP, zcol=zcol, ...)
            } else callGeneric(obj.SP, zcol=zcol, ...)
          }
)
#'
#'
#' Method merge.
#' @rdname merge
#' @title Merge an AQuadtree object with a data.frame
#' @details Merges the AQuadtree object data with the data.frame on the columns "cellCode" and cellNum"
#' @aliases merge
#' @importFrom sp merge
#' @export
#'
#' @param x an object of class AQuadtree.
#' @param y an object of class data.frame
#' @return An AQuadtree object where the data is extended with the input data.frame
#'
setMethod("merge", signature = c(x="AQuadtree", y="data.frame"),
          function(x, y){
            if (!(all(c("cellCode","cellNum") %in% names(x))))
              stop("first object does not contain 'cellcode' and 'cellNum' attributes", call.=FALSE)
            if (!(all(c("cellCode","cellNum") %in% names(y))))
              stop("second object does not contain 'cellcode' and 'cellNum' attributes", call.=FALSE)
            return (
              new("AQuadtree",
                  sp::merge(as(x, "SpatialPolygonsDataFrame"), y, by=c("cellCode","cellNum"), sort=F),
                  dim = x@dim,
                  layers = x@layers,
                  threshold = as.numeric(NA),
                  thresholdField = as.character(NA),
                  colnames = c(x@colnames, names(y)[!names(y) %in% c('cellCode', 'cellNum')]),
                  loss = as.numeric(NA)
              )
            )
          }
)
#' Method area.QT
#' @rdname area.QT
#' @title area.QT AQuadtree-method
#' @details Get the areas of the Quadtree grid cells in square meters
#'
#' @aliases area.QT
#' @export
#'
#' @param obj an object of class AQuadtree.
#' @param residual logical; if TRUE cells marked as residual cells are included
#' @param ... passed through.
#'
setGeneric("area.QT", function(obj, residual=TRUE, ...) standardGeneric("area.QT"))

#' @return area of Quadtree grid cells in square meters
#' @rdname area.QT
#' @export
setMethod("area.QT", "AQuadtree",
          function(obj, residual=FALSE, ...){
            obj.SP<-as(obj, "SpatialPolygonsDataFrame")
            if (!residual) obj.SP<-obj.SP[!obj.SP$residual,]
            return((obj@dim/2^(obj.SP@data$level-1))^2)
          }
)

