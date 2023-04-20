#' @title Create a Grid grid covering a given geographic zone.
#' @description
#' \code{createGrid} returns a SpatialPolygons object representing a grid
#' covering a given geographic zone following the INSPIRE Specification on
#' Geographical Grid Systems. Each polygon will be identified with it's CellCode
#' code.
#' @importFrom sp bbox GridTopology is.projected spChFIDs CRS proj4string
#'   as.SpatialPolygons.GridTopology
#' @details
#' INSPIRE Specification on Geographical Grid Systems\cr
#' The objective of the coding system is to generate unique identifiers for each
#' point, for any of the recommended resolutions.\cr
#' The cellCode is a text string, composed of cell size and cell coordinates.
#' Cell codes start with the cell's size prefix. The cell size is denoted in meter (m)
#' for cell sizes below 1000m and kilometre (km) for cell sizes from 1000m and
#' above.\cr
#' Examples: a 100 meter cell has an identifier starting with “100m”, the
#' identifier  of a 10000 meter cell starts with “10km”.\cr
#' The coordinate part of the cell code reflects the distance of the lower left
#' grid cell corner from the false origin of the CRS. In order to reduce the
#' length of the string, Easting (E) and Northing (N) values are divided by
#' 10^n (n is the number of zeros in the cell size value). Example for a cell
#' size of 10000 meters: The number of zeros in the cell size value is 4.
#' The resulting divider for Easting and Northing values is 10^4 = 10000.\cr
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
#' @param zone object of class "SpatialPoints", "SpatialPointsDataFrame",
#' "SpatialPolygons" or "SpatialPolygonsDataFrame" specifying the zone to
#' be covered by the grid.
#' @param dim a single integer specifying the initial cell sizes in meters, defaults to 1000.
#' @param intersect, logical, if TRUE the resulting grid will be
#' intersected with the given zone. If zone is of class SpatialPoints, only cells
#' containing points will be kept on the resulting grid. If zone is of
#' class SpatialPolygons, only cells inside or partially inside polygons
#' in zone will be kept on the resulting grid.
#' Defaults to TRUE
#' @param outline, logical, if TRUE the resulting grid will be
#' clipped with the outlines of the given zone. Only applicable if zone is of
#' class SpatialPolygons.
#' Defaults to FALSE
#' @return SpatialPolygons dataset representing a grid with squared cells of
#' the given size.
#' @export
#' @examples
#'
#' data("BarcelonaPop")
#' BarcelonaPop.INSPIRE_GRID<-createGrid(BarcelonaPop)
#' plot(BarcelonaPop.INSPIRE_GRID)
#'
#' \dontrun{
#' BarcelonaPop.INSPIRE_GRID.10km<-createGrid(BarcelonaPop, 10000, intersect=FALSE)
#' plot(BarcelonaPop.INSPIRE_GRID.10km)
#'
#' data("BarcelonaCensusTracts")
#' Barcelona.INSPIRE_GRID<-createGrid(BarcelonaCensusTracts, outline=TRUE)
#' plot(Barcelona.INSPIRE_GRID)
#' }
#'
createGrid <- function(zone, dim=1000, intersect=TRUE, outline=FALSE) {
    ## aux function to determine the number of trailing zeros
    trailingZeros<-function(x){
      i<-0
      while (x %% 10 == 0 ) {
        x<-x%/%10
        i<-i+1
      }
      return(i)
    }
    # requires sf
    if (outline) requireNamespace("sf", quietly = TRUE)
    if (missing(zone)) stop("argument 'zone' is missing, with no default", call.="FALSE")
    stopifnot(dim>0, inherits(zone, c("SpatialPolygons", "SpatialPoints")))
    zoneBbox<-bbox(zone)
    if (zoneBbox[1,"min"]<0 || zoneBbox[2,"min"]<0){
        stop("zone outside limits", call.="FALSE")
    }
    gridTop<-GridTopology(
                c(zoneBbox[1,"min"]%/%dim*dim+dim/2,zoneBbox[2,"min"]%/%dim*dim+dim/2),
                c(dim,dim),
                c((zoneBbox[1,"max"]-zoneBbox[1,"min"])%/%dim+2,(zoneBbox[2,"max"]-zoneBbox[2,"min"])%/%dim+2))
    if (!is.na(is.projected(zone)) && is.projected(zone)) {
      SPGrid.polygons<-as.SpatialPolygons.GridTopology(gridTop, slot(zone, "proj4string"))
    } else {
      SPGrid.polygons<-as.SpatialPolygons.GridTopology(gridTop)
    }

    SPGrid.sf <- sf::st_as_sf(SPGrid.polygons)

    if (intersect) {
      message("intersecting...\n")
      zone.sf<-sf::st_as_sf(zone)
      SPGrid.sf<-SPGrid.sf[sapply(sf::st_intersects(SPGrid.sf, zone.sf, sparse = T), any),]
    }
    if (outline & class(zone) %in% c("SpatialPolygons", "SpatialPolygonsDataFrame")){
      message("creating outline...\n")
      SPGrid.sf<-sf::st_intersection(sf::st_union(zone.sf), SPGrid.sf)
    }

    return(as(SPGrid.sf, "Spatial"))
}
