#' @title Add SpatialPoints to an AQuadtree obtject.
#' @description
#' Given an object of class AQuadtree and an object of class SpatialPoints or
#' SpatialPointsDataFrame for the same area, \code{pointsToAQuadtree}
#' returns a new object of class AQuadtree aggregating the data from the points
#' to the cells where each point fall.
#' @importFrom methods new as
#' @importFrom sp CRS proj4string identicalCRS
#' @importFrom dplyr %>% summarise_ group_by summarise_at funs
#' @details
#' The function \code{pointsToAQuadtree} returns a new AQuadtree object with
#' the input set of points aggregated to the input AQuadtree object. The function
#' creates a “p.total” attribute to compute the total
#' number of points aggregated to each cell of the input AQuadtree.
#' If points is an object of class SpatialPointsDataFrame, the function
#' summarises numeric attributes in the dataframe using the \code{mean}
#' function, and deploys factor attributes creating a new attribute for each label of the
#' factor to calculate the count. The attributes added to the resulting
#' AQuadtree object are prefixed with “p.”.
#' @param qt object of class "AQuadtree".
#' @param points object of class "SpatialPoints" or "SpatialPointsDataFrame".
#' @return AQuadtree with the information of the given set of points aggregated
#' at each corresponding cell of the given AQuadtree.
#' @export
#' @examples
#' data("BarcelonaPop")
#' Barcelona.QT<-AQuadtree(BarcelonaPop)
#' BcnWomen75yPop<-BarcelonaPop[BarcelonaPop$sex=='woman' & BarcelonaPop$age>=75, 'age']
#' Barcelona.extended.QT<-pointsToAQuadtree(Barcelona.QT, BcnWomen75yPop)
#'
#' \dontrun{
#' ## not an AQuadtree object
#' pointsToAQuadtree(CharlestonCensusTracts, CharlestonPop)
#'
#' ## spatial object not projected
#' sp.not.projected<-spTransform(CharlestonPop,CRS("+proj=longlat +datum=NAD27"))
#' is.projected(sp.not.projected)
#' pointsToAQuadtree(AQuadtree(CharlestonPop), sp.not.projected)
#'
#' }
pointsToAQuadtree<-function(qt, points){
  cellCode<-NULL
  if (missing(qt)) stop("argument 'qt' is missing, with no default", call.="FALSE")
  stopifnot(class(qt)=="AQuadtree")
  if (length(qt)==0) stop("argument 'qt' has length 0", call.="FALSE")
  if (missing(points)) stop("argument 'points' is missing, with no default", call.="FALSE")
  if (length(points)==0) stop("argument 'points' has length 0", call.="FALSE")
  if (!inherits(points, "SpatialPoints")) stop("argument 'points' is not a 'SpatialPoints' or 'SpatialPointsDataFrame' object", call.="FALSE")
  stopifnot(is.projected(qt), is.projected(points), identicalCRS(qt, points))
  if (any(bbox(points)<0)) stop("negative bbox not permited, use a different projection", call.="FALSE")


  #create summarising expression
  summariseExpr<-c("total"="n()")
  colnames<-names(points)

  if (!is.null(colnames)) {
    # treat possible factors within colnames
    summariseCols<-
      unlist(mapply(function(col){
                      if (is.numeric(points@data[,col])){
                        return(setNames(paste0('mean', "(", col, ")"), col))
                      } else if (is.factor(points@data[,col])){
                        lev<-levels(points@data[,col])
                        newCols<-paste0(col, ".", lev)
                        if (length(lev)>5) stop(sprintf("factor column %s has more than 5 levels", col), call.="FALSE")
                        # decompose factor creating a new column for factor level
                        points@data[newCols]<<-1*(points@data[rep(col, length(lev))]==as.list(lev))
                        colnames<<-c(colnames, newCols)   # add new created columns to colnames
                        return (setNames(paste0('sum', "(", newCols, ")"), newCols))
                      } else return(NULL)  # only aggregate factor and numeric attributes
                    }, colnames, USE.NAMES = FALSE, SIMPLIFY = FALSE))

    summariseExpr<-c(summariseExpr, summariseCols)
    summariseExpr<-summariseExpr[!duplicated(summariseExpr)]

    names(summariseExpr)<-paste0('p.',  names(summariseExpr))
  }
  maxLayers<-max(qt$level)
  points_ID<-as.data.frame(spatialPointsCellCodes(points, dim=slot(qt,'dim'),layers=maxLayers))
  points.agg <- data.frame(
    'cellCode' = character(),
    'cellNum' = character()
  )

  len <- 0
  for (i in 2:maxLayers) {
    len <- len + ceiling(log10(2^(2*(i-1))))
    aux <- points_ID %>%
      group_by(cellCode, 'cellNum'=substr(points_ID$cellNum, 1, len)) %>%
      summarise_(.dots=summariseExpr)

    aux<-do.call(data.frame,aux)
    aux<-aux[paste(aux$cellCode,aux$cellNum) %in% paste(qt$cellCode,qt$cellNum),]

    points_ID <- points_ID[!paste(points_ID$cellCode, substr(points_ID$cellNum,1,len)) %in% paste(aux$cellCode,aux$cellNum),]

    points.agg<-rbind(points.agg,aux)
  }

  # the rest of the point should aggregate to the cells of first level (residual or not)
  aux <- points_ID %>%
        group_by(cellCode, 'cellNum'=substr(points_ID$cellNum, 1, 0)) %>%
        summarise_(.dots=summariseExpr)
  aux<-do.call(data.frame,aux)
  points.agg<-rbind(points.agg,aux)

  merge(qt, points.agg)
  return (
    new("AQuadtree",
        sp::merge(as(qt, "SpatialPolygonsDataFrame"), points.agg, by=c("cellCode","cellNum"), sort=F),
        dim = qt@dim,
        layers = qt@layers,
        threshold = as.numeric(NA),
        thresholdField = as.character(NA),
        colnames = c(qt@colnames, names(points.agg)[!names(points.agg) %in% c('cellCode', 'cellNum')]),
        loss = as.numeric(NA)
    )
  )
}
