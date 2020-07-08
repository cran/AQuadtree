#' @title Create a Quadtree grid to anonymise spatial point data
#' @description
#' \code{createAQuadtree} returns a SpatialPointsDataFrame representing a Quadtree
#' hierarchical geographic dataset. The resulting grid contains varying size cells
#' depending on a given threshold and column.
#' with identifiers
#' A \code{cellCode} and \code{cellNum} is created for each cell as in INSPIRE
#' Specification on Geographical Grid Systems.
#' @importFrom stats setNames
#' @importFrom methods new
#' @importFrom sp coordinates SpatialPolygons Polygons Polygon CRS proj4string
#'  SpatialPolygonsDataFrame
#' @importFrom dplyr %>% summarise_ group_by summarise_at
#' @details
#' Given a set of points a varying size Quadtree grid is created performing a
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
#' cellCode and cellNum for each cell in the Quadtree.
#' The objective of the coding system is to generate unique
#' identifiers for each cell, for any of the resolutions.\cr
#' The cellCode is a text string, composed of cell size and cell coordinates.
#' Cell codes start with a cell size prefix. The cell size is denoted in meter (m)
#' for cell sizes below 1000 m and kilometre (km) for cell sizes from 1000 m and
#' above.\cr
#' Examples: a 100 meter cell has an identifier starting with “100m”, the
#' identifier  of a 10000 meter cell starts with “10km”.\cr
#' The coordinate part of the cell code reflects the distance of the lower left
#' grid cell corner from the false origin of the CRS. In order to reduce the
#' length of the string, Easting (E) and Northing (N) values are divided by
#' 10^n (n is the number of zeros in the cell size value). Example for a cell
#' size of 10000 meters: The number of zeros in the cell size value is 4.
#' The resulting divider for Easting and Northing values is 10^4 = 10000.\cr
#' The cellNum is a sequence of concatenated integers identifying all the
#' hierarchical partitions of the main cell in which the point resides.
#' For instance, the cellNum of the top right cell would be 416 (fourth
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
#'
#' @param points object of class "SpatialPoints" or "SpatialPointsDataFrame".
#' @param dim a single integer specifying the initial cell sizes in meters, defaults to 1000.
#' @param layers a single integer specifying the number of divisions of the
#' initial cells, defaults to 5.
#' @param colnames character or character vector specifying the
#' columns to summarise in the resulting quadtree. For columns of class factor,
#' a column for each factor level cill be created.
#' @param threshold number. The threshold minimum value each cell must have
#' in the column \code{thresholdField}.
#' @param thresholdField character or character vector specifying the
#' columns to which the \code{threshold} value will apply. If not specified,
#' threshold value will be applied over the total cell points number.
#' ThresholdField must be one of the colnames.
#' @param funs character or character vector specifying the summary
#' functions for each of the \code{colnames}. If vector, the size must be the
#' same as colnames.
#' @param as character indicating return type, if "AQuadtree" a quadtree
#' class element will be returned, otherwise a SpatialPolygonsDataFrame
#' will ber returned. Defaults to "Spatial".
#' @param ineq.threshold inequality threshold value to be considered on the
#' disaggregation process. Forces disaggregation under the given inequality
#' threshold.
#' @param loss.threshold loss threshold value to be considered on the
#' disaggregation process. Stops disaggregation when there's much loss
#' (i.e loss rate > ineq.threshold ).
#' @return SpatialPolygonsDataFrame representing a varying size Quadtree
#' aggregation for the given points.
#' @export
#' @examples
#' data("CharlestonPop")
#' aQuadtree.Charleston<-createAQuadtree(CharlestonPop, threshold=10,
#'   colnames="sex", thresholdField=c("sex.male", "sex.female"))
#'
createAQuadtree <- function(points, dim=1000, layers=5, colnames=NULL, threshold=100, thresholdField=NULL, funs=NULL, as="Spatial", ineq.threshold=0.25, loss.threshold=0.4) {
  cellCode<-NULL
  ContainerID<-NULL
  ## aux function to determine the number of trailing zeros
  trailingZeros<-function(x){
    i<-0
    while (x %% 10 == 0 ) {
      x<-x%/%10
      i<-i+1
    }
    return(i)
  }
  ## f.Ineq inequality function
  f.Ineq<-function(x){
    x <- x[!(x == 0)]
    Th <- x/mean(x)
    Th <- sum(x * log(Th))
    Th <- Th/sum(x)
  }
  ## f.Loss cell loss percentage
  f.Loss<-function(x){
    L <- sum(x[x<threshold])/sum(x)
  }
  #stopifnot(require("sp"), require("dplyr"))
  if (missing(points)) stop("argument 'points' is missing, with no default", call.="FALSE")
  if (length(points)==0) stop("argument 'points' has length 0", call.="FALSE")
  stopifnot(dim>0, layers>=2)
  if (!inherits(points, "SpatialPoints")) stop("argument 'points' is not a 'SpatialPoints' or 'SpatialPointsDataFrame' object", call.="FALSE")
  if (!is.projected(points)) stop("spatial data must be projected", call.="FALSE")
  if (layers>10) stop("maximum 10 layers allowed", call.="FALSE")
  if (any(bbox(points)<0)) stop("negative bbox not permited, use a different projection", call.="FALSE")
  if (!(all(colnames %in% names(points)))) {
    stop(sprintf("some colnames (%s) not in object names (%s)", paste(colnames[!(colnames %in% names(points))] , collapse=", "), paste(names(points), collapse=", ")), call.="FALSE")
  }

  if (is.null(thresholdField)) {
    thresholdField<-"total"
  } else {
    # add thresholdField fields to selection (colnames)
    colnamesToAdd<-NULL
    for (f in thresholdField) {
      if (!(f %in% names(points))) {
        # control if fieldsToAdd come from factors
        f_<-unlist(strsplit(f, ".", fixed = TRUE))[1]
        if ((f_ %in% names(points))){
          colnamesToAdd<-c(colnamesToAdd, f_)
        } else {
          stop(sprintf("thresholdField (%s) not in object names (%s)", f_, paste(names(points), collapse=", ")), call.="FALSE")
        }
      } else {
        if (is.factor(points@data[,f])) stop(sprintf("thresholdField (%s) is a factor", f), call.="FALSE")
        if (!(f %in% colnames)) colnamesToAdd<-c(colnamesToAdd, f)
      }
    }
    colnamesToAdd<-colnamesToAdd[!duplicated(colnamesToAdd)]
    if (!is.null(funs) & length(funs)>1) funs<-c(funs, rep("sum", length(colnamesToAdd)))
    colnames<-c(colnames, colnamesToAdd)

  }
  if (length(funs)>1 & length(funs)!=length(colnames)) {
    stop("colnames and funs parameters do not have same length", call.="FALSE")
  }

  #create summarising expression
  summariseExpr<-c("total"="n()")

  if (!is.null(colnames)) {
    if (is.null(funs)) funs <- "sum"
    # treat possible factors within colnames
    summariseCols<-
        unlist(mapply(function(col, f){
            if (is.factor(points@data[,col])){
              lev<-levels(points@data[,col])
              newCols<-paste0(col, ".", lev)
              if (length(lev)>5) stop(sprintf("factor column %s has more than 5 levels", col), call.="FALSE")
              # decompose factor creating a new column for factor level
              points@data[newCols]<<-1*(points@data[rep(col, length(lev))]==as.list(lev))
              colnames<<-c(colnames, newCols)   # add new created columns to colnames
              return (setNames(paste0(f, "(", newCols, ")"), newCols))
            }else{
              return(setNames(paste0(f, "(", col, ")"), col))
            }
          }, colnames, funs, USE.NAMES = FALSE, SIMPLIFY = FALSE))
    colnames<-colnames[sapply(colnames, function(c) !is.factor(points@data[,c]))]
    colnames<-colnames[!duplicated(colnames)]

    summariseExpr<-c(summariseExpr, summariseCols)
    summariseExpr<-summariseExpr[!duplicated(summariseExpr)]
  }
  sizePrefix<-ifelse(dim>=1000, paste0(dim/1000, "km"), paste0(dim, "m"))

  # create points data.frame from input points' coordinates
  if (is.null(colnames)) {
    pts<-data.frame(x=coordinates(points)[,1], y=coordinates(points)[,2])   # there's no extra columns to keep
  } else {
    pts<-data.frame(x=coordinates(points)[,1], y=coordinates(points)[,2], points[colnames]@data)   #  keep extra columns given by colnames
  }
  # add x, y cell origin to each point
  pts$CellOrigin.x<-as.integer(pts$x%/%dim*dim)
  pts$CellOrigin.y<-as.integer(pts$y%/%dim*dim)
  # calculate string CellCode of the form "1kmNyyyyExxxx"
  zerosToRemove<-as.integer(10^trailingZeros(dim))
  cellCodeE<-pts$CellOrigin.x/zerosToRemove
  cellCodeN<-pts$CellOrigin.y/zerosToRemove
  lenCellCode<-nchar(max(cellCodeE, cellCodeN))
  pts$cellCode<-paste0(sizePrefix, "N",formatC(cellCodeN, width=lenCellCode, flag=0, mode = 'integer'), "E",formatC(cellCodeE, width=lenCellCode, flag=0, mode = 'integer'))
  # cellNumPos stores the number of digits for each subpart of the cell codes
  #   cellNumPosStart<-c( 1,14,15,17,19,22)
  #   cellNumPosStop <-c(13,14,16,18,21,25)
  cellNumPosStart<-c(1)
  cellNumPosStop<-c(nchar(pts$cellCode[1]))
  # calculate cellNum of the form "nmmooopppp..." with n 1:4   mm 01:16  ooo 001:256   pppp 0001:1024 ...
  pts$cellNum<-''
  for (i in 2:layers){
    zeros<-ceiling(log10(2^(2*(i-1))))
    cellNumPosStop[i]<-cellNumPosStop[i-1]+zeros
    cellNumPosStart[i]<-cellNumPosStop[i-1]+1
    size<-dim/2^(i-1)
    pts$cellNum<-paste0(pts$cellNum, formatC((pts$x-pts$CellOrigin.x)%/%size + (2^(i-1))*(pts$y-pts$CellOrigin.y)%/%size+1, width=zeros, flag=0, mode = 'integer'))
  }
  pts[c("x", "y", "CellOrigin.x", "CellOrigin.y")]<-NULL
  pts$cellCodesStr<-paste0(pts$cellCode, pts$cellNum)
  ## message("layer:1\n")
  prevGrid <- pts %>% group_by(cellCode) %>% summarise_(.dots=summariseExpr)

  prevGrid<-prevGrid[apply(prevGrid[,thresholdField]>=threshold,1, all),]   # remove high level cells with underthreshold population
  if (nrow(prevGrid)==0) stop("empty set, try a smaller threshold", call.="FALSE")
  prevGrid$level<-1

  #removed.Points<-pts[!(pts$cellCode %in% prevGrid$cellCode),]   # keep points before removing them
  removed.Points<-data.frame()

  pts<-pts[pts$cellCode %in% prevGrid$cellCode,]  # remove points corresponding to cells with underthreshold population
  maxLayer<-1
  # create elements of the quadtree aggregating at each level
  quadtree.Elements<-data.frame()
  for (i in 2:layers){
    ## message("layer:", i, "\n")
    actualcellNumPos<-cellNumPosStop[i]
    pts$cellCode<-substr(pts$cellCodesStr, 1, cellNumPosStop[i])

    actualGrid <- pts %>% group_by(cellCode) %>% summarise_(.dots=summariseExpr)
    actualGrid$level<-i
    maxLayer<-i
    actualGrid$ContainerID<-substr(actualGrid$cellCode, 1, cellNumPosStop[i-1])

    actualGrid<-actualGrid[actualGrid$ContainerID  %in% actualGrid[apply(actualGrid[,thresholdField]>=threshold,1, all),]$ContainerID,]
    if (nrow(actualGrid)==0) {
      break
    }
    codes.Ineq<-as.data.frame(actualGrid %>% group_by(ContainerID) %>% summarise_at(thresholdField,f.Ineq))
    ## do not aggregate when there's much inequality between cell population (i.e inequality index > ineq.threshold )
    if (length(thresholdField)>1)
      codes.Ineq.high<-codes.Ineq[apply(codes.Ineq[,thresholdField]>ineq.threshold,1, any), 1]
    else
      codes.Ineq.high<-codes.Ineq[codes.Ineq[,thresholdField]>ineq.threshold, 1]

    codes.Loss<-as.data.frame(actualGrid %>% group_by(ContainerID) %>% summarise_at(thresholdField,f.Loss))
    ## do not disaggregate when there's much Loss
    if (length(thresholdField)>1)
      codes.Loss.high<-codes.Loss[apply(codes.Loss[,thresholdField]>loss.threshold,1, any), 1]
    else
      codes.Loss.high<-codes.Loss[codes.Loss[,thresholdField]>loss.threshold, 1]

    codesToAggregate<-actualGrid[apply(actualGrid[, thresholdField] < threshold, 1, any),]$ContainerID

    codesToAggregate<-codesToAggregate[!(codesToAggregate %in% codes.Ineq.high) | (codesToAggregate %in% codes.Loss.high)]

    #keep cells discarded before removing them
    removed.Points<-rbind(removed.Points, pts[pts$cellCode %in% actualGrid[ !(actualGrid$ContainerID %in% codesToAggregate) & apply(actualGrid[, thresholdField] < threshold, 1, any),]$cellCode,])

    actualGrid<-actualGrid[ !(actualGrid$ContainerID %in% codesToAggregate) & apply(actualGrid[, thresholdField]>=threshold, 1, all), ]

    prevGrid<-prevGrid[!(prevGrid$cellCode %in% actualGrid$ContainerID),]

    actualGrid$ContainerID<-NULL

    quadtree.Elements<-rbind(quadtree.Elements, prevGrid)

    #remove all the points already in quadtree.Elements
    pts<-pts[pts$cellCode %in% actualGrid$cellCode,]

    prevGrid<-actualGrid

    if (nrow(pts)==0) break  #  finish when there's no points left
  }

  # add remaining elements from last layer
  quadtree.Elements<-rbind(quadtree.Elements, prevGrid)
  quadtree.Elements$residual<-FALSE

  if (nrow(removed.Points)>0) {
    # aggregate removed cells and add them to the quadtree
    removed.Points$cellCode<-substr(removed.Points$cellCode, 1, cellNumPosStop[1])
    removed.Points$ContainerID<-NULL

    removed.Points<-removed.Points %>% group_by(cellCode) %>% summarise_(.dots=summariseExpr)
    ## keep only aggregation of removed point over the threshold value
    removed.Points<-removed.Points[apply(removed.Points[thresholdField]>=threshold, 1, all),]
    if (nrow(removed.Points)>0) {
      removed.Points$level<-1
      removed.Points$residual<-TRUE
      quadtree.Elements<-rbind(quadtree.Elements, removed.Points)
    }
  }
  # convert quadtree elements to SpatialPolygons
  i<--1
  quadtree.SP<-SpatialPolygons(
    sapply(quadtree.Elements$cellCode,
           function(cellCode){
             i<<-i+1
             IDs<-substring(cellCode, cellNumPosStart, cellNumPosStop)
             IDs<-IDs[IDs!=""]
             IDs.length<-length(IDs)
             actualDim<-dim / 2^(IDs.length-1)
             elemNum<-as.integer(dim/actualDim)
             if (IDs.length==1){
               coords.x<-as.integer(strsplit(strsplit(IDs, "N")[[1]][2], "E")[[1]][2])*zerosToRemove
               coords.y<-as.integer(strsplit(strsplit(IDs, "N")[[1]][2], "E")[[1]][1])*zerosToRemove
             } else {
               coords.x<-as.integer(strsplit(strsplit(IDs, "N")[[1]][2], "E")[[1]][2]) * zerosToRemove + ((as.integer(IDs[IDs.length]) - 1) %% elemNum) * actualDim
               coords.y<-as.integer(strsplit(strsplit(IDs, "N")[[1]][2], "E")[[1]][1]) * zerosToRemove + ((as.integer(IDs[IDs.length]) - 1) %/% elemNum) * actualDim
             }
             pol.points<-cbind(
               x=c(coords.x, coords.x+actualDim, coords.x+actualDim, coords.x, coords.x),
               y=c(coords.y, coords.y, coords.y+actualDim, coords.y+actualDim, coords.y)
             )
             return (list( Polygons(list(Polygon(pol.points)), i)))
           }, USE.NAMES = FALSE), proj4string=CRS(proj4string(points), doCheckCRSArgs=TRUE))
  quadtree.Elements$cellNum<-substr(quadtree.Elements$cellCode, cellNumPosStart[2], cellNumPosStop[length(cellNumPosStop)])
  quadtree.Elements$cellCode<-substr(quadtree.Elements$cellCode, 1, cellNumPosStop[1])
  colOrder<-c("cellCode", "cellNum", "level", "residual", setdiff(colnames(quadtree.Elements), c("cellCode", "cellNum", "level", "residual")))
  quadtree.Elements<-quadtree.Elements[colOrder]
  quadtree<-SpatialPolygonsDataFrame(quadtree.SP, as.data.frame(quadtree.Elements), match.ID = FALSE)

  quadtree@bbox<-bbox(quadtree)
  if (as=="AQuadtree") {
    return(
      new("AQuadtree",
          quadtree,
          dim=dim,
          layers=maxLayer,
          colnames=c("total", colnames),
          threshold=threshold,
          thresholdField=thresholdField,
          loss=nrow(points)-sum(quadtree$total)
      )
    )
  } else {
    return(quadtree)
  }
}
