#' @title Join two AQuadtree objects from the same area, to compare their data
#' @description
#' Given two objects of class AQuadtree for the same area, wich, for instance,
#' may contain data from two different periods, \code{joinAQuadtrees}
#' returns a new object of class AQuadtree with the common zones at the lowest
#' shared level, summarising the data from both AQuadtrees.
#' @importFrom methods new as
#' @importFrom stats weighted.mean
#' @importFrom sp SpatialPolygons Polygons Polygon CRS proj4string
#' proj4string<- SpatialPolygonsDataFrame spChFIDs
#' @importFrom dplyr summarise_at funs
#' @details
#' The function \code{joinAQuadtrees} creates a new AQuadtree object from two
#' given AQuadtree objects with data from the same area. The data of the
#' two given objects is summarised at the smallest possible cells shared by
#' both given objects. All the input data is maintained on the new created
#' object. This function can be used to join the different attributes from
#' the same area or information from different periods.
#' @param qt1 object of class "AQuadtree" containing the first object to join.
#' @param qt2 object of class "AQuadtree" containing the second object to join.
#' @param withResiduals logical indicating if \code{residual} cells should be
#' maintained (TRUE) or not (FALSE, default).
#' @param mean.1 character or character vector specifying the columns in the
#' first AQuadtreeto which a weighted mean should be computed. By default
#' the aggregation function used is \code{sum}.
#' @param mean.2 character or character vector specifying the columns in the
#' first AQuadtreeto which a weighted mean should be computed. By default
#' the aggregation function used is \code{sum}.
#' @return AQuadtree with the information of the two given objects summarised
#' at the lowest level shared by both objects.
#' @export
#' @examples
#' data("CharlestonPop")
#' CharlestonPop.AQT_1<-AQuadtree(CharlestonPop, layers = 2)
#' CharlestonPop.AQT_2<-AQuadtree(CharlestonPop, colnames="sex",
#'        thresholdField=c("sex.male", "sex.female"), layers = 2)
#' CharlestonPop.AQT_1_2<-joinAQuadtrees(CharlestonPop.AQT_1, CharlestonPop.AQT_2)
#'
#' \dontrun{
#' ## non AQuadtree objects
#' joinAQuadtrees(CharlestonPop, CharlestonCensusTracts)
#' }
joinAQuadtrees<-function(qt1, qt2, withResiduals=FALSE, mean.1=NULL, mean.2=NULL){
  .=NULL
  if (missing(qt1)) stop("argument 'qt1' is missing, with no default", call.="FALSE")
  if (missing(qt2)) stop("argument 'qt2' is missing, with no default", call.="FALSE")
  stopifnot(class(qt1)=="AQuadtree", class(qt2)=="AQuadtree", class(withResiduals)=="logical")
  stopifnot(is.projected(qt1), is.projected(qt2), proj4string(qt1)==proj4string(qt2))
  if (qt1@dim != qt2@dim) stop("initial dimensions of 'qt1' and 'qt2' differ", call.="FALSE")
  if (!(all(mean.1 %in% qt1@colnames))) {
    stop(sprintf("some 'mean.1' vars (%s) not in object names (%s)", paste(mean.1[!(mean.1 %in% qt1@colnames)] , collapse=", "), paste(qt1@colnames, collapse=", ")), call.="FALSE")
  }
  if (!(all(mean.2 %in% qt2@colnames))) {
    stop(sprintf("some 'mean.2' vars (%s) not in object names (%s)", paste(mean.2[!(mean.2 %in% qt2@colnames)] , collapse=", "), paste(qt2@colnames, collapse=", ")), call.="FALSE")
  }

  if (length(intersect(qt1$cellCode, qt2$cellCode))==0) stop("no common cells found", call.="FALSE")

  names(qt1)<-sapply(names(qt1), function(n){if (n %in% qt1@colnames) paste0(n,".1") else n}, simplify = TRUE, USE.NAMES = FALSE)
  qt1@colnames<-paste0(qt1@colnames,".1")
  names(qt2)<-sapply(names(qt2), function(n){if (n %in% qt2@colnames) paste0(n,".2") else n}, simplify = TRUE, USE.NAMES = FALSE)
  qt2@colnames<-paste0(qt2@colnames,".2")

  if (!is.null(mean.1)) {
    mean.1<-paste0(mean.1, ".1")
    sum.1<-qt1@colnames[!(qt1@colnames %in% mean.1)]
  } else sum.1<-qt1@colnames
  if (!is.null(mean.2)) {
    mean.2<-paste0(mean.2, ".2")
    sum.2<-qt2@colnames[!(qt2@colnames %in% mean.2)]
  } else  sum.2<-qt2@colnames


  qt.act<-SpatialPolygonsDataFrame(SpatialPolygons(list()), data=data.frame())
  proj4string(qt.act)<-proj4string(qt1)

  layerNumber<-max(qt1@layers, qt2@layers)
  cellCodes<-unique(union(qt1$cellCode, qt2$cellCode))
  for (mainCell in cellCodes) {
    qt1.act<-as(qt1[qt1$cellCode==mainCell,], "SpatialPolygonsDataFrame")
    qt2.act<-as(qt2[qt2$cellCode==mainCell,], "SpatialPolygonsDataFrame")
    if (length(qt1.act)==0 || length(qt2.act)==0) next

    if (length(qt1.act[qt1.act$cellNum=="" & !qt1.act$residual,])>0) {
      currentCell.sp<-as(qt1.act[qt1.act$cellNum=="" & !qt1.act$residual,], "SpatialPolygons")
      df1<-qt1.act[qt1.act$cellNum=="" & !qt1.act$residual,]@data
      if (is.null(mean.2))
        df2<-summarise_at(qt2.act@data, sum.2, funs(sum))
      else
        df2<-cbind(summarise_at(qt2.act@data, sum.2, funs(sum)), summarise_at(qt2.act@data, mean.2, funs(weighted.mean(., w=qt2.act$total.2))))

      qt.act<-rbind(qt.act, SpatialPolygonsDataFrame(currentCell.sp, data.frame(df1, df2, stringsAsFactors=FALSE), match.ID = FALSE))
      next
    } else if (length(qt2.act[qt2.act$cellNum=="" & !qt2.act$residual,])>0) {
        currentCell.sp<-as(qt2.act[qt2.act$cellNum=="" & !qt2.act$residual,], "SpatialPolygons")
        if (is.null(mean.1))
          df1<-summarise_at(qt1.act@data, sum.1, funs(sum))
        else
          df1<-cbind(summarise_at(qt1.act@data, sum.1, funs(sum)), summarise_at(qt1.act@data, mean.1, funs(weighted.mean(., w=qt1.act$total.1))))
        df2<-qt2.act[qt2.act$cellNum=="" & !qt2.act$residual,]@data
        qt.act<-rbind(qt.act, SpatialPolygonsDataFrame(currentCell.sp, data.frame(df1, df2, stringsAsFactors=FALSE), match.ID = FALSE))
        next
    }

    if (withResiduals) {
      if (length(qt1.act[qt1.act$residual,])>0) {
        if(length(qt2.act[qt2.act$residual,])>0) {
          qt.act<-rbind(qt.act,
                        SpatialPolygonsDataFrame(
                          as(qt1.act[qt1.act$residual,], "SpatialPolygons"),
                          data.frame(qt1.act[qt1.act$residual,]@data, qt2.act[qt2.act$residual,qt2@colnames], stringsAsFactors=FALSE), match.ID = FALSE))
        } else {
          auxdf<-as.data.frame(t(rep(0, length(qt2@colnames))))
          names(auxdf)<-qt2@colnames
          qt.act<-rbind(qt.act,
                        SpatialPolygonsDataFrame(
                          as(qt1.act[qt1.act$residual,], "SpatialPolygons"),
                          data.frame(qt1.act[qt1.act$residual,]@data, auxdf, stringsAsFactors=FALSE), match.ID = FALSE))
        }
      } else if (length(qt2.act[qt2.act$residual,])>0) {
        auxdf<-as.data.frame(t(rep(0, length(qt1@colnames))))
        names(auxdf)<-qt1@colnames
        qt.act<-rbind(qt.act,
                      SpatialPolygonsDataFrame(
                        as(qt2.act[qt2.act$residual,], "SpatialPolygons"),
                        data.frame(auxdf, qt2.act[qt2.act$residual,]@data, stringsAsFactors=FALSE), match.ID = FALSE))
      }
    }

    pos<-0
    for (i in 2:layerNumber) {
      pos<-pos+ceiling(log10(2^(2*(i-1))))

      for (currentCell in qt1.act[qt1.act$level==i,]$cellNum) {
        currentCell.sp<-as(qt1.act[qt1.act$cellNum==currentCell,], "SpatialPolygons")
        df2<-qt2.act@data[substr(qt2.act$cellNum,1,pos)==currentCell, ]
        if (nrow(df2)==0) next

        df1<-qt1.act[qt1.act$cellNum==currentCell,]@data

        if (is.null(mean.2))
          df2<-summarise_at(df2, sum.2, funs(sum))
        else
          df2<-cbind(
                summarise_at(df2, sum.2, funs(sum)),
                summarise_at(df2, mean.2, funs(weighted.mean(., w=df2[, 'total.2']))))
        qt.act<-rbind(qt.act, SpatialPolygonsDataFrame(currentCell.sp, data.frame(df1, df2, stringsAsFactors=FALSE), match.ID = FALSE))
      }

      qt2.act<-qt2.act[!(substr(qt2.act$cellNum, 1, pos)%in%qt1.act[qt1.act$level==i, ]$cellNum),]
      qt1.act<-qt1.act[qt1.act$level!=i, ]

      for (currentCell in qt2.act[qt2.act$level==i,]$cellNum) {
        currentCell.sp<-as(qt2.act[qt2.act$cellNum==currentCell,], "SpatialPolygons")
        df1<-qt1.act@data[substr(qt1.act$cellNum,1,pos)==currentCell, ]
        if (nrow(df1)==0) next

        df2<-qt2.act[qt2.act$cellNum==currentCell,]@data

        if (is.null(mean.1))
          df1<-summarise_at(df1, sum.1, funs(sum))
        else
          df1<-cbind(
                summarise_at(df1, sum.1, funs(sum)),
                summarise_at(df1, mean.1, funs(weighted.mean(., w=df1[,'total.1']))))
        qt.act<-rbind(qt.act, SpatialPolygonsDataFrame(currentCell.sp, data.frame(df1, df2, stringsAsFactors=FALSE), match.ID = FALSE))
      }

      qt1.act<-qt1.act[!(substr(qt1.act$cellNum, 1, pos)%in%qt2.act[qt2.act$level==i, ]$cellNum), ]
      qt2.act<-qt2.act[qt2.act$level!=i, ]

    }
  }
  qt.act<-spChFIDs(qt.act, as.character(1:length(qt.act)))
  qt.act@data<-qt.act@data[c("cellCode", "cellNum", "level", "residual", qt1@colnames, qt2@colnames)]

  return(
    new("AQuadtree",
        qt.act,
        dim=qt1@dim,
        layers=layerNumber,
        threshold=as.numeric(NA),
        thresholdField = as.character(NA),
        colnames= c(qt1@colnames, qt2@colnames),
        loss=as.numeric(NA)
    )
  )
}
