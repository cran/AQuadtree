% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createAQuadtree.R
\name{createAQuadtree}
\alias{createAQuadtree}
\title{Create a Quadtree grid to anonymise spatial point data}
\usage{
createAQuadtree(
  points,
  dim = 1000,
  layers = 5,
  colnames = NULL,
  threshold = 100,
  thresholdField = NULL,
  funs = NULL,
  as = "Spatial",
  ineq.threshold = 0.25,
  loss.threshold = 0.4
)
}
\arguments{
\item{points}{object of class "SpatialPoints" or "SpatialPointsDataFrame".}

\item{dim}{a single integer specifying the initial cell sizes in meters, defaults to 1000.}

\item{layers}{a single integer specifying the number of divisions of the
initial cells, defaults to 5.}

\item{colnames}{character or character vector specifying the
columns to summarise in the resulting quadtree. For columns of class factor,
a column for each factor level cill be created.}

\item{threshold}{number. The threshold minimum value each cell must have
in the column \code{thresholdField}.}

\item{thresholdField}{character or character vector specifying the
columns to which the \code{threshold} value will apply. If not specified,
threshold value will be applied over the total cell points number.
ThresholdField must be one of the colnames.}

\item{funs}{character or character vector specifying the summary
functions for each of the \code{colnames}. If vector, the size must be the
same as colnames.}

\item{as}{character indicating return type, if "AQuadtree" a quadtree
class element will be returned, otherwise a SpatialPolygonsDataFrame
will ber returned. Defaults to "Spatial".}

\item{ineq.threshold}{inequality threshold value to be considered on the
disaggregation process. Forces disaggregation under the given inequality
threshold.}

\item{loss.threshold}{loss threshold value to be considered on the
disaggregation process. Stops disaggregation when there's much loss
(i.e loss rate > ineq.threshold ).}
}
\value{
SpatialPolygonsDataFrame representing a varying size Quadtree
aggregation for the given points.
}
\description{
\code{createAQuadtree} returns a SpatialPointsDataFrame representing a Quadtree
hierarchical geographic dataset. The resulting grid contains varying size cells
depending on a given threshold and column.
with identifiers
A \code{cellCode} and \code{cellNum} is created for each cell as in INSPIRE
Specification on Geographical Grid Systems.
}
\details{
Given a set of points a varying size Quadtree grid is created performing a
bottom-up aggregation considering a minimum threshold for each cell.
Cells with a value under the threshold for the \code{thresholdField} are
aggregated to the upper level in a quadtree manner.\cr
When no \code{thresholdField} is given, total number  of points in the cell
will be used, and so, given a threshold of k, none of the cells in the
resulting grid have a value less than k individuals as in a k-anonymity model.\cr
The Quadtree produced balances information loss and accuracy. For instance,
for the set of cells in the left image, where numbers in the cells represent
the values in the \code{thresholdField}, using a \code{threshold} value of 100,
the resulting Quadtree will be the one on the right. As we can see, some cells
will be discarded, and some aggregated to maintain as much information as
possible, keeping at the same time as much disaggregation as possible\cr
\if{html}{\figure{QTexampleA.png}{options: width=260 alt="62.5m2 cells"}}
\if{latex}{\figure{QTexampleA.png}{options: width=4.5cm}}
\if{html}{\figure{QTexampleB.png}{options: width=250 alt="resulting Quadtree"}}
\if{latex}{\figure{QTexampleB.png}{options: width=4.4cm}}\cr
The INSPIRE coding system for cell identifiers will be used to generate a
cellCode and cellNum for each cell in the Quadtree.
The objective of the coding system is to generate unique
identifiers for each cell, for any of the resolutions.\cr
The cellCode is a text string, composed of cell size and cell coordinates.
Cell codes start with a cell size prefix. The cell size is denoted in meter (m)
for cell sizes below 1000 m and kilometre (km) for cell sizes from 1000 m and
above.\cr
Examples: a 100 meter cell has an identifier starting with “100m”, the
identifier  of a 10000 meter cell starts with “10km”.\cr
The coordinate part of the cell code reflects the distance of the lower left
grid cell corner from the false origin of the CRS. In order to reduce the
length of the string, Easting (E) and Northing (N) values are divided by
10^n (n is the number of zeros in the cell size value). Example for a cell
size of 10000 meters: The number of zeros in the cell size value is 4.
The resulting divider for Easting and Northing values is 10^4 = 10000.\cr
The cellNum is a sequence of concatenated integers identifying all the
hierarchical partitions of the main cell in which the point resides.
For instance, the cellNum of the top right cell would be 416 (fourth
in first partition, sixteenth in second partition)\cr
The input object must be projected and units should be in 'meters'
because the system uses the INSPIRE coding system.
}
\examples{
data("CharlestonPop")
aQuadtree.Charleston<-createAQuadtree(CharlestonPop, threshold=10,
  colnames="sex", thresholdField=c("sex.male", "sex.female"))

}
\seealso{
\itemize{
  \item{
   D2.8.I.2 INSPIRE Specification on Geographical Grid Systems – Guidelines
   \url{https://inspire.ec.europa.eu/documents/Data_Specifications/INSPIRE_Specification_GGS_v3.0.1.pdf}
  }
  \item{
   EEA reference grid dataset
   \url{https://data.europa.eu/euodp/data/dataset/data_eea-reference-grids-2}
  }
}
}
