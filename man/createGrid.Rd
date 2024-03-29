% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createGrid.R
\name{createGrid}
\alias{createGrid}
\title{Create a Grid grid covering a given geographic zone.}
\usage{
createGrid(zone, dim = 1000, intersect = TRUE, outline = FALSE)
}
\arguments{
\item{zone}{object of class "SpatialPoints", "SpatialPointsDataFrame",
"SpatialPolygons" or "SpatialPolygonsDataFrame" specifying the zone to
be covered by the grid.}

\item{dim}{a single integer specifying the initial cell sizes in meters, defaults to 1000.}

\item{intersect, }{logical, if TRUE the resulting grid will be
intersected with the given zone. If zone is of class SpatialPoints, only cells
containing points will be kept on the resulting grid. If zone is of
class SpatialPolygons, only cells inside or partially inside polygons
in zone will be kept on the resulting grid.
Defaults to TRUE}

\item{outline, }{logical, if TRUE the resulting grid will be
clipped with the outlines of the given zone. Only applicable if zone is of
class SpatialPolygons.
Defaults to FALSE}
}
\value{
SpatialPolygons dataset representing a grid with squared cells of
the given size.
}
\description{
\code{createGrid} returns a SpatialPolygons object representing a grid
covering a given geographic zone following the INSPIRE Specification on
Geographical Grid Systems. Each polygon will be identified with it's CellCode
code.
}
\details{
INSPIRE Specification on Geographical Grid Systems\cr
The objective of the coding system is to generate unique identifiers for each
point, for any of the recommended resolutions.\cr
The cellCode is a text string, composed of cell size and cell coordinates.
Cell codes start with the cell's size prefix. The cell size is denoted in meter (m)
for cell sizes below 1000m and kilometre (km) for cell sizes from 1000m and
above.\cr
Examples: a 100 meter cell has an identifier starting with “100m”, the
identifier  of a 10000 meter cell starts with “10km”.\cr
The coordinate part of the cell code reflects the distance of the lower left
grid cell corner from the false origin of the CRS. In order to reduce the
length of the string, Easting (E) and Northing (N) values are divided by
10^n (n is the number of zeros in the cell size value). Example for a cell
size of 10000 meters: The number of zeros in the cell size value is 4.
The resulting divider for Easting and Northing values is 10^4 = 10000.\cr
}
\examples{

data("BarcelonaPop")
BarcelonaPop.INSPIRE_GRID<-createGrid(BarcelonaPop)
plot(BarcelonaPop.INSPIRE_GRID)

\dontrun{
BarcelonaPop.INSPIRE_GRID.10km<-createGrid(BarcelonaPop, 10000, intersect=FALSE)
plot(BarcelonaPop.INSPIRE_GRID.10km)

data("BarcelonaCensusTracts")
Barcelona.INSPIRE_GRID<-createGrid(BarcelonaCensusTracts, outline=TRUE)
plot(Barcelona.INSPIRE_GRID)
}

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
