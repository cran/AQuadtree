#'
#'
#'
#' Radomly created population points for Charleston, SC MSA, USA.
#'
#' A dataset containing randomly created population for Charleston, SC MSA,
#' USA.
#' Population has been created randomly with distributions of census tracts
#' from the dataset Charleston1,	2000 Census Tract Data for Charleston, SC
#' MSA and counties
#' @seealso https://spatial.uchicago.edu/sample-data
#'
#' @format A SpatialPointsDataFrame with 54619 rows and 3 attributes:
#' \describe{
#'   \item{age}{group age the individual as a factor with
#'    levels: "under16", "16_65", "over65"}
#'   \item{sex}{sex of the individual as a factor with levels: "male", "female"}
#'   \item{origin}{origin of the individual as a factor with
#'    levels: "asian", "black", "hisp", "multi_ra", "white"}
#' }
#'
"CharlestonPop"
#'
#'
#'
#' Census tract borders of Charleston, SC MSA, USA.
#'
#' A SpatiaPolygons object containing the Census tract borders of
#' Charleston, SC MSA, USA.
#'
#' @format A SpatialPolygons object with 117 polygons
#'
"CharlestonCensusTracts"
#'
#'
#'
#' Radomly created population points for Barcelona city in Catalonia.
#'
#' A dataset containing randomly created population for the Barcelona city
#' in Catalonia for the year 2018.
#' Population has been created randomly with the real distributions of census
#' tracts from the dataset dividing the total population by 20
#' (\url{https://www.bcn.cat/estadistica/catala/dades/tpob/pad/padro/a2018/edat/index.htm}).
#'
#' @format A SpatialPointsDataFrame with 81359 rows and 2 attributes:
#' \describe{
#'   \item{age}{age the individual}
#'   \item{sex}{sex of the individual as a factor with levels: "man", "woman"}
#' }
#'
"BarcelonaPop"
#'
#'
#'
#' Census tract borders of Barcelona city in Catalonia.
#'
#' A SpatiaPolygons object containing the Census tract borders of
#' Barcelona city in Catalonia.
#'
#' @format A SpatialPolygons object
#'
"BarcelonaCensusTracts"
