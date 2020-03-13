## ----setup, include = FALSE---------------------------------------------------
options(width=80)
knitr::opts_chunk$set(
  collapse = TRUE,
  warning=FALSE, 
  message=FALSE,
  fig.show='hold',
  tidy.opts=list(width.cutoff=80),
  tidy=TRUE,
  comment = "##"
)

library(knitr)
hook_output = knit_hooks$get('output')
knit_hooks$set(output = function(x, options) {
  # this hook is used only when the linewidth option is not NULL
  if (!is.null(n <- options$linewidth)) {
    x = knitr:::split_lines(x)
    # any lines wider than n should be wrapped
    x = unlist(sapply(x, function(x){
      if (nchar(x) > n) {
        paste(strwrap(x, width = n), collapse = paste0('\n', options$comment, ' '))
      } else {
        x
      }
    }, simplify = T, USE.NAMES = FALSE))
  }
  hook_output(x, options)
})


library(AQuadtree)

## ----echo=FALSE, fig.align='center', out.width="60%", fig.cap="\\label{fig:Figure 1}Three level quadtree splitting cell numbering example. Initial cell on the (left);  first quadtree subdivision (center); second quadtree subdivision (right)", fig.show='hold'----
knitr::include_graphics('images/Fig1.png')

## ----echo=FALSE, fig.align='center', out.width="25%", fig.cap="\\label{fig:Figure 2}Set of spatial points (a) and the corresponding 62.5m grid with no threshold restrictions (b) (the numbers indicate the points aggregated in each cell).", fig.subcap=rep("", 4), fig.show='hold'----
knitr::include_graphics(c('images/Fig2a.png','images/Fig2b.png'))

## ----echo=FALSE, fig.align='center', out.width="24%", fig.cap="\\label{fig:Figure 3}Disaggregation examples with threshold value 17. No disaggregation and no loss (a); disaggregation with suppression of 4 points (b) ; more disaggregation with suppression of 12 points (c); maximum disaggregation with suppression of 29 points (d).", fig.subcap=rep("", 4), fig.show='hold'----

knitr::include_graphics(c('images/Fig3a.png','images/Fig3b.png','images/Fig3c.png','images/Fig3d.png'))

## ----echo=FALSE, fig.align='center', out.width="28%", fig.cap="\\label{fig:Figure 4}Example of a residual cell.", fig.show='hold'----
knitr::include_graphics('images/Fig4.png')

## -----------------------------------------------------------------------------
example.QT<-AQuadtree(CharlestonPop)
class(example.QT)

## ----echo=2:4, fig.align='center', out.width="40%", fig.cap="AQuadtree plot and spplot"----
oldpar<-par(mar = c(0,0,0,0))
bcn.QT<-AQuadtree(BarcelonaPop)
plot(bcn.QT)
spplot(bcn.QT, by.density=TRUE)
par(oldpar)

## ---- linewidth=90------------------------------------------------------------
charleston.QT<-AQuadtree(CharlestonPop, dim = 10000, layers = 4)
summary(charleston.QT)

## ---- linewidth=90------------------------------------------------------------
class(BarcelonaPop$sex)
levels(BarcelonaPop$sex)
bcn.QT<-AQuadtree(BarcelonaPop, colnames = names(BarcelonaPop), funs = c('mean', 'sum'))
summary(bcn.QT)

## ---- linewidth=90------------------------------------------------------------
bcn.QT<-AQuadtree(BarcelonaPop, colnames = c('age','sex'), 
	funs = c('mean', 'sum'), threshold=17, 
	thresholdField=c("sex.man", "sex.woman"))
summary(bcn.QT)

## ----echo=2:5, fig.align='center', out.width="40%", fig.cap="\\label{fig:Figure 6}Examples of the effect of the ineq.threshold parameter."----
oldpar<-par(mar = c(0,0,0,0))
bcn.QT <- AQuadtree(BarcelonaPop, threshold = 5, ineq.threshold = 0.01)
plot(bcn.QT)
bcn.QT <- AQuadtree(BarcelonaPop, threshold = 5, ineq.threshold = 0.5)
plot(bcn.QT)
par(oldpar)

## ---- linewidth=90------------------------------------------------------------
bcn.QT<-AQuadtree(BarcelonaPop)
slotNames(bcn.QT)

## ---- linewidth=90------------------------------------------------------------
names(bcn.QT)
head(bcn.QT)

## ---- linewidth=90------------------------------------------------------------
data("BarcelonaPop", package = "AQuadtree")
summary(BarcelonaPop)

## ---- linewidth=90------------------------------------------------------------
data("CharlestonPop", package = "AQuadtree")
summary(CharlestonPop)

## -----------------------------------------------------------------------------
devtools::session_info("AQuadtree")

