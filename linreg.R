rm(list = ls())

linreg <- setRefClass("linreg",
                      fields = list(
                        formula <-"formula",
                        input <- "data.frame",
                        x <- "matrix",
                        y <- "vector",
                        coef_hat <- "matrix",
                        y_hat <- "vector",
                        res_hat <- "vector",
                        df <- "numeric"
                        ), 
                      methods = list(
                        print <- function(){},
                        coef <- function(){},
                        plot <- function(){},
                        summary <- function(){},
                        resid <- function(){},
                        pred <- function(){}
                        )
                      )
#linreg <- function(x,y){
  #setting up the elements of matrix
  #x <- as.matrix(x)
  #y<- as.vector(y)
  #nc <- NCOL(x)
  #nr <- length(y)
  #QRdec <- qr(x)
  # compute (x'x)^(-1) x'y
  #coeff <- solve.qr(QRdec,y)
  # compute the degrees of freedom
  #DegFr <- nrow(x)-nrow(y)
  # compute standard deviation
  #StDev <- sum((y- x%*%coeff)**2)/DegFr}
# compute 