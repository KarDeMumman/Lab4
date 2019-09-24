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
