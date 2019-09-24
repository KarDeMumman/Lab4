rm(list = ls())

linreg <- setRefClass("linreg",
                      contains = "input", #fill the characteristics of the data here
                      #what are the inheritence features from the parent calss.
                      # what types of input the class can inherit. (I think)
                      fields  = list(
                        formula <-"formula",
                        data <- "data.frame",
                        x <- "matrix",
                        y <- "vector",
                        coef_hat <- "matrix",
                        y_hat <- "vector",
                        res_hat <- "vector",
                        df <- "numeric"
                        ), 
                      methods = list(
                        initialize <- function(formula,data){
                          x <<- as.matrix(data)
                          y <<- as.vector()
                          #solving for the estimation coefficients.
                          #should be changed to QRd
                          coef_hat <<- as.matrix(solve(t(x)%*%x)%*%t(x)%*%y)
                          y_hat <<- as.vector()
                          res_hat <<- as.vector(y-y_hat)
                          df <<- as.numeric()
                          },
                        #use inline (\n) for printing. 
                        #cat does not do it by default
                        ### 'fill' and label lines:
                        #cat(paste(letters, 100* 1:26),"\n", fill = TRUE, labels = paste0("{", 1:10, "}:"))
                        
                        print <- function(){cat()},
                        coef <- function(){
                          coef <- as.vector(coef_hat)
                        names(coef) <- rownames(coef_hat)
                        return(coef)
                        },
                        plot <- function(){return()},
                        summary <- function(){return()},
                        resid <- function(){return(res_hat)},
                        pred <- function(){return(y_hat)}
                        )
                      )
