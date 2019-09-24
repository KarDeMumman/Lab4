rm(list = ls())

linreg <- setRefClass("linreg",
                      fields  = list(
                        formula <-"formula",
                        data <- "data.frame",
                        X <- "matrix",
                        y <- "vector",
                        coef_hat <- "matrix",
                        y_hat <- "vector",
                        resids <- "vector",
                        df <- "numeric",
                        rsqrd <- "numeric",
                        corrMatrix <- "matrix",
                        fstat <- "vector",
                        sigma2 <- "numeric"
                      ), 
                      methods = list(
                        initialize <- function(formula,data){
                          #model.matrix
                          stopifnot(is.data.frame(data))
                          stopifnot(all.vars(formula %in% colnames(data)))
                          X <<- model.matrix(formula, data)
                          y <<- data[all.vars(forumla)[1]]
                          #solving for the estimation coefficients.
                          #should be changed to QRd
                          #for calculating the summary stat and the coefficients 
                          coef_hat <<- as.matrix(solve(t(x)%*%x)%*%t(x)%*%y) 
                          y_hat <<- as.vector(x%*%coef_hat)
                          resids <<- as.vector(y-y_hat) #for calculating the summary stat
                          df <<- as.numeric(nrow(x)-ncol(x)) #for calculating the summary stat
                          sigma2 <<- as.numeric() #for calculating the summary stat
                          fstat <<- as.vector() #for calculating the summary stat
                          rsqrd <<- as.numeric() #for calculating the summary stat
                          corrMatrix <<- as.matrix() #for calculating the summary stat
                        },
                        #use inline (\n) for printing. 
                        #cat does not do it by default
                        ### 'fill' and label lines:
                        #cat(paste(letters, 100* 1:26),"\n", fill = TRUE, labels = paste0("{", 1:10, "}:"))
                        print <- function(){cat()},
                        coef <- function(){
                          coef <- as.vector(coef_hat)
                          names(coef) <- colnames(X)
                          return(coef)
                        },
                        plot <- function(){
                          require("ggplot2")
                          return()},
                        summary <- function(){return()},
                        resid <- function(){return(resids)},
                        pred <- function(){return(y_hat)}
                      )
)
