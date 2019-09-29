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
                          stopifnot(all.vars(formula) %in% colnames(data))
                          #Setting up the model
                          X <<- model.matrix(formula, data)
                          #Setting up the dependent variable
                          y <<- data[all.vars(formula)[1]]
                          #QR decomposition modified GramS-chmidt
                          #Setting up the Q-matrix
                        
                         
                          nR <- length(y)
                          pC <- NCOL(X)
                          #Creating the matrix of zeros (for all the data and variables)
                          Q <- matrix(0, nR,pC)
                          #creating the R matrix of size pC*pC
                          R <- matrix(0, pC, pC)
                          #Creating the matrix of residuals
                          E <- matrix(0, nR,pC)
                          #matrix of beta0s
                          E[,1] = X[,1]
                          
                          for (iter in 2:ncol(X)) {
                            E[, iter] = X[, iter]
                            for (num in seq(1, (iter - 1), 1)) {
                              E[, iter] = E[, iter] - ((sum(E[, num]*X[, iter]) / 
                                                          sum(E[, num]*E[, num])) * (E[, num]))
                              }
                          }
                          #l2 normalization
                          Q = apply(E, 2, function(X) { X / sqrt(sum(X*X)) })
                          R = t(Q) %*% X
                          
                          
                          coef_hat <<- as.matrix(backsolve(R,qb)) 
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
