rm(list = ls())

linreg <- setRefClass("linreg",
                      fields  = list(
                        formula ="formula",
                        data = "data.frame",
                        X = "matrix",
                        y = "vector",
                        coef_hat = "matrix",
                        y_hat = "vector",
                        resids = "vector",
                        degfree = "numeric",
                        #rsqrd <- "numeric",
                        sigma2 = "numeric",
                        variance = "numeric"
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
                          #QR decomposition GramS-chmidt
                          #Setting up the Q-matrix
                          n <- nrow(X)
                          m <- ncol(X)
                          #Creating the matrix of zeros (for all the data and variables)
                          Q <- matrix(0, n,m)
                          #creating the R matrix of size m*m
                          R <- matrix(0, m, m)
                          #Creating the projection matrix
                          V <- matrix(0, n,m)
                          V[,1] <- X[,1]
                          for (i in 2:m) {
                            V[,i]<- X[,i]
                            for (j in seq(1, (i-1),1)) {
                              #l2 normalization
                              V[,i] <- V[,i] - ((sum(V[,j]*V[,i]) /sum(V[,j]*V[,j])) * (V[,j]))
                              }
                          }
                          Q <- apply(V, 2, function(X) { X / sqrt(sum(X*X)) })
                          #upper triangle(values are approximately zero)
                          R <- t(Q) %*% X
                          yqt <- t(Q)%*%as.matrix(y)
                          #finding the betahats
                          beta_hat <<- as.vector(backsolve(R,yqt)) 
                          #using the round(R) changes the values of beta_hat
                          names(beta_hat) <<- colnames(X) 
                          y_hat <<- as.vector(X%*%beta_hat)
                          resids <<- as.vector(y-y_hat) #for calculating the summary stat
                          degfree <<- as.numeric(nrow(X)-ncol(X)) #for calculating the summary stat
                          sigma2 <<- as.numeric((t(resids)%*%as.matrix(resids))/degfree) #for calculating the summary stat
                          variance <<- as.numeric(sqrt(sigma2),2) #finding the variance
                          #tstat <<- as.vector()
                          # rsqrd <<- as.numeric() #for calculating the summary stat
                          #corrMatrix <<- as.matrix() #for calculating the summary stat
                        },
                        coef <- function(){return(beta_hat)},
                        print <- function(){
                          print.default(format(coef()), print.gap = 2L,quote = FALSE)
                          cat("\n")
                          },
                        plot <- function(){
                          require("ggplot2")
                          require("ggThemeAssist")
                          return()
                          },
                        summary <- function(){
                          return()
                          },
                        resid <- function(){
                          return(resids)
                          },
                        pred <- function(){
                          return(y_hat)
                          }
                      )
)
