rm(list = ls())

linreg <- setRefClass("linreg",
                      fields  = list(
                        formula ="formula",
                        data = "data.frame",
                        X = "matrix",
                        y = "vector",
                        beta_hat = "matrix",
                        y_hat = "vector",
                        resids = "vector",
                        std_res = "vector",
                        degfree = "numeric",
                        sigma2 = "numeric",
                        variance = "numeric",
                        var_beta = "numeric",
                        t_value = "numeric",
                        p_value = "numeric"
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
                          Q <<- matrix(0, n,m)
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
                          std_res <<- as.vector(abs(scale(resids))) #standardized residuals
                          degfree <<- as.numeric(nrow(X)-ncol(X)) #for calculating the summary stat
                          sigma2 <<- as.numeric((t(resids)%*%as.matrix(resids))/degfree) #for calculating the summary stat
                          variance <<- as.numeric(sqrt(sigma2))
                          var_beta <<- as.numeric(sigma2 * solve((t(X)%*%X))) #finding the variance
                          t_value <<- as.numeric(beta_hat/(sqrt(diag(var_beta))))
                          p_value <<- as.numeric(2*pt(-abs(t_value),degfree))#getting the p_values
                        },
                        coef <- function(){return(beta_hat)},
                        print <- function(){
                          cat("\n","Call:","\n",paste("linreg(formula = ", format(formula), ", data = ", parse , ")\n\n ", sep = ""))
                          print.default(format(coef()), print.gap = 2L,quote = FALSE)
                          cat("\n")
                          },
                        plot <- function(){
                          #The required packages for plotting and setting up the 
                          #background image
                          LinReg::require("ggplot2")
                          LinReg::require("png")
                          plt_theme <- theme(
                            
                            plot.background = element_rect(fill = '#b8f5f9'),    # Background of the entire plot
                            
                            panel.background = element_rect(fill = '#54f4ff'),   # Background of plotting area
                            panel.border = element_rect(colour = "#c2c8d6", fill = NA),       # Border around plotting area.
                            # fill argument should be NA
                            panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                            colour = "white"),   # Major grid lines
                            panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                            colour = "white"),   # Minor grid lines
                            axis.line = element_line(colour = "black") 
                          )
                          
                          liu_logo <- png::readPNG("LiU_primary_white.png")
                          dat1 <- data.frame(y_hat, resids)
                          pl1<-ggplot(dat1, aes(y_hat, resids))+
                          geom_point()+
                          geom_smooth(method="lm", na.rm=TRUE, color="blue")+
                          xlab(paste("Fitted Values"))+
                          ylab("Residuals")+
                          ggtitle("Residuals vs Fitted Plot")+
                            plt_theme
                          
                          dat2 <- data.frame(y_hat, std_res)
                          pl2<-ggplot(dat, aes(y_hat, std_res))+
                            geom_point()+
                            geom_smooth(method="lm", na.rm=TRUE, color="blue")+
                            xlab(paste("Fitted Values"))+
                            ylab("Residuals")+
                            ggtitle("Standardized Residuals vs Fitted Plot")+
                            plt_theme
                          return(list(pl1,pl2))
                          },
                        summary <- function(){
                          cat("linreg(formula = ", 
                              format(formula), ", data = ", parse,
                              ") :\n\n ", sep = "")
                          summ_stat <- matrix(NA,nrow = 6)
                          summ_stat[,1] <- beta_hat
                          summ_stat[,2] <- var_beta
                          summ_stat[,3] <- t_value
                          summ_stat[,4] <- p_value
                          summ_stat[,5] <- deg_free
                          summ_stat[,6] <-for (i in 1:length(p_value)){
                            if (p_value[i]<0.001) significance[i]<-"***"
                            else if (p_value[i]<0.01) significance[i]<-"**"
                            else if (p_value[i]<0.05) significance[i]<-"*"
                            else if (p_value[i]<0.1) significance[i]<-"."
                            else significance[i]<-" "
                          }
                          colnames(summ_stat) <- c("Estimated values", "Standard Error of Estimates", "t_value", "Pr(>|t|)", "degrees of freedom", "significance")
                          cat(paste(summ_stat),"\n",
                              sep="\n")
                          },
                        resid <- function(){
                          return(resids)
                          },
                        pred <- function(){
                          return(y_hat)
                          }
                      )
)
