rm(list = ls())
#' @title Linear regression
#' Linear regression function.
#' @author c("Nastaran Meftahi", "Sarah Walid Alsaadi")
#' @param formula of linear regression as objects in linreg calss
#' @param data a dataset
#' @param formula = Petal.Length~Sepal.Width+Sepal.Length
#' @param data = iris
#' @return func<-linreg$new(formula,data)
#' @exportClass linreg
#' @export linreg
#' @return coef()
#' @return print()
#' @return resid()
#' @return plot()
#' @return summary()
#' @return pred()
#' @export
# export
# exportClass 


linreg <- setRefClass(Class = "linreg",
                      
                      fields  = list(
                        formula ="formula",
                        data = "data.frame",
                        X = "matrix",
                        y = "matrix",
                        beta_hat = "matrix",
                        y_hat = "matrix",
                        u_hat = "matrix",
                        std_u_hat = "vector",
                        df = "numeric",
                        sigma2 = "numeric",
                        variance = "numeric",
                        variance.beta = "matrix",
                        t_value = "numeric",
                        p_value = "numeric",
                        par = "character"
                      ),
                      methods = list(
                        
                        initialize = function(formula,data){
                          #model.matrix
                          stopifnot(is.data.frame(data))
                          stopifnot(all.vars(formula) %in% colnames(data))
                          # create informative labels for data sets and plots
                          data <<- data
                          formula <<- formula
                          #Setting up the model
                          X <<- model.matrix(formula, data)
                          #Setting up the dependent variable
                          ycols <-all.vars(formula)[1]
                          y <<- as.matrix(data[ycols])
                          par <<- deparse(substitute(data))
                          beta_hat <<- as.matrix(solve(t(X)%*%as.matrix(X))%*%as.matrix(t(X))%*%as.matrix(y))
                          y_hat <<- as.matrix(as.matrix(X)%*%as.matrix(beta_hat))
                          u_hat <<- as.matrix(y-y_hat) #for calculating the summary stat
                          std_u_hat <<- as.vector(-abs(scale(u_hat))) #standardized residuals
                          df <<- as.numeric(nrow(X)-ncol(X)) #for calculating the summary stat
                          sigma2 <<- as.numeric(as.matrix(t(u_hat)%*%as.matrix(u_hat))/df) #for calculating the summary stat
                          variance <<- as.numeric(sqrt(sigma2))
                          variance.beta <<- as.matrix(sigma2 * solve(as.matrix(t(X)%*%as.matrix((X))))) #finding the variance
                          t_value <<- as.numeric(beta_hat/(sqrt(diag(variance.beta))))
                          p_value <<- as.numeric(2*pt(-abs(t_value),df))#getting the p_values
                        },
                        coef = function(){
                          coef<- as.vector(beta_hat)
                          names(coef) <- rownames(beta_hat)
                          return(beta_hat)
                        },
                        
                        print = function(){
                          coef <- numeric()
                          cat("Call:\nlinreg(formula = ",
                              Reduce(paste,deparse(formula)),
                              ", data = ",
                              par,
                              ")\n\nCoefficients:\n ",sep="")
                          
                          for (i in 1:length(beta_hat)){
                            coef[i]<-beta_hat[i]
                          }
                          coef.names<-paste(" ",rownames(beta_hat)[1],"",sep="")
                          for (i in (2:length(coef))){
                            coef.names<-paste(coef.names," ",rownames(beta_hat)[i],sep="")
                          }
                          cat(coef.names)
                          cat("\n")
                          cat(coef,sep="")
                          #cat("\n","Call:","\n",paste("linreg(formula = ", format(formula), ", data = ", parse , ")\n\n ", sep = ""))
                          #print.default(format(coef()), print.gap = 2L,quote = FALSE)
                          #cat("\n")
                        },
                        plot = function(){
                          #The packages for plotting
                          library("ggplot2")
                          library("png")
                          library("magick")
                          #background image
                          plt_theme <- theme(
                            plot.background = element_rect(fill = '#b8f5f9',linetype = 2),    # Background of the entire plot
                            
                            panel.background = element_rect(fill = '#54f4ff'),   # Background of plotting area
                            panel.border = element_rect(colour = "#c2c8d6", fill = NA),       # Border around plotting area.
                            # fill argument should be NA
                            panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                            colour = "white"),   # Major grid lines
                            panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                            colour = "white"),   # Minor grid lines
                            axis.line = element_line(colour = "black")
                          )
                          liu_logo <- image_read_svg("https://upload.wikimedia.org/wikipedia/commons/5/57/Linkoping_University_Logo.svg")
                          dat1 <- data.frame(y_hat, u_hat)
                          pl1<-ggplot(dat1, aes(y_hat, u_hat))+
                            geom_point()+
                            geom_smooth(method="lm", na.rm=TRUE, color="blue")+
                            xlab(paste("Fitted Values"))+
                            ylab("Residuals")+
                            ggtitle("Residuals vs Fitted Plot")+
                            plt_theme
                          
                          dat2 <- data.frame(y_hat, std_u_hat)
                          pl2<-ggplot(dat2, aes(y_hat, std_u_hat))+
                            geom_point()+
                            geom_smooth(method="lm", na.rm=TRUE, color="blue")+
                            xlab(paste("Fitted Values"))+
                            ylab("Residuals")+
                            ggtitle("Standardized Residuals vs Fitted Plot")+
                            plt_theme
                          
                          return(list(pl1,pl2))
                        },
                        resid = function(){
                          return(u_hat)
                        },
                        pred = function(){
                          return(y_hat)
                        },
                        
                        summary = function(){
                          cat("linreg(formula = ", format(formula), ", data = ", par, ") :\n\n ", sep = "")
                          #result<- setNames(object = data.frame(cbind(beta_hat,sqrt(diag(variance.beta)),t_value,round(p_value,2),significance(p_value))),
                          #                nm =c("Estimates", "Std. Error","t value", "p value","Significance"))
                          result<-data.frame(cbind(round(beta_hat,2),round(sqrt(diag(variance.beta)),2),round(t_value,2),round(p_value,2),significance(p_value)))
                          colnames(result) <-c("Estimates", "Std. Error","t value", "p value","Significance")
                          tabl(result)
                          cat("\n Residual standard error: ", variance, " on ", df, " degrees of freedom", sep = "")
                        }
                      )
)

#' Significance
#' Helper function to encode significance based on a given p_value
#' @param p_value The given p_value
significance <- function(p_value) {
  significance <- ifelse(p_value < 0.001, "***",
                         ifelse(p_value < 0.01, "**",
                                (ifelse(p_value < 0.05, "*",
                                        (ifelse(p_value < 0.1, "."," "))))))
  return(significance)
}

#' Tabl
#' Helper function to print a given result
#' @param result The result to print
tabl<-function(result){
  print(result)
}