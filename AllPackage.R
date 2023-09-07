# data
SimpleLinearRegression<-function(feature,target){

  # b0 , b1
  # sum squares
  sx2 = sum((feature)^2)
  sy2 = sum((target)^2)
  k = 1

  # for sxx , sxy
  sum_xy = sum((feature) * (target))
  x_bar = mean(feature)
  y_bar = mean(target)
  n = length(feature)

  # sxx , syy , sxy
  sxy = sum_xy - n * x_bar * y_bar
  sxx = sx2 - n * (x_bar^2)
  syy = sy2 - n * (y_bar^2)

  # b0 , b1
  b1 = sxy / sxx
  b0 = y_bar - b1 * x_bar

  # ssr , sse
  SSR = (b1^2) * sxx
  SST = syy
  SSE = SST - SSR

  # for anova table
  MSE = SSE / (n-2)
  MSR = SSR

  cat("\nY =" , b0 , "+(" , b1, ")X")


  anova_table <- data.frame(
    Source = c("Regression" , " Residuals" ,"TOTAL"  ),
    SE     = c(SSR , SSE , SST),
    DOF    = c(k , n-2 , n-1),
    MS     = c(MSR , MSE , "" ),
    f      = c(MSR / MSE  , "" , "")
  )

  cat("\nSXX =" , sxx , " , SYY =" , syy , " , SXY =" , sxy )

  cat("\nB0 = " , b0  , " , B1 = " , b1)

  cat("\nSSR =" , SSR , " , SSE =" , SSE ,"\n"  )

  print ("#################Anova Table : ######################")

  print(anova_table)

  # scatter plot
  plot(feature, target,col='black',  main = "Scatter Plot with Fitted Line",
       col.main='black',xlab="input",ylab="output")
  # To fitted line to the scatter plot
  abline(a = b0, b = b1 , col='red',lwd=3)

  # equation
  cat("\nY =" , b0 , "+" , b1, "X\n")

  # 95% C.I for b0 estimator
  cat("Enter Confidence interval: ")

  # Read the input as a numeric value
  CI <- scan(n=1)
  alpha <- 1-CI

  margin0 <- qt(1-(alpha/2),df=n-2)*sqrt(MSE*((1/n)+((x_bar^2)/sxx)))

  lowerinterval <- b0 - margin0

  upperinterval <- b0 + margin0

  cat("\n95% confidence interval for 'b0' = [ ", lowerinterval, " , ", upperinterval, " ]")

  # 95% C.I for b1 estimator
  margin1 <- qt(0.975,df=n-2)*sqrt(MSE/sxx)

  lowerinterval <- b1 - margin1

  upperinterval <- b1 + margin1

  cat("\n95% confidence interval for 'b1 estimator' = [ ", lowerinterval, " , ", upperinterval, " ]")

  # 95% C.I for mean response
  x0= 10
  y0=b0 + b1*x0
  margin_mean_response <-  qt(0.975,df=n-2)*sqrt(MSE*((1/n)+((x0-x_bar)^2/sxx)))

  lowerinterval <- y0 - margin_mean_response

  upperinterval <- y0 + margin_mean_response

  cat("\n95% confidence interval for 'mean response' = [ ", lowerinterval, " , ", upperinterval, " ]")

  # 95% C.I for new observations
  margin_new_observations <-  qt(0.975,df=n-2)*sqrt(1+MSE*((1/n)+((x0-x_bar)^2/sxx)))

  lowerinterval <- y0 - margin_new_observations

  upperinterval <- y0 + margin_new_observations

  cat("\n95% confidence interval for 'new observations' = [ ", lowerinterval, " , ", upperinterval, " ]")

}
MultipleLinearRegression <- function(data) {
  #num of observations
  n <-  nrow(data)

  # b0 with ones
  b0_coef <- rep(1, n)

  # features
  features <- subset(data, select = c("bedrooms", "bathrooms","sqft_living","sqft_above"))
  y <- data$price
  # p
  k = ncol(features)
  p = k+1
  
  # Design Matrix
  X <-matrix(cbind(b0_coef, data$bedrooms , data$bathrooms,data$sqft_living,data$sqft_above), nrow = n)
  
  xtx <- t(X) %*% X
  xty <- t(X) %*% y
  yty <- t(y) %*% y
  
  #calculate beta matrix#
  beta <- solve(xtx) %*%  xty
  

  # predicted response values
  yhat <- X %*% beta         # n*1
  
  # SSE , SST ,SSR
  SSE <- yty -(t(beta) %*% xty)
  SST <- yty  - n * (sum(y)/n)^2
  SSR <- SST - SSE # sum of squared regression

  # variance-covariance matrix for beta
  vcov <-  sigma2*diag(solve(xtx))
  sigma2 <- c(SSE / (n - p)) # estimate of the error variance
  cat("\nvariance covariance matrix :\n" )
  print(vcov)
  conf_level <- 95 / 100 # confidence level
  alpha <- 1 - conf_level # significance level
  tval <- qt(1 - alpha/2, n - p) # t-value for confidence interval
  
  
  # CI for b
  beta_ci <- cbind(beta - tval*sqrt(vcov), beta + tval*sqrt(vcov)) # confidence interval for beta
  beta_ci <- data.frame(beta_ci)
  cat("\nconfidence Interval for B :\n")
  print(beta_ci)


  # anova
  MSR=SSR/k
  MSE=SSE/(n-p)
  anova_table <- data.frame(
    Source = c("Regression" , " Residuals" ,"TOTAL"  ),
    SE     = c(SSR , SSE , SST),
    DOF    = c(k , n-p , n-1),
    MS     = c(MSR , MSE , "" ),
    f      = c(MSR / MSE  , "" , "")
  )

  print ("#################  Anova Table  ######################")
  print(anova_table)



  F <- MSR / MSE  # calculate the F-statistic
  cat("\nF-statistic :\n", F)

  p_value <- pf(F, k, n-p, lower.tail = FALSE)   # calculate the p-value

  R_squared <- SSR / SST  # calculate the coefficient of determination
  cat("\nR_squared :\n",R_squared)

  #Determine the level of significance
  alpha <- 0.05

  #"lower.tail" is a logical value indicating whether to calculate the probability
  #of values greater than "F" (if FALSE)
  # hypothesis test for multiple linear regression
  if (p_value < alpha) {
    cat("\nReject H0\n")
  } else {
    cat("\nDo not Reject H0\n")
  }
  
  ### mean response ##
  x_node <- matrix(X[1,])
  
  mean_response<- X[1,] %*% beta # mean response
  
  
  margin_mean_response <-  tval* sqrt(MSE* t(x_node) %*% solve(t(X) %*% X) %*% x_node)
  

  lowerinterval <- mean_response - margin_mean_response
  upperinterval <- mean_response + margin_mean_response

  cat("\nMean Response Interval : (" , lowerinterval , "," , upperinterval , ")\n")


  ### new observation ##
  ##Random sample##
  new_obs <- c(1, 35, 70,12,12)
  y_new_obs <- new_obs%*% beta

  margin_new_observations <-  tval*sqrt(MSE*(1+t(x_node) %*% solve(t(X) %*% X) %*% x_node))

  n_lowerinterval <- y_new_obs - margin_mean_response
  n_upperinterval <- y_new_obs + margin_mean_response

  cat("\nNew Observation Interval: (" , n_lowerinterval , "," , n_upperinterval , ")\n")


  residuals <- data$price - yhat
  standard_errors <- c(residuals) / as.vector(sqrt(MSE))


  # Plot the standard errors against the fitted values of price
  plot(yhat, standard_errors, xlab = "Fitted values of price", ylab = "Standard errors")
  
  # Add horizontal lines at y = -3 and y = 3
  abline(h = c(-3, 3), lty = 2)
  
  # Add a title to the plot
  title("Standard errors vs. fitted values of price")





  # Plot the standard errors against the fitted values of price
  plot(data$bedrooms, standard_errors, xlab = "Fitted values of price", ylab = "Standard errors")
  
  # Add horizontal lines at y = -3 and y = 3
  abline(h = c(-3, 3), lty = 2)
  
  # Add a title to the plot
  title("Standard errors vs. bedrooms feature")




  ############################
  # Plot the standard errors against the fitted values of price
  plot(data$bathrooms, standard_errors, xlab = "Fitted values of price", ylab = "Standard errors")
  
  # Add horizontal lines at y = -3 and y = 3
  abline(h = c(-3, 3), lty = 2)
  
  # Add a title to the plot
  title("Standard errors vs. bathrooms feature")





  # Plot the standard errors against the fitted values of price
  plot(data$sqft_living, standard_errors, xlab = "Fitted values of price", ylab = "Standard errors")
  
  # Add horizontal lines at y = -3 and y = 3
  abline(h = c(-3, 3), lty = 2)
  
  # Add a title to the plot
  title("Standard errors vs. sqft_living feature")




  # Plot the standard errors against the fitted values of price
  plot(data$sqft_above, standard_errors, xlab = "Fitted values of price", ylab = "Standard errors")
  
  # Add horizontal lines at y = -3 and y = 3
  abline(h = c(-3, 3), lty = 2)
  
  # Add a title to the plot
  title("Standard errors vs. sqft_above feature")

}



# data
data_type <- readline(prompt = "Enter data type (csv, txt, xls, xlsx, sav, rds, json): ")
data_path <- readline(prompt = "Enter path to data file: ")
# C:/Users/Blu-Ray/Desktop/kc_house_data.csv
if (data_type == "csv" || data_type == "txt") {
  data <- read.csv(data_path)
} else if (data_type == "xls" || data_type == "xlsx") {
  data <- readxl::read_excel(data_path)
} else if (data_type == "sav") {
  data <- haven::read_sav(data_path)
} else if (data_type == "rds") {
  data <- readRDS(data_path)
} else if (data_type == "json") {
  data <- jsonlite::fromJSON(data_path)
}
head(data)

# for simple linear regression
feature = data$sqft_living
target = data$price
SimpleLinearRegression( feature ,target )


