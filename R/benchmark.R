#==================================================================================================
#                     Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
#==================================================================================================

suppressMessages(library(RcppURT))

run <- function()
{
  sizes = c(100,150,200,250,300,350,400,450,500,1000,1500,2000,2500,3000,3500,4000,4500,5000)

  for (i in 1:length(sizes)) {

    # generating Wiener process
    data = cumsum(rnorm(n=sizes[i]))

    if (sizes[i] < 1000) niter = 10000
    else niter = 1000
        
    # with R6 classes
    start1 = Sys.time()
    for (k in 1:niter) {
      test = ADF_d$new(data, method='AIC')
      # uncomment this line and comment the one above to switch to single precision
      #test = ADF_f$new(data, method='AIC')
      test$statistic()
    }
    end1 = Sys.time()

    # with Rcpp functions
    start2 = Sys.time()
    for (k in 1:niter) {
      test = ADFtest_d(data, method='AIC')
      # uncomment this line and comment the one above to switch to single precision
      #test = ADFtest_f(data, method='AIC')
    }
    end2 = Sys.time()

    cat(sprintf("%8d", sizes[i]))
    cat(sprintf("%8.1f", end1 - start1))
    cat(sprintf("%8.1f\n", end2 - start2))
  }
}
