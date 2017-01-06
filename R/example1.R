#==================================================================================================
#                     Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
#==================================================================================================

suppressMessages(library(RcppURT))

run <- function()
{
  x = as.matrix(read.table("../data/x.csv", sep=","))
  y = as.matrix(read.table("../data/y.csv", sep=","))

  # running OLS regression as in ./examples/example1.cpp using double precision type
  fit = OLS_d$new(y, x, stats=TRUE)
  fit$show()

  # running OLS regression as in ./examples/example1.cpp using single precision type
  fit = OLS_f$new(y, x, stats=TRUE)
  fit$show()

  # running first ADF test as in ./examples/example2.cpp using double precision type
  test = ADF_d$new(y, lags=10, trend='ct')
  test$show()

  # running second ADF test as in ./examples/example2.cpp using double precision type
  test$method = 'AIC'
  test$bootstrap = TRUE
  test$niter = 10000
  test$show()
}
