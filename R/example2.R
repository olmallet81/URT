#==================================================================================================
#                     Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
#==================================================================================================

suppressMessages(library(RcppURT))

run <- function()
{
  x = as.matrix(read.table("../data/x.csv", sep=","))
  y = as.matrix(read.table("../data/y.csv", sep=","))

  # running OLS regression as in ./examples/example1.cpp using double precision type
  fit = OLSreg_d(y, x, stats=TRUE, output=TRUE)

  # running OLS regression as in ./examples/example1.cpp using single precision type
  fit = OLSreg_f(y, x, stats=TRUE, output=TRUE)

  # running first ADF test as in ./examples/example2.cpp using double precision type
  test = ADFtest_d(y, lags=10, trend='ct', output=TRUE)

  # running second ADF test as in ./examples/example2.cpp using double precision type
  test = ADFtest_d(y, method='AIC', trend='ct', output=TRUE, bootstrap=TRUE, niter=10000)
}
