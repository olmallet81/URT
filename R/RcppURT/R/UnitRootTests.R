#==================================================================================================
#                     Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
#==================================================================================================

#--------------------------------------------------------------------------------------------------

# C++ class ADF<double> as a function
ADFtest_d <- function(data, lags = "", method = "AIC", trend = "c", output = FALSE, bootstrap = FALSE, niter = 1000) {
  if (lags == "") {
    .Call('ADFtest_d_2', PACKAGE = 'RcppURT', data_ = data, method_ = method, trend_ = trend, output_ = output, bootstrap_ = bootstrap, niter_ = niter)
  } else {
    .Call('ADFtest_d_1', PACKAGE = 'RcppURT', data_ = data, lags_ = lags, trend_ = trend, output_ = output, bootstrap_ = bootstrap, niter_ = niter)
  }
}

#--------------------------------------------------------------------------------------------------

# C++ class DFGLS<double> as a function
DFGLStest_d <- function(data, lags = "", method = "AIC", trend = "c", output = FALSE, bootstrap = FALSE, niter = 1000) {
  if (lags == "") {
    .Call('DFGLStest_d_2', PACKAGE = 'RcppURT', data_ = data, method_ = method, trend_ = trend, output_ = output, bootstrap_ = bootstrap, niter_ = niter)
  } else {
    .Call('DFGLStest_d_1', PACKAGE = 'RcppURT', data_ = data, lags_ = lags, trend_ = trend, output_ = output, bootstrap_ = bootstrap, niter_ = niter)
  }
}

#--------------------------------------------------------------------------------------------------

# C++ class PP<double> as a function
PPtest_d <- function(data, lags = "", lags_type = "long", trend = "c", test_type = "tau", output = FALSE, bootstrap = FALSE, niter = 1000) {
  if (lags == "") {
    .Call('PPtest_d_2', PACKAGE = 'RcppURT', data_ = data, lags_type_ = lags_type, trend_ = trend, test_type_ = test_type, output_ = output, bootstrap_ = bootstrap, niter_ = niter)
  } else {
    .Call('PPtest_d_1', PACKAGE = 'RcppURT', data_ = data, lags_ = lags, trend_ = trend, test_type_ = test_type, output_ = output, bootstrap_ = bootstrap, niter_ = niter)
  }
}

#--------------------------------------------------------------------------------------------------

# C++ class KPSS<double> as a function
KPSStest_d <- function(data, lags = "", lags_type = "long", trend = "c", output = FALSE) {
  if (lags == "") {
    .Call('KPSStest_d_2', PACKAGE = 'RcppURT', data_ = data, lags_type_ = lags_type, trend_ = trend, output_ = output)
  } else {
    .Call('KPSStest_d_1', PACKAGE = 'RcppURT', data_ = data, lags_ = lags, trend_ = trend, output_ = output)
  }
}

#--------------------------------------------------------------------------------------------------
###################################################################################################
#--------------------------------------------------------------------------------------------------

# C++ class ADF<float> as a function
ADFtest_f <- function(data, lags = "", method = "AIC", trend = "c", output = FALSE, bootstrap = FALSE, niter = 1000) {
  if (lags == "") {
    .Call('ADFtest_f_2', PACKAGE = 'RcppURT', data_ = data, method_ = method, trend_ = trend, output_ = output, bootstrap_ = bootstrap, niter_ = niter)
  } else {
    .Call('ADFtest_f_1', PACKAGE = 'RcppURT', data_ = data, lags_ = lags, trend_ = trend, output_ = output, bootstrap_ = bootstrap, niter_ = niter)
  }
}

#--------------------------------------------------------------------------------------------------

# C++ class DFGLS<float> as a function
DFGLStest_f <- function(data, lags = "", method = "AIC", trend = "c", output = FALSE, bootstrap = FALSE, niter = 1000) {
  if (lags == "") {
    .Call('DFGLStest_f_2', PACKAGE = 'RcppURT', data_ = data, method_ = method, trend_ = trend, output_ = output, bootstrap_ = bootstrap, niter_ = niter)
  } else {
    .Call('DFGLStest_f_1', PACKAGE = 'RcppURT', data_ = data, lags_ = lags, trend_ = trend, output_ = output, bootstrap_ = bootstrap, niter_ = niter)
  }
}

#--------------------------------------------------------------------------------------------------

# C++ class PP<float> as a function
PPtest_f <- function(data, lags = "", lags_type = "long", trend = "c", test_type = "tau", output = FALSE, bootstrap = FALSE, niter = 1000) {
  if (lags == "") {
    .Call('PPtest_f_2', PACKAGE = 'RcppURT', data_ = data, lags_type_ = lags_type, trend_ = trend, test_type_ = test_type, output_ = output, bootstrap_ = bootstrap, niter_ = niter)
  } else {
    .Call('PPtest_f_1', PACKAGE = 'RcppURT', data_ = data, lags_ = lags, trend_ = trend, test_type_ = test_type, output_ = output, bootstrap_ = bootstrap, niter_ = niter)
  }
}

#--------------------------------------------------------------------------------------------------

# C++ class KPSS<float> as a function
KPSStest_f <- function(data, lags = "", lags_type = "long", trend = "c", output = FALSE) {
  if (lags == "") {
    .Call('KPSStest_f_2', PACKAGE = 'RcppURT', data_ = data, lags_type_ = lags_type, trend_ = trend, output_ = output)
  } else {
    .Call('KPSStest_f_1', PACKAGE = 'RcppURT', data_ = data, lags_ = lags, trend_ = trend, output_ = output)
  }
}

#--------------------------------------------------------------------------------------------------
