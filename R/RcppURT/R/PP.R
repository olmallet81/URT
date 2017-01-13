#==================================================================================================
#                     Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
#==================================================================================================

suppressMessages(library(R6))

#--------------------------------------------------------------------------------------------------

# C++ class PP<double> exposed to R
PP_d <- R6Class("PP_d",
  private = list(
    ptr = NULL
  ),
  public = list(
    initialize = function(data = NA, lags = "", lags_type = "long", trend = "c", test_type = "tau", regression = FALSE) {
      if (lags == "") {
        private$ptr <- .Call('PP_d__new2', PACKAGE = 'RcppURT', data_ = data, lags_type_ = lags_type, trend_ = trend, test_type_ = test_type, regression_ = regression)
      } else {
        private$ptr <- .Call('PP_d__new1', PACKAGE = 'RcppURT', data_ = data, lags_ = lags, trend_ = trend, test_type_ = test_type, regression_ = regression)
      }
    },
    statistic = function() {
      .Call('UnitRoot_d_statistic', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    pvalue = function() {
      .Call('UnitRoot_d_pvalue', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    show = function() {
      .Call('UnitRoot_d_show', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    }
  ),
  active = list(
    stat = function() {
      .Call('UnitRoot_d_get_stat', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    pval = function() {
      .Call('UnitRoot_d_get_pval', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    ols = function() {
      return(OLSlist_d(private$ptr))
    },
    trends = function() {
      .Call('UnitRoot_d_get_trends', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    bootstrap = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_d_get_bootstrap', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        .Call('UnitRoot_d_bootstrap', PACKAGE = 'RcppURT', ptr_ = private$ptr, bootstrap_ = value)
      }
    },
    lags = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_d_get_lags', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        .Call('UnitRoot_d_lags', PACKAGE = 'RcppURT', ptr_ = private$ptr, lags_ = value)
      }
    },
    lags_type = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_d_get_lags_type', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        Call('UnitRoot_d_lags_type', PACKAGE = 'RcppURT', ptr_ = private$ptr, lags_type_ = value)
      }
    },
    niter = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_d_get_niter', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        .Call('UnitRoot_d_niter', PACKAGE = 'RcppURT', ptr_ = private$ptr, niter_ = value)
      }
    },
    regression = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_d_get_regression', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        .Call('UnitRoot_d_regression', PACKAGE = 'RcppURT', ptr_ = private$ptr, regression_ = value)
      }
    },
    test_type = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_d_get_test_type', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        .Call('UnitRoot_d_test_type', PACKAGE = 'RcppURT', ptr_ = private$ptr, test_type_ = value)
      }
    },
    trend = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_d_get_trend', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        .Call('UnitRoot_d_trend', PACKAGE = 'RcppURT', ptr_ = private$ptr, trend_ = value)
      }
    }
  )
)

#--------------------------------------------------------------------------------------------------
###################################################################################################
#--------------------------------------------------------------------------------------------------

# C++ class PP<float> exposed to R
PP_f <- R6Class("PP_f",
  private = list(
    ptr = NULL
  ),
  public = list(
    initialize = function(data = NA, lags = "", lags_type = "long", trend = "c", test_type = "tau", regression = FALSE) {
      if (lags == "") {
        private$ptr <- .Call('PP_f__new2', PACKAGE = 'RcppURT', data_ = data, lags_type_ = lags_type, trend_ = trend, test_type_ = test_type, regression_ = regression)
      } 
      else {
        private$ptr <- .Call('PP_f__new1', PACKAGE = 'RcppURT', data_ = data, lags_ = lags, trend_ = trend, test_type_ = test_type, regression_ = regression)
      }
    },
    statistic = function() {
      .Call('UnitRoot_f_statistic', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    pvalue = function() {
      .Call('UnitRoot_f_pvalue', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    show = function() {
      .Call('UnitRoot_f_show', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    }
  ),
  active = list(
    stat = function() {
      .Call('UnitRoot_f_get_stat', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    pval = function() {
      .Call('UnitRoot_f_get_pval', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    ols = function() {
      return(OLSlist_f(private$ptr))
    },
    trends = function() {
      .Call('UnitRoot_f_get_trends', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    bootstrap = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_f_get_bootstrap', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        .Call('UnitRoot_f_bootstrap', PACKAGE = 'RcppURT', ptr_ = private$ptr, bootstrap_ = value)
      }
    },
    lags = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_f_get_lags', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        .Call('UnitRoot_f_lags', PACKAGE = 'RcppURT', ptr_ = private$ptr, lags_ = value)
      }
    },
    lags_type = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_f_get_lags_type', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } 
      else {
        Call('UnitRoot_f_lags_type', PACKAGE = 'RcppURT', ptr_ = private$ptr, lags_type_ = value)
      }
    },
    niter = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_f_get_niter', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        .Call('UnitRoot_f_niter', PACKAGE = 'RcppURT', ptr_ = private$ptr, niter_ = value)
      }
    },
    regression = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_f_get_regression', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        .Call('UnitRoot_f_regression', PACKAGE = 'RcppURT', ptr_ = private$ptr, regression_ = value)
      }
    },
    test_type = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_f_get_test_type', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        .Call('UnitRoot_f_test_type', PACKAGE = 'RcppURT', ptr_ = private$ptr, test_type_ = value)
      }
    },
    trend = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_f_get_trend', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        .Call('UnitRoot_f_trend', PACKAGE = 'RcppURT', ptr_ = private$ptr, trend_ = value)
      }
    }
  )
)

#--------------------------------------------------------------------------------------------------
