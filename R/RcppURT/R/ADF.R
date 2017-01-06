#==================================================================================================
#                     Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
#==================================================================================================

suppressMessages(library(R6))

#--------------------------------------------------------------------------------------------------

# C++ class ADF<double> exposed to R
ADF_d <- R6Class('ADF_d',
  private = list(
    ptr = NULL
  ),
  public = list(
    initialize = function(data = NA, lags = '', method = 'AIC', trend = 'c', regression = FALSE) {
      if (lags == '') {
        private$ptr <- .Call('ADF_d__new2', PACKAGE = 'RcppURT', data_ = data, method_ = method, trend_ = trend, regression_ = regression)
      } else {
        private$ptr <- .Call('ADF_d__new1', PACKAGE = 'RcppURT', data_ = data, lags_ = lags, trend_ = trend, regression_ = regression)
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
    level = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_d_get_level', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        .Call('UnitRoot_d_level', PACKAGE = 'RcppURT', ptr_ = private$ptr, level_ = value)
      }
    },
    max_lags = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_d_get_max_lags', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        .Call('UnitRoot_d_max_lags', PACKAGE = 'RcppURT', ptr_ = private$ptr, max_lags_ = value)
      }
    },
    method = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_d_get_method', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        .Call('UnitRoot_d_method', PACKAGE = 'RcppURT', ptr_ = private$ptr, method_ = value)
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

# C++ class ADF<float> exposed to R
ADF_f <- R6Class('ADF_f',
  private = list(
    ptr = NULL
  ),
  public = list(
    initialize = function(data = NA, lags = '', method = 'AIC', trend = 'c', regression = FALSE) {
      if (lags == '') {
        private$ptr <- .Call('ADF_f__new2', PACKAGE = 'RcppURT', data_ = data, method_ = method, trend_ = trend, regression_ = regression)
      } else {
        private$ptr <- .Call('ADF_f__new1', PACKAGE = 'RcppURT', data_ = data, lags_ = lags, trend_ = trend, regression_ = regression)
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
    level = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_f_get_level', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        .Call('UnitRoot_f_level', PACKAGE = 'RcppURT', ptr_ = private$ptr, level_ = value)
      }
    },
    max_lags = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_f_get_max_lags', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        .Call('UnitRoot_f_max_lags', PACKAGE = 'RcppURT', ptr_ = private$ptr, max_lags_ = value)
      }
    },
    method = function(value) {
      if (missing(value)) {
        .Call('UnitRoot_f_get_method', PACKAGE = 'RcppURT', ptr_ = private$ptr)
      } else {
        .Call('UnitRoot_f_method', PACKAGE = 'RcppURT', ptr_ = private$ptr, method_ = value)
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

