#==================================================================================================
#                     Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
#==================================================================================================

suppressMessages(library(R6))

#--------------------------------------------------------------------------------------------------

# C++ class KPSS<double> exposed to R
KPSS_d <- R6Class("KPSS_d",
  private = list(
    ptr = NULL
  ),
  public = list(
    initialize = function(data = NA, lags = "", lags_type = "long", trend = "c") {
      if (lags == "") {
        private$ptr <- .Call('KPSS_d__new2', PACKAGE = 'RcppURT', data_ = data, lags_type_ = lags_type, trend_ = trend)
      } else {
        private$ptr <- .Call('KPSS_d__new1', PACKAGE = 'RcppURT', data_ = data, lags_ = lags, trend_ = trend)
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
    trends = function() {
      .Call('UnitRoot_d_get_trends', PACKAGE = 'RcppURT', ptr_ = private$ptr)
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

# C++ class KPSS<float> exposed to R
KPSS_f <- R6Class("KPSS_f",
  private = list(
    ptr = NULL
  ),
  public = list(
    initialize = function(data = NA, lags = "", lags_type = "long", trend = "c") {
      if (lags == "") {
        private$ptr <- .Call('KPSS_f__new2', PACKAGE = 'RcppURT', data_ = data, lags_type_ = lags_type, trend_ = trend)
      } 
      else {
        private$ptr <- .Call('KPSS_f__new1', PACKAGE = 'RcppURT', data_ = data, lags_ = lags, trend_ = trend)
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
    trends = function() {
      .Call('UnitRoot_f_get_trends', PACKAGE = 'RcppURT', ptr_ = private$ptr)
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
