#==================================================================================================
#                     Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
#==================================================================================================

suppressMessages(library(R6))

#--------------------------------------------------------------------------------------------------

# C++ class OLS<double> exposed to R 
OLS_d <- R6Class("OLS_d",
  private = list(
    ptr = NULL
  ),
  public = list(
    initialize = function(y = NA, x = NA, stats = FALSE) {
      private$ptr <- .Call('OLS_d__new', PACKAGE = 'RcppURT', y_ = y, x_ = x, stats_ = stats)
    },
    get_stats = function(y = NA, x = NA) {
      .Call('OLS_d_get_stats', PACKAGE = 'RcppURT', ptr_ = private$ptr, y_ = y, x_ = x)
    },
    show = function() {
      .Call('OLS_d_show', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    }
  ),
  active = list(
    param = function() {
      .Call('OLS_d_get_param', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    resid = function() {
      .Call('OLS_d_get_resid', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    t_stat = function() {
      .Call('OLS_d_get_t_stat', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    var = function() {
      .Call('OLS_d_get_var', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    MSE = function() {
      .Call('OLS_d_get_MSE', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    ndf = function() {
      .Call('OLS_d_get_ndf', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    R2 = function() {
      .Call('OLS_d_get_R2', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    adj_R2 = function() {
      .Call('OLS_d_get_adj_R2', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    F_stat = function() {
      .Call('OLS_d_get_F_stat', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    DW_stat = function() {
      .Call('OLS_d_get_DW_stat', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    }
  )
)

#--------------------------------------------------------------------------------------------------
###################################################################################################
#--------------------------------------------------------------------------------------------------

# C++ class OLS<float> exposed to R
OLS_f <- R6Class("OLS_f",
  private = list(
    ptr = NULL
  ),
  public = list(
    initialize = function(y = NA, x = NA, stats = FALSE) {
      private$ptr <- .Call('OLS_f__new', PACKAGE = 'RcppURT', y_ = y, x_ = x, stats_ = stats)
    },
    get_stats = function(y = NA, x = NA) {
      .Call('OLS_f_get_stats', PACKAGE = 'RcppURT', ptr_ = private$ptr, y_ = y, x_ = x)
    },
    show = function() {
      .Call('OLS_f_show', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    }
  ),
  active = list(
    param = function() {
      .Call('OLS_f_get_param', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    resid = function() {
      .Call('OLS_f_get_resid', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    t_stat = function() {
      .Call('OLS_f_get_t_stat', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    var = function() {
      .Call('OLS_f_get_var', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    MSE = function() {
      .Call('OLS_f_get_MSE', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    ndf = function() {
      .Call('OLS_f_get_ndf', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    R2 = function() {
      .Call('OLS_f_get_R2', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    adj_R2 = function() {
      .Call('OLS_f_get_adj_R2', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    F_stat = function() {
      .Call('OLS_f_get_F_stat', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    },
    DW_stat = function() {
      .Call('OLS_f_get_DW_stat', PACKAGE = 'RcppURT', ptr_ = private$ptr)
    }
  )
)

#--------------------------------------------------------------------------------------------------

