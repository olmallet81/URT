#==================================================================================================
#                     Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
#==================================================================================================

#--------------------------------------------------------------------------------------------------

# OLS<double> as a list for unit-root tests
OLSlist_d <- function(ptr) {

    ols_ptr <- .Call('UnitRoot_d_get_ols', PACKAGE = 'RcppURT', ptr_ = ptr)

    return(list("param" = .Call('OLS_d_get_param', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
                "resid" = .Call('OLS_d_get_resid', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
                "t_stat" = .Call('OLS_d_get_t_stat', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
                "var" = .Call('OLS_d_get_var', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
                "MSE" = .Call('OLS_d_get_MSE', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
                "ndf" = .Call('OLS_d_get_ndf', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
                "R2" = .Call('OLS_d_get_R2', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
                "adj_R2" = .Call('OLS_d_get_adj_R2', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
                "F_stat" = .Call('OLS_d_get_F_stat', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
                "DW_stat" = .Call('OLS_d_get_DW_stat', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
				"IC" = .Call('OLS_d_get_IC', PACKAGE = 'RcppURT', ptr_ = ols_ptr)))
}

#--------------------------------------------------------------------------------------------------
###################################################################################################
#--------------------------------------------------------------------------------------------------

# OLS<float> as a list for unit-root tests
OLSlist_f <- function(ptr) {

    ols_ptr <- .Call('UnitRoot_f_get_ols', PACKAGE = 'RcppURT', ptr_ = ptr)

    return(list("param" = .Call('OLS_f_get_param', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
                "resid" = .Call('OLS_f_get_resid', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
                "t_stat" = .Call('OLS_f_get_t_stat', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
                "var" = .Call('OLS_f_get_var', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
                "MSE" = .Call('OLS_f_get_MSE', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
                "ndf" = .Call('OLS_f_get_ndf', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
                "R2" = .Call('OLS_f_get_R2', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
                "adj_R2" = .Call('OLS_f_get_adj_R2', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
                "F_stat" = .Call('OLS_f_get_F_stat', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
                "DW_stat" = .Call('OLS_f_get_DW_stat', PACKAGE = 'RcppURT', ptr_ = ols_ptr),
				"IC" = .Call('OLS_f_get_IC', PACKAGE = 'RcppURT', ptr_ = ols_ptr)))
}


#--------------------------------------------------------------------------------------------------
