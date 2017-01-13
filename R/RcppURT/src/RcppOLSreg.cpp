//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#include "../include/URT.hpp"

#include <RcppArmadillo.h>

using namespace urt;

//-------------------------------------------------------------------------------------------------

// C++ class OLS<double> as a function
RcppExport SEXP OLSreg_d(SEXP y_, SEXP x_, SEXP stats_, SEXP output_)
{
   Rcpp::NumericVector yr(y_);
   arma::vec y(yr.begin(), yr.size(), false);
   Rcpp::NumericMatrix xr(x_);
   arma::mat x(xr.begin(), xr.nrow(), xr.ncol(), false);
   bool stats = Rcpp::as<bool>(stats_);
   bool output = Rcpp::as<bool>(output_);

   OLS<double> fit;

   try {
      fit = OLS<double>(y, x, stats);
      if (output) fit.show();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::List::create(Rcpp::Named("param") = fit.param,
                             Rcpp::Named("resid") = fit.resid,
                             Rcpp::Named("t_stat") = fit.resid,
                             Rcpp::Named("var") = fit.resid,
                             Rcpp::Named("MSE") = fit.MSE,
                             Rcpp::Named("ndf") = fit.ndf,
                             Rcpp::Named("R2") = fit.R2,
                             Rcpp::Named("adj_R2") = fit.adj_R2,
                             Rcpp::Named("F_stat") = fit.F_stat,
                             Rcpp::Named("DW_stat") = fit.DW_stat);
}

//-------------------------------------------------------------------------------------------------
//#################################################################################################
//-------------------------------------------------------------------------------------------------

// C++ class OLS<float> as a function
RcppExport SEXP OLSreg_f(SEXP y_, SEXP x_, SEXP stats_, SEXP output_)
{
   const arma::fvec& y = Rcpp::as<arma::fvec>(y_);
   const arma::fmat& x = Rcpp::as<arma::fmat>(x_);

   // NB: unfortunately for single precision type we cannot use the same method as with double precision type for converting SEXP objects into Armadillo arrays as we would have to cast NumericVector and NumericMatrix into float and then generate another copy, moreover R do not support float type so we would get NaN values

   bool stats = Rcpp::as<bool>(stats_);
   bool output = Rcpp::as<bool>(output_);

   OLS<float> fit;

   try {
      fit = OLS<float>(y, x, stats);
      if (output) fit.show();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::List::create(Rcpp::Named("param") = fit.param,
                             Rcpp::Named("resid") = fit.resid,
                             Rcpp::Named("t_stat") = fit.resid,
                             Rcpp::Named("var") = fit.resid,
                             Rcpp::Named("MSE") = fit.MSE,
                             Rcpp::Named("ndf") = fit.ndf,
                             Rcpp::Named("R2") = fit.R2,
                             Rcpp::Named("adj_R2") = fit.adj_R2,
                             Rcpp::Named("F_stat") = fit.F_stat,
                             Rcpp::Named("DW_stat") = fit.DW_stat);
 
}

//-------------------------------------------------------------------------------------------------
