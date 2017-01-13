//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#include "../include/URT.hpp"

#include <RcppArmadillo.h>

using namespace urt;

//-------------------------------------------------------------------------------------------------

// C++ class OLS<double> exposed to Rcpp

// create an external pointer to an OLS<double> object
RcppExport SEXP OLS_d__new(SEXP y_, SEXP x_, SEXP stats_) 
{
   // convert inputs to appropriate C++ types
   Rcpp::NumericVector yr(y_);
   arma::vec y(yr.begin(), yr.size(), false);
   Rcpp::NumericMatrix xr(x_);
   arma::mat x(xr.begin(), xr.nrow(), xr.ncol(), false); 
   bool stats = Rcpp::as<bool>(stats_); 

   OLS<double>* p = nullptr;

   try {
      p = new OLS<double>(y, x, stats);
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   // create a pointer to an OLS<double> object and wrap it as an external pointer
   Rcpp::XPtr<OLS<double>> ptr(p, true);

   // return the external pointer to the R side
   return ptr;
}

RcppExport SEXP OLS_d_get_stats(SEXP ptr_, SEXP y_, SEXP x_)
{
   Rcpp::XPtr<OLS<double>> ptr(ptr_);

   Rcpp::NumericVector yr(y_);
   arma::vec y(yr.begin(), yr.size(), false);
   Rcpp::NumericMatrix xr(x_);
   arma::mat x(xr.begin(), xr.nrow(), xr.ncol(), false); 

   ptr->get_stats(y, x);

   return R_NilValue;
}

RcppExport SEXP OLS_d_show(SEXP ptr_)
{
   Rcpp::XPtr<OLS<double>> ptr(ptr_);

   try {
      ptr->show();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return R_NilValue;
}

RcppExport SEXP OLS_d_get_param(SEXP ptr_)
{
   Rcpp::XPtr<OLS<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->param);
}

RcppExport SEXP OLS_d_get_resid(SEXP ptr_)
{
   Rcpp::XPtr<OLS<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->resid);
}

RcppExport SEXP OLS_d_get_t_stat(SEXP ptr_)
{
   Rcpp::XPtr<OLS<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->t_stat);
}

RcppExport SEXP OLS_d_get_var(SEXP ptr_)
{
   Rcpp::XPtr<OLS<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->var);
}

RcppExport SEXP OLS_d_get_MSE(SEXP ptr_)
{
   Rcpp::XPtr<OLS<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->MSE);
}

RcppExport SEXP OLS_d_get_ndf(SEXP ptr_)
{
   Rcpp::XPtr<OLS<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->ndf);
}

RcppExport SEXP OLS_d_get_R2(SEXP ptr_)
{
   Rcpp::XPtr<OLS<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->R2);
}

RcppExport SEXP OLS_d_get_adj_R2(SEXP ptr_)
{
   Rcpp::XPtr<OLS<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->adj_R2);
}

RcppExport SEXP OLS_d_get_F_stat(SEXP ptr_)
{
   Rcpp::XPtr<OLS<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->F_stat);
}

RcppExport SEXP OLS_d_get_DW_stat(SEXP ptr_)
{
   Rcpp::XPtr<OLS<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->DW_stat);
}

RcppExport SEXP OLS_d_get_IC(SEXP ptr_)
{
   Rcpp::XPtr<OLS<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->IC);
}

//-------------------------------------------------------------------------------------------------
//#################################################################################################
//-------------------------------------------------------------------------------------------------

// C++ class OLS<float> exposed to Rcpp

// create an external pointer to an OLS<float> object
RcppExport SEXP OLS_f__new(SEXP y_, SEXP x_, SEXP stats_) 
{
   const arma::fvec& y = Rcpp::as<arma::fvec>(y_);
   const arma::fmat& x = Rcpp::as<arma::fmat>(x_);
   bool stats = Rcpp::as<bool>(stats_);    

   OLS<float>* p = nullptr;

   try {
      p = new OLS<float>(y, x, stats);
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   Rcpp::XPtr<OLS<float>> ptr(p, true);

   return ptr;
}

RcppExport SEXP OLS_f_get_stats(SEXP ptr_, SEXP y_, SEXP x_)
{
   Rcpp::XPtr<OLS<float>> ptr(ptr_);

   const arma::fvec& y = Rcpp::as<arma::fvec>(y_);
   const arma::fmat& x = Rcpp::as<arma::fmat>(x_); 

   ptr->get_stats(y, x);

   return R_NilValue;
}

RcppExport SEXP OLS_f_show(SEXP ptr_)
{
   Rcpp::XPtr<OLS<float>> ptr(ptr_);

   try {
      ptr->show();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return R_NilValue;
}

RcppExport SEXP OLS_f_get_param(SEXP ptr_)
{
   Rcpp::XPtr<OLS<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->param);
}

RcppExport SEXP OLS_f_get_resid(SEXP ptr_)
{
   Rcpp::XPtr<OLS<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->resid);
}

RcppExport SEXP OLS_f_get_t_stat(SEXP ptr_)
{
   Rcpp::XPtr<OLS<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->t_stat);
}

RcppExport SEXP OLS_f_get_var(SEXP ptr_)
{
   Rcpp::XPtr<OLS<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->var);
}

RcppExport SEXP OLS_f_get_MSE(SEXP ptr_)
{
   Rcpp::XPtr<OLS<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->MSE);
}

RcppExport SEXP OLS_f_get_ndf(SEXP ptr_)
{
   Rcpp::XPtr<OLS<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->ndf);
}

RcppExport SEXP OLS_f_get_R2(SEXP ptr_)
{
   Rcpp::XPtr<OLS<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->R2);
}

RcppExport SEXP OLS_f_get_adj_R2(SEXP ptr_)
{
   Rcpp::XPtr<OLS<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->adj_R2);
}

RcppExport SEXP OLS_f_get_F_stat(SEXP ptr_)
{
   Rcpp::XPtr<OLS<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->F_stat);
}

RcppExport SEXP OLS_f_get_DW_stat(SEXP ptr_)
{
   Rcpp::XPtr<OLS<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->DW_stat);
}

RcppExport SEXP OLS_f_get_IC(SEXP ptr_)
{
   Rcpp::XPtr<OLS<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->IC);
}

//-------------------------------------------------------------------------------------------------





