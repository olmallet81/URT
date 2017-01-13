//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#include "../include/URT.hpp"

#include <RcppArmadillo.h>

using namespace urt;

//-------------------------------------------------------------------------------------------------

// C++ class UnitRoot<double> exposed to Rcpp

// access variable bootsrap
RcppExport SEXP UnitRoot_d_get_bootstrap(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->bootstrap);
}

// access and modify variable bootstrap
RcppExport SEXP UnitRoot_d_bootstrap(SEXP ptr_, SEXP bootstrap_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   ptr->bootstrap = Rcpp::as<bool>(bootstrap_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_d_get_lags(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->lags);
}

RcppExport SEXP UnitRoot_d_lags(SEXP ptr_, SEXP lags_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   ptr->lags = Rcpp::as<int>(lags_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_d_get_lags_type(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->lags_type);
}

RcppExport SEXP UnitRoot_d_lags_type(SEXP ptr_, SEXP lags_type_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   ptr->lags_type = Rcpp::as<const char*>(lags_type_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_d_get_level(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->level);
}

RcppExport SEXP UnitRoot_d_level(SEXP ptr_, SEXP level_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   ptr->level = Rcpp::as<float>(level_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_d_get_max_lags(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->max_lags);
}

RcppExport SEXP UnitRoot_d_max_lags(SEXP ptr_, SEXP max_lags_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   ptr->max_lags = Rcpp::as<int>(max_lags_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_d_get_method(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->method);
}

RcppExport SEXP UnitRoot_d_method(SEXP ptr_, SEXP method_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   ptr->method = Rcpp::as<const char*>(method_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_d_get_niter(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->niter);
}

RcppExport SEXP UnitRoot_d_niter(SEXP ptr_, SEXP niter_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   ptr->niter = Rcpp::as<int>(niter_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_d_get_regression(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->regression);
}

RcppExport SEXP UnitRoot_d_regression(SEXP ptr_, SEXP regression_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   ptr->regression = Rcpp::as<bool>(regression_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_d_get_test_type(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->test_type);
}

RcppExport SEXP UnitRoot_d_test_type(SEXP ptr_, SEXP test_type_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   ptr->test_type = Rcpp::as<const char*>(test_type_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_d_get_trend(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->trend);
}

RcppExport SEXP UnitRoot_d_trend(SEXP ptr_, SEXP trend_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   ptr->trend = Rcpp::as<const char*>(trend_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_d_get_stat(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->get_stat());
}

RcppExport SEXP UnitRoot_d_get_pval(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->get_pval());
}

RcppExport SEXP UnitRoot_d_get_ols(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   Rcpp::XPtr<OLS<double>> ols_ptr(new OLS<double>(ptr->get_ols()), true);
   return ols_ptr;
}

RcppExport SEXP UnitRoot_d_get_trends(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);
   return Rcpp::wrap(ptr->get_trends());
}

RcppExport SEXP UnitRoot_d_statistic(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);

   double stat = 0.0;

   try {
      stat = ptr->statistic();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::wrap(stat);
}

RcppExport SEXP UnitRoot_d_pvalue(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);

   double pval = 0.0;

   try {
      pval = ptr->pvalue();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::wrap(pval);
}

RcppExport SEXP UnitRoot_d_show(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<double>> ptr(ptr_);

   try {
      ptr->show();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return R_NilValue;
}

//-------------------------------------------------------------------------------------------------

// C++ class ADF<double> exposed to Rcpp

RcppExport SEXP ADF_d__new1(SEXP data_, SEXP lags_, SEXP trend_, SEXP regression_) 
{
   Rcpp::NumericVector x(data_);
   arma::vec data(x.begin(), x.size(), false);
   int lags = Rcpp::as<int>(lags_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool regression = Rcpp::as<bool>(regression_);

   Rcpp::XPtr<UnitRoot<double>> ptr(new ADF<double>(data, lags, trend, regression), true);

   return ptr;
}

RcppExport SEXP ADF_d__new2(SEXP data_, SEXP method_, SEXP trend_, SEXP regression_) 
{
   Rcpp::NumericVector x(data_);
   arma::vec data(x.begin(), x.size(), false);
   const char* method = Rcpp::as<const char*>(method_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool regression = Rcpp::as<bool>(regression_);

   Rcpp::XPtr<UnitRoot<double>> ptr(new ADF<double>(data, method, trend, regression), true);

   return ptr;
}

//-------------------------------------------------------------------------------------------------

// C++ class DFGLS<double> exposed to Rcpp

RcppExport SEXP DFGLS_d__new1(SEXP data_, SEXP lags_, SEXP trend_, SEXP regression_) 
{
   Rcpp::NumericVector x(data_);
   arma::vec data(x.begin(), x.size(), false);
   int lags = Rcpp::as<int>(lags_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool regression = Rcpp::as<bool>(regression_);

   Rcpp::XPtr<UnitRoot<double>> ptr(new DFGLS<double>(data, lags, trend, regression), true);

   return ptr;
}

RcppExport SEXP DFGLS_d__new2(SEXP data_, SEXP method_, SEXP trend_, SEXP regression_) 
{
   Rcpp::NumericVector x(data_);
   arma::vec data(x.begin(), x.size(), false);
   const char* method = Rcpp::as<const char*>(method_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool regression = Rcpp::as<bool>(regression_);

   Rcpp::XPtr<UnitRoot<double>> ptr(new DFGLS<double>(data, method, trend, regression), true);

   return ptr;
}

//-------------------------------------------------------------------------------------------------

// C++ class PP<double> exposed to Rcpp

RcppExport SEXP PP_d__new1(SEXP data_, SEXP lags_, SEXP trend_, SEXP test_type_, SEXP regression_) 
{
   Rcpp::NumericVector x(data_);
   arma::vec data(x.begin(), x.size(), false);
   int lags = Rcpp::as<int>(lags_);
   const char* trend = Rcpp::as<const char*>(trend_);
   const char* test_type = Rcpp::as<const char*>(test_type_);
   bool regression = Rcpp::as<bool>(regression_);

   Rcpp::XPtr<UnitRoot<double>> ptr(new PP<double>(data, lags, trend, test_type, regression), true);

   return ptr;
}

RcppExport SEXP PP_d__new2(SEXP data_, SEXP lags_type_, SEXP trend_, SEXP test_type_, SEXP regression_) 
{
   Rcpp::NumericVector x(data_);
   arma::vec data(x.begin(), x.size(), false);
   const char* lags_type = Rcpp::as<const char*>(lags_type_);
   const char* trend = Rcpp::as<const char*>(trend_);
   const char* test_type = Rcpp::as<const char*>(test_type_);
   bool regression = Rcpp::as<bool>(regression_);

   Rcpp::XPtr<UnitRoot<double>> ptr(new PP<double>(data, lags_type, trend, test_type, regression), true);

   return ptr;
}

//-------------------------------------------------------------------------------------------------

// C++ class KPSS<double> exposed to Rcpp

RcppExport SEXP KPSS_d__new1(SEXP data_, SEXP lags_, SEXP trend_) 
{
   Rcpp::NumericVector x(data_);
   arma::vec data(x.begin(), x.size(), false);
   int lags = Rcpp::as<int>(lags_);
   const char* trend = Rcpp::as<const char*>(trend_);

   Rcpp::XPtr<UnitRoot<double>> ptr(new KPSS<double>(data, lags, trend), true);

   return ptr;
}

RcppExport SEXP KPSS_d__new2(SEXP data_, SEXP lags_type_, SEXP trend_) 
{
   Rcpp::NumericVector x(data_);
   arma::vec data(x.begin(), x.size(), false);
   const char* lags_type = Rcpp::as<const char*>(lags_type_);
   const char* trend = Rcpp::as<const char*>(trend_);

   Rcpp::XPtr<UnitRoot<double>> ptr(new KPSS<double>(data, lags_type, trend), true);

   return ptr;
}

//-------------------------------------------------------------------------------------------------
//#################################################################################################
//-------------------------------------------------------------------------------------------------

// C++ class UnitRoot<float> exposed to Rcpp

// access variable bootsrap
RcppExport SEXP UnitRoot_f_get_bootstrap(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->bootstrap);
}

// access and modify variable bootstrap
RcppExport SEXP UnitRoot_f_bootstrap(SEXP ptr_, SEXP bootstrap_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   ptr->bootstrap = Rcpp::as<bool>(bootstrap_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_f_get_lags(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->lags);
}

RcppExport SEXP UnitRoot_f_lags(SEXP ptr_, SEXP lags_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   ptr->lags = Rcpp::as<int>(lags_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_f_get_lags_type(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->lags_type);
}

RcppExport SEXP UnitRoot_f_lags_type(SEXP ptr_, SEXP lags_type_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   ptr->lags_type = Rcpp::as<const char*>(lags_type_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_f_get_level(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->level);
}

RcppExport SEXP UnitRoot_f_level(SEXP ptr_, SEXP level_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   ptr->level = Rcpp::as<float>(level_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_f_get_max_lags(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->max_lags);
}

RcppExport SEXP UnitRoot_f_max_lags(SEXP ptr_, SEXP max_lags_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   ptr->max_lags = Rcpp::as<int>(max_lags_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_f_get_method(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->method);
}

RcppExport SEXP UnitRoot_f_method(SEXP ptr_, SEXP method_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   ptr->method = Rcpp::as<const char*>(method_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_f_get_niter(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->niter);
}

RcppExport SEXP UnitRoot_f_niter(SEXP ptr_, SEXP niter_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   ptr->niter = Rcpp::as<int>(niter_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_f_get_regression(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->regression);
}

RcppExport SEXP UnitRoot_f_regression(SEXP ptr_, SEXP regression_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   ptr->regression = Rcpp::as<bool>(regression_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_f_get_test_type(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->test_type);
}

RcppExport SEXP UnitRoot_f_test_type(SEXP ptr_, SEXP test_type_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   ptr->test_type = Rcpp::as<const char*>(test_type_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_f_get_trend(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->trend);
}

RcppExport SEXP UnitRoot_f_trend(SEXP ptr_, SEXP trend_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   ptr->trend = Rcpp::as<const char*>(trend_);
   return R_NilValue;
}

RcppExport SEXP UnitRoot_f_get_stat(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->get_stat());
}

RcppExport SEXP UnitRoot_f_get_pval(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->get_pval());
}

RcppExport SEXP UnitRoot_f_get_ols(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   Rcpp::XPtr<OLS<float>> ols_ptr(new OLS<float>(ptr->get_ols()), true);
   return ols_ptr;
}

RcppExport SEXP UnitRoot_f_get_trends(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);
   return Rcpp::wrap(ptr->get_trends());
}

RcppExport SEXP UnitRoot_f_statistic(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);

   float stat = 0.0;

   try {
      stat = ptr->statistic();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::wrap(stat);
}

RcppExport SEXP UnitRoot_f_pvalue(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);

   float pval = 0.0;

   try {
      pval = ptr->pvalue();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::wrap(pval);
}

RcppExport SEXP UnitRoot_f_show(SEXP ptr_)
{
   Rcpp::XPtr<UnitRoot<float>> ptr(ptr_);

   try {
      ptr->show();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return R_NilValue;
}

//-------------------------------------------------------------------------------------------------

// C++ class ADF<float> exposed to Rcpp

RcppExport SEXP ADF_f__new1(SEXP data_, SEXP lags_, SEXP trend_, SEXP regression_) 
{
   const arma::fvec& data = Rcpp::as<arma::fvec>(data_);
   int lags = Rcpp::as<int>(lags_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool regression = Rcpp::as<bool>(regression_);

   Rcpp::XPtr<UnitRoot<float>> ptr(new ADF<float>(data, lags, trend, regression), true);

   return ptr;
}

RcppExport SEXP ADF_f__new2(SEXP data_, SEXP method_, SEXP trend_, SEXP regression_) 
{
   const arma::fvec& data = Rcpp::as<arma::fvec>(data_);
   const char* method = Rcpp::as<const char*>(method_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool regression = Rcpp::as<bool>(regression_);

   Rcpp::XPtr<UnitRoot<float>> ptr(new ADF<float>(data, method, trend, regression), true);

   return ptr;
}

//-------------------------------------------------------------------------------------------------

// C++ class DFGLS<float> exposed to Rcpp

RcppExport SEXP DFGLS_f__new1(SEXP data_, SEXP lags_, SEXP trend_, SEXP regression_) 
{
   const arma::fvec& data = Rcpp::as<arma::fvec>(data_);
   int lags = Rcpp::as<int>(lags_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool regression = Rcpp::as<bool>(regression_);

   Rcpp::XPtr<UnitRoot<float>> ptr(new DFGLS<float>(data, lags, trend, regression), true);

   return ptr;
}

RcppExport SEXP DFGLS_f__new2(SEXP data_, SEXP method_, SEXP trend_, SEXP regression_) 
{
   const arma::fvec& data = Rcpp::as<arma::fvec>(data_);
   const char* method = Rcpp::as<const char*>(method_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool regression = Rcpp::as<bool>(regression_);

   Rcpp::XPtr<UnitRoot<float>> ptr(new DFGLS<float>(data, method, trend, regression), true);

   return ptr;
}

//-------------------------------------------------------------------------------------------------

// C++ class PP<float> exposed to Rcpp

RcppExport SEXP PP_f__new1(SEXP data_, SEXP lags_, SEXP trend_, SEXP test_type_, SEXP regression_) 
{
   const arma::fvec& data = Rcpp::as<arma::fvec>(data_);
   int lags = Rcpp::as<int>(lags_);
   const char* trend = Rcpp::as<const char*>(trend_);
   const char* test_type = Rcpp::as<const char*>(test_type_);
   bool regression = Rcpp::as<bool>(regression_);

   Rcpp::XPtr<UnitRoot<float>> ptr(new PP<float>(data, lags, trend, test_type, regression), true);

   return ptr;
}

RcppExport SEXP PP_f__new2(SEXP data_, SEXP lags_type_, SEXP trend_, SEXP test_type_, SEXP regression_) 
{
   const arma::fvec& data = Rcpp::as<arma::fvec>(data_);
   const char* lags_type = Rcpp::as<const char*>(lags_type_);
   const char* trend = Rcpp::as<const char*>(trend_);
   const char* test_type = Rcpp::as<const char*>(test_type_);
   bool regression = Rcpp::as<bool>(regression_);

   Rcpp::XPtr<UnitRoot<float>> ptr(new PP<float>(data, lags_type, trend, test_type, regression), true);

   return ptr;
}

//-------------------------------------------------------------------------------------------------

// C++ class KPSS<float> exposed to Rcpp

RcppExport SEXP KPSS_f__new1(SEXP data_, SEXP lags_, SEXP trend_) 
{
   const arma::fvec& data = Rcpp::as<arma::fvec>(data_);
   int lags = Rcpp::as<int>(lags_);
   const char* trend = Rcpp::as<const char*>(trend_);

   Rcpp::XPtr<UnitRoot<float>> ptr(new KPSS<float>(data, lags, trend), true);

   return ptr;
}

RcppExport SEXP KPSS_f__new2(SEXP data_, SEXP lags_type_, SEXP trend_) 
{
   const arma::fvec& data = Rcpp::as<arma::fvec>(data_);
   const char* lags_type = Rcpp::as<const char*>(lags_type_);
   const char* trend = Rcpp::as<const char*>(trend_);

   Rcpp::XPtr<UnitRoot<float>> ptr(new KPSS<float>(data, lags_type, trend), true);

   return ptr;
}

//-------------------------------------------------------------------------------------------------
