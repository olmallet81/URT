//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#include "../include/URT.hpp"

#include <RcppArmadillo.h>

using namespace urt;

//-------------------------------------------------------------------------------------------------

// ADF<double> test as a function for test for a given number of lags
RcppExport SEXP ADFtest_d_1(SEXP data_, SEXP lags_, SEXP trend_, SEXP output_, SEXP bootstrap_, SEXP niter_)
{
   Rcpp::NumericVector x(data_);
   arma::vec data(x.begin(), x.size(), false);
   int lags = Rcpp::as<int>(lags_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool output = Rcpp::as<bool>(output_);
   bool bootstrap = Rcpp::as<bool>(bootstrap_);
  
   ADF<double> test(data, lags, trend);

   if (bootstrap) {
      test.bootstrap = true;
      test.niter = Rcpp::as<int>(niter_);
   }

   try {
      if (output) test.show(); else test.pvalue();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::List::create(Rcpp::Named("stat") = test.get_stat(),
                             Rcpp::Named("pval") = test.get_pval(),
                             Rcpp::Named("lags") = test.lags); 
}

//-------------------------------------------------------------------------------------------------

// ADF<double> test as a function for test with lag length optimization
RcppExport SEXP ADFtest_d_2(SEXP data_, SEXP method_, SEXP trend_, SEXP output_, SEXP bootstrap_, SEXP niter_)
{
   Rcpp::NumericVector x(data_);
   arma::vec data(x.begin(), x.size(), false);
   const char* method = Rcpp::as<const char*>(method_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool output = Rcpp::as<bool>(output_);
   bool bootstrap = Rcpp::as<bool>(bootstrap_);

   ADF<double> test(data, method, trend);

   if (bootstrap) {
      test.bootstrap = true;
      test.niter = Rcpp::as<int>(niter_);
   }

   try {
      if (output) test.show(); else test.pvalue();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::List::create(Rcpp::Named("stat") = test.get_stat(),
                             Rcpp::Named("pval") = test.get_pval(),
                             Rcpp::Named("lags") = test.lags); 
}

//-------------------------------------------------------------------------------------------------

// DFGLS<double> test as a function for test for a given number of lags
RcppExport SEXP DFGLStest_d_1(SEXP data_, SEXP lags_, SEXP trend_, SEXP output_, SEXP bootstrap_, SEXP niter_)
{
   Rcpp::NumericVector x(data_);
   arma::vec data(x.begin(), x.size(), false);
   int lags = Rcpp::as<int>(lags_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool output = Rcpp::as<bool>(output_);
   bool bootstrap = Rcpp::as<bool>(bootstrap_);

   DFGLS<double> test(data, lags, trend);

   if (bootstrap) {
      test.bootstrap = true;
      test.niter = Rcpp::as<int>(niter_);
   }

   try {
      if (output) test.show(); else test.pvalue();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::List::create(Rcpp::Named("stat") = test.get_stat(),
                             Rcpp::Named("pval") = test.get_pval(),
                             Rcpp::Named("lags") = test.lags); 
}

//-------------------------------------------------------------------------------------------------

// DFGLS<double> test as a function for test with lag length optimization
RcppExport SEXP DFGLStest_d_2(SEXP data_, SEXP method_, SEXP trend_, SEXP output_, SEXP bootstrap_, SEXP niter_)
{
   Rcpp::NumericVector x(data_);
   arma::vec data(x.begin(), x.size(), false);
   const char* method = Rcpp::as<const char*>(method_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool output = Rcpp::as<bool>(output_);
   bool bootstrap = Rcpp::as<bool>(bootstrap_);

   DFGLS<double> test(data, method, trend);

   if (bootstrap) {
      test.bootstrap = true;
      test.niter = Rcpp::as<int>(niter_);
   }

   try {
      if (output) test.show(); else test.pvalue();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::List::create(Rcpp::Named("stat") = test.get_stat(),
                             Rcpp::Named("pval") = test.get_pval(),
                             Rcpp::Named("lags") = test.lags); 
}

//-------------------------------------------------------------------------------------------------

// PP<double> test as a function for test for a given number of lags
RcppExport SEXP PPtest_d_1(SEXP data_, SEXP lags_, SEXP trend_, SEXP test_type_, SEXP output_, SEXP bootstrap_, SEXP niter_)
{
   Rcpp::NumericVector x(data_);
   arma::vec data(x.begin(), x.size(), false);
   int lags = Rcpp::as<int>(lags_);
   const char* trend = Rcpp::as<const char*>(trend_);
   const char* test_type = Rcpp::as<const char*>(test_type_);
   bool output = Rcpp::as<bool>(output_);
   bool bootstrap = Rcpp::as<bool>(bootstrap_);

   PP<double> test(data, lags, trend, test_type);

   if (bootstrap) {
      test.bootstrap = true;
      test.niter = Rcpp::as<int>(niter_);
   }

   try {
      if (output) test.show(); else test.pvalue();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::List::create(Rcpp::Named("stat") = test.get_stat(),
                             Rcpp::Named("pval") = test.get_pval(),
                             Rcpp::Named("lags") = test.lags); 
}

//-------------------------------------------------------------------------------------------------

// PP<double> test as a function for test with default lag value
RcppExport SEXP PPtest_d_2(SEXP data_, SEXP lags_type_, SEXP trend_, SEXP test_type_, SEXP output_, SEXP bootstrap_, SEXP niter_)
{
   Rcpp::NumericVector x(data_);
   arma::vec data(x.begin(), x.size(), false);
   const char* lags_type = Rcpp::as<const char*>(lags_type_);
   const char* trend = Rcpp::as<const char*>(trend_);
   const char* test_type = Rcpp::as<const char*>(test_type_);
   bool output = Rcpp::as<bool>(output_);
   bool bootstrap = Rcpp::as<bool>(bootstrap_);

   PP<double> test(data, lags_type, trend, test_type);

   if (bootstrap) {
      test.bootstrap = true;
      test.niter = Rcpp::as<int>(niter_);
   }

   try {
      if (output) test.show(); else test.pvalue();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::List::create(Rcpp::Named("stat") = test.get_stat(),
                             Rcpp::Named("pval") = test.get_pval(),
                             Rcpp::Named("lags") = test.lags); 
}

//-------------------------------------------------------------------------------------------------

// KPSS<double> test as a function for test for a given number of lags
RcppExport SEXP KPSStest_d_1(SEXP data_, SEXP lags_, SEXP trend_, SEXP output_)
{
   Rcpp::NumericVector x(data_);
   arma::vec data(x.begin(), x.size(), false);
   int lags = Rcpp::as<int>(lags_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool output = Rcpp::as<bool>(output_);

   KPSS<double> test(data, lags, trend);

   try {
      if (output) test.show(); else test.pvalue();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::List::create(Rcpp::Named("stat") = test.get_stat(),
                             Rcpp::Named("pval") = test.get_pval(),
                             Rcpp::Named("lags") = test.lags); 
}

//-------------------------------------------------------------------------------------------------

// KPSS<double> test as a function for test with default lag value
RcppExport SEXP KPSStest_d_2(SEXP data_, SEXP lags_type_, SEXP trend_, SEXP output_)
{
   Rcpp::NumericVector x(data_);
   arma::vec data(x.begin(), x.size(), false);
   const char* lags_type = Rcpp::as<const char*>(lags_type_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool output = Rcpp::as<bool>(output_);

   KPSS<double> test(data, lags_type, trend);

   try {
      if (output) test.show(); else test.pvalue();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::List::create(Rcpp::Named("stat") = test.get_stat(),
                             Rcpp::Named("pval") = test.get_pval(),
                             Rcpp::Named("lags") = test.lags); 
}

//-------------------------------------------------------------------------------------------------
//#################################################################################################
//-------------------------------------------------------------------------------------------------

// ADF<float> test as a function for test for a given number of lags
RcppExport SEXP ADFtest_f_1(SEXP data_, SEXP lags_, SEXP trend_, SEXP output_, SEXP bootstrap_, SEXP niter_)
{
   const arma::fvec& data = Rcpp::as<arma::fvec>(data_);
   int lags = Rcpp::as<int>(lags_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool output = Rcpp::as<bool>(output_);
   bool bootstrap = Rcpp::as<bool>(bootstrap_);

   ADF<float> test(data, lags, trend);

   if (bootstrap) {
      test.bootstrap = true;
      test.niter = Rcpp::as<int>(niter_);
   }

   try {
      if (output) test.show(); else test.pvalue();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::List::create(Rcpp::Named("stat") = test.get_stat(),
                             Rcpp::Named("pval") = test.get_pval(),
                             Rcpp::Named("lags") = test.lags);  
}

//-------------------------------------------------------------------------------------------------

// ADF<float> test as a function for test with lag length optimization
RcppExport SEXP ADFtest_f_2(SEXP data_, SEXP method_, SEXP trend_, SEXP output_, SEXP bootstrap_, SEXP niter_)
{
   const arma::fvec& data = Rcpp::as<arma::fvec>(data_);
   const char* method = Rcpp::as<const char*>(method_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool output = Rcpp::as<bool>(output_);
   bool bootstrap = Rcpp::as<bool>(bootstrap_);

   ADF<float> test(data, method, trend);

   if (bootstrap) {
      test.bootstrap = true;
      test.niter = Rcpp::as<int>(niter_);
   }

   try {
      if (output) test.show(); else test.pvalue();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::List::create(Rcpp::Named("stat") = test.get_stat(),
                             Rcpp::Named("pval") = test.get_pval(),
                             Rcpp::Named("lags") = test.lags); 
}

//-------------------------------------------------------------------------------------------------

// DFGLS<double> test as a function for test for a given number of lags
RcppExport SEXP DFGLStest_f_1(SEXP data_, SEXP lags_, SEXP trend_, SEXP output_, SEXP bootstrap_, SEXP niter_)
{
   const arma::fvec& data = Rcpp::as<arma::fvec>(data_);
   int lags = Rcpp::as<int>(lags_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool output = Rcpp::as<bool>(output_);
   bool bootstrap = Rcpp::as<bool>(bootstrap_);

   DFGLS<float> test(data, lags, trend);

   if (bootstrap) {
      test.bootstrap = true;
      test.niter = Rcpp::as<int>(niter_);
   }

   try {
      if (output) test.show(); else test.pvalue();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::List::create(Rcpp::Named("stat") = test.get_stat(),
                             Rcpp::Named("pval") = test.get_pval(),
                             Rcpp::Named("lags") = test.lags);   
}

//-------------------------------------------------------------------------------------------------

// DFGLS<double> test as a function for test with lag length optimization
RcppExport SEXP DFGLStest_f_2(SEXP data_, SEXP method_, SEXP trend_, SEXP output_, SEXP bootstrap_, SEXP niter_)
{
   const arma::fvec& data = Rcpp::as<arma::fvec>(data_);
   const char* method = Rcpp::as<const char*>(method_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool output = Rcpp::as<bool>(output_);
   bool bootstrap = Rcpp::as<bool>(bootstrap_);

   DFGLS<float> test(data, method, trend);

   if (bootstrap) {
      test.bootstrap = true;
      test.niter = Rcpp::as<int>(niter_);
   }

   try {
      if (output) test.show(); else test.pvalue();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::List::create(Rcpp::Named("stat") = test.get_stat(),
                             Rcpp::Named("pval") = test.get_pval(),
                             Rcpp::Named("lags") = test.lags);  
}

//-------------------------------------------------------------------------------------------------

// PP<double> test as a function for test for a given number of lags
RcppExport SEXP PPtest_f_1(SEXP data_, SEXP lags_, SEXP trend_, SEXP test_type_, SEXP output_, SEXP bootstrap_, SEXP niter_)
{
   const arma::fvec& data = Rcpp::as<arma::fvec>(data_);
   int lags = Rcpp::as<int>(lags_);
   const char* trend = Rcpp::as<const char*>(trend_);
   const char* test_type = Rcpp::as<const char*>(test_type_);
   bool output = Rcpp::as<bool>(output_);
   bool bootstrap = Rcpp::as<bool>(bootstrap_);

   PP<float> test(data, lags, trend, test_type);

   if (bootstrap) {
      test.bootstrap = true;
      test.niter = Rcpp::as<int>(niter_);
   }

   try {
      if (output) test.show(); else test.pvalue();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::List::create(Rcpp::Named("stat") = test.get_stat(),
                             Rcpp::Named("pval") = test.get_pval(),
                             Rcpp::Named("lags") = test.lags); 
}

//-------------------------------------------------------------------------------------------------

// PP<double> test as a function for test with default lag value
RcppExport SEXP PPtest_f_2(SEXP data_, SEXP lags_type_, SEXP trend_, SEXP test_type_, SEXP output_, SEXP bootstrap_, SEXP niter_)
{
   const arma::fvec& data = Rcpp::as<arma::fvec>(data_);
   const char* lags_type = Rcpp::as<const char*>(lags_type_);
   const char* trend = Rcpp::as<const char*>(trend_);
   const char* test_type = Rcpp::as<const char*>(test_type_);
   bool output = Rcpp::as<bool>(output_);
   bool bootstrap = Rcpp::as<bool>(bootstrap_);

   PP<float> test(data, lags_type, trend, test_type);

   if (bootstrap) {
      test.bootstrap = true;
      test.niter = Rcpp::as<int>(niter_);
   }

   try {
      if (output) test.show(); else test.pvalue();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::List::create(Rcpp::Named("stat") = test.get_stat(),
                             Rcpp::Named("pval") = test.get_pval(),
                             Rcpp::Named("lags") = test.lags); 
}

//-------------------------------------------------------------------------------------------------

// KPSS<double> test as a function for test for a given number of lags
RcppExport SEXP KPSStest_f_1(SEXP data_, SEXP lags_, SEXP trend_, SEXP output_)
{
   const arma::fvec& data = Rcpp::as<arma::fvec>(data_);
   int lags = Rcpp::as<int>(lags_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool output = Rcpp::as<bool>(output_);

   KPSS<float> test(data, lags, trend);

   try {
      if (output) test.show(); else test.pvalue();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::List::create(Rcpp::Named("stat") = test.get_stat(),
                             Rcpp::Named("pval") = test.get_pval(),
                             Rcpp::Named("lags") = test.lags); 
}

//-------------------------------------------------------------------------------------------------

// KPSS<double> test as a function for test with default lag value
RcppExport SEXP KPSStest_f_2(SEXP data_, SEXP lags_type_, SEXP trend_, SEXP output_)
{
   const arma::fvec& data = Rcpp::as<arma::fvec>(data_);
   const char* lags_type = Rcpp::as<const char*>(lags_type_);
   const char* trend = Rcpp::as<const char*>(trend_);
   bool output = Rcpp::as<bool>(output_);

   KPSS<float> test(data, lags_type, trend);

   try {
      if (output) test.show(); else test.pvalue();
   } catch(std::exception &ex) {	
	  forward_exception_to_r(ex);
   }

   return Rcpp::List::create(Rcpp::Named("stat") = test.get_stat(),
                             Rcpp::Named("pval") = test.get_pval(),
                             Rcpp::Named("lags") = test.lags);  
}

//-------------------------------------------------------------------------------------------------

