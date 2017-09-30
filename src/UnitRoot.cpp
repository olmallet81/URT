//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#include "../include/URT.hpp"

#ifdef _OPENMP 
  #include <omp.h>
  // getting maximum number of threads available
  static const int MAX_THREADS = omp_get_max_threads();
#endif

namespace urt {

//=================================================================================================

// constructor for running test for a given number of lags
template <typename T>
UnitRoot<T>::UnitRoot(const Vector<T>& data, int lags, const std::string& trend, bool regression)
{
   // copying data 
   this->data = data; 
   // setting pointer to data      
   set_data();
   #ifdef USE_ARMA
   nobs = data.n_rows;
   #elif defined(USE_BLAZE) || defined(USE_EIGEN) 
   nobs = data.size();
   #endif        
   this->lags = lags;
   this->trend = trend;
   this->regression = regression;
}

/*-------------------------------------------------------------------------------------------------*/
    
// constructor for running test with lag length optmization
template <typename T>
UnitRoot<T>::UnitRoot(const Vector<T>& data, const std::string& method, const std::string& trend, bool regression)
{
   // copying data 
   this->data = data;
   // setting pointer to data
   set_data();
   #ifdef USE_ARMA
   nobs = data.n_rows;
   #elif defined(USE_BLAZE) || defined(USE_EIGEN) 
   nobs = data.size();
   #endif   
   this->method = method;
   this->trend = trend;
   this->regression = regression;
}

//*************************************************************************************************

// get test statistic
template <typename T>
const T& UnitRoot<T>::get_stat() const
{
   return stat;
}

//*************************************************************************************************

// get test pvalue
template <typename T>
const T& UnitRoot<T>::get_pval() const
{
   return pval;
}

//*************************************************************************************************

// get test OLS regression results
template <typename T>
const OLS<T>& UnitRoot<T>::get_ols() const
{
   return *result;
}

//*************************************************************************************************

// get test valid trends
template <typename T>
const std::vector<std::string>& UnitRoot<T>::get_trends() const
{
   return valid_trends;
}

//*************************************************************************************************

// set pointer to the original data
template <typename T>
void UnitRoot<T>::set_data()
{
    ptr = &data;
}

//*************************************************************************************************

// set number of lags, checking input validity
// this function needs to know optim so it needs to be run after set_method()
template <typename T>
void UnitRoot<T>::set_lags()
{
   if (!optim) {
      // number of lags cannot be strictly negative
      if (lags < 0) {
         lags = 0;
         std::cout << "\n  WARNING: number of lags cannot be negative, it has been set to 0 by default.\n\n";
      }
      // if user has switched from a default lags value to a value of his choice (for all tests)
      if (!lags_type.empty() && lags != prev_lags) {
         lags_type = std::string();
      }
   } else {
      // max_lags cannot be strictly negative
      if (max_lags < 0) {
         max_lags = 0;
         std::cout << "\n  WARNING: maximum number of lags cannot be negative, it has been set to a default value (L12-rule).\n\n";
      }
      // if user has switched from a default max_lags value to a value of his choice or to a given number of lags (for ADF and DFGLS tests only) 
      if (!lags_type.empty() && (lags != prev_lags || max_lags != prev_max_lags)) {
         lags_type = std::string();
      }
      // if user has switched from a test with optimized lag length to a test for a given lag (for ADF and DFGLS tests)
      if (lags != prev_lags) {
         method = std::string();
         max_lags = 0;
         prev_max_lags = 0;
         optim = false;
      }
   }

   // computing max_lags value for lag length optimization for ADF and DFGLS tests if max_lags has not been set or if a negative value has been chosen by mistake, if the user has defined max_lags but lags_type too, lags_type will have the last word as its role is to reinitialize max_lags default value
   // or computing default lags value for KPSS and PP tests and (ADF and DFGLS too as an option)
   if ((optim && !max_lags) || !lags_type.empty()) {
      // short => L4-rule (Schwert 1989)
      if (lags_type == "short") {
         max_lags = int(4 * pow(0.01 * nobs, 0.25));
      }
      // long => L12-rule (Schwert 1989)
      else {
         max_lags = int(12 * pow(0.01 * nobs, 0.25));
      }
      // updating lags only for PP and KPSS tests, for ADF and DFGLS tests lags will be updated at the next optimization or set back to prev_lags if max_lags, trend, method and level are the same as before 
      lags = max_lags; 
   }

   // for all tests if the number of lags is still different than its previous value 
   if (!optim && lags != prev_lags) {
      new_test = true;
      new_lags = true;
      prev_lags = lags;
   }
   // only for ADF and DFGLS tests when chosing a max_lags value
   else if (optim && max_lags != prev_max_lags) {
      new_test = true;
      new_lags = true;
      prev_max_lags = max_lags;
   }
}

//*************************************************************************************************

// set lags type long or short for PP and KPSS default lags value or ADF and DFGLS default maxlags value
template <typename T>
void UnitRoot<T>::set_lags_type()
{
   // skipping this function if lags_type is empty or did not change
   if (lags_type.empty() || lags_type == prev_lags_type) {
      return;
   }
   if (lags_type != "long" && lags_type != "short") {
      std::cout << "\n  WARNING: unknown default type of lags, long has been selected by default.\n\n";
      // default lags type is long
      lags_type = "long";
   }
   prev_lags_type = lags_type;
}

//*************************************************************************************************

// set statistic level for optimal lags selection (for method T-STAT only)
template <typename T>
void UnitRoot<T>::set_level()
{
   // skipping this function if method is not T-STAT or level did not change
   if (level == prev_level) {
       return;
   }
   // the level cannot be negative
   if (level < 0) {
      level *= -1;
   }
   if (level != prev_level) {
      new_test = true;
      new_level = true;
      prev_level = level;
   }
}

//*************************************************************************************************

// set method for lag length optimization, checking input validity
template <typename T>
void UnitRoot<T>::set_method()
{
   // skipping this function if method is empty or did not change
   if (method.empty() || method == prev_method) {
      return;
   }
   // checking input method validity
   auto p = std::find(valid_methods.begin(), valid_methods.end(), method);
   // setting default method to AIC
   if (p == valid_methods.end()) {
      method = "AIC";
      std::cout << "\n  WARNING: unknown method for lag length optimization, AIC has been selected by default.\n\n";
      std::cout << "  Possible methods for this test are: " << valid_methods << "\n";
   }
   if (method != prev_method) {
      new_method = true;
      new_test = true;
   }
   // optim is a shorcut to know whether we are in presence of test with lag length optimization or not, it avoids calling everytime method.empty()
   optim = true;
}

//*************************************************************************************************

// set Information Criterion function and correction coefficient
// this function needs to know max_lags so it must be run after set_lags()
template <typename T>
void UnitRoot<T>::set_IC()
{   
   // skipping this function if method did not change
   if (!new_method) { 
      return;
   }
   // setting functor to Information Criterion function
   if (method[0] != 'M') {
      // classic IC
      ICfunc = &UnitRoot<T>::IC;
   } else {
      // modified IC
      ICfunc = &UnitRoot<T>::MIC;
   }
   // setting Information Criterion correction coefficient
   // Bayesian Information Criterion (Schwartz-Bayes)
   if (method == "BIC" || method == "MBIC") {
      ICcc = log(nobs - 1 - max_lags);
   }
   // Hannan-Quinn Information Criterion
   else if (method == "HQC" || method == "MHQC") {
      ICcc = 2.0 * log(log(nobs - 1 - max_lags));
   }
   // default: Akaike Information Criterion
   else {
      ICcc = 2.0;
   }
}

//*************************************************************************************************

// set numnber of iterations for bootstrap
template <typename T>
void UnitRoot<T>::set_niter()
{
   // checking number of iterations niter
   if (niter <= 0) {
      niter = 1000;
      std::cout << "\n  WARNING: in UnitRoot<T>::bootstrap(), number of iterations cannot be null or negative it has been set to 1000 by default.\n\n";
   } else if (niter < 1000) {
      std::cout << "\n  WARNING: in UnitRoot<T>::bootstrap(), number of iterations might be too small for computing p-value.\n\n";
   }
   if (niter != prev_niter) {
      new_niter = true;
      prev_niter = niter;
   }
}

//*************************************************************************************************

// set statistic type, checking input validity
template <typename T>
void UnitRoot<T>::set_test_type()
{
   // skipping this function if test_type did not change
   if (test_type == prev_test_type) {
      return;
   }
   // setting default test_type to tau
   if (test_type != "rho" && test_type != "tau") {
      test_type = "tau";
      std::cout << "\n  WARNING: unknown test type, tau has been selected by default.\n\n";
   }
   if (test_type != prev_test_type) {
      new_test = true;
      new_test_type = true;
      prev_test_type = test_type;
   }
}

//*************************************************************************************************

// set regression trend, checking input validity
template <typename T>
void UnitRoot<T>::set_trend()
{
   // skipping this function if trend did not change
   if (trend == prev_trend) {
      return;
   }
   // checking input trend validity
   auto p = std::find(valid_trends.begin(), valid_trends.end(), trend);
   // setting default trend to constant
   if (p == valid_trends.end()) {
      trend = "c";
      trend_type = "constant";
      npar = 2;
      std::cout << "\n  WARNING: unknown regression trend selected, regression with constant term has been selected by default.\n\n";
      std::cout << "  Possible trends for this test are " << valid_trends << "\n"; 
   } else {
      if (trend == "c") {
         trend_type = "constant";
         npar = 2;
      }
      else if (trend == "nc") {
         trend_type = "no constant";
         npar = 1;
      }
      else if (trend == "ct") {
         trend_type = "constant trend";
         npar = 3;
      }
      else if (trend == "ctt") {
         trend_type = "quadratic trend";
         npar = 4;
      }
   }
   if (trend != prev_trend) {
      new_test = true;
      new_trend = true;
      prev_trend = trend;
   }
}

//*************************************************************************************************

// OLS demeaning or detrending
template <typename T>
void UnitRoot<T>::ols_detrend()
{  
   // initializing w => constant case
   #ifdef USE_ARMA
   Matrix<T> w(nobs, 1, arma::fill::ones); 
   #elif USE_BLAZE
   Matrix<T> w(nobs, 1);
   w = forEach(w, [](T val){ return 1; });
   #elif USE_EIGEN
   Matrix<T> w = Matrix<T>::Ones(nobs, 1);
   #endif

   // constant trend case
   if (trend == "ct") {
      // appending new column to w
      #ifdef USE_ARMA
      w.insert_cols(1, arma::cumsum(w));
      #elif USE_BLAZE
      w.resize(nobs, w.columns() + 1);
      std::partial_sum(&w(0, 0), &w(nobs, 0), &w(0, 1));
      #elif USE_EIGEN
      w.conservativeResize(Eigen::NoChange, w.cols() + 1);
      w.col(w.cols() - 1) = Vector<T>::LinSpaced(w.rows(), 0, w.rows() - 1);
      #endif
   }
   else if (trend == "nc") {
      throw std::invalid_argument("\n  ERROR: in UnitRoot<T>::ols_detrend(), no detrending possible when regression trend set to no constant.\n\n");
   }

   // running OLS
   OLS<T> fit(*ptr, w);
   // getting detrended data
   z = *ptr - w * fit.param;
   // setting pointer to detrended data
   ptr = &z;
}

//*************************************************************************************************

// GLS demeaning or detrending
template <typename T>
void UnitRoot<T>::gls_detrend()
{
   T alpha;
   int nc;
   // detrend case
   if (trend == "ct") {
      alpha = 1.0 - 13.5 / nobs;
      nc = 2;
   }
   // demean case
   else if (trend == "c") { 
      alpha = 1.0 - 7.0 / nobs;
      nc = 1;
   } else {
      throw std::invalid_argument("\n  ERROR: in UnitRoot<T>::gls_detrend(), no detrending possible when regression trend set to no constant.\n\n");
   }

   Matrix<T> u(nobs, nc); 
   Vector<T> v(nobs);
   Matrix<T> w(nobs, nc);

   // initializing u, v and w
   u(0, 0) = 1;
   v[0] = ptr->operator[](0);
   w(0, 0) = 1;

   if (nc == 2) {
      w(0, 1) = 1;
      u(0, 1) = 1;
   }   
    
   for (int i = 1; i < nobs; ++i) {
      w(i, 0) = 1;
      u(i, 0) = 1 - alpha;

      if (nc == 2) {
         w(i, 1) = i + 1;
         u(i, 1) = w(i, 1) - alpha * i;
      }
      v[i] = ptr->operator[](i) - alpha * ptr->operator[](i - 1);
   } 

   // running OLS
   OLS<T> fit(v, u);
   // getting detrended data
   z = *ptr - w * fit.param;
   // setting pointer to detrended data
   ptr = &z;
}

//*************************************************************************************************

// run ADF test for a given lag
template <typename T>
void UnitRoot<T>::adf_regression()
{
   // number of lines
   int nr = nobs - lags - 1;
   // number of columns
   int nc = npar + lags;

   // checking if enough data are available to evaluate ADF model coefficients by OLS
   if (nr  < nc) {
      throw std::invalid_argument("\n  ERROR: in UnitRoot<T>::adf_regression(), more data required to compute ADF test for " + std::to_string(lags) + " lags, at least " + std::to_string(nc - nr) + " element(s) need(s) to be added or the number of lags to be reduced.\n\n");
   }

   // NB: the user might see this exception thrown again even after having increased the data dimension in the case of lag length optimization as max_lags is function of the data dimension, a second adjustment will be then necessary before the data dimension is accepted

   // allocating memory
   #ifdef USE_ARMA
   x.set_size(nr, nc);
   y.set_size(nr);
   #elif defined(USE_BLAZE) || defined(USE_EIGEN)
   x.resize(nr, nc);
   y.resize(nr);
   #endif

   for (int i = 0; i < nr; ++i) {
      // filling vector of dependent variable
      y[i] = ptr->operator[](i + lags + 1) - ptr->operator[](i + lags);
      // filling matrix of independent variables
      x(i, 0) = ptr->operator[](i + lags);
      // testing if model contains constant term or more
      if (npar >= 2) { 
         x(i, 1) = 1; 
         if (npar >= 3) {
            x(i, 2) = i + 1;
            if (npar == 4) {
               x(i, 3) = x(i, 2) * x(i, 2);
            }
         }
      }   
      // computing lag difference terms
      for (int j = npar; j < nc ; ++j) {
         x(i, j) = ptr->operator[](i - j + nc) - ptr->operator[](i - j + nc - 1);
      }  
   }

   // if OLS regression statistics are required
   if (!regression) {
      result = std::make_shared<OLS<T>>(y, x);
   } else {
      result = std::make_shared<OLS<T>>(y, x, true);
      prev_regression = true;
   }

   // updating test t-statistic
   stat = result->t_stat[0];
}

//*************************************************************************************************

// initialize dependent and independent variables in ADF test regression for lag length optimization
template <typename T>
void UnitRoot<T>::initialize_adf()
{
   // number of lines in y vector and x matrix
   this->nrows = nobs - 1;
   // number of columns
   int ncols = npar + max_lags;

   // checking if enough data are available to evaluate ADF model coefficients by OLS
   if (nrows - max_lags < ncols) {
      throw std::invalid_argument("\n  ERROR: in UnitRoot<T>::adf_regression(), more data required to compute ADF test for " + std::to_string(max_lags) + " lags, at least " + std::to_string(ncols - nrows + max_lags) + " element(s) need(s) to be added or the number of lags to be reduced.\n\n");
   }
    // NB: the user might see this exception thrown again even after having increased the data dimension in the case of lag length optimization as max_lags is function of the data dimension, a second adjustment will be then necessary before the data dimension is accepted

   // allocating memory
   #ifdef USE_ARMA
   x.set_size(nrows, ncols);
   y.set_size(nrows);
   #elif defined(USE_EIGEN) || defined(USE_BLAZE)
   x.resize(nrows, ncols);
   y.resize(nrows);
   #endif

   #ifdef _OPENMP
   #pragma omp parallel for num_threads(MAX_THREADS)
   #endif
   for (int i = 0; i < nrows; ++i) {
      // filling vector of dependent variable
      y[i] = ptr->operator[](i + 1) - ptr->operator[](i);
      // filling matrix of independent variables
      x(i, 0) = ptr->operator[](i);
      // testing if model contains constant term or more
      if (npar >= 2) { 
         x(i, 1) = 1; 
         if (npar >= 3) {
            x(i, 2) = i + 1;
            if (npar == 4) {
               x(i, 3) = x(i, 2) * x(i, 2);
            }
         }
      }
      // computing lag difference terms
      for (int j = 0; j < max_lags ; ++j) {
         if (j < i) {
            x(i, j + npar) = ptr->operator[](i - j) - ptr->operator[](i - j - 1);
         } else {
            x(i, j + npar) = 0;
         }
      }   
   }
}

//*************************************************************************************************

// optimize ADF test lag length using Information Criterion
template <typename T>
void UnitRoot<T>::optimize_lag()
{
   // if new method only (with an optimization previously been run and not being T-STAT) we use the same results as before without recomputing the tests for all lags from max_lags to 0, we only compute the new information criterion value
   if (new_method && !prev_method.empty() && prev_method[0] != 'T') {
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(MAX_THREADS)
      #endif
      for (int i = max_lags; i > -1; --i) {
         results[i]->IC = (this->*ICfunc)(results[i]);
      }
   }
   // if not we need to recompute all tests from max_lags to 0
   else {
      // clearing results vector of previous results
      results.clear();
      // allocating memory to results vector
      results.resize(max_lags + 1);
 
      // This for loop will be executed by a number of threads equal to MAX_THREADS or will ignore preprocessor directives if compiler does not support OpenMP or if OpenMP not called
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(MAX_THREADS)
      #endif
      // searching for numnber of lags minimizing IC
      for (int i = max_lags; i > -1; --i) {
         // selecting y sub-vector and x sub-matrix
         #ifdef USE_ARMA
       	Vector<T> ysub(y.memptr() + i, nrows - i, false);
         Matrix<T> xsub(x.submat(i, 0, nrows - 1, npar + i - 1));
         #elif USE_BLAZE
         Vector<T> ysub(subvector(y, i, nrows - i));
         Matrix<T> xsub(submatrix(x, i, 0, nrows - i, npar + i));
         #elif USE_EIGEN
         Vector<T> ysub(y.segment(i, nrows - i));
         Matrix<T> xsub(x.block(i, 0, nrows - i, npar + i));
         #endif
         // adding test result for current number of lags to results vector      
         results[i] = std::make_shared<OLS<T>>(ysub, xsub);         
         // recording number of lags for current test
         results[i]->lags = i;
         // computing new Information Criterion for other methods
         results[i]->IC = (this->*ICfunc)(results[i]);
      }
   }

   // getting ADF test OLS result with smallest Information Criterion
   auto best_ols = std::min_element(results.begin(), results.end(), [](const std::shared_ptr<OLS<T>>& ols1, const std::shared_ptr<OLS<T>>& ols2){return T(ols1->IC) < T(ols2->IC);});

   result = *best_ols;
   // updating number of lags
   lags = result->lags;

   // if OLS regression statistics are required
   if (regression) {
      // selecting y sub-vector and x sub-matrix
      #ifdef USE_ARMA
      Vector<T> ysub(y.memptr() + lags, nrows - lags, false);
      Matrix<T> xsub(x.submat(lags, 0, nrows - 1, npar + lags - 1));
      #elif USE_BLAZE
      Vector<T> ysub(subvector(y, lags, nrows - lags));
      Matrix<T> xsub(submatrix(x, lags, 0, nrows - lags, npar + lags)); 
      #elif USE_EIGEN
      Vector<T> ysub(y.segment(lags, nrows - lags));
      Matrix<T> xsub(x.block(lags, 0, nrows - lags, npar + lags));
      #endif
      // computing OLS regression statistics
      result->get_stats(ysub, xsub);
      prev_regression = true;
   }

   // updating test t-statistic
   stat = result->t_stat[0];
   // updating previous lags value;
   prev_lags = lags;
   // updating previous method
   prev_method = method;
}

//*************************************************************************************************

// NB: for better comparison when optimizing lags using Information Criterion, all calculations are made with the same dimension: nobs - 1 - maxlags

// compute Information Criterion
template <typename T>
T UnitRoot<T>::IC(const std::shared_ptr<OLS<T>>& res)
{
    // number of regressors is:
    // number of lags difference terms + 1
    // plus 1 for the model with constant
    // plus 2 for the model with constant trend
    // plus 3 for the model with quadratic trend
  
   #ifdef USE_ARMA
   // getting number of regressors
   int k = res->param.n_elem;
   // getting number of residuals
   int n = res->resid.n_elem;
   // getting residuals sub-vector 
   Vector<T> z(res->resid.memptr() + max_lags - res->lags, n - max_lags + res->lags, false);
   #elif defined(USE_BLAZE) || defined(USE_EIGEN) 
   // getting number of regressors
   int k = res->param.size();
   // getting number of residuals
   int n = res->resid.size();
   // getting residuals sub-vector
   #ifdef USE_BLAZE
   Vector<T> z(subvector(res->resid, max_lags - res->lags, n - max_lags + res->lags));
   #elif USE_EIGEN
   Vector<T> z(res->resid.segment(max_lags - res->lags, n - max_lags + res->lags));
   #endif
   #endif
   // computing factor
   T factor = 1.0 / (n - max_lags + res->lags);
   // computing residuals variance
   #ifdef USE_ARMA
   T sigma2 = arma::as_scalar(z.t() * z) * factor;
   #elif USE_BLAZE
   T sigma2 = (blaze::trans(z) * z) * factor;
   #elif USE_EIGEN
   T sigma2 = z.dot(z) * factor;
   #endif
   // computing criterion
   return log(sigma2) + ICcc * k * factor;
}

//*************************************************************************************************

// compute Modified Information Criterion
template <typename T>
T UnitRoot<T>::MIC(const std::shared_ptr<OLS<T>>& res)
{ 
   #ifdef USE_ARMA
   // getting number of residuals
   int n = res->resid.n_elem;
   // getting residuals sub-vector
   Vector<T> z(res->resid.memptr() + max_lags - res->lags, n - max_lags + res->lags, false);
   #elif defined(USE_BLAZE) || defined(USE_EIGEN)
   // getting number of residuals
   int n = res->resid.size();
   // getting residuals sub-vector
   #elif USE_BLAZE
   Vector<T> z(subvector(res->resid, max_lags - res->lags, n - max_lags + res->lags));
   #ifdef USE_EIGEN
   Vector<T> z(res->resid.segment(max_lags - res->lags, n - max_lags + res->lags));
   #endif
   #endif
   // computing factor
   T factor = 1.0 / (n - max_lags + res->lags);
   #ifdef USE_ARMA
   // computing residuals variance
   T sigma2 = arma::as_scalar(z.t() * z) * factor;
   // getting data lagged 1 term sub-vector
   Vector<T> y(ptr->memptr() + max_lags, nobs - max_lags - 1, false);
   // computing y^2
   T y2 = arma::as_scalar(y.t() * y);
   #elif USE_BLAZE
   // computing residuals variance
   T sigma2 = (blaze::trans(z) * z) * factor;
   // getting data lagged 1 term sub-vector
   Vector<T> y(subvector(*ptr, max_lags, nobs - max_lags - 1));
   // computing y^2
   T y2 = blaze::trans(y) * y; 
   #elif USE_EIGEN
   // computing residuals variance
   T sigma2 = z.dot(z) * factor;
   // getting data lagged 1 term sub-vector
   Vector<T> y(ptr->segment(max_lags, nobs - max_lags - 1));
   // computing y^2
   T y2 = y.transpose() * y;
   #endif
   // computing tau
   T tau = (res->param[0] * res->param[0] * y2) / sigma2; 
   // computing criterion
   return log(sigma2) + ICcc * (tau + res->lags) * factor;
}

//*************************************************************************************************

// select optimal ADF test lag using lags difference terms t-statistics
template <typename T>
void UnitRoot<T>::select_lag()
{
   bool opt_lag_found = false;

   // clearing results vector of previous results
   results.clear();
   // allocating memory to results vector
   results.resize(max_lags + 1);

   // searching for optimal lags
   for (int i = max_lags; i > -1; --i) {
      // selecting y sub-vector and x sub-matrix
      #ifdef USE_ARMA
      Vector<T> ysub(y.memptr() + i, nrows - i, false);
      Matrix<T> xsub(x.submat(i, 0, nrows - 1, npar + i - 1));
      #elif USE_BLAZE
      Vector<T> ysub(subvector(y, i, nrows - i));
      Matrix<T> xsub(submatrix(x, i, 0, nrows - i, npar + i)); 
      #elif USE_EIGEN
      Vector<T> ysub(y.segment(i, nrows - i));
      Matrix<T> xsub(x.block(i, 0, nrows - i, npar + i));
      #endif
      // adding test result for current number of lags to results vector      
      results[i] = std::make_shared<OLS<T>>(ysub, xsub);
      // stopping at the first significant lags difference term for T-STAT method
      if (fabs(results[i]->t_stat[i]) > level) {
         // updating result
         result = results[i];
         // if OLS regression statistics are required
         if (regression) {
            result->get_stats(ysub, xsub);
            prev_regression = true;
         }
         // recording number of lags 
         result->lags = i;
         opt_lag_found = true;
         // breaking the loop
         i = -1; 
      }
   }

   if (!opt_lag_found) {
      std::cout << "\n  WARNING: no optimal number of lags found with method T-STAT and level " << level << " please try again with a (positive) lower level.\n\n";
      // updating result
      result = results[0];
      // if OLS regression statistics are required
      if (regression) {
         // selecting y sub-vector and x sub-matrix
         #ifdef USE_ARMA
         Vector<T> ysub(y.memptr(), nrows, false);
         Matrix<T> xsub(x.submat(0, 0, nrows - 1, npar - 1));
         #elif USE_BLAZE
         Vector<T> ysub(subvector(y, 0, nrows - 0));
         Matrix<T> xsub(submatrix(x, 0, 0, nrows - 0, npar + 0));
         #elif USE_EIGEN
         Vector<T> ysub(y.segment(0, nrows));
         Matrix<T> xsub(x.block(0, 0, nrows, npar));
         #endif
         result->get_stats(ysub, xsub);
         prev_regression = true;
      }
      // recording number of lags 
      result->lags = 0;
   }

   // updating number of lags
   lags = result->lags;
   // updating test t-statistic
   stat = result->t_stat[0];
   // updating previous lags value;
   prev_lags = lags;   
   // updating previous method
   prev_method = method; 
}

//*************************************************************************************************

// compute ADF tests 
template <typename T>
void UnitRoot<T>::compute_adf()
{
   // ADF test with lag length optimization
   if (optim) {
      // declaring functor to lag optimization function
      void (UnitRoot<T>::*optim_func)();

      if (method[0] != 'T') {
         // setting information criterion (for all methods excepted T-STAT)
         set_IC();
         // optimizing lag length
         optim_func = &UnitRoot<T>::optimize_lag;
      } else {
         // setting statistic level for optimal lag selection (for T-STAT method only)
         set_level();
         // selecting optimal lag length
         optim_func = &UnitRoot<T>::select_lag;
      }
      // if new trend or new max_lags (for all methods)
      if (new_trend || new_lags) {
         // computing vector and matrix for OLS 
         initialize_adf();
         // turning new_method to false to force a new optimization as trend or lags is new
         new_method = false;
         // optimizing lag
         (this->*optim_func)();             
         new_trend = false;
         new_lags = false;
      }
      // if new method only (for all methods)
      else if (new_method) { 
         // optimizing lag
         (this->*optim_func)();
         new_method = false;
      }
      // if new level (for T-STAT method only)
      else if (method[0] == 'T' && new_level) {
         // selecting optimal lag length
         (this->*optim_func)();
         new_level = false;
      }
      // if the only change is regression set to true 
      else if (regression && !prev_regression) {
         // selecting y sub-vector and x sub-matrix
         #ifdef USE_ARMA
         Vector<T> ysub(y.memptr() + lags, nrows - lags, false);
         Matrix<T> xsub(x.submat(lags, 0, nrows - 1, npar + lags - 1));
         #elif USE_BLAZE
         Vector<T> ysub(subvector(y, lags, nrows - lags));
         Matrix<T> xsub(submatrix(x, lags, 0, nrows - lags, npar + lags)); 
         #elif USE_EIGEN
         Vector<T> ysub(y.segment(lags, nrows - lags));
         Matrix<T> xsub(x.block(lags, 0, nrows - lags, npar + lags));
         #endif
         result->get_stats(ysub, xsub);
         prev_regression = true;
      }
      // if user did not modify any parameter we set lags back to prev_lags as lags has become max_lags in set_lags()
      else {
         lags = prev_lags;
      }
   }
   // ADF test for a given lag
   else {
      // if new trend or new lags or new data from boostrap
      if (new_trend || new_lags || new_data) {
         adf_regression();
         new_trend = false;
         new_lags = false;
      }
      // if the only change is regression set to true 
      else if (regression && !prev_regression) {
         result->get_stats(y, x);
         prev_regression = true;
      } 
   }
}

//*************************************************************************************************

template <typename T>
void UnitRoot<T>::run_bootstrap()
{  
   // optim has to be set to false for avoiding optimization in case it was set to true
   bool saved_optim = optim;
   optim = false;
   // saving result from OLS and test stat value (they will be modified during bootstrap)
   T saved_stat = stat;
   OLS<T> saved_result = *result;   
 
   // computing centred residuals from original test
   #ifdef USE_ARMA 
   Vector<T> eps = result->resid - arma::mean(result->resid);
   #elif USE_BLAZE
   Vector<T> eps(result->resid);
   T m = std::accumulate(&eps[0], &eps[eps.size()], 0.0) / eps.size();
   eps = forEach(eps, [&m](const T& val){ return val - m; });  
   #elif USE_EIGEN
   T m = result->resid.mean();
   Vector<T> eps = result->resid - Vector<T>::Constant(result->resid.size(), 1, m);   
   #endif 

   // computing vector of lags difference terms u dimension
   int nu = nobs - 1;

   Vector<T> u, delta;

   if (lags > 0) {
      #ifdef USE_ARMA 
      // initializing vector of lags difference terms
      u = ptr->subvec(1, lags - 1) - ptr->subvec(0, lags - 2);
      // taking lags difference terms coefficients from original test in reverse order
      delta = arma::flipud(result->param.subvec(npar, result->param.size() - 1));
      #elif USE_BLAZE
      // initializing vector of lags difference terms
      u = subvector(*ptr, 1, lags - 1) - subvector(*ptr, 0, lags - 1); 
      // taking lags difference terms coefficients from original test in reverse order
      delta = subvector(result->param, npar, result->param.size() - npar);  
      std::reverse(&delta[0], &delta[delta.size()]);
      #elif USE_EIGEN
      // initializing vector of lags difference terms
      u = ptr->segment(1, lags - 1) - ptr->segment(0, lags - 1); 
      // taking lags difference terms coefficients from original test in reverse order
      delta = result->param.segment(npar, result->param.size() - npar).reverse();
      #endif    
   }
   #if defined(USE_ARMA) || defined(USE_BLAZE)
   u.resize(nu);
   #elif USE_EIGEN
   u.conservativeResize(nu, Eigen::NoChange);
   #endif

   // allowing test statistic recalculation even if all parameters remain identical
   new_data = true;

   // initializing new path 
   #ifdef USE_ARMA 
   Vector<T> new_path(ptr->memptr(), lags + 1, false);
   new_path.resize(nobs);
   #elif USE_BLAZE
   Vector<T> new_path(subvector(*ptr, 0, lags + 1));
   new_path.resize(nobs);
   #elif USE_EIGEN
   Vector<T> new_path(ptr->segment(0, lags + 1));
   new_path.conservativeResize(nobs, Eigen::NoChange);
   #endif
    
   // vector of new statistics
   std::vector<T> stats(niter);

   // initializing Boost random uniform (integer) numbers generator
   boost::uniform_int<> udistrib(0, eps.size() - 1); 
   boost::variate_generator<boost::mt19937&, boost::uniform_int<>> runif(rng, udistrib);

   T term = 0;

   // bootstrapping Unit-Root test, computing new paths
   for(int k = niter; k--; ) { 
      for (int i = lags; i < nu; ++i) {
         switch (lags) {
            case 0:
               break;
            default:
               #ifdef USE_ARMA 
               term = arma::as_scalar(arma::Col<T>(u.memptr() + i - lags, lags, false).t() * delta);
               #elif USE_BLAZE
               term = blaze::trans(subvector(u, i - lags, lags)) * delta;
               #elif USE_EIGEN 
               term = u.segment(i - lags, lags).transpose() * delta;          
               #endif 
         }
         u[i] = term + eps[runif()];
         new_path[i + 1] = new_path[i] + u[i];
      }

      // setting ptr to new_path
      ptr = &new_path;
      // computing new test statistic
      stats[k] = statistic();
   }

   // sorting stats
   std::sort(stats.data(), stats.data() + niter);

   // computing critical values from bootstrap
   for (int i = 0; i < 15; ++i) {
      critical_values[i] = .5 * (stats[floor((niter - 1) * probas[i])] + stats[floor((niter - 1) * probas[i]) + 1]);
   }

   // locking back test statistic recalculation if all parameters remain identical
   new_data = false;
   // setting pointer back to original data
   set_data();
   // getting back original test parameters and results
   optim = saved_optim;
   stat = saved_stat;
   result = std::make_shared<OLS<T>>(saved_result);

    // NB: if we want to obtain the same critical values as shown in most of tables we will use gaussian random numbers, however this is wrong as we do not know the residual distribution, then the best way of obtaining the true critical values is to use the model residuals and shuffle them to obtain new paths. However this method is valid when the sample size is large enough, indeed when too small the residuals diversity will be small and the critical values estimations will poor when using lags.
}

//*************************************************************************************************

// compute critical values from probabilities
template <typename T>
void UnitRoot<T>::compute_cv()
{
   // computing adjusted number of observations
   int n = nobs - lags - 1;

   for (int i = 0; i < 15; ++i) {
      // computing critical value
      critical_values[i] = 0;

      int n0 = coeff_ptr->at(probas[i]).at(0).size();

      for (int j = 0; j < n0; ++j) {
         critical_values[i] += coeff_ptr->at(probas[i]).at(0).at(j) / pow(n, j);
      }

      int n1 = coeff_ptr->at(probas[i]).at(1).size();

      for (int j = 0; j < n1; ++j) {
         critical_values[i] += coeff_ptr->at(probas[i]).at(1).at(j) * pow( T(lags) / n, j + 1);
      }
   }
}

//*************************************************************************************************

// compute p-value by linear interpolation from critical values
template <typename T>
void UnitRoot<T>::compute_pval()
{
   // if stat is smaller than critical value for first probability (in absolute value)
   if (stat <= critical_values[0]) {
      pval = probas[0];
   } else {
      for (int i = 1; i < 15; ++i) {
         if (stat <= critical_values[i]) {
            pval = probas[i - 1] + (stat - critical_values[i - 1]) * (probas[i] - probas[i - 1]) / (critical_values[i] - critical_values[i - 1]);
            break;
         }
      }
   }

   // if stat is greater than critical value for last probability in absolute value
   if (stat > critical_values[14]) {
      pval = probas[14];
   }
}

//*************************************************************************************************

// compute p-value
template <typename T>
const T& UnitRoot<T>::pvalue()
{
   if (std::isnan(stat)) {
      pval =  -std::nan("");
   } 
   else if (!bootstrap) {
      // if a new test has been run or p-value was previously computed by bootstrap
      if (new_test || prev_bootstrap) {
         // computing critical values
         compute_cv();
         // computing p-value
         compute_pval();

         new_test = false;
      }
      // NB: if test parameters remain identical but p-value was previously computed using bootstrap, p-value will be recomputed using the critical values coefficients
   } 
   else {
      // setting number of iterations
      set_niter();

      // if a new test has been run or new number of iterations or p-value was not previously computed by bootstrap
      if (new_test || new_niter || !prev_bootstrap) {
         // running bootstrap (computing critical values)
         run_bootstrap();
         // computing p-value
         compute_pval();

         new_test = false;
         new_niter = false;
      }
   }
   // NB: new_test avoids p-value recomputation when test parameters have remained identical
   prev_bootstrap = bootstrap;

   return pval;
}

//*************************************************************************************************

// output test results
template <typename T>
void UnitRoot<T>::show()
{
   // setting precision for floating numbers output
   std::cout << std::fixed << std::setprecision(3);

   // for PP test only
   std::string stat_type;

   if (test_type == "tau")  {      
      stat_type = " (Z-tau)";
   }
   else if (test_type == "rho") { 
     stat_type = " (Z-rho)";
   }

   // outputting test name
   std::cout << "\n  " + test_name + " Test Results" + stat_type + "\n";
   std::cout << "  ====================================\n";

   // outputting test statistic
   std::cout << "  Statistic" << std::setw(27) << stat << "\n";

   // outputting p-value
   std::string s;
   (bootstrap) ? s = " (*)\n" : s = "\n";

   std::cout << "  P-value";

   if (pval <= probas[0]) {
      std::cout << std::setw(29) << "< 0.001";
   }
   else if (pval >= probas[14]) {
      std::cout << std::setw(29) << "> 0.999";
   } else {
      std::cout << std::setw(29) << pval;
   }
   std::cout << s;

   // if lag length has been optimized (for ADF or DFGLS)
   if (optim) {
      std::cout << "  Optimal Lags" << std::setw(24) << lags << "\n";

      if (method == "T-STAT") {
         std::cout << "  Method" << std::setw(30) << method << "\n";
      } else {
         std::cout << "  Criterion" << std::setw(27) << method << "\n";
      }
   } else {    
      std::cout << "  Lags" << std::setw(32) << lags << "\n";
   }

   std::cout << "  Trend" << std::setw(31) << trend_type << "\n"; 
   std::cout << "  ------------------------------------\n\n";

   // outputting test hypothesis
   std::cout << "  Test Hypothesis\n";
   std::cout << "  ------------------------------------\n";

   if (test_name == "KPSS") {
      std::cout << "  H0: The process is weakly stationary" << "\n";
      std::cout << "  H1: The process contains a unit root" << "\n";
   } else {
      std::cout << "  H0: The process contains a unit root" << "\n";
      std::cout << "  H1: The process is weakly stationary" << "\n";
   }
   std::cout << "\n";

   // outputting critical values
   std::cout << "  Critical Values" << s;
   std::cout << "  ---------------\n";

   // declaring vector of critical value indexes
   std::vector<int> idx;

   // KPSS is a one-sided right-tailed test so we need to select different critical values than ADF, DFGLS and PP which are one-sided left-tailed tests
   if (test_name == "KPSS") {
      idx = {12, 10, 9};
   } 
   else {
      idx = {2, 4, 5};
   }
  
   if (std::isnan(stat)) {
      std::cout << "   1% " << std::setw(11) << -std::nan("") << "\n";
      std::cout << "   5% " << std::setw(11) << -std::nan("") << "\n";
      std::cout << "  10% " << std::setw(11) << -std::nan("") << "\n";
   } 
   else {
      std::cout << "   1% " << std::setw(11) << critical_values[idx[0]] << "\n";
      std::cout << "   5% " << std::setw(11) << critical_values[idx[1]] << "\n";
      std::cout << "  10% " << std::setw(11) << critical_values[idx[2]] << "\n";
   }

   if (bootstrap) {
      std::cout << "\n  (*) computed by bootstrap\n";
   }
   std::cout << "\n";

   // outputting test conclusion
   std::cout << "  Test Conclusion\n";
   std::cout << "  ---------------\n";

   if (pval <= 0.01) {
      std::cout << "  We can reject H0 at the 1% significance level\n";
   }
   else if (pval <= 0.05) {
      std::cout << "  We can reject H0 at the 5% significance level\n";
   }
   else if (pval <= 0.10) {
      std::cout << "  We can reject H0 at the 10% significance level\n";
   } 
   else if (!std::isnan(pval)) {
      std::cout << "  We cannot reject H0\n";
   } 
   else {
      std::cout << "  We cannot conclude, nan produced\n";
   }
   std::cout << "\n";

   // outputting OLS results
   if (regression) {
       result->show();    
   }
}

//*************************************************************************************************

// << operator overload for outputting results
template <typename T>
std::ostream& operator<<(std::ostream& out, urt::UnitRoot<T>& test)
{
   test.show();

   return out;
}

//=================================================================================================

}
