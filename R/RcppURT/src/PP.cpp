//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#include "../include/URT.hpp"

namespace urt { 

//=================================================================================================

// parameter constructor for computing PP test for a given number of lags
template <typename T>
PP<T>::PP(const Vector<T>& data, int lags, const std::string& trend, const std::string& test_type, bool regression) : ur(data, lags, trend, regression)
{   
   ur::test_name = _test_name;
   this->lags = lags;
   ur::test_type = test_type;
   ur::valid_trends = _valid_trends;
}

//*************************************************************************************************

// parameter constructor for computing PP test for a default number of lags (long or short)
template <typename T>
PP<T>::PP(const Vector<T>& data, const std::string& lags_type, const std::string& trend, const std::string& test_type, bool regression) : ur(data, 0, trend, regression)
{
   ur::test_name = _test_name;
   ur::test_type = test_type;
   ur::valid_trends = _valid_trends;
   ur::lags_type = lags_type;
}

//*************************************************************************************************

// compute test statistic
template <class T>
void PP<T>::compute_stat()
{
   // number of residuals
   int n = ur::nobs - 1;
   // computing gamma vector for PP test statistics computation
   Vector<T> G(lags + 1);

#ifdef USE_ARMA
   // pointer to test residual elements memory
   T* ptr = ur::result->resid.memptr();
   for (int j = 0; j <= lags; ++j) {
      G[j] = arma::as_scalar(Vector<T>(ptr + j, n - j, false).t() * Vector<T>(ptr, n - j, false));
   }
#elif USE_BLAZE
   for (int j = 0; j <= lags; ++j) {
      G[j] = blaze::trans(subvector(ur::result->resid, j, n - j)) * subvector(ur::result->resid, 0, n - j);
   }
#elif USE_EIGEN
   for (int j = 0; j <= lags; ++j) {
      G[j] = ur::result->resid.segment(j, n - j).transpose() * ur::result->resid.segment(0, n - j);
   }
#endif

   G /= n; 

   T term = 0;
   // computing long-run variance estimator (Newey-West)
   for (int j = 1; j <= lags; ++j) {
      term += (1.0 - j / (lags + 1.0)) * G[j]; 
   }
   T L2 = G[0] + 2 * term;
   T rho = ur::result->param[0];
   // OLS mean square error
   T S2 = ur::result->MSE;
   // OLS variance of rho
   T sigma2 = ur::result->var[0];
   // computing PP test statistic
   if (ur::test_type == "rho") {
      // normalized statistic
      ur::stat = n * rho - .5 * n * n * sigma2 * (L2 - G[0]) / S2;
   } else {
      // t-statistic
      ur::stat = sqrt(G[0] / L2) * rho / sqrt(sigma2) - .5 * (L2 - G[0]) * n * sqrt(sigma2 / (L2 * S2));
   }
}

//*************************************************************************************************

// compute PP test statistic
template <typename T>
const T& PP<T>::statistic()
{
   // setting type of lags (if a default type of lags value has been chosen)
   ur::set_lags_type();
   // setting number of lags
   ur::lags = this->lags;
   ur::set_lags();
   this->lags = ur::lags;
   // this test is not augmented, ur::lags must always be 0 in ur::adf_regression() and ur::bootstrap()
   ur::lags = 0;
   // setting regression trend 
   ur::set_trend();
   // setting test type tau or rho
   ur::set_test_type();
   // if new trend or bootstrap
   if (ur::new_trend || ur::new_data) { 
      // computing ADF regression
      ur::adf_regression();
      // computing test statistic
      this->compute_stat();
      ur::new_trend = false;
   }
   // if new lags, new maxlags or new test type we only compute the statistic
   else if (ur::new_lags || ur::new_test_type) {
      // computing PP test statistic
      this->compute_stat();
      ur::new_lags = false;
      ur::new_test_type = false;
   }
   return ur::stat;
}

//*************************************************************************************************

// compute PP test p-value
template <typename T>
const T& PP<T>::pvalue()
{
   // computing test statistic
   this->statistic();
   // setting critical values coeff pointer
   if (ur::test_type == "tau") {
      ur::coeff_ptr = &coeff_pptau.at(ur::trend);
   } else {
      ur::coeff_ptr = &coeff_pprho.at(ur::trend);
   }
   // computing p-value
   ur::pvalue();
   return ur::pval;
}

//*************************************************************************************************

// output PP test results
template <typename T>
void PP<T>::show()
{
   // in case user modified method for PP it should always be empty
   ur::method = std::string();
   // computing p-value
   this->pvalue();
   // setting ur::lags from 0 to back to lags
   ur::lags = this->lags;
   // outputting results
   ur::show();
   // setting ur::lags back to 0 (in case bootstrap is run next)
   ur::lags = 0;
}

//=================================================================================================

}


