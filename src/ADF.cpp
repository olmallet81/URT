//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#include "../include/URT.hpp"

namespace urt {

//=================================================================================================

// parameter constructor for computing ADF test for a given number of lags
template <typename T>
ADF<T>::ADF(const Vector<T>& data, int lags, const std::string& trend, bool regression) : ur(data, lags, trend, regression)
{
   ur::test_name = _test_name;
   ur::valid_trends = _valid_trends;
}

//*************************************************************************************************

// parameter constructor for computing ADF test with lag length optimization
template <typename T>
ADF<T>::ADF(const Vector<T>& data, const std::string& method, const std::string& trend, bool regression) : ur(data, method, trend, regression)
{
   ur::test_name = _test_name;
   ur::valid_trends = _valid_trends;
}

//*************************************************************************************************

// compute test statistic
template <typename T>
const T& ADF<T>::statistic()
{
   // setting type of lags (if a default type of lags value has been chosen)
   ur::set_lags_type();
   // setting optimization method
   ur::set_method();
   // setting number of lags
   ur::set_lags();
   // setting regression trend
   ur::set_trend();
   // computing ADF test
   ur::compute_adf();

   return ur::stat;
}

//*************************************************************************************************

// compute test p-value
template <class T>
const T& ADF<T>::pvalue()
{
   // computing test statistic
   this->statistic();
   // setting critical values coefficients pointer
   ur::coeff_ptr = &coeff_adf.at(ur::trend);
   // computing p-value
   ur::pvalue();

   return ur::pval;
} 

//*************************************************************************************************

// output test results
template <class T>
void ADF<T>::show()
{
   // in case user modified test type, for ADF it should always be empty
   ur::test_type = std::string();  
   // computing p-value
   this->pvalue();
   // outputting results
   ur::show();
}

//=================================================================================================

}

